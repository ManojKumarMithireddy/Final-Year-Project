"""
Quantum router — Grover POC (local simulator) endpoints.
"""

import time
import datetime
from typing import Optional

from fastapi import APIRouter, HTTPException, Depends
from qiskit import transpile
from qiskit_aer import AerSimulator
from qiskit_aer.noise import NoiseModel, depolarizing_error

from db import get_db
from auth import get_current_user, get_optional_user
from models.schemas import QuantumToyRequest, BioQuantumRequest
from services.grover import build_grover_circuit, build_grover_step_circuits
from services.bio_grover import (
    fetch_patient_dna, apply_brca1_mutation, get_marker_seq,
    build_patient_nodes, run_bio_grover_local,
    build_bio_step_circuits,
    encode_dna, decode_bits,
    MARKER_GENE_NAME, MARKER_VARIANT, MARKER_REGION_DESC, PATIENT_ACCESSION,
)

router = APIRouter(prefix="/api/search", tags=["quantum"])


def _get_timestamp(client_ts: Optional[str]) -> datetime.datetime:
    """Use browser-supplied ISO timestamp if valid, else fall back to server UTC."""
    if client_ts:
        try:
            return datetime.datetime.fromisoformat(client_ts.replace("Z", "+00:00"))
        except (ValueError, TypeError):
            pass
    return datetime.datetime.now(datetime.timezone.utc)


@router.post("/quantum-simulation-poc")
async def quantum_toy(request: QuantumToyRequest, current_user: Optional[dict] = Depends(get_optional_user)):
    """
    Executes a Grover search circuit on the local Aer simulator (Proof of Concept).

    The 'target_bits' field is expected to be OTP-XOR-encrypted by the client
    as part of the PIR demonstration. An optional depolarizing noise model simulates
    real NISQ hardware gate errors.

    Scale note: A 4-qubit circuit searches 2^4=16 states. Real genomic-scale Grover's
    (for ~100,000 windows) would require ~17 qubits with fault-tolerant error correction —
    well beyond current NISQ-era hardware capabilities.
    """
    target_bitstring = request.target_bits
    if not target_bitstring or not all(c in '01' for c in target_bitstring):
        target_bitstring = "101"

    qc, iterations = build_grover_circuit(target_bitstring)
    n_qubits = len(target_bitstring)

    simulator = AerSimulator()
    transpiled_qc = transpile(qc, simulator)

    # Apply optional depolarizing noise model to simulate real hardware
    noise_model = None
    if request.noise_level > 0:
        noise_model = NoiseModel()
        err_1q = depolarizing_error(request.noise_level, 1)
        err_2q = depolarizing_error(min(request.noise_level * 10, 0.9), 2)
        err_3q = depolarizing_error(min(request.noise_level * 15, 0.9), 3)
        noise_model.add_all_qubit_quantum_error(err_1q, ['h', 'x', 'id', 's', 'sdg', 't', 'tdg'])
        noise_model.add_all_qubit_quantum_error(err_2q, ['cx', 'cz'])   # 2-qubit gates
        noise_model.add_all_qubit_quantum_error(err_3q, ['ccx'])         # 3-qubit Toffoli gate

    start_time = time.perf_counter()
    result = simulator.run(transpiled_qc, shots=1024, noise_model=noise_model).result()
    exec_time_ms = (time.perf_counter() - start_time) * 1000

    counts = result.get_counts()
    measured_state = max(counts, key=counts.get)
    circuit_diagram = str(qc.draw(output="text"))

    result_data = {
        "measured_state": measured_state,
        "counts": counts,
        "execution_time_ms": exec_time_ms,
        "iterations": iterations,
        "qubits": n_qubits,
        "circuit_diagram": circuit_diagram,
        "noise_level": request.noise_level,
        "step_circuits": build_grover_step_circuits(target_bitstring),
    }

    if current_user:
        db = get_db()
        await db.history.insert_one({
            "email": current_user["sub"],
            "type": "quantum_local_sim",
            "target_bits": request.target_bits,
            "measured_state": measured_state,
            "execution_time_ms": exec_time_ms,
            "noise_level": request.noise_level,
            "timestamp": _get_timestamp(request.client_timestamp),
        })

    return result_data


@router.get("/circuit-info")
async def circuit_info(target_bits: str = "0000"):
    """
    Returns the Grover circuit diagram and step-by-step data for a given
    target bitstring WITHOUT running a simulation. Used by the History page
    to visualize past results on demand.
    """
    if not target_bits or not all(c in '01' for c in target_bits) or len(target_bits) > 18:
        raise HTTPException(status_code=400, detail="target_bits must be 1-18 characters of 0/1.")

    qc, iterations = build_grover_circuit(target_bits)
    circuit_diagram = str(qc.draw(output="text", fold=-1))
    step_circuits = build_grover_step_circuits(target_bits)

    return {
        "circuit_diagram": circuit_diagram,
        "step_circuits": step_circuits,
        "iterations": iterations,
        "qubits": len(target_bits),
    }


@router.get("/quantum-poc/bio-marker")
async def get_bio_marker(n_codons: int = 2, force: bool = False, offset: int = 0):
    """
    Lightweight endpoint — returns the current disease marker info without running a simulation.
    force=true clears the NCBI sequence cache and re-fetches from NCBI.
    Also returns nearby_windows (offsets -5..+5) so the frontend can cycle locally.
    """
    reference_seq = fetch_patient_dna(force=force)
    if not reference_seq:
        raise HTTPException(status_code=502, detail="Failed to fetch BRCA1 reference from NCBI. The NCBI service may be temporarily unavailable.")
    mutant_seq  = apply_brca1_mutation(reference_seq)
    marker_seq  = get_marker_seq(mutant_seq, n_codons, offset=offset)
    target_bits = encode_dna(marker_seq)
    # Pre-compute windows -5..+5; flag windows whose bits appear in the healthy reference
    # (in_reference=True → non-specific marker, found in both patients → false positive).
    ref_bits = {nd["bits"] for nd in build_patient_nodes(reference_seq, n_codons)}
    nearby_windows = []
    for o in range(-5, 6):
        wdna  = get_marker_seq(mutant_seq, n_codons, offset=o)
        wbits = encode_dna(wdna)
        nearby_windows.append({
            "offset":       o,
            "dna":          wdna,
            "bits":         wbits,
            "in_reference": wbits in ref_bits,
        })
    return {
        "marker_dna":        marker_seq,
        "marker_bits":       target_bits,
        "marker_variant":    MARKER_VARIANT,
        "marker_region":     MARKER_REGION_DESC,
        "marker_gene":       MARKER_GENE_NAME,
        "patient_accession": PATIENT_ACCESSION,
        "n_codons":          n_codons,
        "n_qubits":          n_codons * 6,
        "offset":            offset,
        "nearby_windows":    nearby_windows,
    }


# ---------------------------------------------------------------------------
# Bio-Grover Constrained Search PoC (BRCA1)
# ---------------------------------------------------------------------------

@router.post("/quantum-poc/bio-local")
async def bio_grover_local(
    request: BioQuantumRequest,
    current_user: Optional[dict] = Depends(get_optional_user),
):
    """
    Constrained Grover search for BRCA1 disease detection.

    1. Fetches the BRCA1 mRNA reference sequence (NM_007294.4) from NCBI.
    2. Splits it into non-overlapping nodes of (n_codons * 3) nucleotides.
    3. Fetches the c.5266dupC disease marker node from NCBI.
    4. Encodes the marker as the Grover search target (2 bits/nucleotide).
    5. Runs constrained Grover on the Aer simulator — search space is ONLY
       the loaded DNA nodes, not all 2^n states.

    n_codons=1 -> 6 qubits,  n_codons=2 -> 12 qubits,  n_codons=3 -> 18 qubits
    """
    n_codons = request.n_codons
    n_qubits = n_codons * 6

    # Fetch BRCA1 reference from NCBI
    reference_seq = fetch_patient_dna()
    if not reference_seq:
        raise HTTPException(status_code=502, detail="Failed to fetch patient DNA from NCBI (NM_007294.4).")

    # Mutant sequence = reference + c.5266dupC insertion (always used for marker extraction)
    mutant_seq       = apply_brca1_mutation(reference_seq)
    marker_seq_clean = get_marker_seq(mutant_seq, n_codons, offset=request.marker_offset)  # disease marker from mutant

    # Patient sequence depends on whether this patient carries the mutation
    patient_seq  = mutant_seq if request.has_mutation else reference_seq
    patient_label = "Carrier (c.5266dupC present)" if request.has_mutation else "Healthy Control"

    # Build patient node table — scale with qubit capacity
    # We use enough sequence to fill 2^n_qubits slots at most (preventing huge tables).
    # Exception: n_codons=1 has ≤4^3=64 unique 3-nt windows (the full 2^6 state space),
    # so using the full sequence is safe and necessary — the 192-nt cap (64×3) would
    # exclude INSERTION_POS=5266, making the carrier's mutation codon invisible.
    max_nodes    = min(1 << n_qubits, 4096)
    use_len      = len(patient_seq) if n_codons == 1 else max_nodes * n_codons * 3
    patient_used = patient_seq[:use_len]
    patient_nodes = build_patient_nodes(patient_used, n_codons)

    if not patient_nodes:
        raise HTTPException(status_code=500, detail="Patient node table is empty after filtering.")

    # Use custom marker if provided (frontend override), otherwise auto-compute from mutant sequence
    if request.custom_marker_dna:
        if len(request.custom_marker_dna) != n_codons * 3:
            raise HTTPException(
                status_code=422,
                detail=f"custom_marker_dna must be {n_codons * 3} characters for n_codons={n_codons}.",
            )
        marker_seq_clean = request.custom_marker_dna

    target_bits = encode_dna(marker_seq_clean)

    # Check marker uniqueness: does the mutant marker appear in the healthy reference?
    # If yes, the marker is too short/non-specific — Grover will "find" it in both patients.
    ref_nodes_bits  = {nd["bits"] for nd in build_patient_nodes(reference_seq[:use_len], n_codons)}
    marker_in_reference = target_bits in ref_nodes_bits

    # Run constrained Grover on Aer
    try:
        sim_result = run_bio_grover_local(patient_nodes, target_bits, shots=1024)
    except Exception as exc:
        raise HTTPException(status_code=500, detail=f"Quantum simulation error: {exc}")

    # Build a node table preview (up to 64 nodes) — the target node is ALWAYS
    # included even when it falls beyond the first 64 positions (e.g. position
    # ~877 for n_codons=2 at INSERTION_POS 5266), otherwise the frontend table
    # would show "DETECTED" with no highlighted row.
    PREVIEW_SIZE = 64
    target_node = next((nd for nd in patient_nodes if nd["bits"] == target_bits), None)
    if target_node:
        early = [nd for nd in patient_nodes[:PREVIEW_SIZE - 1] if nd["bits"] != target_bits]
        preview_src = sorted(early + [target_node], key=lambda nd: nd["position"])
    else:
        preview_src = patient_nodes[:PREVIEW_SIZE]
    nodes_preview = [
        {"position": nd["position"], "dna": nd["dna"], "bits": nd["bits"],
         "is_target": nd["bits"] == target_bits}
        for nd in preview_src
    ]

    response = {
        **sim_result,
        "n_codons":             n_codons,
        "n_qubits":             n_qubits,
        "has_mutation":         request.has_mutation,
        "patient_label":        patient_label,
        "patient_accession":    PATIENT_ACCESSION,
        "marker_gene":          MARKER_GENE_NAME,
        "marker_variant":       MARKER_VARIANT,
        "marker_region":        MARKER_REGION_DESC,
        "marker_dna":           marker_seq_clean,
        "marker_bits":          target_bits,
        "marker_in_reference":  marker_in_reference,
        "total_nodes":          len(patient_nodes),
        "nodes_preview":        nodes_preview,
    }

    if current_user:
        db = get_db()
        await db.history.insert_one({
            "email":            current_user["sub"],
            "type":             "bio_grover_local",
            "run_id":           request.run_id,
            "target_bits":      target_bits,
            "measured_state":   sim_result["measured_state"],
            "n_codons":         n_codons,
            "n_qubits":         n_qubits,
            "has_mutation":     request.has_mutation,
            "patient_label":    patient_label,
            "marker_variant":   MARKER_VARIANT,
            "marker_dna":       marker_seq_clean,
            "custom_marker_dna": request.custom_marker_dna or None,
            "detection_result": sim_result["detection_result"],
            "confidence":       sim_result["confidence"],
            "execution_time_ms": sim_result["execution_time_ms"],
            "timestamp":        _get_timestamp(request.client_timestamp),
        })

    return response


@router.get("/quantum-poc/bio-circuit-info")
async def bio_circuit_info(
    target_bits: str,
    n_codons: int = 1,
    has_mutation: bool = True,
):
    """
    Rebuild the BRCA1 Grover step-circuits for a history entry.
    Uses the module-level NCBI sequence cache — fast after the first simulation run.
    target_bits / n_codons / has_mutation are read directly from the stored history record.
    """
    if not (1 <= n_codons <= 3):
        raise HTTPException(status_code=422, detail="n_codons must be 1, 2, or 3.")
    if not target_bits or not all(c in "01" for c in target_bits):
        raise HTTPException(status_code=422, detail="target_bits must be a binary string.")
    if len(target_bits) != n_codons * 6:
        raise HTTPException(status_code=422, detail=f"target_bits length must be {n_codons * 6} for n_codons={n_codons}.")

    n_qubits = n_codons * 6

    reference_seq = fetch_patient_dna()
    if not reference_seq:
        raise HTTPException(status_code=502, detail="Failed to fetch BRCA1 reference from NCBI.")

    patient_seq = apply_brca1_mutation(reference_seq) if has_mutation else reference_seq
    max_nodes   = min(1 << n_qubits, 4096)
    use_len     = len(patient_seq) if n_codons == 1 else max_nodes * n_codons * 3
    patient_nodes = build_patient_nodes(patient_seq[:use_len], n_codons)

    if not patient_nodes:
        raise HTTPException(status_code=500, detail="Patient node table is empty.")

    try:
        step_circuits = build_bio_step_circuits(patient_nodes, target_bits)
    except Exception as exc:
        raise HTTPException(status_code=500, detail=f"Circuit build error: {exc}")

    return {"step_circuits": step_circuits, "n_qubits": n_qubits}
