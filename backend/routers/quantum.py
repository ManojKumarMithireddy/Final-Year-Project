"""
Quantum router — Grover POC (local simulator) and IBM Cloud QPU
submit / status endpoints.
"""

import time
import datetime
from typing import Optional

from fastapi import APIRouter, HTTPException, Depends
from qiskit import transpile
from qiskit_aer import AerSimulator
from qiskit_aer.noise import NoiseModel, depolarizing_error
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager

from db import get_db
from auth import get_current_user, get_optional_user
from crypto import decrypt_value
from models.schemas import QuantumToyRequest, IBMSubmitRequest, IBMStatusRequest, BioQuantumRequest, BioIBMSubmitRequest
from services.grover import build_grover_circuit, build_grover_step_circuits
from services.bio_grover import (
    fetch_patient_dna, fetch_disease_marker,
    build_patient_nodes, build_bio_grover_circuit, run_bio_grover_local,
    encode_dna, decode_bits,
    MARKER_GENE_NAME, MARKER_VARIANT, MARKER_REGION_DESC, PATIENT_ACCESSION,
)

# IBM Cloud Imports (optional)
try:
    from qiskit_ibm_runtime import QiskitRuntimeService, SamplerV2 as Sampler
    IBM_RUNTIME_AVAILABLE = True
except ImportError:
    IBM_RUNTIME_AVAILABLE = False

router = APIRouter(prefix="/api/search", tags=["quantum"])


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
            "timestamp": datetime.datetime.now(datetime.timezone.utc),
        })

    return result_data


@router.get("/circuit-info")
async def circuit_info(target_bits: str = "0000"):
    """
    Returns the Grover circuit diagram and step-by-step data for a given
    target bitstring WITHOUT running a simulation. Used by the History page
    to visualize past results on demand.
    """
    if not target_bits or not all(c in '01' for c in target_bits) or len(target_bits) > 8:
        raise HTTPException(status_code=400, detail="target_bits must be 1-8 characters of 0/1.")

    qc, iterations = build_grover_circuit(target_bits)
    circuit_diagram = str(qc.draw(output="text", fold=-1))
    step_circuits = build_grover_step_circuits(target_bits)

    return {
        "circuit_diagram": circuit_diagram,
        "step_circuits": step_circuits,
        "iterations": iterations,
        "qubits": len(target_bits),
    }


@router.post("/quantum-poc/ibm-submit")
async def ibm_submit_job(request: IBMSubmitRequest, current_user: dict = Depends(get_current_user)):
    """
    Submits a Grover search circuit to an IBM Cloud QPU.

    SECURITY: IBM credentials (api_key, crn) are fetched server-side from the
    authenticated user's saved credentials — they are NOT transmitted in the request body.
    """
    if not IBM_RUNTIME_AVAILABLE:
        raise HTTPException(status_code=500, detail="qiskit-ibm-runtime not installed")

    db = get_db()
    creds = await db.credentials.find_one({"email": current_user["sub"]})
    if not creds or not creds.get("ibm_api_key"):
        raise HTTPException(
            status_code=400,
            detail="No IBM Cloud credentials saved. Please save your API key and CRN in the credentials panel first.",
        )

    try:
        service = QiskitRuntimeService(
            channel="ibm_cloud",
            token=decrypt_value(creds["ibm_api_key"]),
            instance=creds["ibm_crn"],
        )
        backend = service.least_busy(operational=True, simulator=False)

        target_bitstring = request.target_bits
        if not target_bitstring or not all(c in '01' for c in target_bitstring):
            target_bitstring = "101"

        qc, iterations = build_grover_circuit(target_bitstring)

        pm = generate_preset_pass_manager(backend=backend, optimization_level=3)
        isa_circuit = pm.run(qc)

        sampler = Sampler(mode=backend)
        job = sampler.run([isa_circuit])

        job_data = {
            "job_id": job.job_id(),
            "backend": backend.name,
            "circuit_diagram": str(qc.draw(output="text")),
            "status": job.status(),
            "iterations": iterations,
            "execution_time_ms": 0,
        }

        await db.history.insert_one({
            "email": current_user["sub"],
            "type": "quantum_ibm_submit",
            "target_bits": request.target_bits,
            "job_id": job.job_id(),
            "backend": backend.name,
            "timestamp": datetime.datetime.now(datetime.timezone.utc),
        })

        return job_data
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/quantum-poc/ibm-status")
async def ibm_job_status(request: IBMStatusRequest, current_user: dict = Depends(get_current_user)):
    """
    Checks the status of an IBM Cloud job.
    IBM credentials are fetched server-side from the authenticated user's saved credentials.
    """
    if not IBM_RUNTIME_AVAILABLE:
        raise HTTPException(status_code=500, detail="qiskit-ibm-runtime not installed")

    db = get_db()
    creds = await db.credentials.find_one({"email": current_user["sub"]})
    if not creds or not creds.get("ibm_api_key"):
        raise HTTPException(status_code=400, detail="No IBM Cloud credentials found. Please save your credentials first.")

    try:
        service = QiskitRuntimeService(
            channel="ibm_cloud",
            token=decrypt_value(creds["ibm_api_key"]),
            instance=creds["ibm_crn"],
        )
        job = service.job(request.job_id)
        status = job.status()

        if status == "DONE":
            result = job.result()
            counts = result[0].data.meas.get_counts()
            measured_state = max(counts, key=counts.get)
            return {
                "status": status,
                "counts": counts,
                "measured_state": measured_state,
                "execution_time_ms": 0,
            }
        return {"status": status}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


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

    # Fetch patient DNA from NCBI
    patient_seq = fetch_patient_dna()
    if not patient_seq:
        raise HTTPException(status_code=502, detail="Failed to fetch patient DNA from NCBI (NM_007294.4).")

    # Fetch disease marker from NCBI
    marker_seq = fetch_disease_marker(n_codons)
    if not marker_seq or len(marker_seq) < n_codons * 3:
        raise HTTPException(status_code=502, detail="Failed to fetch BRCA1 disease marker from NCBI.")

    marker_seq_clean = marker_seq[: n_codons * 3]

    # Build patient node table — scale with qubit capacity
    # We use enough sequence to fill 2^n_qubits slots at most (preventing huge tables)
    max_nodes    = min(1 << n_qubits, 4096)
    use_len      = max_nodes * n_codons * 3
    patient_used = patient_seq[:use_len]
    patient_nodes = build_patient_nodes(patient_used, n_codons)

    if not patient_nodes:
        raise HTTPException(status_code=500, detail="Patient node table is empty after filtering.")

    target_bits = encode_dna(marker_seq_clean)

    # Run constrained Grover on Aer
    try:
        sim_result = run_bio_grover_local(patient_nodes, target_bits, shots=1024)
    except Exception as exc:
        raise HTTPException(status_code=500, detail=f"Quantum simulation error: {exc}")

    # Build a node table preview (first 64 nodes for UI)
    nodes_preview = [
        {"position": nd["position"], "dna": nd["dna"], "bits": nd["bits"],
         "is_target": nd["bits"] == target_bits}
        for nd in patient_nodes[:64]
    ]

    response = {
        **sim_result,
        "n_codons":          n_codons,
        "n_qubits":          n_qubits,
        "patient_accession": PATIENT_ACCESSION,
        "marker_gene":       MARKER_GENE_NAME,
        "marker_variant":    MARKER_VARIANT,
        "marker_region":     MARKER_REGION_DESC,
        "marker_dna":        marker_seq_clean,
        "marker_bits":       target_bits,
        "total_nodes":       len(patient_nodes),
        "nodes_preview":     nodes_preview,
    }

    if current_user:
        db = get_db()
        await db.history.insert_one({
            "email":            current_user["sub"],
            "type":             "bio_grover_local",
            "n_codons":         n_codons,
            "marker_variant":   MARKER_VARIANT,
            "detection_result": sim_result["detection_result"],
            "confidence":       sim_result["confidence"],
            "timestamp":        datetime.datetime.now(datetime.timezone.utc),
        })

    return response


@router.post("/quantum-poc/bio-ibm-submit")
async def bio_grover_ibm_submit(
    request: BioIBMSubmitRequest,
    current_user: dict = Depends(get_current_user),
):
    """
    Submit a constrained BRCA1 Grover circuit to IBM Cloud QPU.
    IBM credentials fetched server-side (same as toy PoC).
    n_codons capped at 2 (12 qubits) to fit real QPU topology.
    """
    if not IBM_RUNTIME_AVAILABLE:
        raise HTTPException(status_code=500, detail="qiskit-ibm-runtime not installed")

    db    = get_db()
    creds = await db.credentials.find_one({"email": current_user["sub"]})
    if not creds or not creds.get("ibm_api_key"):
        raise HTTPException(
            status_code=400,
            detail="No IBM Cloud credentials saved. Please save your API key and CRN first.",
        )

    n_codons = request.n_codons
    n_qubits = n_codons * 6

    # Fetch sequences
    patient_seq = fetch_patient_dna()
    if not patient_seq:
        raise HTTPException(status_code=502, detail="Failed to fetch patient DNA from NCBI.")
    marker_seq = fetch_disease_marker(n_codons)
    if not marker_seq or len(marker_seq) < n_codons * 3:
        raise HTTPException(status_code=502, detail="Failed to fetch BRCA1 marker from NCBI.")

    marker_seq_clean = marker_seq[: n_codons * 3]
    max_nodes         = min(1 << n_qubits, 512)   # smaller cap for QPU
    patient_nodes     = build_patient_nodes(patient_seq[: max_nodes * n_codons * 3], n_codons)
    target_bits       = encode_dna(marker_seq_clean)

    try:
        qc, k, n_unique, n_unconstrained, target_found, _ = \
            build_bio_grover_circuit(patient_nodes, target_bits)

        service = QiskitRuntimeService(
            channel="ibm_cloud",
            token=decrypt_value(creds["ibm_api_key"]),
            instance=creds["ibm_crn"],
        )
        backend    = service.least_busy(operational=True, simulator=False)
        pm         = generate_preset_pass_manager(backend=backend, optimization_level=3)
        isa_circuit = pm.run(qc)
        sampler    = Sampler(mode=backend)
        job        = sampler.run([isa_circuit])

        await db.history.insert_one({
            "email":          current_user["sub"],
            "type":           "bio_grover_ibm_submit",
            "n_codons":       n_codons,
            "marker_variant": MARKER_VARIANT,
            "job_id":         job.job_id(),
            "backend":        backend.name,
            "timestamp":      datetime.datetime.now(datetime.timezone.utc),
        })

        return {
            "job_id":          job.job_id(),
            "backend":         backend.name,
            "status":          job.status(),
            "iterations":      k,
            "n_qubits":        n_qubits,
            "n_unique":        n_unique,
            "n_unconstrained": n_unconstrained,
            "target_found":    target_found,
            "marker_dna":      marker_seq_clean,
            "marker_bits":     target_bits,
            "circuit_diagram": str(qc.draw(output="text")),
        }
    except Exception as exc:
        raise HTTPException(status_code=500, detail=str(exc))
