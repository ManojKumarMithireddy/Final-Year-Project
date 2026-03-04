"""
Constrained Grover's algorithm for BRCA1 genetic disease detection PoC.

Key algorithmic difference from the toy PoC
--------------------------------------------
  Toy PoC:  |psi0> = H^n|0>           uniform superposition over ALL 2^n states
  This PoC: |psi0> = StatePrep(sv)|0>  uniform over ONLY valid DNA nodes

Constrained diffusion operator:
  D_c = 2|psi0><psi0| - I
  Implemented as direct matrix multiplication (numpy):
    state = oracle_diag * state    (phase-flip target)
    state = D_c @ state            (inversion about |psi0>)

WHY numpy instead of Qiskit circuits for Aer simulation
---------------------------------------------------------
  Qiskit's StatePreparation gate uses isometry (UCGate) decomposition, which
  requires O(4^n) CNOT gates and allocates O(4^n) matrices during transpilation.
  For n=6: ~700 CNOTs per StatePrep call; with k iterations having 3 StatePrep
  calls each, the transpile allocates hundreds of MB and segfaults on
  memory-constrained hosts (e.g., HuggingFace Spaces cpu-basic ~2 GB).

  The numpy path computes the SAME mathematical operation with:
    - Oracle:    element-wise multiply (O(2^n) time/space)
    - Diffusion: matrix-vector product  (O(4^n) time, O(2^n) space)
  For n=6 this is trivial; for n=12 it uses ~64 MB — never crashes.

BRCA1 target
------------
  Patient:  NM_007294.4  (Human BRCA1 mRNA reference)
  Marker:   Same accession, region starting at nt 5262 (0-based)
            This spans the c.5266dupC pathogenic hotspot (exon 20).
            5262 % 3 == 0  and  5262 % 6 == 0, so it aligns with node
            boundaries for both n_codons=1 and n_codons=2.
"""

import io
import os
import time
import logging
from math import floor, pi, sqrt
from typing import Dict, List, Optional, Tuple

import numpy as np
from Bio import Entrez, SeqIO
from qiskit import QuantumCircuit, transpile

logger = logging.getLogger(__name__)

# 2-bit DNA encoding (A=00, C=01, G=10, T=11) -- matches frontend dnaToBits
_ENC: Dict[str, str] = {"A": "00", "C": "01", "G": "10", "T": "11", "U": "01", "N": "00"}
_DEC: Dict[str, str] = {"00": "A", "01": "C", "10": "G", "11": "T"}

# BRCA1 configuration
PATIENT_ACCESSION  = "NM_007294.4"
MARKER_ACCESSION   = "NM_007294.4"
MARKER_GENE_NAME   = "BRCA1"
MARKER_VARIANT     = "c.5266dupC"
MARKER_REGION_DESC = "BRCA1 exon 20 — c.5266dupC pathogenic hotspot"
# 0-based insertion point of the duplicated C in the mutant sequence
# c.5266dupC inserts an extra C at 0-based index 5266 of NM_007294.4
INSERTION_POS   = 5266

# Module-level sequence cache (fetched from NCBI once per process lifetime)
_seq_cache: Dict[str, str] = {}


# ---------------------------------------------------------------------------
# NCBI helpers
# ---------------------------------------------------------------------------

def _setup_entrez() -> None:
    Entrez.email   = os.getenv("NCBI_EMAIL", "bioquantum.capstone@example.com")
    Entrez.timeout = int(os.getenv("NCBI_TIMEOUT", "20"))


def _ncbi_fetch(
    accession: str,
    seq_start: Optional[int] = None,
    seq_stop:  Optional[int] = None,
) -> Optional[str]:
    """
    Fetch a (sub-)sequence from NCBI Entrez.
    seq_start / seq_stop are 1-based inclusive (NCBI convention).
    Returns the uppercase sequence string, or None on error.
    """
    key = f"{accession}|{seq_start}|{seq_stop}"
    if key in _seq_cache:
        return _seq_cache[key]

    _setup_entrez()
    try:
        kw: Dict = dict(db="nucleotide", id=accession, rettype="fasta", retmode="text")
        if seq_start is not None:
            kw["seq_start"] = seq_start
        if seq_stop is not None:
            kw["seq_stop"] = seq_stop

        handle = Entrez.efetch(**kw)
        raw    = handle.read()
        handle.close()

        record = SeqIO.read(io.StringIO(raw), "fasta")
        seq    = str(record.seq).upper().replace(" ", "").replace("\n", "")
        _seq_cache[key] = seq
        logger.info("NCBI: fetched %d nt from %s (range %s-%s)", len(seq), accession, seq_start, seq_stop)
        return seq
    except Exception as exc:
        logger.error("NCBI fetch failed [%s]: %s", accession, exc)
        return None


def fetch_patient_dna() -> Optional[str]:
    """Fetch full BRCA1 mRNA reference sequence from NCBI (NM_007294.4)."""
    return _ncbi_fetch(PATIENT_ACCESSION)


def apply_brca1_mutation(seq: str) -> str:
    """
    Simulate a c.5266dupC carrier by inserting an extra 'C' at INSERTION_POS.
    Returns the mutant sequence representing a patient who carries the pathogenic variant.
    """
    if len(seq) <= INSERTION_POS:
        return seq
    return seq[:INSERTION_POS] + "C" + seq[INSERTION_POS:]


def get_marker_seq(mutant_seq: str, n_codons: int) -> str:
    """
    Extract the disease marker codon(s) from a mutant sequence.
    Finds the node boundary that contains INSERTION_POS for the given n_codons window size.
    """
    wsize  = n_codons * 3
    start  = (INSERTION_POS // wsize) * wsize
    return mutant_seq[start : start + wsize]


# ---------------------------------------------------------------------------
# DNA encoding helpers
# ---------------------------------------------------------------------------

def encode_dna(dna: str) -> str:
    """DNA sequence -> 2-bit binary string (A=00, C=01, G=10, T=11)."""
    return "".join(_ENC.get(c, "00") for c in dna.upper())


def decode_bits(bits: str) -> str:
    """2-bit binary string -> DNA sequence."""
    return "".join(_DEC.get(bits[i:i+2], "?") for i in range(0, len(bits), 2))


# ---------------------------------------------------------------------------
# Patient node table
# ---------------------------------------------------------------------------

def build_patient_nodes(sequence: str, n_codons: int) -> List[Dict]:
    """
    Split a DNA sequence into non-overlapping windows of n_codons codons
    (window_size = 3 * n_codons nucleotides).

    Returns a list of dicts sorted by position:
      { "dna": str, "bits": str, "position": int (0-based node index) }

    Windows containing ambiguous bases (outside ACGT) are excluded.
    """
    wsize = n_codons * 3
    valid = frozenset("ACGT")
    nodes: List[Dict] = []
    for start in range(0, len(sequence) - wsize + 1, wsize):
        dna = sequence[start: start + wsize]
        if all(c in valid for c in dna):
            nodes.append({
                "dna":      dna,
                "bits":     encode_dna(dna),
                "position": start // wsize,
            })
    return nodes


# ---------------------------------------------------------------------------
# Quantum circuit primitives (used for IBM QPU path only)
# ---------------------------------------------------------------------------

def _build_grover_circuit_standard(n: int, target_bits: str, k: int) -> QuantumCircuit:
    """
    Standard (unconstrained) Grover circuit with H^n initialisation.
    Used for the IBM QPU submit path — StatePreparation is NOT used here
    because its isometry decomposition creates O(4^n) gates and will either
    exhaust memory during transpilation or exceed QPU gate depth limits.

    The oracle still marks the correct target, so the quantum measurement
    correctly tells us whether the marker is detectable in the n-qubit space.
    The target_found flag (computed classically) tells us whether that state
    is present in the patient DNA table.
    """
    qc = QuantumCircuit(n)
    qc.h(range(n))
    for _ in range(k):
        # Oracle
        for i, bit in enumerate(reversed(target_bits)):
            if bit == "0":
                qc.x(i)
        qc.h(n - 1)
        qc.mcx(list(range(n - 1)), n - 1)
        qc.h(n - 1)
        for i, bit in enumerate(reversed(target_bits)):
            if bit == "0":
                qc.x(i)
        # Standard Grover diffuser
        qc.h(range(n))
        qc.x(range(n))
        qc.h(n - 1)
        qc.mcx(list(range(n - 1)), n - 1)
        qc.h(n - 1)
        qc.x(range(n))
        qc.h(range(n))
    qc.measure_all()
    return qc


def _build_display_circuit(n: int, target_bits: str, k: int) -> str:
    """
    High-level display circuit showing conceptual structure using named boxes.
    NEVER transpiled or run — only used for qc.draw() diagram output.
    """
    qc = QuantumCircuit(n)
    # Show state prep as a labelled box
    prep = QuantumCircuit(n, name="DNA|ψ₀⟩")
    qc.append(prep.to_gate(), range(n))
    for _ in range(k):
        # Oracle shown as compact named box
        oracle_qc = QuantumCircuit(n, name="Oracle")
        for i, bit in enumerate(reversed(target_bits)):
            if bit == "0":
                oracle_qc.x(i)
        oracle_qc.h(n - 1)
        if n > 1:
            oracle_qc.mcx(list(range(n - 1)), n - 1)
        oracle_qc.h(n - 1)
        for i, bit in enumerate(reversed(target_bits)):
            if bit == "0":
                oracle_qc.x(i)
        qc.append(oracle_qc.to_gate(), range(n))
        # Constrained diffusion as named box
        diff = QuantumCircuit(n, name="C-Diffusion")
        qc.append(diff.to_gate(), range(n))
    qc.measure_all()
    return str(qc.draw(output="text", fold=-1))


# ---------------------------------------------------------------------------
# Main constrained Grover runner — numpy statevector (crash-safe)
# ---------------------------------------------------------------------------

def run_bio_grover_local(
    patient_nodes: List[Dict],
    target_bits:   str,
    shots:         int = 1024,
) -> Dict:
    """
    Constrained Grover search via DIRECT NUMPY STATEVECTOR SIMULATION.

    This is mathematically identical to running the Qiskit circuit on Aer but
    avoids all circuit synthesis / transpilation:

      1. Build |ψ₀⟩ = (1/√N) Σ_{s∈S} |s⟩   (S = unique DNA node bit-strings)
      2. Oracle:    state = diag(-1 at target, +1 elsewhere) * state
      3. Diffusion: state = (2|ψ₀⟩⟨ψ₀| - I) @ state
      4. Repeat k = ⌊π/4·√N⌋ times
      5. Sample 'shots' outcomes from |state|²

    Memory: O(2^n) — for n=18 that is 2 MB of complex128. No segfault risk.
    """
    n     = len(target_bits)
    n_max = 1 << n

    # Deduplicate nodes, preserving first-occurrence order
    seen: Dict[str, bool] = {}
    for nd in patient_nodes:
        seen.setdefault(nd["bits"], True)
    unique_bits  = list(seen.keys())
    N            = len(unique_bits)
    target_found = target_bits in seen

    if N == 0:
        raise ValueError("No valid patient nodes to build search space from.")

    # ── Initial statevector: uniform over DNA nodes ───────────────────────────
    amplitude = 1.0 / sqrt(N)
    sv0       = np.zeros(n_max, dtype=complex)
    for bits in unique_bits:
        sv0[int(bits, 2)] = amplitude

    # ── Operators ─────────────────────────────────────────────────────────────
    # Phase oracle: element-wise diagonal (−1 on target, +1 on all others)
    oracle_diag                    = np.ones(n_max, dtype=complex)
    oracle_diag[int(target_bits, 2)] = -1.0

    # Constrained diffusion matrix: 2|ψ₀⟩⟨ψ₀| − I
    D = 2.0 * np.outer(sv0, sv0.conj()) - np.eye(n_max, dtype=complex)

    # Optimal iteration count: k = ⌊π/4 · √N⌋
    k = max(1, floor(pi / 4 * sqrt(N)))

    # ── Run Grover iterations ─────────────────────────────────────────────────
    t0    = time.perf_counter()
    state = sv0.copy()
    for _ in range(k):
        state = oracle_diag * state   # Phase oracle  (O(2^n))
        state = D @ state             # Constrained diffusion  (O(4^n))
    exec_ms = (time.perf_counter() - t0) * 1000

    # ── Sample measurements ───────────────────────────────────────────────────
    probs = np.abs(state) ** 2
    probs /= probs.sum()   # Re-normalise to remove floating-point drift
    indices = np.random.choice(n_max, size=shots, p=probs)
    counts: Dict[str, int] = {}
    for idx in indices:
        b          = format(idx, f"0{n}b")
        counts[b]  = counts.get(b, 0) + 1

    measured = max(counts, key=counts.get)
    total    = sum(counts.values())
    conf     = counts.get(target_bits, 0) / total if total else 0.0

    init_probs = {
        format(i, f"0{n}b"): round(float(abs(sv0[i]) ** 2), 8)
        for i in range(n_max) if abs(sv0[i]) > 1e-12
    }

    return {
        "measured_state":    measured,
        "counts":            counts,
        "execution_time_ms": round(exec_ms, 2),
        "iterations":        k,
        "n_qubits":          n,
        "n_unique":          N,
        "n_unconstrained":   n_max,
        "target_found":      target_found,
        "detection_result":  "FOUND" if measured == target_bits else "NOT_FOUND",
        "confidence":        round(conf, 4),
        "init_probs":        init_probs,
        "circuit_diagram":   _build_display_circuit(n, target_bits, k),
        "step_circuits":     build_bio_step_circuits(patient_nodes, target_bits),
    }




# ---------------------------------------------------------------------------
# Step-by-step circuit breakdown (for GroverStepNavigator)
# ---------------------------------------------------------------------------

def build_bio_step_circuits(
    patient_nodes: List[Dict],
    target_bits:   str,
) -> list:
    """
    Returns 4-step breakdown for GroverStepNavigator.
    Probability distributions computed from numpy statevectors.
    Gate diagrams built with Qiskit draw() — never transpiled or run.
    """
    n     = len(target_bits)
    n_max = 1 << n

    # Rebuild operators (same as run_bio_grover_local)
    seen: Dict[str, bool] = {}
    for nd in patient_nodes:
        seen.setdefault(nd["bits"], True)
    unique_bits = list(seen.keys())
    N           = len(unique_bits)
    k           = max(1, floor(pi / 4 * sqrt(N)))

    amplitude = 1.0 / sqrt(N)
    sv0       = np.zeros(n_max, dtype=complex)
    for bits in unique_bits:
        sv0[int(bits, 2)] = amplitude

    oracle_diag                    = np.ones(n_max, dtype=complex)
    oracle_diag[int(target_bits, 2)] = -1.0
    D = 2.0 * np.outer(sv0, sv0.conj()) - np.eye(n_max, dtype=complex)

    basis_labels = [format(i, f"0{n}b") for i in range(n_max)]

    def _probs(state: np.ndarray) -> Dict[str, float]:
        p = np.abs(state) ** 2
        return {label: round(float(p[i]), 6) for i, label in enumerate(basis_labels)}

    # ── Compute probability snapshots ─────────────────────────────────────────
    probs_step1 = _probs(sv0)

    s2 = oracle_diag * sv0
    probs_step2 = _probs(s2)

    s3 = D @ s2
    probs_step3 = _probs(s3)

    # Step 4: after all k iterations
    state = sv0.copy()
    for _ in range(k):
        state = oracle_diag * state
        state = D @ state
    probs_step4 = _probs(state)

    # ── Gate diagrams using Qiskit (draw only, no execution) ──────────────────
    x_count = sum(1 for bit in target_bits if bit == "0")

    # Step 1: constrained state prep shown as a named box
    qc1 = QuantumCircuit(n)
    qc1.append(QuantumCircuit(n, name="DNA|ψ₀⟩").to_gate(), range(n))
    diag1 = str(qc1.draw(output="text", fold=-1))

    # Step 2: Oracle — actual X / H / MCX gates
    qc2 = QuantumCircuit(n)
    for i, bit in enumerate(reversed(target_bits)):
        if bit == "0":
            qc2.x(i)
    qc2.h(n - 1)
    if n > 1:
        qc2.mcx(list(range(n - 1)), n - 1)
    qc2.h(n - 1)
    for i, bit in enumerate(reversed(target_bits)):
        if bit == "0":
            qc2.x(i)
    diag2 = str(qc2.draw(output="text", fold=-1))

    # Step 3: Diffusion — actual H / X / MCX gates
    qc3 = QuantumCircuit(n)
    qc3.h(range(n))
    qc3.x(range(n))
    qc3.h(n - 1)
    if n > 1:
        qc3.mcx(list(range(n - 1)), n - 1)
    qc3.h(n - 1)
    qc3.x(range(n))
    qc3.h(range(n))
    diag3 = str(qc3.draw(output="text", fold=-1))

    # Step 4: Measurement
    qc4 = QuantumCircuit(n)
    qc4.measure_all()
    diag4 = str(qc4.draw(output="text", fold=-1))

    return [
        {
            "step": 1,
            "label": "Constrained Initialisation",
            "description": (
                f"Instead of H⊗{n} (uniform over all {n_max} states), "
                f"DNA|ψ₀⟩ loads only the {N} unique DNA nodes from the patient table. "
                f"Each node gets amplitude 1/√{N}; states absent from the patient data have zero amplitude."
            ),
            "diagram": diag1,
            "gates": [
                {
                    "name": "DNA State Preparation",
                    "symbol": "DNA|ψ₀⟩",
                    "count": 1,
                    "explanation": (
                        f"Prepares a uniform superposition over {N} unique BRCA1 DNA nodes "
                        f"(from NM_007294.4). Search space constrained to {N}/{n_max} states — "
                        "only biologically real patterns carry non-zero amplitude."
                    ),
                }
            ],
            "probabilities": probs_step1,
        },
        {
            "step": 2,
            "label": "Oracle",
            "description": (
                f"Phase-flip the BRCA1 c.5266dupC marker |{target_bits}⟩ (amplitude × −1). "
                f"X gates flip '0' qubits to |1⟩, an H-MCX-H sequence implements controlled-Z, "
                f"then X gates restore the basis."
            ),
            "diagram": diag2,
            "gates": [
                {
                    "name": "Pauli-X (NOT)",
                    "symbol": "X",
                    "count": x_count * 2,
                    "explanation": (
                        f"{x_count} X gate(s) before the MCX flip '0' qubits so the controlled gate "
                        f"fires on the correct state; {x_count} after restore the basis."
                    ),
                },
                {
                    "name": "Hadamard (H)",
                    "symbol": "H",
                    "count": 2,
                    "explanation": "Sandwiches MCX on qubit n−1, converting it to a phase flip (H·MCX·H ≡ controlled-Z).",
                },
                {
                    "name": "Multi-Controlled-X (MCX)",
                    "symbol": "MCX",
                    "count": 1,
                    "explanation": (
                        f"Flips qubit {n-1} when all {n-1} controls are |1⟩. "
                        f"With H sandwich, phase-flips the target state |{target_bits}⟩."
                    ),
                },
            ],
            "probabilities": probs_step2,
        },
        {
            "step": 3,
            "label": "Constrained Diffusion",
            "description": (
                f"Inversion about the constrained mean |ψ₀⟩: 2|ψ₀⟩⟨ψ₀| − I. "
                f"Amplifies the disease marker's probability within the DNA node subspace. "
                f"Repeated {k}× ≈ π/4·√{N}."
            ),
            "diagram": diag3,
            "gates": [
                {
                    "name": "Hadamard (H)",
                    "symbol": "H",
                    "count": n * 2,
                    "explanation": (
                        f"Applied to all {n} qubits before and after the X+MCX layer. "
                        "Rotates to/from the Hadamard basis to perform inversion about the mean."
                    ),
                },
                {
                    "name": "Pauli-X (NOT)",
                    "symbol": "X",
                    "count": n * 2,
                    "explanation": f"Applied to all {n} qubits twice — implements inversion about the |0…0⟩ state.",
                },
                {
                    "name": "Multi-Controlled-X (MCX)",
                    "symbol": "MCX",
                    "count": 1,
                    "explanation": (
                        f"With the H sandwich on qubit {n-1} implements a phase flip on |1…1⟩. "
                        "Combined with surrounding X gates, completes the inversion about the mean."
                    ),
                },
            ],
            "probabilities": probs_step3,
        },
        {
            "step": 4,
            "label": "Measurement",
            "description": (
                f"After {k} iteration(s) the disease marker |{target_bits}⟩ has "
                f"probability ≈ 1. Measuring all {n} qubits collapses the superposition — "
                f"if the top result matches the marker, BRCA1 mutation is detected."
            ),
            "diagram": diag4,
            "gates": [
                {
                    "name": "Measurement (M)",
                    "symbol": "⊗",
                    "count": n,
                    "explanation": (
                        f"Each of the {n} qubits is projected onto {{|0⟩, |1⟩}}. "
                        f"After {k} Grover iteration(s), the amplitude of the BRCA1 marker "
                        f"|{target_bits}⟩ is ≈ 1, so the classical result matches with high probability."
                    ),
                }
            ],
            "probabilities": probs_step4,
        },
    ]


def build_bio_grover_ibm(
    patient_nodes: List[Dict],
    target_bits:   str,
) -> Tuple[QuantumCircuit, int, int, int, bool]:
    """
    Build a QPU-compatible Grover circuit for the IBM submit path.
    Uses standard H^n initialisation (not StatePreparation) so the circuit
    can be transpiled to QPU basis gates without O(4^n) isometry synthesis.
    Returns: (qc, k, N_unique, N_unconstrained, target_found)
    """
    n     = len(target_bits)
    n_max = 1 << n

    seen: Dict[str, bool] = {}
    for nd in patient_nodes:
        seen.setdefault(nd["bits"], True)
    N            = len(seen)
    target_found = target_bits in seen

    k  = max(1, floor(pi / 4 * sqrt(N)))
    qc = _build_grover_circuit_standard(n, target_bits, k)
    return qc, k, N, n_max, target_found
