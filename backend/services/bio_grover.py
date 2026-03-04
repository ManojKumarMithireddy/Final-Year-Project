"""
Constrained Grover's algorithm for BRCA1 genetic disease detection PoC.

Key algorithmic difference from the toy PoC
--------------------------------------------
  Toy PoC:  |psi0> = H^n|0>           uniform superposition over ALL 2^n states
  This PoC: |psi0> = StatePrep(sv)|0>  uniform over ONLY valid DNA nodes

Constrained diffusion operator:
  D_c = 2|psi0><psi0| - I
  Implemented as: StatePrep . (2|0><0|-I) . StatePrep_dagger

This restricts amplitude flow entirely to the loaded DNA nodes,
so the algorithm cannot amplify states that have no real patient data.

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
from qiskit.circuit.library import StatePreparation
from qiskit_aer import AerSimulator

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
# 0-based nt start; 5262 % 3 == 0 and 5262 % 6 == 0 (node-boundary-aligned)
MARKER_NT_START = 5262

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


def fetch_disease_marker(n_codons: int) -> Optional[str]:
    """
    Fetch the disease marker node from NCBI.
    Returns n_codons * 3 nucleotides from the BRCA1 c.5266dupC hotspot region.
    """
    start_1b = MARKER_NT_START + 1          # NCBI is 1-based
    stop_1b  = MARKER_NT_START + n_codons * 3
    return _ncbi_fetch(MARKER_ACCESSION, seq_start=start_1b, seq_stop=stop_1b)


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
# Quantum circuit primitives
# ---------------------------------------------------------------------------

def _oracle(n: int, target_bits: str) -> QuantumCircuit:
    """
    Phase oracle: applies -1 phase to |target_bits> and +1 to all others.
    Uses reversed(target_bits) -> qubit-index mapping matching Qiskit's
    measurement output (leftmost char = q_{n-1}, rightmost = q_0).
    """
    qc = QuantumCircuit(n, name="Oracle")
    for i, bit in enumerate(reversed(target_bits)):
        if bit == "0":
            qc.x(i)
    qc.h(n - 1)
    qc.mcx(list(range(n - 1)), n - 1)
    qc.h(n - 1)
    for i, bit in enumerate(reversed(target_bits)):
        if bit == "0":
            qc.x(i)
    return qc


def _zero_reflection(n: int) -> QuantumCircuit:
    """
    Implements 2|0...0><0...0| - I.
    Inner core of the constrained diffuser: StatePrep . ZeroReflect . StatePrep_dag
    """
    qc = QuantumCircuit(n, name="ZeroReflect")
    qc.x(range(n))
    qc.h(n - 1)
    qc.mcx(list(range(n - 1)), n - 1)
    qc.h(n - 1)
    qc.x(range(n))
    return qc


# ---------------------------------------------------------------------------
# Main constrained Grover circuit builder
# ---------------------------------------------------------------------------

def build_bio_grover_circuit(
    patient_nodes: List[Dict],
    target_bits:   str,
) -> Tuple[QuantumCircuit, int, int, int, bool, Dict[str, float]]:
    """
    Build the constrained Grover circuit for BRCA1 disease detection.

    Circuit structure:
      1. DNA|psi0> = StatePreparation(sv)   [uniform over valid DNA nodes only]
      2. k iterations of:
           Oracle(target)                   [phase-flip |target_bits>]
           ConstrainedDiffusion             [2|psi0><psi0| - I]
      3. Measure all

    where k = floor(pi/4 * sqrt(|S|))  and  |S| = number of unique DNA nodes.

    Returns: (qc, k, n_unique, n_unconstrained, target_found, init_probs)
    """
    n       = len(target_bits)
    n_max   = 1 << n   # 2^n

    # Deduplicate while preserving insertion order
    seen: Dict[str, bool] = {}
    for nd in patient_nodes:
        seen.setdefault(nd["bits"], True)
    unique_bits  = list(seen.keys())
    N            = len(unique_bits)
    target_found = target_bits in seen

    if N == 0:
        raise ValueError("No valid patient nodes to build search space from.")

    # Uniform statevector over unique DNA node indices
    amplitude = 1.0 / sqrt(N)
    sv        = np.zeros(n_max, dtype=complex)
    for bits in unique_bits:
        sv[int(bits, 2)] = amplitude

    state_prep = StatePreparation(sv, label="DNA|psi0>")

    # Constrained diffusion: StatePrep_dag . ZeroReflect . StatePrep
    diffusion = QuantumCircuit(n, name="Constrained Diffusion")
    diffusion.append(state_prep.inverse(), range(n))
    diffusion.compose(_zero_reflection(n), inplace=True)
    diffusion.append(state_prep, range(n))

    # Optimal iteration count k = floor(pi/4 * sqrt(N))
    k = max(1, floor(pi / 4 * sqrt(N)))

    # Assemble full circuit
    qc = QuantumCircuit(n)
    qc.append(state_prep, range(n))
    for _ in range(k):
        qc.append(_oracle(n, target_bits), range(n))
        qc.append(diffusion, range(n))
    qc.measure_all()

    # Initial probability distribution (directly from sv, no simulation needed)
    init_probs: Dict[str, float] = {
        format(i, f"0{n}b"): round(float(abs(sv[i]) ** 2), 8)
        for i in range(n_max)
        if abs(sv[i]) > 1e-12
    }

    return qc, k, N, n_max, target_found, init_probs


# ---------------------------------------------------------------------------
# Local Aer simulation entry point
# ---------------------------------------------------------------------------

def run_bio_grover_local(
    patient_nodes: List[Dict],
    target_bits:   str,
    shots:         int = 1024,
) -> Dict:
    """
    Execute the constrained Grover circuit on the local Aer simulator.
    Returns a dict ready to merge into the FastAPI response payload.
    """
    qc, k, n_unique, n_unconstrained, target_found, init_probs = \
        build_bio_grover_circuit(patient_nodes, target_bits)

    n = len(target_bits)
    simulator  = AerSimulator()
    transpiled = transpile(qc, simulator)

    t0      = time.perf_counter()
    result  = simulator.run(transpiled, shots=shots).result()
    exec_ms = (time.perf_counter() - t0) * 1000

    counts   = result.get_counts()
    measured = max(counts, key=counts.get)
    total    = sum(counts.values())
    conf     = counts.get(target_bits, 0) / total if total else 0.0

    return {
        "measured_state":     measured,
        "counts":             counts,
        "execution_time_ms":  round(exec_ms, 2),
        "iterations":         k,
        "n_qubits":           n,
        "n_unique":           n_unique,
        "n_unconstrained":    n_unconstrained,
        "target_found":       target_found,
        "detection_result":   "FOUND" if measured == target_bits else "NOT_FOUND",
        "confidence":         round(conf, 4),
        "init_probs":         init_probs,
        "circuit_diagram":    str(qc.draw(output="text", fold=-1)),
    }
