"""
Grover's algorithm circuit construction utilities.

Provides factory functions for building complete and step-by-step
Grover search QuantumCircuits for a given target bitstring.
"""

import numpy as np
from functools import lru_cache
from qiskit import QuantumCircuit


def build_grover_circuit(target_bitstring: str):
    """
    Factory function: builds a complete Grover search QuantumCircuit for the given target bitstring.

    Eliminates code duplication between local simulator and IBM Cloud execution paths.

    NOTE on Privacy / PIR Demonstration:
    The 'target_bitstring' received here may be XOR-encrypted by the client using a
    One-Time Pad (OTP) before transmission. The server executes Grover's algorithm
    without knowing the original (plaintext) query — it only sees the encrypted bits.
    The client recovers the true result by XOR-ing the returned state with the same OTP key.

    This approximates a Privacy-Preserving Information Retrieval (PIR) protocol in intent:
    the server cannot distinguish an encrypted query from a random bitstring, so it does
    not learn which genomic pattern the user is searching for.

    ⚠️ Caveat: A cryptographically rigorous PIR would require an oblivious oracle — i.e.,
    the oracle circuit itself should not reveal which state is being marked. The current
    implementation is a pedagogical demonstration, not a production-secure PIR scheme.

    Returns: (QuantumCircuit, iterations)
    """
    n_qubits = len(target_bitstring)
    N = 2 ** n_qubits
    iterations = int(np.floor(np.pi / 4 * np.sqrt(N)))

    qc = QuantumCircuit(n_qubits)

    # Initialize superposition
    qc.h(range(n_qubits))

    for _ in range(iterations):
        # --- Oracle: phase-flip the target state ---
        for i, bit in enumerate(reversed(target_bitstring)):
            if bit == '0':
                qc.x(i)
        qc.h(n_qubits - 1)
        qc.mcx(list(range(n_qubits - 1)), n_qubits - 1)
        qc.h(n_qubits - 1)
        for i, bit in enumerate(reversed(target_bitstring)):
            if bit == '0':
                qc.x(i)

        # --- Grover Diffuser: inversion about the mean ---
        qc.h(range(n_qubits))
        qc.x(range(n_qubits))
        qc.h(n_qubits - 1)
        qc.mcx(list(range(n_qubits - 1)), n_qubits - 1)
        qc.h(n_qubits - 1)
        qc.x(range(n_qubits))
        qc.h(range(n_qubits))

    qc.measure_all()
    return qc, iterations


@lru_cache(maxsize=64)
def build_grover_step_circuits(target_bitstring: str) -> list:
    """
    Builds 4 isolated QuantumCircuit objects — one per Grover phase — and
    returns their text-art diagrams along with a label and description.

    Steps:
      1. Initialization  – Hadamard superposition
      2. Oracle          – Phase-flip the target state (one iteration)
      3. Diffusion       – Inversion-about-the-mean (one iteration)
      4. Measurement     – Collapse all qubits to classical bits
    """
    n = len(target_bitstring)
    N = 2 ** n
    iters = int(np.floor(np.pi / 4 * np.sqrt(N)))
    x_count = sum(1 for bit in target_bitstring if bit == '0')

    # ── Build the full circuit with barriers for step segmentation ──────
    qc_full = QuantumCircuit(n)
    
    # Step 1: Initialization
    qc_full.h(range(n))
    qc_full.barrier()
    
    # Step 2: Oracle
    for i, bit in enumerate(reversed(target_bitstring)):
        if bit == '0':
            qc_full.x(i)
    qc_full.h(n - 1)
    qc_full.mcx(list(range(n - 1)), n - 1)
    qc_full.h(n - 1)
    for i, bit in enumerate(reversed(target_bitstring)):
        if bit == '0':
            qc_full.x(i)
    qc_full.barrier()
    
    # Step 3: Diffusion
    qc_full.h(range(n))
    qc_full.x(range(n))
    qc_full.h(n - 1)
    qc_full.mcx(list(range(n - 1)), n - 1)
    qc_full.h(n - 1)
    qc_full.x(range(n))
    qc_full.h(range(n))
    qc_full.barrier()
    
    # Step 4: Measurement
    qc_full.measure_all()
    
    # Draw the full circuit to text without line wrapping
    full_text = str(qc_full.draw(output='text', fold=-1))
    lines = full_text.split('\n')
    
    max_len = max(len(line) for line in lines)
    padded_lines = [line.ljust(max_len) for line in lines]
    
    # Analyze all lines to find the columns where barriers '░' are drawn.
    barrier_cols = set()
    for line in padded_lines:
        for i, char in enumerate(line):
            if char == '░':
                barrier_cols.add(i)
                
    barrier_cols = sorted(list(barrier_cols))
    
    # Group adjacent barrier columns into clusters
    final_barrier_splits = []
    if barrier_cols:
        cluster_end = barrier_cols[0]
        for c in barrier_cols[1:]:
            if c > cluster_end + 3:
                final_barrier_splits.append(cluster_end)
            cluster_end = c
        final_barrier_splits.append(cluster_end)
                
    partial_diagrams = []
    if len(final_barrier_splits) >= 3:
        col_starts = [0] + [c + 1 for c in final_barrier_splits]
        col_ends = final_barrier_splits + [max_len]
        
        for step_idx in range(4):
            if step_idx < len(col_starts):
                s_start = col_starts[step_idx]
                s_end = col_ends[step_idx]
            else:
                s_start, s_end = 0, max_len
                
            step_lines = []
            for line in padded_lines:
                if not line.strip():
                    step_lines.append(line)
                    continue
                step_lines.append(line[s_start:s_end])
                
            partial_diagrams.append("\n".join(step_lines))
    else:
        # Fallback if barriers not found properly
        partial_diagrams = [full_text] * 4

    # ── Statevector probability distributions at each cumulative step ─────
    # We build cumulative circuits to compute the actual quantum state after
    # each phase, then extract the probability for every basis state.
    from qiskit.quantum_info import Statevector

    basis_labels = [format(i, f'0{n}b') for i in range(N)]

    def _probs_dict(circuit):
        """Run statevector simulation and return {bitstring: probability}."""
        sv = Statevector.from_instruction(circuit)
        probs = sv.probabilities_dict()
        # Ensure every basis state is present (some may have 0 probability)
        return {label: round(probs.get(label, 0.0), 6) for label in basis_labels}

    # Cumulative circuit: Step 1 — Initialization (Hadamard on all qubits)
    qc_cumul = QuantumCircuit(n)
    qc_cumul.h(range(n))
    probs_step1 = _probs_dict(qc_cumul)

    # Cumulative circuit: Step 2 — + Oracle
    for i, bit in enumerate(reversed(target_bitstring)):
        if bit == '0':
            qc_cumul.x(i)
    qc_cumul.h(n - 1)
    qc_cumul.mcx(list(range(n - 1)), n - 1)
    qc_cumul.h(n - 1)
    for i, bit in enumerate(reversed(target_bitstring)):
        if bit == '0':
            qc_cumul.x(i)
    probs_step2 = _probs_dict(qc_cumul)

    # Cumulative circuit: Step 3 — + Diffusion
    qc_cumul.h(range(n))
    qc_cumul.x(range(n))
    qc_cumul.h(n - 1)
    qc_cumul.mcx(list(range(n - 1)), n - 1)
    qc_cumul.h(n - 1)
    qc_cumul.x(range(n))
    qc_cumul.h(range(n))
    probs_step3 = _probs_dict(qc_cumul)

    # Step 4 (Measurement) — probabilities are identical to post-diffusion
    probs_step4 = probs_step3

    return [
        {
            "step": 1,
            "label": "Initialization",
            "description": (
                f"Hadamard (H) gates applied to all {n} qubits create an equal superposition "
                f"of all N = 2\u207f = {N} basis states. Each state starts with amplitude 1/\u221a{N}."
            ),
            "diagram": partial_diagrams[0],
            "gates": [
                {
                    "name": "Hadamard (H)",
                    "symbol": "H",
                    "count": n,
                    "explanation": (
                        f"Applied to all {n} qubits simultaneously. "
                        f"Maps |0\u27e9 \u2192 (|0\u27e9+|1\u27e9)/\u221a2, placing each qubit in a 50/50 superposition. "
                        f"Together they create a uniform superposition of all {N} basis states — "
                        f"amplitude 1/\u221a{N} on each."
                    ),
                },
            ],
            "probabilities": probs_step1,
        },
        {
            "step": 2,
            "label": "Oracle",
            "description": (
                f"The oracle marks the target state |{target_bitstring}\u27e9 by flipping its phase "
                f"(amplitude × \u22121). X-gates convert |0\u27e9 qubits to |1\u27e9, a multi-controlled-Z "
                f"applies the phase flip, then X-gates restore the basis."
            ),
            "diagram": partial_diagrams[1],
            "gates": [
                {
                    "name": "Pauli-X (NOT)",
                    "symbol": "X",
                    "count": x_count * 2,
                    "explanation": (
                        f"{x_count} X gate(s) applied before and {x_count} after the controlled phase flip. "
                        "Flips '0' qubits to |1\u27e9 so the multi-controlled gate fires on the correct state. "
                        "Acts like a classical NOT: |0\u27e9 \u2192 |1\u27e9, |1\u27e9 \u2192 |0\u27e9."
                    ),
                },
                {
                    "name": "Hadamard (H)",
                    "symbol": "H",
                    "count": 2,
                    "explanation": (
                        "Sandwiches the MCX on the top qubit to convert it into a phase-flip (controlled-Z). "
                        "H\u00b7MCX\u00b7H on the target qubit implements a multi-controlled-Z without a native CZ gate."
                    ),
                },
                {
                    "name": "Multi-Controlled-X (MCX)",
                    "symbol": "MCX",
                    "count": 1,
                    "explanation": (
                        f"Flips qubit {n-1} only when all {n-1} control qubit(s) are |1\u27e9. "
                        f"With the H sandwich this realises a phase flip — the marked state gains amplitude \u00d7\u22121."
                    ),
                },
            ],
            "probabilities": probs_step2,
        },
        {
            "step": 3,
            "label": "Diffusion Operator",
            "description": (
                f"Inversion about the mean: H gates transform to the computational basis, "
                f"X gates flip all qubits, a multi-controlled-Z inverts the all-|1\u27e9 state, "
                f"then X and H restore. This amplifies the marked state's amplitude each iteration. "
                f"Repeated {iters}\u00d7 \u2248 \u03c0/4 \u00b7 \u221a{N}."
            ),
            "diagram": partial_diagrams[2],
            "gates": [
                {
                    "name": "Hadamard (H)",
                    "symbol": "H",
                    "count": n * 2,
                    "explanation": (
                        f"Applied to all {n} qubits twice (before and after the X+MCX layer). "
                        "The first H\u2297\u207f rotates into the computational basis; the final H\u2297\u207f undoes that, "
                        "completing the 2|\u03c8\u27e9\u27e8\u03c8|\u22121 reflection about the average amplitude."
                    ),
                },
                {
                    "name": "Pauli-X (NOT)",
                    "symbol": "X",
                    "count": n * 2,
                    "explanation": (
                        f"Applied to all {n} qubits before and after the MCX. "
                        "Converts the all-zeros state to all-ones so the MCX can target it, "
                        "implementing 'inversion about zero' as the inner core of the diffuser."
                    ),
                },
                {
                    "name": "Multi-Controlled-X (MCX)",
                    "symbol": "MCX",
                    "count": 1,
                    "explanation": (
                        f"With the H sandwich on qubit {n-1} implements a multi-controlled-Z on |1\u22ef1\u27e9. "
                        "After surrounding X gates this phase-flips the |0\u22ef0\u27e9 state, completing inversion about the mean."
                    ),
                },
            ],
            "probabilities": probs_step3,
        },
        {
            "step": 4,
            "label": "Measurement",
            "description": (
                f"After {iters} iteration(s) the target state |{target_bitstring}\u27e9 has "
                f"probability \u2248 1. Measuring all {n} qubits collapses the superposition "
                f"and returns the target bitstring with high probability."
            ),
            "diagram": partial_diagrams[3],
            "gates": [
                {
                    "name": "Measurement (M)",
                    "symbol": "\u29fa",
                    "count": n,
                    "explanation": (
                        f"Each of the {n} qubits is projected onto {{|0\u27e9, |1\u27e9}}. "
                        f"After {iters} Grover iteration(s) the amplitude of |{target_bitstring}\u27e9 is \u2248 1, "
                        "so the classical outcome matches the target with high probability. "
                        "This irreversible step converts quantum superposition into a classical bitstring."
                    ),
                },
            ],
            "probabilities": probs_step4,
        },
    ]
