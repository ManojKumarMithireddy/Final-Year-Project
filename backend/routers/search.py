"""
Search router — classical sliding window, theoretical quantum simulation,
and NCBI sequence fetching endpoints.
"""

import time
import math
from fastapi import APIRouter, HTTPException

from services.ncbi import fetch_sequence, ACCESSION_RE
from models.schemas import SearchRequest, QuantumSimulationRequest

router = APIRouter(prefix="/api", tags=["search"])


@router.get("/sequence/{accession}")
async def get_sequence(accession: str):
    """Fetches a DNA sequence from NCBI by accession ID."""
    if not ACCESSION_RE.match(accession):
        raise HTTPException(
            status_code=400,
            detail="Invalid accession format. Use only letters, numbers, dots, dashes, or underscores (max 30 chars).",
        )
    seq = fetch_sequence(accession)
    if not seq:
        raise HTTPException(status_code=404, detail="Sequence not found or NCBI error.")
    return {"accession": accession, "length": len(seq), "sequence": seq}


@router.post("/search/classical")
async def classical_search(request: SearchRequest):
    """Performs a classical sliding window search and measures execution time."""
    dataset = request.dataset.upper()
    target = request.target.upper()
    k = len(target)
    n_windows = len(dataset) - k + 1

    if n_windows <= 0:
        raise HTTPException(status_code=400, detail="Target longer than dataset.")

    start_time = time.perf_counter()
    matches = []
    for i in range(n_windows):
        if dataset[i:i + k] == target:
            matches.append(i)
    end_time = time.perf_counter()

    return {
        "n_windows": n_windows,
        "matches": matches,
        "execution_time_ms": (end_time - start_time) * 1000,
        "classical_queries": n_windows,
    }


@router.post("/search/quantum-simulation")
async def quantum_simulation(request: QuantumSimulationRequest):
    """
    Calculates THEORETICAL Grover metrics.
    Note: These are analytically derived values (O(√N) oracle calls), NOT measured
    from an actual quantum circuit of size N. The chart in the frontend visualizes
    this theoretical scaling advantage as a comparison, not an experimental result.
    """
    N = request.n_windows
    if N <= 0:
        return {"error": "Invalid window count"}

    quantum_queries = int(math.floor((math.pi / 4) * math.sqrt(N)))
    return {
        "n_windows": N,
        "quantum_queries": quantum_queries,
        "advantage_ratio": round(N / max(1, quantum_queries), 2),
    }
