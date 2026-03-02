"""
NCBI sequence fetching and DNA encoding utilities.
"""

import os
import re
import io
import logging
from Bio import Entrez, SeqIO

logger = logging.getLogger(__name__)

# --- ACCESSION VALIDATION ---
# Only allow safe NCBI accession formats (e.g. NM_007294, AF123456.1)
ACCESSION_RE = re.compile(r'^[A-Za-z0-9_.\-]{1,30}$')

# --- NCBI CONFIGURATION ---
# NCBI Terms of Service require a real contact email.
# Warn loudly at startup if the placeholder is still in use.
_ncbi_email = os.getenv("NCBI_EMAIL", "")
if not _ncbi_email:
    logger.warning(
        "NCBI_EMAIL is not set. NCBI Terms of Service require a real email address. "
        "Set NCBI_EMAIL in your .env file."
    )
    _ncbi_email = "bioquantum.capstone@example.com"  # fallback for dev only
Entrez.email = _ncbi_email

# 10-second timeout on all NCBI requests to prevent hung workers.
Entrez.timeout = int(os.getenv("NCBI_TIMEOUT", "10"))


def fetch_sequence(accession_id: str):
    """Fetches a DNA sequence from NCBI by accession ID (capped at 100k bases)."""
    try:
        handle = Entrez.efetch(db="nucleotide", id=accession_id, rettype="fasta", retmode="text")
        raw_fasta = handle.read()
        handle.close()
        record = SeqIO.read(io.StringIO(raw_fasta), "fasta")
        return str(record.seq[:100000]).upper()
    except Exception as e:
        logger.error("NCBI fetch failed for accession '%s': %s", accession_id, e)
        return None


def dna_to_bits(dna_seq: str) -> str:
    """Encodes a DNA string to a binary string using 2-bit encoding (A=00, C=01, G=10, T=11)."""
    mapping = {'A': '00', 'C': '01', 'G': '10', 'T': '11'}
    return "".join(mapping.get(base, '00') for base in dna_seq)
