"""
NCBI sequence fetching and DNA encoding utilities.
"""

import os
import re
import io
from Bio import Entrez, SeqIO

# --- ACCESSION VALIDATION ---
# Only allow safe NCBI accession formats (e.g. NM_007294, AF123456.1)
ACCESSION_RE = re.compile(r'^[A-Za-z0-9_.\-]{1,30}$')

# --- NCBI CONFIGURATION ---
# Read from env var to comply with NCBI terms of service.
# See: https://www.ncbi.nlm.nih.gov/books/NBK25497/#chapter2.Usage_Guidelines_and_Requiremen
Entrez.email = os.getenv("NCBI_EMAIL", "bioquantum.capstone@example.com")


def fetch_sequence(accession_id: str):
    """Fetches a DNA sequence from NCBI by accession ID (capped at 100k bases)."""
    try:
        handle = Entrez.efetch(db="nucleotide", id=accession_id, rettype="fasta", retmode="text")
        raw_fasta = handle.read()
        handle.close()
        record = SeqIO.read(io.StringIO(raw_fasta), "fasta")
        return str(record.seq[:100000]).upper()
    except Exception as e:
        print(f"NCBI Exception: {e}")
        return None


def dna_to_bits(dna_seq: str) -> str:
    """Encodes a DNA string to a binary string using 2-bit encoding (A=00, C=01, G=10, T=11)."""
    mapping = {'A': '00', 'C': '01', 'G': '10', 'T': '11'}
    return "".join(mapping.get(base, '00') for base in dna_seq)
