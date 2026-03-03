---
title: BioQuantum API
emoji: 🧬
colorFrom: blue
colorTo: indigo
sdk: docker
app_port: 7860
pinned: false
---

# BioQuantum — Hybrid Quantum Search API

FastAPI backend for the BioQuantum hybrid classical-quantum genomic search platform.

- Grover's algorithm proof-of-concept (Qiskit + Qiskit-Aer local simulator)
- NCBI Entrez DNA sequence fetch (Biopython)
- JWT + Google OAuth authentication
- MongoDB (Motor async driver)
- Optional IBM Cloud QPU submission

**Frontend:** deployed on Vercel  
**Stack:** Python 3.11 · FastAPI · Qiskit 2.x · Qiskit-Aer · Motor · PyJWT · Cryptography
