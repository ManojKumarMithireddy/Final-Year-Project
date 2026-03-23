## 🚀 Live Deployment

| Service  | URL |
|----------|-----|
| Frontend | https://bioquantum-frontend.vercel.app |
| Backend API | https://manojkumarmithireddy-bioquantum-api.hf.space/docs |

> **Backend secrets still needed:** Replace placeholder values in  
> [HF Space Settings → Secrets](https://huggingface.co/spaces/ManojKumarMithireddy/bioquantum-api/settings):  
> MONGODB_URI, GOOGLE_CLIENT_ID, BREVO_API_KEY, NCBI_EMAIL.

---
# BioQuantum Hybrid Search

A full-stack, production-patterned web application demonstrating a **hybrid classical-quantum genomic sliding-window search**. The system fetches real DNA sequences from NCBI, runs a classical O(N) search, and compares it against the theoretical O(√N) Grover's algorithm advantage — with a working Grover circuit Proof of Concept running on a local Qiskit simulator or IBM Cloud QPU.

---

## Architecture

```text
bioquantum_hybrid_search/
├── backend/              # FastAPI (Python)
│   ├── main.py           # App factory & Uvicorn entrypoint (run `python main.py`)
│   ├── auth.py           # JWT and password hashing logic
│   ├── crypto.py         # PIR OTP-XOR encryption
│   ├── db.py             # MongoDB async client (Motor)
│   ├── email_service.py  # Brevo API email integration
│   ├── rate_limit.py     # SlowAPI rate limiting configuration
│   ├── models/           # Data validation and database schemas
│   │   └── schemas.py    # Pydantic models
│   ├── routers/          # Modular API endpoints
│   │   ├── auth.py       # Login, register, Google OAuth
│   │   ├── credentials.py# IBM credential management
│   │   ├── quantum.py    # Grover endpoints & IBM submission
│   │   ├── search.py     # Classical search endpoints
│   │   └── user.py       # History and user profile endpoints
│   ├── services/         # Core business logic
│   │   ├── bio_grover.py # DNA sequence to quantum circuit mapping
│   │   ├── grover.py     # Underlying Qiskit circuit building
│   │   └── ncbi.py       # Biopython Entrez API fetch
│   ├── tests/            # Unit testing
│   │   └── test_core.py  # Pytest core logic tests
│   └── venv/             # Python virtual environment (not committed)
└── frontend/             # React + Vite + TailwindCSS
    └── src/
        ├── pages/        # Main route views
        │   ├── Dashboard.jsx     # Tabbed router interface
        │   ├── History.jsx       # Quantum job history
        │   ├── Login.jsx         # Auth page (Google + Email)
        │   ├── ResetPassword.jsx # Password recovery
        │   └── VerifyEmail.jsx   # Email verification logic
        └── components/   # UI Building Blocks
            ├── AmplitudeChart.jsx      # Dynamic Statevector probability visualization
            ├── BackendSleepBanner.jsx  # Cold start notification for free-tier hosting
            ├── ComplexityChart.jsx     # Recharts O(N) vs O(√N) line chart
            ├── ErrorBoundary.jsx       # Global React error boundaries
            ├── GroverPOC.jsx           # Main Interactive Grover circuit & PIR pipeline
            ├── GroverStepNavigator.jsx # Step-by-step interactive circuit composition
            ├── HybridSearchPanel.jsx   # NCBI fetch + classical/quantum comparison
            └── SessionTimeoutModal.jsx # Secure session auto-logout management
```

---

## Prerequisites

| Tool | Version |
|------|---------|
| Python | 3.10+ |
| Node.js | 18+ |
| MongoDB | 6+ (local or Atlas) |

---

## Setup & Running

### 1. Backend

```powershell
cd backend

# Create and activate virtual environment
python -m venv venv
.\venv\Scripts\activate

# Install dependencies
pip install fastapi uvicorn[standard] motor python-jose[cryptography] passlib[bcrypt] \
            google-auth biopython qiskit qiskit-aer qiskit-ibm-runtime \
            python-dotenv ntplib pytest

# Configure environment
copy .env.example .env   # then edit with your values
```

**`backend/.env`** (create manually):
```env
MONGODB_URI=mongodb://localhost:27017
JWT_SECRET=your_strong_secret_here
JWT_ALGORITHM=HS256
GOOGLE_CLIENT_ID=your_google_client_id.apps.googleusercontent.com
```

```powershell
# Start the backend API server locally (Uvicorn starts automatically)
python main.py
```

### 2. Frontend

```powershell
cd frontend
npm install

# Configure environment (already set up)
# Edit VITE_GOOGLE_CLIENT_ID in frontend/.env
```

**`frontend/.env`**:
```env
VITE_GOOGLE_CLIENT_ID=your_google_client_id.apps.googleusercontent.com
VITE_API_BASE=http://localhost:8000/api
```

```powershell
npm run dev
# Open http://localhost:5173
```

---

## Running Tests

```powershell
cd backend
.\venv\Scripts\activate
python -m pytest tests/ -v
```

Expected output:
```
tests/test_core.py::test_dna_to_bits                  PASSED
tests/test_core.py::test_grover_iteration_formula[...] PASSED  (×4 parametrized)
tests/test_core.py::test_build_grover_circuit_structure PASSED
tests/test_core.py::test_password_hash_roundtrip       PASSED
tests/test_core.py::test_accession_sanitization        PASSED
```

---

## Key Features

| Feature | Description |
|---------|-------------|
| **NCBI Integration** | Live DNA sequences via Biopython `Entrez.efetch` |
| **Classical Search** | Timed O(N) sliding-window search in Python |
| **Modular Backend** | Clean FastAPI architecture (split into routers, models, services) |
| **Grover POC & PIR** | Executable Grover circuit with OTP-XOR client-side blind querying |
| **Circuit UI** | Step-by-step interactive composition with ASCII visualization & explanations |
| **Amplitude Viz** | Dynamic statevector simulation displaying basis state probabilities |
| **Noise Model** | Optional depolarizing noise (0–5%) on local simulator to model NISQ hardware |
| **UX & Security** | Modern UI with session timeouts, loaders, and dual auth (Google OAuth + JWT) |
| **IBM Security** | IBM API credentials stored safely server-side |

---

## Quantum Scale Limitations

> ⚠️ **Important for academic honesty:**

This project demonstrates Grover's algorithm at a **small scale (2–6 qubits, searching 4–64 states)**. The "quantum advantage" chart is **theoretically derived**, not experimentally measured.

Real genomic-scale Grover's search over ~100,000 windows would require:
- **~17 qubits** (log₂(100,000))
- **Fault-tolerant error correction** (not available on current NISQ hardware)
- Multi-controlled gates that decompose super-polynomially without error correction

The Proof of Concept tab demonstrates the **algorithm structure and privacy protocol** — not full quantum speedup on real genomic data.

---

## PIR Privacy Model

The "blind query" pipeline (Tab 2) simulates a Privacy-Preserving Information Retrieval protocol:

1. **Client** XOR-encrypts query bits with a random One-Time Pad key
2. **Server** runs Grover's on the encrypted bits (cannot determine original query)
3. **Client** XOR-decrypts the returned state to recover the result

**Caveat:** A cryptographically rigorous PIR protocol requires an *oblivious oracle* — the circuit itself should not reveal which state is being amplified. The current implementation is a pedagogical approximation of PIR intent.

---

## API Reference

| Method | Endpoint | Description |
|--------|----------|-------------|
| POST | `/api/auth/register` | Email/password registration |
| POST | `/api/auth/login` | Email/password login |
| POST | `/api/auth/google` | Google OAuth login |
| GET | `/api/sequence/{accession}` | Fetch DNA from NCBI |
| POST | `/api/search/classical` | Sliding-window classical search |
| POST | `/api/search/quantum-simulation` | Theoretical Grover metrics |
| POST | `/api/search/quantum-simulation-poc` | Execute local Grover circuit |
| POST | `/api/search/quantum-poc/ibm-submit` | Submit job to IBM Cloud QPU |
| POST | `/api/search/quantum-poc/ibm-status` | Check IBM job status |
| GET | `/api/history` | User's simulation history |

---

## References

- Grover, L. K. (1996). A fast quantum mechanical algorithm for database search. *STOC '96*.
- Nielsen & Chuang — *Quantum Computation and Quantum Information*, Chapter 6.
- [Qiskit Textbook](https://learning.quantum.ibm.com/)
- [NCBI Entrez API](https://www.ncbi.nlm.nih.gov/books/NBK25497/)
