from fastapi import FastAPI, HTTPException, Depends
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel, Field
from typing import Optional
import time
import math
import re
import os
import numpy as np
from Bio import Entrez, SeqIO
import io
import datetime
from itsdangerous import URLSafeTimedSerializer, SignatureExpired, BadSignature

from db import get_db
from auth import verify_google_token, create_access_token, get_current_user, get_optional_user, verify_password, get_password_hash
from email_service import send_verification_email, send_resend_verification_email, send_password_reset_email

# ── Email-verification token (separate secret + 1-hour TTL) ───────────────────
VERIFY_SECRET = os.getenv("JWT_SECRET", "verify_secret_fallback") + "_email_verify"
_verify_serializer = URLSafeTimedSerializer(VERIFY_SECRET)

def _make_verify_token(email: str) -> str:
    return _verify_serializer.dumps(email, salt="email-verify")

def _check_verify_token(token: str, max_age: int = 3600) -> str:
    """Returns the email embedded in the token, raises HTTPException on failure."""
    try:
        email = _verify_serializer.loads(token, salt="email-verify", max_age=max_age)
        return email
    except SignatureExpired:
        raise HTTPException(status_code=400, detail="Verification link has expired. Please request a new one.")
    except BadSignature:
        raise HTTPException(status_code=400, detail="Invalid verification token.")

#venv\Scripts\uvicorn main:app --reload --port 8000

# Qiskit Imports
from qiskit import QuantumCircuit, transpile
from qiskit_aer import AerSimulator
from qiskit_aer.noise import NoiseModel, depolarizing_error
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager

# IBM Cloud Imports
try:
    from qiskit_ibm_runtime import QiskitRuntimeService, SamplerV2 as Sampler
    IBM_RUNTIME_AVAILABLE = True
except ImportError:
    IBM_RUNTIME_AVAILABLE = False

# NTP clock check (optional dependency)
try:
    import ntplib
    NTPLIB_AVAILABLE = True
except ImportError:
    NTPLIB_AVAILABLE = False

app = FastAPI(title="BioQuantum Hybrid Quantum Search API")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:5173", "http://127.0.0.1:5173", "http://localhost:5174", "http://127.0.0.1:5174"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# --- ACCESSION VALIDATION ---
# Only allow safe NCBI accession formats (e.g. NM_007294, AF123456.1)
ACCESSION_RE = re.compile(r'^[A-Za-z0-9_.\-]{1,30}$')

# --- NCBI CONFIGURATION ---
# NOTE: Replace with your real email to comply with NCBI terms of service.
# See: https://www.ncbi.nlm.nih.gov/books/NBK25497/#chapter2.Usage_Guidelines_and_Requiremen
Entrez.email = "bioquantum.capstone@example.com"


@app.on_event("startup")
async def check_system_clock():
    """
    Pre-flight check: verify system clock is close to NTP time.
    A clock skew > 5s can cause Google OAuth token validation to fail
    with 'Token used too early' errors even with clock_skew_in_seconds set.
    """
    if NTPLIB_AVAILABLE:
        try:
            c = ntplib.NTPClient()
            response = c.request('pool.ntp.org', version=3)
            from datetime import timezone
            ntp_time = datetime.datetime.fromtimestamp(response.tx_time, tz=timezone.utc)
            local_time = datetime.datetime.now(tz=timezone.utc)
            delta = abs((ntp_time - local_time).total_seconds())
            if delta > 5:
                print(f"⚠️  WARNING: System clock is {delta:.1f}s off NTP. Google token auth may fail.")
            else:
                print(f"✅  Clock sync OK: {delta:.2f}s from NTP (pool.ntp.org)")
        except Exception as e:
            print(f"⚠️  NTP check failed (non-critical): {e}")
    else:
        print("ℹ️  ntplib not installed — clock skew check skipped. Install with: pip install ntplib")


# --- MODELS ---

class GoogleLoginRequest(BaseModel):
    token: str

class StandardAuthRequest(BaseModel):
    email: str
    password: str
    name: Optional[str] = None

class ResendVerificationRequest(BaseModel):
    email: str

class ForgotPasswordRequest(BaseModel):
    email: str

class ResetPasswordRequest(BaseModel):
    token: str
    new_password: str

class SaveCredentialsRequest(BaseModel):
    api_key: str
    crn: str

class SearchRequest(BaseModel):
    dataset: str
    target: str

class QuantumSimulationRequest(BaseModel):
    n_windows: int

class QuantumToyRequest(BaseModel):
    target_bits: str
    noise_level: float = Field(default=0.0, ge=0.0, le=0.05,
        description="Depolarizing noise probability per gate (0.0 = ideal, 0.05 = 5% error rate).")

class IBMSubmitRequest(BaseModel):
    """IBM credentials are fetched server-side from the authenticated user's saved credentials."""
    target_bits: str

class IBMStatusRequest(BaseModel):
    """IBM credentials are fetched server-side from the authenticated user's saved credentials."""
    job_id: str


# --- CORE QUANTUM HELPER ---

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


# --- DNA HELPER ---
def fetch_sequence(accession_id: str):
    try:
        handle = Entrez.efetch(db="nucleotide", id=accession_id, rettype="fasta", retmode="text")
        raw_fasta = handle.read()
        handle.close()
        record = SeqIO.read(io.StringIO(raw_fasta), "fasta")
        return str(record.seq[:100000]).upper()
    except Exception as e:
        print(f"NCBI Exception: {e}")
        return None

def dna_to_bits(dna_seq: str):
    mapping = {'A': '00', 'C': '01', 'G': '10', 'T': '11'}
    return "".join(mapping.get(base, '00') for base in dna_seq)


# --- AUTH ENDPOINTS ---

@app.post("/api/auth/register", status_code=202)
async def register_user(request: StandardAuthRequest):
    """
    Register a new local user.
    - Saves the user with is_verified=False.
    - Sends a verification email via Brevo.
    - Does NOT return a JWT — the user must verify their email first.
    """
    db = get_db()
    existing_user = await db.users.find_one({"email": request.email})
    if existing_user:
        raise HTTPException(status_code=400, detail="Email already registered")

    name = request.name or request.email.split("@")[0]
    user = {
        "email": request.email,
        "name": name,
        "hashed_password": get_password_hash(request.password),
        "auth_provider": "local",
        "is_verified": False,
        "created_at": datetime.datetime.utcnow()
    }
    await db.users.insert_one(user)

    # Build signed token and send verification email
    token = _make_verify_token(request.email)
    app_base = os.getenv("APP_BASE_URL", "http://localhost:5173")
    verify_url = f"{app_base}/verify-email?token={token}"
    try:
        send_verification_email(request.email, name, verify_url)
    except RuntimeError as e:
        # Don't block registration if email fails — surface a clear warning
        print(f"⚠️  Email send failed: {e}")
        return {
            "message": "Account created but email delivery failed. Check BREVO_API_KEY / BREVO_SENDER_EMAIL in .env.",
            "verify_url_debug": verify_url   # handy for local dev without a real key
        }

    return {"message": f"Verification email sent to {request.email}. Please check your inbox."}


@app.get("/api/auth/verify-email")
async def verify_email(token: str):
    """
    Validates the signed token from the verification link.
    On success: sets is_verified=True and returns a JWT so the user is
    immediately signed in without a separate login step.
    """
    email = _check_verify_token(token)   # raises HTTPException on bad/expired token

    db = get_db()
    user = await db.users.find_one({"email": email})
    if not user:
        raise HTTPException(status_code=404, detail="User not found.")
    if user.get("is_verified"):
        # Already verified — still hand back a JWT (idempotent, user may click link twice)
        access_token = create_access_token({"sub": user["email"], "name": user.get("name")})
        return {"access_token": access_token, "user": {"email": user["email"], "name": user.get("name")}}

    await db.users.update_one({"email": email}, {"$set": {"is_verified": True}})
    access_token = create_access_token({"sub": user["email"], "name": user.get("name")})
    return {"access_token": access_token, "user": {"email": user["email"], "name": user.get("name")}}


@app.post("/api/auth/resend-verification")
async def resend_verification(request: ResendVerificationRequest):
    """
    Regenerates and resends the verification email.
    Returns a generic 200 regardless of whether the email exists
    to prevent user enumeration.
    """
    db = get_db()
    user = await db.users.find_one({"email": request.email})
    if user and not user.get("is_verified") and user.get("auth_provider") == "local":
        token = _make_verify_token(request.email)
        app_base = os.getenv("APP_BASE_URL", "http://localhost:5173")
        verify_url = f"{app_base}/verify-email?token={token}"
        try:
            send_resend_verification_email(request.email, user.get("name", "User"), verify_url)
        except RuntimeError as e:
            raise HTTPException(status_code=500, detail=str(e))
    return {"message": "If that email is registered and unverified, a new link has been sent."}


# ── FORGOT / RESET PASSWORD ────────────────────────────────────────────────────

# Separate serializer salt keeps reset tokens distinct from email-verify tokens
_RESET_SALT = "password-reset"


@app.post("/api/auth/forgot-password")
async def forgot_password(request: ForgotPasswordRequest):
    """
    Generates a signed, time-limited (1 h) password-reset link and emails it.
    Always returns 200 to prevent user enumeration.
    """
    db = get_db()
    user = await db.users.find_one({"email": request.email})
    # Only send if the user exists, is a local (password) account, and is verified
    if user and user.get("auth_provider") == "local" and user.get("hashed_password"):
        token = _verify_serializer.dumps(request.email, salt=_RESET_SALT)
        app_base = os.getenv("APP_BASE_URL", "http://localhost:5173")
        reset_url = f"{app_base}/reset-password?token={token}"
        try:
            send_password_reset_email(request.email, user.get("name", "User"), reset_url)
        except RuntimeError as e:
            print(f"⚠️  Password reset email failed: {e}")
            # Surface URL in dev so users can still reset without Brevo configured
            return {
                "message": "If that email is registered, a reset link has been sent.",
                "reset_url_debug": reset_url,
            }
    return {"message": "If that email is registered, a reset link has been sent."}


@app.post("/api/auth/reset-password")
async def reset_password(request: ResetPasswordRequest):
    """
    Validates the signed token and updates the user's password.
    """
    try:
        email = _verify_serializer.loads(request.token, salt=_RESET_SALT, max_age=3600)
    except SignatureExpired:
        raise HTTPException(status_code=400, detail="Reset link has expired. Please request a new one.")
    except BadSignature:
        raise HTTPException(status_code=400, detail="Invalid or tampered reset token.")

    # Basic server-side password strength check
    import re as _re
    pwd = request.new_password
    if len(pwd) < 8 or not _re.search(r'[A-Z]', pwd) or not _re.search(r'[0-9!@#$%^&*]', pwd):
        raise HTTPException(
            status_code=422,
            detail="Password must be at least 8 characters and include an uppercase letter and a number or special character."
        )

    db = get_db()
    user = await db.users.find_one({"email": email})
    if not user or user.get("auth_provider") != "local":
        raise HTTPException(status_code=404, detail="User not found.")

    new_hash = get_password_hash(pwd)
    await db.users.update_one({"email": email}, {"$set": {"hashed_password": new_hash}})
    return {"message": "Password updated successfully. You can now sign in with your new password."}


@app.post("/api/auth/login")
async def login_user(request: StandardAuthRequest):
    db = get_db()
    user = await db.users.find_one({"email": request.email})
    if not user or "hashed_password" not in user:
        raise HTTPException(status_code=400, detail="Invalid email or password")
    if not verify_password(request.password, user["hashed_password"]):
        raise HTTPException(status_code=400, detail="Invalid email or password")
    # Block login until email is verified
    if not user.get("is_verified", False):
        raise HTTPException(status_code=403, detail="EMAIL_NOT_VERIFIED")
    access_token = create_access_token({"sub": user["email"], "name": user.get("name")})
    return {"access_token": access_token, "user": {"email": user["email"], "name": user.get("name")}}

@app.post("/api/auth/google")
async def google_login(request: GoogleLoginRequest):
    try:
        idinfo = verify_google_token(request.token)
        email = idinfo.get("email")
        name = idinfo.get("name")
        picture = idinfo.get("picture")

        db = get_db()
        user = await db.users.find_one({"email": email})
        if not user:
            user = {
                "email": email,
                "name": name,
                "picture": picture,
                "auth_provider": "google",
                "created_at": datetime.datetime.utcnow()
            }
            await db.users.insert_one(user)

        access_token = create_access_token({"sub": email, "name": name})
        return {"access_token": access_token, "user": {"email": email, "name": name, "picture": picture}}
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))

@app.get("/api/user/me")
async def get_me(current_user: dict = Depends(get_current_user)):
    db = get_db()
    user = await db.users.find_one({"email": current_user["sub"]}, {"_id": 0})
    return user

@app.post("/api/credentials")
async def save_credentials(request: SaveCredentialsRequest, current_user: dict = Depends(get_current_user)):
    db = get_db()
    await db.credentials.update_one(
        {"email": current_user["sub"]},
        {"$set": {"ibm_api_key": request.api_key, "ibm_crn": request.crn, "updated_at": datetime.datetime.utcnow()}},
        upsert=True
    )
    return {"status": "success"}

@app.get("/api/credentials")
async def get_credentials(current_user: dict = Depends(get_current_user)):
    db = get_db()
    creds = await db.credentials.find_one({"email": current_user["sub"]}, {"_id": 0})
    if not creds:
        return {"ibm_api_key": "", "ibm_crn": ""}
    return creds

@app.get("/api/history")
async def get_history(current_user: dict = Depends(get_current_user)):
    db = get_db()
    cursor = db.history.find({"email": current_user["sub"]}, {"_id": 0}).sort("timestamp", -1).limit(50)
    history = await cursor.to_list(length=50)
    return history


# --- SEQUENCE & SEARCH ENDPOINTS ---

@app.get("/api/sequence/{accession}")
async def get_sequence(accession: str):
    """Fetches a DNA sequence from NCBI by accession ID."""
    if not ACCESSION_RE.match(accession):
        raise HTTPException(
            status_code=400,
            detail="Invalid accession format. Use only letters, numbers, dots, dashes, or underscores (max 30 chars)."
        )
    seq = fetch_sequence(accession)
    if not seq:
        raise HTTPException(status_code=404, detail="Sequence not found or NCBI error.")
    return {"accession": accession, "length": len(seq), "sequence": seq}

@app.post("/api/search/classical")
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
        "classical_queries": n_windows
    }

@app.post("/api/search/quantum-simulation")
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
        "advantage_ratio": round(N / max(1, quantum_queries), 2)
    }

@app.post("/api/search/quantum-simulation-poc")
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
        "noise_level": request.noise_level
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
            "timestamp": datetime.datetime.utcnow()
        })

    return result_data

@app.post("/api/search/quantum-poc/ibm-submit")
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
            detail="No IBM Cloud credentials saved. Please save your API key and CRN in the credentials panel first."
        )

    try:
        service = QiskitRuntimeService(
            channel="ibm_cloud",
            token=creds["ibm_api_key"],
            instance=creds["ibm_crn"]
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
            "execution_time_ms": 0
        }

        await db.history.insert_one({
            "email": current_user["sub"],
            "type": "quantum_ibm_submit",
            "target_bits": request.target_bits,
            "job_id": job.job_id(),
            "backend": backend.name,
            "timestamp": datetime.datetime.utcnow()
        })

        return job_data
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.post("/api/search/quantum-poc/ibm-status")
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
            token=creds["ibm_api_key"],
            instance=creds["ibm_crn"]
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
                "execution_time_ms": 0
            }
        return {"status": status}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
