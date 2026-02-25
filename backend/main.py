from fastapi import FastAPI, HTTPException, Depends
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from typing import Optional, List, Dict, Any
import time
import math
import numpy as np
import secrets
from Bio import Entrez, SeqIO
import io
import datetime

from db import get_db
from auth import verify_google_token, create_access_token, get_current_user, get_optional_user, verify_password, get_password_hash

#venv\Scripts\uvicorn main:app --reload --port 8000
# Qiskit Imports
from qiskit import QuantumCircuit, transpile
from qiskit_aer import AerSimulator
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager

# IBM Cloud Imports
try:
    from qiskit_ibm_runtime import QiskitRuntimeService, SamplerV2 as Sampler
    IBM_RUNTIME_AVAILABLE = True
except ImportError:
    IBM_RUNTIME_AVAILABLE = False

app = FastAPI(title="BioQuantum Hybrid Quantum Search API")

# Allow CORS for local React development
app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:5173", "http://127.0.0.1:5173", "http://localhost:5174", "http://127.0.0.1:5174"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# --- NCBI CONFIGURATION ---
Entrez.email = "bioquantum.capstone@example.com" 

# --- MODELS ---
class GoogleLoginRequest(BaseModel):
    token: str

class StandardAuthRequest(BaseModel):
    email: str
    password: str
    name: Optional[str] = None

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

class IBMSubmitRequest(BaseModel):
    target_bits: str
    api_key: str
    crn: str

class IBMStatusRequest(BaseModel):
    job_id: str
    api_key: str
    crn: str

# --- HELPER FUNCTIONS ---
def fetch_sequence(accession_id: str):
    try:
        handle = Entrez.efetch(db="nucleotide", id=accession_id, rettype="fasta", retmode="text")
        # Read the raw text
        raw_fasta = handle.read()
        handle.close()
        
        # Parse it
        record = SeqIO.read(io.StringIO(raw_fasta), "fasta")
        # Return up to 100,000 bases to keep payloads reasonable but large enough for the demo
        return str(record.seq[:100000]).upper()
    except Exception as e:
        print(f"NCBI Exception: {e}")
        return None

def dna_to_bits(dna_seq: str):
    mapping = {'A': '00', 'C': '01', 'G': '10', 'T': '11'}
    return "".join(mapping.get(base, '00') for base in dna_seq)


# --- AUTH ENDPOINTS ---

@app.post("/api/auth/register")
async def register_user(request: StandardAuthRequest):
    db = get_db()
    existing_user = await db.users.find_one({"email": request.email})
    if existing_user:
        raise HTTPException(status_code=400, detail="Email already registered")
        
    user = {
        "email": request.email,
        "name": request.name or request.email.split("@")[0],
        "hashed_password": get_password_hash(request.password),
        "auth_provider": "local",
        "created_at": datetime.datetime.utcnow()
    }
    await db.users.insert_one(user)
    
    access_token = create_access_token({"sub": user["email"], "name": user["name"]})
    return {"access_token": access_token, "user": {"email": user["email"], "name": user["name"]}}

@app.post("/api/auth/login")
async def login_user(request: StandardAuthRequest):
    db = get_db()
    user = await db.users.find_one({"email": request.email})
    
    if not user or "hashed_password" not in user:
        raise HTTPException(status_code=400, detail="Invalid email or password")
        
    if not verify_password(request.password, user["hashed_password"]):
        raise HTTPException(status_code=400, detail="Invalid email or password")
        
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

# --- ENDPOINTS ---

@app.get("/api/sequence/{accession}")
async def get_sequence(accession: str):
    """Fetches a sequence from NCBI."""
    seq = fetch_sequence(accession)
    if not seq:
        raise HTTPException(status_code=404, detail="Sequence not found or NCBI error.")
    return {"accession": accession, "length": len(seq), "sequence": seq}

@app.post("/api/search/classical")
async def classical_search(request: SearchRequest):
    """Performs a classical sliding window search and measures time."""
    dataset = request.dataset.upper()
    target = request.target.upper()
    k = len(target)
    n_windows = len(dataset) - k + 1
    
    if n_windows <= 0:
        raise HTTPException(status_code=400, detail="Target longer than dataset.")

    start_time = time.perf_counter()
    matches = []
    
    # Simple sliding window
    for i in range(n_windows):
        window = dataset[i:i+k]
        if window == target:
            matches.append(i)
            
    end_time = time.perf_counter()
    execution_time_ms = (end_time - start_time) * 1000
    
    # O(N) queries: comparing the window to the target n_windows times
    classical_queries = n_windows 
    
    return {
        "n_windows": n_windows,
        "matches": matches,
        "execution_time_ms": execution_time_ms,
        "classical_queries": classical_queries
    }

@app.post("/api/search/quantum-simulation")
async def quantum_simulation(request: QuantumSimulationRequest):
    """Calculates theoretical Grover metrics."""
    N = request.n_windows
    if N <= 0:
         return {"error": "Invalid window count"}
         
    # Theoretical Grover iterations (calls to oracle)
    quantum_queries = int(math.floor((math.pi / 4) * math.sqrt(N)))
    
    return {
        "n_windows": N,
        "quantum_queries": quantum_queries,
        "advantage_ratio": round(N / max(1, quantum_queries), 2) if quantum_queries > 0 else 0
    }

@app.post("/api/search/quantum-simulation-poc")
async def quantum_toy(request: QuantumToyRequest, current_user: Optional[dict] = Depends(get_optional_user)):
    """Executes a Grover search for proof of concept."""
    target_bitstring = request.target_bits
    
    if not target_bitstring or not all(c in '01' for c in target_bitstring):
        target_bitstring = "101" # fallback default
        
    n_qubits = len(target_bitstring)
    qc = QuantumCircuit(n_qubits)
    
    # Initialization
    qc.h(range(n_qubits))
    
    # Optimal Grover iterations for N=8
    N = 2 ** n_qubits
    iterations = int(np.floor(np.pi / 4 * np.sqrt(N)))
    
    start_time = time.perf_counter()
    
    for _ in range(iterations):
        # Oracle (Marks the state)
        for i, bit in enumerate(reversed(target_bitstring)):
            if bit == '0': qc.x(i)
        
        qc.h(n_qubits-1)
        qc.mcx(list(range(n_qubits-1)), n_qubits-1)
        qc.h(n_qubits-1)
        
        for i, bit in enumerate(reversed(target_bitstring)):
            if bit == '0': qc.x(i)
                
        # Diffuser
        qc.h(range(n_qubits))
        qc.x(range(n_qubits))
        qc.h(n_qubits-1)
        qc.mcx(list(range(n_qubits-1)), n_qubits-1)
        qc.h(n_qubits-1)
        qc.x(range(n_qubits))
        qc.h(range(n_qubits))

    qc.measure_all()
    
    # Simulation
    simulator = AerSimulator()
    transpiled_qc = transpile(qc, simulator)
    result = simulator.run(transpiled_qc, shots=1024).result()
    counts = result.get_counts()
    
    end_time = time.perf_counter()
    exec_time_ms = (end_time - start_time) * 1000
    
    measured_state = max(counts, key=counts.get)
    
    # Generate terminal-style text output for frontend
    circuit_diagram = str(qc.draw(output="text"))
    
    result_data = {
        "measured_state": measured_state,
        "counts": counts,
        "execution_time_ms": exec_time_ms,
        "iterations": iterations,
        "qubits": n_qubits,
        "circuit_diagram": circuit_diagram
    }
    
    if current_user:
        db = get_db()
        await db.history.insert_one({
            "email": current_user["sub"],
            "type": "quantum_local_sim",
            "target_bits": request.target_bits,
            "measured_state": measured_state,
            "execution_time_ms": exec_time_ms,
            "timestamp": datetime.datetime.utcnow()
        })
        
    return result_data

@app.post("/api/search/quantum-poc/ibm-submit")
async def ibm_submit_job(request: IBMSubmitRequest, current_user: Optional[dict] = Depends(get_optional_user)):
    """Submits a Grover search to IBM Cloud."""
    if not IBM_RUNTIME_AVAILABLE:
        raise HTTPException(status_code=500, detail="qiskit-ibm-runtime not installed")
        
    try:
        service = QiskitRuntimeService(
            channel="ibm_cloud", 
            token=request.api_key, 
            instance=request.crn
        )
        backend = service.least_busy(operational=True, simulator=False)
        
        target_bitstring = request.target_bits
        if not target_bitstring or not all(c in '01' for c in target_bitstring):
            target_bitstring = "101"
            
        n_qubits = len(target_bitstring)
        qc = QuantumCircuit(n_qubits)
        
        qc.h(range(n_qubits))
        N = 2 ** n_qubits
        iterations = int(np.floor(np.pi / 4 * np.sqrt(N)))
        
        for _ in range(iterations):
            for i, bit in enumerate(reversed(target_bitstring)):
                if bit == '0': qc.x(i)
            qc.h(n_qubits-1)
            qc.mcx(list(range(n_qubits-1)), n_qubits-1)
            qc.h(n_qubits-1)
            for i, bit in enumerate(reversed(target_bitstring)):
                if bit == '0': qc.x(i)
                    
            qc.h(range(n_qubits))
            qc.x(range(n_qubits))
            qc.h(n_qubits-1)
            qc.mcx(list(range(n_qubits-1)), n_qubits-1)
            qc.h(n_qubits-1)
            qc.x(range(n_qubits))
            qc.h(range(n_qubits))

        qc.measure_all()
        
        pm = generate_preset_pass_manager(backend=backend, optimization_level=3)
        isa_circuit = pm.run(qc)
        
        sampler = Sampler(mode=backend)
        job = sampler.run([isa_circuit])
        
        circuit_diagram = str(qc.draw(output="text"))
        
        job_data = {
            "job_id": job.job_id(),
            "backend": backend.name,
            "circuit_diagram": circuit_diagram,
            "status": job.status(),
            "iterations": iterations,
            "execution_time_ms": 0 # to match interface
        }
        
        if current_user:
            db = get_db()
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
async def ibm_job_status(request: IBMStatusRequest):
    """Checks the status of an IBM Cloud job."""
    if not IBM_RUNTIME_AVAILABLE:
        raise HTTPException(status_code=500, detail="qiskit-ibm-runtime not installed")
        
    try:
        service = QiskitRuntimeService(
            channel="ibm_cloud", 
            token=request.api_key, 
            instance=request.crn
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
                "execution_time_ms": 0 # to match interface
            }
        else:
            return {"status": status}
            
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
