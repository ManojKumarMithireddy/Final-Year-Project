"""
Pydantic request models for the BioQuantum API.
"""

import re
from pydantic import BaseModel, Field, EmailStr, field_validator
from typing import Optional

# Compiled once at module level for reuse
_UPPERCASE_RE = re.compile(r'[A-Z]')
_SPECIAL_RE   = re.compile(r'[0-9!@#$%^&*]')


def _validate_password_strength(value: str) -> str:
    """Shared password-strength rule: ≥8 chars, ≥1 uppercase, ≥1 digit or special char."""
    if len(value) < 8:
        raise ValueError("Password must be at least 8 characters long.")
    if not _UPPERCASE_RE.search(value):
        raise ValueError("Password must contain at least one uppercase letter.")
    if not _SPECIAL_RE.search(value):
        raise ValueError("Password must contain at least one digit or special character (!@#$%^&*).")
    return value


class GoogleLoginRequest(BaseModel):
    token: str

class StandardAuthRequest(BaseModel):
    """Used for login only — no password-strength enforcement so existing users aren't locked out."""
    email: EmailStr
    password: str
    name: Optional[str] = None

class RegisterRequest(BaseModel):
    """Used for new-account registration — enforces email format and password strength at schema level."""
    email: EmailStr
    password: str
    name: Optional[str] = None

    @field_validator("password")
    @classmethod
    def password_strength(cls, v: str) -> str:
        return _validate_password_strength(v)

class ResendVerificationRequest(BaseModel):
    email: EmailStr

class ForgotPasswordRequest(BaseModel):
    email: EmailStr

class ResetPasswordRequest(BaseModel):
    token: str
    new_password: str

    @field_validator("new_password")
    @classmethod
    def new_password_strength(cls, v: str) -> str:
        return _validate_password_strength(v)

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
    client_timestamp: Optional[str] = Field(default=None,
        description="ISO 8601 timestamp from the browser — stored in history if provided.")

class IBMSubmitRequest(BaseModel):
    """IBM credentials are fetched server-side from the authenticated user's saved credentials."""
    target_bits: str
    client_timestamp: Optional[str] = None

class IBMStatusRequest(BaseModel):
    """IBM credentials are fetched server-side from the authenticated user's saved credentials."""
    job_id: str

class BioQuantumRequest(BaseModel):
    """
    Request for the constrained BRCA1 Grover search PoC.
    n_codons controls both the NCBI fetch volume and the qubit count (n_qubits = n_codons * 6).
    has_mutation: True = carrier (c.5266dupC present); False = healthy control (no mutation).
    """
    n_codons: int = Field(default=1, ge=1, le=3,
        description="Number of codons per node (1=6 qubits, 2=12 qubits, 3=18 qubits).")
    has_mutation: bool = Field(default=True,
        description="True = patient carries c.5266dupC; False = healthy reference patient.")
    client_timestamp: Optional[str] = None

class BioIBMSubmitRequest(BaseModel):
    """IBM QPU variant — n_codons fixed to 1 (6 qubits) for free-tier compatibility."""
    n_codons: int = Field(default=1, ge=1, le=1,
        description="Fixed at 1 codon (6 qubits) for IBM free-tier QPU.")
    client_timestamp: Optional[str] = None
