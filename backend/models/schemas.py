"""
Pydantic request models for the BioQuantum API.
"""

from pydantic import BaseModel, Field
from typing import Optional


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
