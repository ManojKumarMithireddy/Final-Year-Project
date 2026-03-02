"""
Credentials router — save / retrieve IBM Cloud credentials.

Security notes:
- IBM API keys are encrypted at rest using Fernet symmetric encryption before
  being written to MongoDB. Set FERNET_KEY in .env to enable encryption.
- GET /credentials returns only a redacted sentinel so the raw key never
  leaves the server after the initial save.
"""

import datetime
from fastapi import APIRouter, Depends

from db import get_db
from auth import get_current_user
from models.schemas import SaveCredentialsRequest
from crypto import encrypt_value, decrypt_value

router = APIRouter(prefix="/api/credentials", tags=["credentials"])


@router.post("")
async def save_credentials(request: SaveCredentialsRequest, current_user: dict = Depends(get_current_user)):
    db = get_db()
    await db.credentials.update_one(
        {"email": current_user["sub"]},
        {"$set": {
            "ibm_api_key": encrypt_value(request.api_key),
            "ibm_crn": request.crn,
            "updated_at": datetime.datetime.now(datetime.timezone.utc),
        }},
        upsert=True,
    )
    return {"status": "success"}


@router.get("")
async def get_credentials(current_user: dict = Depends(get_current_user)):
    db = get_db()
    creds = await db.credentials.find_one({"email": current_user["sub"]}, {"_id": 0})
    if not creds:
        return {"ibm_api_key": "", "ibm_crn": ""}
    # Return a redacted sentinel so the raw key never leaves the server.
    # The frontend only needs to know whether a key is saved, not the key itself.
    has_key = bool(creds.get("ibm_api_key"))
    return {
        "ibm_api_key": "••••••••" if has_key else "",
        "ibm_crn": creds.get("ibm_crn", ""),
    }
