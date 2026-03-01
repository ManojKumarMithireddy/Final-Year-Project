"""
Credentials router — save / retrieve IBM Cloud credentials.
"""

import datetime
from fastapi import APIRouter, Depends

from db import get_db
from auth import get_current_user
from models.schemas import SaveCredentialsRequest

router = APIRouter(prefix="/api/credentials", tags=["credentials"])


@router.post("")
async def save_credentials(request: SaveCredentialsRequest, current_user: dict = Depends(get_current_user)):
    db = get_db()
    await db.credentials.update_one(
        {"email": current_user["sub"]},
        {"$set": {
            "ibm_api_key": request.api_key,
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
    return creds
