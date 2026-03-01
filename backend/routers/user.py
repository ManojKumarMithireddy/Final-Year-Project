"""
User router — profile (/me) and search history endpoints.

SECURITY FIX: /me uses an explicit allowlist projection to prevent
leaking hashed_password and other internal fields to the client.
"""

from fastapi import APIRouter, Depends

from db import get_db
from auth import get_current_user

router = APIRouter(prefix="/api", tags=["user"])


@router.get("/user/me")
async def get_me(current_user: dict = Depends(get_current_user)):
    db = get_db()
    # Explicit allowlist — never return hashed_password, auth_provider internals, etc.
    user = await db.users.find_one(
        {"email": current_user["sub"]},
        {"_id": 0, "email": 1, "name": 1, "picture": 1, "created_at": 1, "is_verified": 1},
    )
    return user


@router.get("/history")
async def get_history(current_user: dict = Depends(get_current_user)):
    db = get_db()
    cursor = db.history.find({"email": current_user["sub"]}, {"_id": 0}).sort("timestamp", -1).limit(50)
    history = await cursor.to_list(length=50)
    return history
