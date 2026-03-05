"""
User router -- profile (/me) and search history endpoints.

Security:
- /me uses an explicit allowlist projection to prevent leaking hashed_password
  and other internal fields to the client.
- /history excludes the email field from returned records (privacy) and
  includes _id serialised as a string so the frontend can use it as a stable
  React key.
"""

import logging
from fastapi import APIRouter, Depends, Query

from db import get_db
from auth import get_current_user

logger = logging.getLogger(__name__)

router = APIRouter(prefix="/api", tags=["user"])


@router.get("/user/me")
async def get_me(current_user: dict = Depends(get_current_user)):
    db = get_db()
    # Explicit allowlist -- never return hashed_password, auth_provider internals, etc.
    user = await db.users.find_one(
        {"email": current_user["sub"]},
        {"_id": 0, "email": 1, "name": 1, "picture": 1, "created_at": 1, "is_verified": 1},
    )
    return user


@router.get("/history")
async def get_history(
    current_user: dict = Depends(get_current_user),
    skip: int = Query(default=0, ge=0, description="Records to skip (pagination)."),
    limit: int = Query(default=20, ge=1, le=50, description="Max records per page."),
):
    db = get_db()
    # Fetch one extra to determine whether more pages exist.
    cursor = db.history.find(
        {"email": current_user["sub"]},
        {"email": 0},
    ).sort("timestamp", -1).skip(skip).limit(limit + 1)
    raw = await cursor.to_list(length=limit + 1)
    has_more = len(raw) > limit
    history = raw[:limit]
    # Serialise ObjectId -> string so the JSON encoder does not choke.
    for doc in history:
        if "_id" in doc:
            doc["_id"] = str(doc["_id"])
    return {"items": history, "has_more": has_more, "skip": skip}