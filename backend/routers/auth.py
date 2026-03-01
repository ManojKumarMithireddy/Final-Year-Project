"""
Auth router — registration, login, Google OAuth, email verification,
forgot / reset password endpoints.
"""

import os
import re
import datetime
from fastapi import APIRouter, HTTPException

from itsdangerous import URLSafeTimedSerializer, SignatureExpired, BadSignature

from db import get_db
from auth import (
    verify_google_token,
    create_access_token,
    get_password_hash,
    verify_password,
)
from models.schemas import (
    GoogleLoginRequest,
    StandardAuthRequest,
    ResendVerificationRequest,
    ForgotPasswordRequest,
    ResetPasswordRequest,
)
from email_service import (
    send_verification_email,
    send_resend_verification_email,
    send_password_reset_email,
)

router = APIRouter(prefix="/api/auth", tags=["auth"])

# ── Email-verification token (separate secret + 1-hour TTL) ───────────────────
VERIFY_SECRET = os.getenv("JWT_SECRET", "verify_secret_fallback") + "_email_verify"
_verify_serializer = URLSafeTimedSerializer(VERIFY_SECRET)

# Separate salt keeps reset tokens distinct from email-verify tokens
_RESET_SALT = "password-reset"

# Password strength regex (compiled once at module level)
_UPPERCASE_RE = re.compile(r'[A-Z]')
_SPECIAL_RE = re.compile(r'[0-9!@#$%^&*]')


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


# ── Registration ──────────────────────────────────────────────────────────────

@router.post("/register", status_code=202)
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
        "created_at": datetime.datetime.now(datetime.timezone.utc),
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
            "verify_url_debug": verify_url,   # handy for local dev without a real key
        }

    return {"message": f"Verification email sent to {request.email}. Please check your inbox."}


# ── Email verification ────────────────────────────────────────────────────────

@router.get("/verify-email")
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


@router.post("/resend-verification")
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


# ── Forgot / Reset Password ──────────────────────────────────────────────────

@router.post("/forgot-password")
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


@router.post("/reset-password")
async def reset_password(request: ResetPasswordRequest):
    """Validates the signed token and updates the user's password."""
    try:
        email = _verify_serializer.loads(request.token, salt=_RESET_SALT, max_age=3600)
    except SignatureExpired:
        raise HTTPException(status_code=400, detail="Reset link has expired. Please request a new one.")
    except BadSignature:
        raise HTTPException(status_code=400, detail="Invalid or tampered reset token.")

    # Basic server-side password strength check
    pwd = request.new_password
    if len(pwd) < 8 or not _UPPERCASE_RE.search(pwd) or not _SPECIAL_RE.search(pwd):
        raise HTTPException(
            status_code=422,
            detail="Password must be at least 8 characters and include an uppercase letter and a number or special character.",
        )

    db = get_db()
    user = await db.users.find_one({"email": email})
    if not user or user.get("auth_provider") != "local":
        raise HTTPException(status_code=404, detail="User not found.")

    new_hash = get_password_hash(pwd)
    await db.users.update_one({"email": email}, {"$set": {"hashed_password": new_hash}})
    return {"message": "Password updated successfully. You can now sign in with your new password."}


# ── Login ─────────────────────────────────────────────────────────────────────

@router.post("/login")
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


@router.post("/google")
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
                "created_at": datetime.datetime.now(datetime.timezone.utc),
            }
            await db.users.insert_one(user)

        access_token = create_access_token({"sub": email, "name": name})
        return {"access_token": access_token, "user": {"email": email, "name": name, "picture": picture}}
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))
