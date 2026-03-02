"""
crypto.py
─────────
Fernet-based symmetric encryption helpers for sensitive values stored in MongoDB
(e.g. IBM Cloud API keys).

Usage
─────
Set FERNET_KEY in your .env file to a URL-safe base64-encoded 32-byte key.
Generate one with:

    python -c "from cryptography.fernet import Fernet; print(Fernet.generate_key().decode())"

If FERNET_KEY is not set, encrypt/decrypt are no-ops with a startup warning
so the app still runs in development without a key configured.
"""

import os
import logging

logger = logging.getLogger(__name__)

_fernet = None

def _get_fernet():
    """Lazily initialise the Fernet instance from the FERNET_KEY env var."""
    global _fernet
    if _fernet is not None:
        return _fernet

    key = os.getenv("FERNET_KEY", "")
    if not key:
        logger.warning(
            "FERNET_KEY is not set — IBM credentials will be stored unencrypted. "
            "Generate a key with: python -c \"from cryptography.fernet import Fernet; "
            "print(Fernet.generate_key().decode())\""
        )
        return None

    try:
        from cryptography.fernet import Fernet
        _fernet = Fernet(key.encode())
        return _fernet
    except Exception as exc:
        logger.error("Invalid FERNET_KEY: %s — credentials stored unencrypted.", exc)
        return None


def encrypt_value(plaintext: str) -> str:
    """
    Encrypt a string value. Returns the ciphertext as a UTF-8 string.
    Returns the plaintext unchanged if FERNET_KEY is not configured.
    """
    f = _get_fernet()
    if f is None:
        return plaintext
    return f.encrypt(plaintext.encode()).decode()


def decrypt_value(ciphertext: str) -> str:
    """
    Decrypt a Fernet-encrypted string. Returns the plaintext.
    Returns the ciphertext unchanged if FERNET_KEY is not configured
    (handles the case where the value was stored before encryption was enabled).
    """
    f = _get_fernet()
    if f is None:
        return ciphertext
    try:
        return f.decrypt(ciphertext.encode()).decode()
    except Exception:
        # Value may have been stored before encryption was enabled — return as-is
        logger.warning("decrypt_value: could not decrypt — returning raw value.")
        return ciphertext
