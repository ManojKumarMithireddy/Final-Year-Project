"""
BioQuantum Hybrid Quantum Search API
─────────────────────────────────────
Thin app factory: configures CORS, lifecycle hooks, and mounts routers.
All endpoint logic lives in the routers/ package.
"""

import os
import logging
import datetime
from contextlib import asynccontextmanager

from fastapi import FastAPI, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse
from slowapi import _rate_limit_exceeded_handler
from slowapi.errors import RateLimitExceeded

from db import close_db
from rate_limit import limiter

# ── Structured logging ────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(name)s  %(message)s",
    datefmt="%Y-%m-%dT%H:%M:%S",
)
logger = logging.getLogger(__name__)

# ── Routers ───────────────────────────────────────────────────────────────────
from routers.auth import router as auth_router
from routers.search import router as search_router
from routers.quantum import router as quantum_router
from routers.credentials import router as credentials_router
from routers.user import router as user_router


# ── NTP clock check (optional dependency) ─────────────────────────────────────
try:
    import ntplib
    NTPLIB_AVAILABLE = True
except ImportError:
    NTPLIB_AVAILABLE = False


# ── Lifespan (replaces deprecated @app.on_event) ─────────────────────────────

@asynccontextmanager
async def lifespan(app: FastAPI):
    """Startup and shutdown lifecycle hooks."""
    # ── Startup ───────────────────────────────────────────────────────────
    if NTPLIB_AVAILABLE:
        try:
            c = ntplib.NTPClient()
            response = c.request('pool.ntp.org', version=3)
            from datetime import timezone
            ntp_time = datetime.datetime.fromtimestamp(response.tx_time, tz=timezone.utc)
            local_time = datetime.datetime.now(tz=timezone.utc)
            delta = abs((ntp_time - local_time).total_seconds())
            if delta > 5:
                logger.warning("System clock is %.1fs off NTP. Google token auth may fail.", delta)
            else:
                logger.info("Clock sync OK: %.2fs from NTP (pool.ntp.org)", delta)
        except Exception as e:
            logger.warning("NTP check failed (non-critical): %s", e)
    else:
        logger.info("ntplib not installed — clock skew check skipped. Install with: pip install ntplib")

    yield  # ── App is running ──

    # ── Shutdown ──────────────────────────────────────────────────────────
    close_db()
    logger.info("MongoDB connection closed.")


# ── App factory ───────────────────────────────────────────────────────────────

app = FastAPI(title="BioQuantum Hybrid Quantum Search API", lifespan=lifespan)

# ── Rate limiting ─────────────────────────────────────────────────────────────
app.state.limiter = limiter
app.add_exception_handler(RateLimitExceeded, _rate_limit_exceeded_handler)

# ── CORS ──────────────────────────────────────────────────────────────────────
# Read from env var for production flexibility; fall back to dev defaults
_default_origins = "http://localhost:5173,http://127.0.0.1:5173,http://localhost:5174,http://127.0.0.1:5174"
_origins = os.getenv("ALLOWED_ORIGINS", _default_origins).split(",")

app.add_middleware(
    CORSMiddleware,
    allow_origins=[o.strip() for o in _origins],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# ── Include routers ───────────────────────────────────────────────────────────
app.include_router(auth_router)
app.include_router(search_router)
app.include_router(quantum_router)
app.include_router(credentials_router)
app.include_router(user_router)


if __name__ == "__main__":
    import uvicorn
    uvicorn.run("main:app", host="127.0.0.1", port=8000, reload=True)
