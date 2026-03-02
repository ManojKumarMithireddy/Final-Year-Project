"""
rate_limit.py
─────────────
Shared SlowAPI limiter instance.
Import `limiter` into main.py to attach it to app.state,
and into any router that needs @limiter.limit decorators.
"""
from slowapi import Limiter
from slowapi.util import get_remote_address

limiter = Limiter(key_func=get_remote_address)
