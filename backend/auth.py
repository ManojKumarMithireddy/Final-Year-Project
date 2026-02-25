import os
import jwt
import datetime
from fastapi import Request, HTTPException, Security, Depends
from fastapi.security import HTTPBearer, HTTPAuthorizationCredentials
from typing import Optional
from google.oauth2 import id_token
from google.auth.transport import requests
from passlib.context import CryptContext

pwd_context = CryptContext(schemes=["bcrypt"], deprecated="auto")

def verify_password(plain_password, hashed_password):
    return pwd_context.verify(plain_password, hashed_password)

def get_password_hash(password):
    return pwd_context.hash(password)

jwt_secret = os.getenv("JWT_SECRET", "super_secret_jwt_key_here_change_in_prod")
jwt_algo = os.getenv("JWT_ALGORITHM", "HS256")
google_client_id = os.getenv("GOOGLE_CLIENT_ID", "")

security = HTTPBearer()

def verify_google_token(token: str):
    try:
        idinfo = id_token.verify_oauth2_token(
            token, 
            requests.Request(), 
            google_client_id,
            clock_skew_in_seconds=30  # Allow 30s leniency for Windows clock drift vs Google NTP
        )
        return idinfo
    except Exception as e:
        print(f"GOOGLE TOKEN VERIFICATION ERROR: {e}")
        raise ValueError(f"Google Token Error: {e}")

def create_access_token(data: dict):
    to_encode = data.copy()
    expire = datetime.datetime.utcnow() + datetime.timedelta(days=7)
    to_encode.update({"exp": expire})
    encoded_jwt = jwt.encode(to_encode, jwt_secret, algorithm=jwt_algo)
    return encoded_jwt

def get_current_user(credentials: HTTPAuthorizationCredentials = Security(security)):
    token = credentials.credentials
    try:
        payload = jwt.decode(token, jwt_secret, algorithms=[jwt_algo])
        return payload
    except jwt.ExpiredSignatureError:
        raise HTTPException(status_code=401, detail="Token has expired")
    except jwt.InvalidTokenError:
        raise HTTPException(status_code=401, detail="Invalid token")

optional_security = HTTPBearer(auto_error=False)

def get_optional_user(credentials: Optional[HTTPAuthorizationCredentials] = Depends(optional_security)):
    if not credentials:
        return None
    try:
        payload = jwt.decode(credentials.credentials, jwt_secret, algorithms=[jwt_algo])
        return payload
    except:
        return None
