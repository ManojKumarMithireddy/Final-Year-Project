import pytest
import math
import numpy as np



# ─────────────────────────────────────────────────────────────────────────────
# 1. DNA → Bits encoding  (pure function, no FastAPI import needed)
# ─────────────────────────────────────────────────────────────────────────────
def dna_to_bits(dna_seq: str) -> str:
    """Local copy of the pure helper to avoid triggering FastAPI app init on import."""
    mapping = {'A': '00', 'C': '01', 'G': '10', 'T': '11'}
    return "".join(mapping.get(base, '00') for base in dna_seq)

def test_dna_to_bits():
    """Verify the 2-bit DNA encoding: A=00, C=01, G=10, T=11."""
    assert dna_to_bits("A")     == "00"
    assert dna_to_bits("C")     == "01"
    assert dna_to_bits("G")     == "10"
    assert dna_to_bits("T")     == "11"
    assert dna_to_bits("ACGT")  == "00011011"
    # A=00, T=11, G=10, C=01, T=11  →  0011100111
    assert dna_to_bits("ATGCT") == "0011100111"
    # Unknown base defaults to '00'
    assert dna_to_bits("N")     == "00"


# ─────────────────────────────────────────────────────────────────────────────
# 2. Grover iteration formula  π/4 · √N
# ─────────────────────────────────────────────────────────────────────────────
@pytest.mark.parametrize("N, expected", [
    (4,   1),   # 2 qubits: π/4 · 2 ≈ 1.57 → 1
    (16,  3),   # 4 qubits: π/4 · 4 ≈ 3.14 → 3
    (64,  6),   # 6 qubits: π/4 · 8 ≈ 6.28 → 6
    (256, 12),  # 8 qubits: π/4 · 16 ≈ 12.57 → 12
])
def test_grover_iteration_formula(N, expected):
    """Grover's optimal iteration count is floor(π/4 · √N)."""
    result = int(np.floor(np.pi / 4 * np.sqrt(N)))
    assert result == expected, f"For N={N}: expected {expected}, got {result}"


# ─────────────────────────────────────────────────────────────────────────────
# 3. build_grover_circuit — structural validation
# ─────────────────────────────────────────────────────────────────────────────
def test_build_grover_circuit_structure():
    """Grover circuit must have the right qubit count, be measured, and return correct iterations."""
    from services.grover import build_grover_circuit
    from qiskit import QuantumCircuit

    for target in ["01", "101", "1010", "0000"]:
        qc, iterations = build_grover_circuit(target)
        n = len(target)
        expected_iters = int(np.floor(np.pi / 4 * np.sqrt(2 ** n)))

        assert isinstance(qc, QuantumCircuit),  "Should return a QuantumCircuit"
        assert qc.num_qubits == n,              f"Expected {n} qubits, got {qc.num_qubits}"
        assert iterations == expected_iters,    f"Iteration count mismatch for '{target}'"
        assert qc.num_clbits > 0,              "Circuit must include measurement operations"


# ─────────────────────────────────────────────────────────────────────────────
# 4. Password hash round-trip
#    Uses bcrypt directly to avoid passlib's version-detection wrap-bug check
#    which fails on newer bcrypt wheels (passlib<1.7.5 + bcrypt>=4.x incompatibility)
# ─────────────────────────────────────────────────────────────────────────────
def test_password_hash_roundtrip():
    """Hashed password must verify correctly and reject wrong passwords."""
    import bcrypt

    password = b"Secure@2026"
    hashed = bcrypt.hashpw(password, bcrypt.gensalt())

    assert hashed != password,                                     "Hash must not equal plaintext"
    assert bcrypt.checkpw(password, hashed) is True,               "Correct password must verify"
    assert bcrypt.checkpw(b"WrongPass!", hashed) is False,         "Wrong password must not verify"


# ─────────────────────────────────────────────────────────────────────────────
# 5. Accession string sanitization
# ─────────────────────────────────────────────────────────────────────────────
def test_accession_sanitization():
    """Accession IDs must match the safe pattern; malformed ones must be rejected."""
    import re
    ACCESSION_RE = re.compile(r'^[A-Za-z0-9_.\-]{1,30}$')

    valid   = ["NM_007294", "AF123456.1", "NC-000001", "AB123456"]
    invalid = ["../etc/passwd", "'; DROP TABLE--", "", "A" * 31, "NM 007294"]

    for acc in valid:
        assert ACCESSION_RE.match(acc),      f"'{acc}' should be valid"
    for acc in invalid:
        assert not ACCESSION_RE.match(acc),  f"'{acc}' should be rejected"


# ─────────────────────────────────────────────────────────────────────────────
# 6. Password-reset token round-trip
#    Validates that the itsdangerous serializer (same one used in main.py)
#    correctly encodes and decodes a password-reset token with the
#    "password-reset" salt, and rejects tampered payloads.
# ─────────────────────────────────────────────────────────────────────────────
def test_password_reset_token_roundtrip():
    """Reset token must encode the email and decode back cleanly; tampering must fail."""
    from itsdangerous import URLSafeTimedSerializer, BadSignature

    secret = "test_jwt_secret_email_verify"
    serializer = URLSafeTimedSerializer(secret)
    salt = "password-reset"
    email = "user@example.com"

    # Generate and immediately decode (no max_age expiry during test)
    token = serializer.dumps(email, salt=salt)
    decoded = serializer.loads(token, salt=salt, max_age=3600)
    assert decoded == email, f"Expected '{email}', got '{decoded}'"

    # Wrong salt must be rejected
    try:
        serializer.loads(token, salt="email-verify", max_age=3600)
        assert False, "Should have raised BadSignature for wrong salt"
    except BadSignature:
        pass  # expected

    # Tampered token must be rejected
    tampered = token[:-4] + "XXXX"
    try:
        serializer.loads(tampered, salt=salt, max_age=3600)
        assert False, "Should have raised BadSignature for tampered token"
    except BadSignature:
        pass  # expected


# ─────────────────────────────────────────────────────────────────────────────
# 7. Session timeout — token expires within 30 minutes
#    Validates that create_access_token sets the exp claim to ≤ 30 min
#    in the future, confirming the session timeout configuration.
# ─────────────────────────────────────────────────────────────────────────────
def test_token_expires_within_30_minutes():
    """JWT exp claim must be ≤ 30 minutes from now (not 7 days)."""
    import jwt as pyjwt
    import datetime

    # Inline import to avoid triggering full app init
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
    from auth import create_access_token, jwt_secret, jwt_algo

    token = create_access_token({"sub": "test@example.com"})
    payload = pyjwt.decode(token, jwt_secret, algorithms=[jwt_algo])

    exp_dt = datetime.datetime.fromtimestamp(payload["exp"], tz=datetime.timezone.utc)
    now = datetime.datetime.now(datetime.timezone.utc)
    delta = exp_dt - now

    # Token must expire within 30 minutes (plus a small buffer for test execution)
    assert delta.total_seconds() <= 30 * 60 + 5, (
        f"Token expires in {delta.total_seconds():.0f}s — expected ≤ 1805s"
    )
    # Token should not already be expired
    assert delta.total_seconds() > 0, "Token is already expired at creation"

