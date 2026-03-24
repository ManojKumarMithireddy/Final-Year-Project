# BioQuantum — Professional Fullstack Code Review

> **Reviewer role:** Senior Fullstack Engineer + BioQuantum Specialist  
> **Review date:** March 4, 2026 (updated — includes BRCA1 dual-patient PoC additions)  
> **Stack:** FastAPI · Motor/MongoDB · Qiskit/Aer · React 19 · Vite · Tailwind CSS

---

## Table of Contents

1. [Executive Summary](#1-executive-summary)
2. [Architecture Overview](#2-architecture-overview)
3. [Security](#3-security)
4. [Backend](#4-backend)
5. [Frontend](#5-frontend)
6. [Quantum Layer](#6-quantum-layer)
7. [BioQuantum Algorithm & Bio Layer](#7-bioquantum-algorithm--bio-layer)
8. [Testing](#8-testing)
9. [DevOps & Deployability](#9-devops--deployability)
10. [Prioritised Action List](#10-prioritised-action-list)
11. [Ratings Summary](#11-ratings-summary)

---

## 1. Executive Summary

BioQuantum is a well-structured full-stack educational project that demonstrates a hybrid classical/quantum approach to genomic sequence search, pairing real NCBI sequence fetching, local Aer simulation, and optional IBM Cloud QPU execution with a polished React SPA.

The codebase shows **good engineering instincts**: proper router/service separation, async MongoDB via Motor, itsdangerous-signed email verification, a custom session-timeout hook, and meaningful inline documentation. The visual layer (Framer Motion animations, Recharts, dark theme) is consistently well-crafted.

The most recent additions — the **BRCA1 c.5266dupC dual-patient constrained Grover PoC** — introduced sophisticated biological framing (carrier vs. healthy control), a numpy statevector simulation to avoid HuggingFace memory crashes, and a marker-specificity diagnostic. These are genuinely impressive for an educational project.

That said, several **security issues** remain, and the BioQuantum layer has a critical memory/correctness flaw that must be addressed. These are documented in the sections below.

---

## 2. Architecture Overview

```
bioquantum_hybrid_search/
├── backend/                  # FastAPI app
│   ├── main.py               # App factory + lifespan + CORS
│   ├── auth.py               # JWT utilities + Google OAuth
│   ├── db.py                 # Motor client singleton
│   ├── email_service.py      # Brevo transactional email
│   ├── models/schemas.py     # Pydantic request models
│   ├── routers/              # auth · search · quantum · credentials · user
│   ├── services/             # grover.py · ncbi.py
│   └── tests/test_core.py    # Pytest unit tests
└── frontend/                 # Vite + React SPA
    ├── src/
    │   ├── App.jsx            # Router, navigation, session timeout
    │   ├── pages/             # Login · Dashboard · History · VerifyEmail · ResetPassword
    │   ├── components/        # GroverPOC · HybridSearchPanel · AmplitudeChart · …
    │   └── lib/api.js         # Axios instance with auth interceptor
    └── package.json
```

**What works well:**

- Clean separation: router → service → DB layer
- Fat-free `main.py` (pure factory/lifecycle)
- `models/schemas.py` centralising all Pydantic models
- `lib/api.js` as a single Axios instance — all callers benefit from the interceptor automatically
- `ErrorBoundary` component preventing white-screen crashes

---

## 3. Security

These are the most important findings and should be addressed first.

### 🔴 Critical

#### 3.1 Hardcoded JWT secret (`auth.py`)

```python
# auth.py
jwt_secret = os.getenv("JWT_SECRET", "super_secret_jwt_key_here_change_in_prod")
```

The fallback string is committed to version control. Any developer who doesn't set the env var ships with a publicly known signing key — all JWTs can be forged. **Remove the fallback entirely** and raise a startup error if `JWT_SECRET` is missing.

```python
jwt_secret = os.environ["JWT_SECRET"]   # fail fast, no default
```

#### 3.2 IBM API key returned in plaintext to the frontend (`routers/credentials.py`)

```python
# GET /api/credentials — returns raw key
return creds   # includes ibm_api_key, ibm_crn
```

The IBM Cloud API key is fetched only to pre-fill the frontend form so the user knows a key is saved. It should never leave the server. Return a redacted sentinel instead:

```python
return {
    "ibm_api_key": "••••••••" if creds.get("ibm_api_key") else "",
    "ibm_crn": creds.get("ibm_crn", ""),
}
```

#### 3.3 IBM credentials encrypted at rest — `FERNET_KEY` silent failure ✅ / ⚠️

Fernet symmetric encryption (`crypto.py`) was added for IBM credentials at rest — a good improvement. However, there is a **silent no-op** path:

```python
# crypto.py
if not FERNET_KEY:
    logger.warning("FERNET_KEY not set — credentials stored unencrypted")
    return plaintext   # silently falls through
```

If `FERNET_KEY` is misconfigured in production (typo, missing env var), credentials are stored in plaintext with only a log warning. This should be a hard failure:

```python
if not FERNET_KEY:
    raise RuntimeError("FERNET_KEY environment variable is required for credential storage.")
```

#### 3.4 Token in `localStorage` — XSS exposure (`App.jsx`, `lib/api.js`)

JWTs in `localStorage` are readable by any JavaScript running on the page, including injected scripts. The preferred pattern for web apps is `httpOnly` cookies, which are inaccessible to JS. This requires a small backend change (`set-cookie` on login) and removing the `Authorization` header injection from `api.js`, but significantly raises the XSS bar.

---

### 🟠 High

#### 3.5 Rate limiting added to auth — not on quantum endpoints ✅ / ⚠️

`slowapi` rate limiting was correctly added to `/api/auth/login` and `/api/auth/register` — good. However, the **quantum endpoints are unprotected**:

- `POST /quantum-poc/bio-local` triggers an NCBI HTTP fetch + heavy numpy matrix operations
- `POST /quantum-poc/bio-ibm-submit` stores credentials and queues IBM QPU jobs

A single user can hammer these endpoints, exhausting NCBI API quotas or queuing many IBM jobs. Add `@limiter.limit("3/minute")` to both quantum routes.

#### 3.6 No email format validation in `StandardAuthRequest`

The `email` field is a plain `str`. Pydantic has a built-in `EmailStr` type (already imported as a dependency via `email-validator`). A malformed email would be registered, then the verification link would go nowhere.

```python
from pydantic import EmailStr
class StandardAuthRequest(BaseModel):
    email: EmailStr
    ...
```

#### 3.7 Password strength enforced in the router, not the schema

Password complexity checks happen inside `routers/auth.py` at runtime, but the `StandardAuthRequest` schema accepts any string. Move enforcement to the Pydantic model with a `@field_validator` so it is consistently applied everywhere the schema is used and appears in the OpenAPI docs.

---

### 🟡 Medium

#### 3.8 `get_history` exposes the user's `email` field in every history record

`routers/user.py` returns history documents with `{"_id": 0}` — email is still included. If the history API is ever shared or logged, it leaks PII. Add `"email": 0` to the projection.

#### 3.9 No HTTPS enforcement / security headers

No `Strict-Transport-Security`, `X-Content-Type-Options`, or `Content-Security-Policy` headers. For production, add `starlette-trustedhost` and a security-header middleware (or sit behind a reverse proxy that adds them).

---

## 4. Backend

### 4.1 Module-level MongoDB client (`db.py`)

```python
# db.py — executed at import time
client = AsyncIOMotorClient(MONGODB_URI)
db = client.bioquantum_hybrid_search
```

This works but is fragile:

- Connection errors at import time crash the whole process with an opaque traceback.
- It makes the DB impossible to mock in unit tests without patching at the module level.

**Recommended approach:** lazy initialisation with a proper FastAPI dependency or a `lifespan`-scoped connection stored on `app.state`.

```python
# In lifespan:
app.state.db = AsyncIOMotorClient(MONGODB_URI).bioquantum_hybrid_search
# In routers:
async def some_route(request: Request):
    db = request.app.state.db
```

### 4.2 `get_db()` called inline in every route handler

Every route calls `get_db()` directly — there is no DI (Dependency Injection) pattern. This makes mocking for integration tests require monkeypatching the module global. Wrapping `get_db` as a FastAPI `Depends` fixture would make test injection clean.

### 4.3 No timeout on NCBI Entrez fetches (`services/ncbi.py`)

`Entrez.efetch` has no timeout. A slow or hung NCBI connection will tie up a worker indefinitely. Wrap with `asyncio.wait_for` or set `Entrez.timeout`.

### 4.4 NCBI placeholder email violates NCBI Terms of Service

```python
Entrez.email = os.getenv("NCBI_EMAIL", "bioquantum.capstone@example.com")
```

NCBI requires a _real_ contact email in case of abuse. The fallback domain `example.com` is explicitly non-operational. Require `NCBI_EMAIL` from the environment with no default, or document prominently in `.env.example` that it must be set.

### 4.5 `requirements.txt` is a full pip freeze

The file contains 55+ packages including transitive/platform dependencies (`colorama`, `rustworkx`, etc.). This makes it fragile on different OS/architecture targets and hides the _direct_ dependencies. Split into:

- `requirements.txt` — direct deps only (fastapi, motor, qiskit, etc.)
- `requirements-lock.txt` or use `pip-tools` / `poetry`

### 4.6 `qiskit-ibm-runtime` not in `requirements.txt`

The IBM Cloud path (`routers/quantum.py`) gracefully handles a missing import, but the package is not listed in requirements — a developer who wants to use IBM hardware has no guidance. Add it as an optional extra with a comment.

### 4.7 Duplicate `dna_to_bits` function

The same function exists in `services/ncbi.py` and is duplicated inline in `tests/test_core.py` (with a comment explaining why). It also exists as `dnaToBits` in the frontend. The backend copy should live in a shared `utils.py` and be imported by both `ncbi.py` and tested directly.

---

## 5. Frontend

### 5.1 Dual Axios usage in `App.jsx`

`App.jsx` imports raw `axios` to register a global 401 interceptor, while all actual API calls go through `lib/api.js`. This is the correct intention (the 401 interceptor needs to catch _all_ axios instances for the window redirect), but it is confusing. Document the reason with a comment, or move the global interceptor into `api.js` and export it so `App.jsx` only imports it without an extra raw `axios` import.

### 5.2 History `key` prop falls back to non-unique `timestamp`

```jsx
// History.jsx
const id = item._id || item.timestamp;
```

`_id` is excluded in the MongoDB projection (`{"_id": 0}`) so `item._id` is always `undefined`. All keys therefore fall back to `item.timestamp`. If two simulations ran within the same second (possible locally), React will warn about duplicate keys and correctly render only one row. The backend should include `_id` in the history response (serialised as a string), or use a generated UUID field.

### 5.3 `GroverPOC.jsx` is 819 lines

This single component handles dataset generation, qubit/noise controls, IBM credential management, session submission, step-by-step circuit navigation, OTP encryption display, and amplitude charting. It should be decomposed:

- `GroverControls.jsx` — qubit, target, noise, backend selectors
- `IBMCredentialsPanel.jsx` — already a logical sub-unit
- `PIRPipeline.jsx` — the three-step encrypt/run/decrypt visualisation
- `CircuitNavigator.jsx` — step tabs + diagram + AmplitudeChart (could be shared with History)

### 5.4 No frontend error handling for circuit-info failures in History

If `GET /api/search/circuit-info` fails (network error, backend down), `loadingCircuit` is cleared but `circuit` remains `null` with no error message rendered — the panel just shows the "No target bits" fallback even when there are target bits. Add an error state to show a retry or failure message.

### 5.5 No input debounce / cancellation on dataset fetch

`HybridSearchPanel.jsx` fires `api.get('/sequence/${accession}')` on button click, but rapid clicks can queue multiple in-flight requests. Disable the button while fetching (already done ✅), but there is no AbortController to cancel a stale request if the component unmounts.

### 5.6 `vite.config.js` has no API proxy

During development, the frontend hits `http://localhost:8000` directly. This means the browser must allow cross-origin requests, and the `VITE_API_BASE` env var must be set. A Vite `server.proxy` entry would eliminate CORS in development entirely and make it work out-of-the-box without an env file:

```js
server: {
  proxy: {
    '/api': 'http://localhost:8000',
  },
}
```

---

## 6. Quantum Layer

### 6.1 Conceptual mismatch (acknowledged, but worth surfacing)

The application combines two parallel tracks that operate at _completely different scales_:

- **Classical search**: runs on real NCBI sequences (up to 100,000 bases) — thousands of search windows.
- **Grover POC**: simulates circuits over at most 2^6 = 64 states (6 qubits).

The `ComplexityChart` plots $O(\sqrt{N})$ vs. $O(N)$ using the _classical_ window count as $N$, but the Grover circuit never actually runs on that many states. The inline warning ("These are analytically derived values") is good, but the chart's x-axis (`~100k windows`) might mislead users into thinking the quantum circuit was run at that scale. Consider adding a second annotation line on the chart at the actual circuit scale.

### 6.2 `build_grover_step_circuits` is called on every POC run and every history open

`build_grover_step_circuits` constructs 4 statevector Qiskit circuits and runs full statevector simulations for each one. For 4–6 qubits this is fast, but for 6 qubits (64 states) it adds meaningful latency. The results are deterministic for a given bitstring, so they could be cached in memory (`@functools.lru_cache`) keyed on the target bitstring:

```python
from functools import lru_cache

@lru_cache(maxsize=64)
def build_grover_step_circuits(target_bitstring: str) -> list:
    ...
```

### 6.3 OTP-XOR "PIR" is a pedagogical approximation (already documented, flagging for completeness)

The code correctly notes in multiple docstrings that the current implementation is _not_ cryptographically rigorous PIR. A genuine PIR scheme requires an oblivious oracle. The documentation is accurate, but the UI labels it as "Private Information Retrieval" in ways that might overstate guarantees to non-expert users. Consider softening to "PIR-inspired Demo".

---

## 7. BioQuantum Algorithm & Bio Layer

This section covers flaws and improvements specific to the BRCA1 constrained Grover PoC introduced in the latest development sprint.

### 🔴 7.1 Out-of-Memory risk for `n_codons=3` (18-qubit diffusion matrix)

The diffusion operator is constructed as a dense 2D numpy array:

```python
# bio_grover.py line ~300
D = 2.0 * np.outer(sv0, sv0.conj()) - np.eye(n_max, dtype=complex)
```

`n_max = 1 << n_qubits`. The memory footprint of `D`:

| `n_codons` | `n_qubits` | `D` matrix size | RAM required |
|------------|------------|-----------------|--------------|
| 1          | 6          | 64 × 64         | ~64 KB ✅   |
| 2          | 12         | 4 096 × 4 096   | ~256 MB ⚠️  |
| 3          | 18         | 262 144 × 262 144| ~1 TB ❌   |

The comment says "Memory: O(2^n)" — this refers only to the statevector, NOT to `D`. The **diffusion matrix is O(4^n)**, which is catastrophic for n=18. HuggingFace cpu-basic has ~2 GB RAM; `n_codons=3` will crash with OOM.

**Fix:** Replace the dense `D @ state` product with the mathematically equivalent three-operation form that is O(2^n) in both time and memory:

```python
# D @ state = 2·⟨ψ₀|state⟩·|ψ₀⟩ − state
def _diffuse(state, sv0):
    overlap = np.dot(sv0.conj(), state)
    return 2.0 * overlap * sv0 - state
```

This avoids ever materialising the 2^n × 2^n matrix. Apply same fix in `build_bio_step_circuits`.

### 🔴 7.2 `build_bio_step_circuits` is redundantly called inside `run_bio_grover_local`

```python
# bio_grover.py line ~344 (inside run_bio_grover_local return dict)
"step_circuits": build_bio_step_circuits(patient_nodes, target_bits),
```

`build_bio_step_circuits` rebuilds `sv0`, `oracle_diag`, `D`, and runs all Grover iterations from scratch — effectively doubling all compute. The diffusion matrix `D` is built **twice** for every `/bio-local` call.

**Fix:** Accept `sv0`, `oracle_diag`, `D` as parameters in `build_bio_step_circuits`, or refactor into a shared helper that `run_bio_grover_local` also calls with pre-built operators.

### 🟠 7.3 Step diagram mismatch — standard H diffuser shown for constrained diffusion

In `build_bio_step_circuits`, Step 3 displays the standard Hadamard diffuser circuit:

```python
# bio_grover.py lines ~429-437
qc3.h(range(n))
qc3.x(range(n))
qc3.h(n-1)
qc3.mcx(...)
...
```

But the actual algorithm uses the **constrained diffusion** `D_c = 2|ψ₀⟩⟨ψ₀| − I`, which cannot be expressed as a simple gate sequence. The displayed circuit is WRONG — it shows unconstrained Grover diffusion. This is educationally misleading.

**Fix:** Replace `qc3` with a named black-box gate (matching how initialisation is shown) and add a label explaining it is the constrained inversion operator:

```python
qc3 = QuantumCircuit(n)
diff = QuantumCircuit(n, name="C-Diff\n2|ψ₀⟩⟨ψ₀|−I")
qc3.append(diff.to_gate(), range(n))
```

### 🟠 7.4 IBM QPU path uses unconstrained Grover (conceptual mismatch with UI)

`_build_grover_circuit_standard` initialises with `qc.h(range(n))` — a standard uniform superposition over all 2^n states, NOT the constrained DNA-node subspace. The UI labels this "Constrained Grover Search" without qualification.

The oracle is correct (marks the target), but the initialisation and diffuser are unconstrained. The quantum circuit is therefore NOT the same algorithm as the local simulation.

**Fix:** Add a clear disclaimer in the IBM UI panel: "IBM QPU uses standard Grover initialisation (H⊗ⁿ) over all 2^n states; the constrained DNA-node subspace is only used in the local simulator."

### 🟠 7.5 `MARKER_ACCESSION` constant is dead code

```python
# bio_grover.py line 56
MARKER_ACCESSION = "NM_007294.4"  # ← never used anywhere
```

`fetch_disease_marker()` was removed when the marker was changed to be derived from the mutant sequence in-memory. The constant is a stale leftover that adds confusion.

### 🟡 7.6 `requirements.txt` has duplicate `slowapi` entry

```
# lines 49 and 51
slowapi==0.1.9
...
slowapi==0.1.9
```

Harmless but sloppy. Remove one entry.

### 🟡 7.7 `qiskit-ibm-runtime` is unpinned

```
qiskit-ibm-runtime>=0.20.0
```

All other dependencies use exact pinned versions. An unpinned `>=` constraint will silently upgrade to a future breaking version. Pin it: `qiskit-ibm-runtime==0.20.0` (or latest tested version).

### 🟡 7.8 No NCBI rate limiting or retry backoff

NCBI Entrez allows maximum 3 requests/second unauthenticated (10/sec with API key). There is no `time.sleep` between calls, no retry/backoff logic, and no `NCBI_API_KEY` support. Rapid test runs or multiple concurrent users could trigger `429` responses from NCBI, which currently propagate as opaque 502 errors.

**Fix:** Add `Entrez.api_key = os.getenv("NCBI_API_KEY")` and a simple exponential backoff wrapper around `_ncbi_fetch`.

### 🟡 7.9 `detection_result` uses most-frequent bitstring, not probability threshold

```python
measured = max(counts, key=counts.get)
detection_result = "FOUND" if measured == target_bits else "NOT_FOUND"
```

If the healthy patient happens to produce target_bits as the most frequent shot (due to randomness or marker collision), the result is `FOUND` with, say, 4% confidence. Currently ANY `measured == target_bits` counts as detection regardless of confidence level.

**Fix:** Add a minimum confidence threshold (e.g., `conf >= 0.15`) before classifying as `FOUND`, or surface both `measured == target_bits` and the confidence value more prominently in the UI so users understand low-confidence results.

### 🟢 7.10 Biological framing is now correct (positive note)

The `apply_brca1_mutation` → `get_marker_seq` pipeline correctly:
- Simulates c.5266dupC by inserting a C at 0-based index 5266
- Extracts the disease marker from the mutant (not from the reference — the original bug)
- Uses the same marker for both carrier and healthy patient comparisons (correct fixed reference point)
- Checks `marker_in_reference` to surface the marker-specificity problem

This is a genuine improvement and the biological logic is sound.

---

## 8. Testing

### 8.1 Coverage

| Area                          | Status    |
| ----------------------------- | --------- |
| DNA encoding                  | ✅ Tested |
| Grover iteration formula      | ✅ Tested |
| Circuit structure             | ✅ Tested |
| Password hashing              | ✅ Tested |
| Accession sanitisation        | ✅ Tested |
| Token round-trip              | ✅ Tested |
| Session expiry                | ✅ Tested |
| BRCA1 mutation helper         | ❌ None   |
| Constrained Grover (numpy)    | ❌ None   |
| Marker specificity check      | ❌ None   |
| API route integration         | ❌ None   |
| Frontend                      | ❌ None   |
| Email service                 | ❌ None   |
| NCBI fetch                    | ❌ None   |

### 8.2 Brittle `sys.path` hack in test file

```python
# test_core.py
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
```

This is repeated twice in the file and is fragile when running pytest from different working directories. Use a `conftest.py` for path setup, or better, make the backend a proper package with a `pyproject.toml` / `setup.py` so pytest resolves imports naturally.

### 8.3 No API integration tests

There are no tests that spin up a `TestClient` and exercise the HTTP layer. FastAPI's built-in `TestClient` (wrapping Starlette's) makes this straightforward. At minimum, tests for the auth flow (register → verify → login → protected route) and the quantum POC endpoint would catch regressions.

---

## 9. DevOps & Deployability

| Concern                     | Status            | Notes                                                                  |
| --------------------------- | ----------------- | ---------------------------------------------------------------------- |
| `.env.example`              | ❌ Missing        | No reference file listing required env vars                            |
| Dockerfile                  | ✅ Present        | HuggingFace Spaces Docker deployment                                   |
| Docker Compose              | ❌ Missing        | No local orchestration of API + MongoDB together                       |
| CI/CD pipeline              | ✅ Present        | GitHub Actions: `ci.yml` (pytest + build), `deploy-backend.yml` (HF), `deploy-frontend.yml` (Vercel) |
| Frontend build optimisation | ⚠️ Default        | No `rollupOptions` code-splitting in `vite.config.js`                  |
| Production CORS             | ⚠️ Partial        | `ALLOWED_ORIGINS` env var present but undocumented                     |
| Logging                     | ✅ Fixed          | `logging.getLogger(__name__)` used in `bio_grover.py` and services     |
| Backend hosted              | ✅ HuggingFace    | `manojkumarmithireddy-bioquantum-api.hf.space`                        |
| Frontend hosted             | ✅ Vercel         | `bioquantum-frontend.vercel.app`                                       |

### Remaining recommended additions:

**`.env.example`** (backend)

```
MONGODB_URI=mongodb://localhost:27017
JWT_SECRET=                    # required — generate with: openssl rand -hex 32
FERNET_KEY=                    # required — generate with: python -c "from cryptography.fernet import Fernet; print(Fernet.generate_key().decode())"
GOOGLE_CLIENT_ID=              # from Google Cloud Console
BREVO_API_KEY=                 # from Brevo dashboard
BREVO_SENDER_EMAIL=            # verified sender address
NCBI_EMAIL=                    # your real email (NCBI ToS)
NCBI_API_KEY=                  # optional — raises rate limit to 10 req/sec
APP_BASE_URL=http://localhost:5173
ALLOWED_ORIGINS=http://localhost:5173
```

**`.env.example`** (frontend)

```
VITE_API_BASE=http://localhost:8000/api
VITE_GOOGLE_CLIENT_ID=         # same as backend GOOGLE_CLIENT_ID
```

---

## 10. Prioritised Action List

### P0 — Critical / Fix immediately

- [x] **[BioQuantum]** Replace dense `D @ state` with `2·⟨ψ₀|state⟩·sv0 − state` in both `run_bio_grover_local` and `build_bio_step_circuits` — eliminates OOM crash for `n_codons=3`
- [ ] **[BioQuantum]** Remove `build_bio_step_circuits` call from inside `run_bio_grover_local` return dict — currently rebuilds all operators twice per request
- [x] Remove hardcoded `JWT_SECRET` default; raise `KeyError` on missing env var
- [x] Redact IBM API key in `GET /api/credentials` response (return `••••••••` if set)
- [ ] Harden `FERNET_KEY` missing: raise `RuntimeError` instead of silent no-op
- [x] Add `.env.example` for both backend and frontend (including `FERNET_KEY`)

### P1 — High / Fix before adding features

- [ ] **[BioQuantum]** Fix Step 3 circuit diagram in `build_bio_step_circuits`: replace standard H diffuser with named black-box `C-Diff 2|ψ₀⟩⟨ψ₀|−I` gate
- [ ] **[BioQuantum]** Add IBM QPU disclaimer in UI: "IBM path uses unconstrained H⊗ⁿ initialisation"
- [ ] **[BioQuantum]** Remove dead `MARKER_ACCESSION` constant from `bio_grover.py`
- [ ] **[BioQuantum]** Remove duplicate `slowapi==0.1.9` from `requirements.txt`
- [ ] **[BioQuantum]** Pin `qiskit-ibm-runtime` to exact version in `requirements.txt`
- [ ] Add `@limiter.limit("3/minute")` to quantum endpoints (`/bio-local`, `/bio-ibm-submit`)
- [x] Switch `email` field to `EmailStr` in `StandardAuthRequest`
- [x] Move password-strength validation into a Pydantic `@field_validator`
- [x] Add `"email": 0` projection in `get_history`

### P2 — Medium / Quality of life

- [ ] **[BioQuantum]** Add `NCBI_API_KEY` support and exponential backoff retry to `_ncbi_fetch`
- [ ] **[BioQuantum]** Add confidence threshold (e.g. ≥15%) before classifying as `FOUND`
- [ ] **[BioQuantum]** Add unit tests for `apply_brca1_mutation`, `get_marker_seq`, `run_bio_grover_local`
- [x] Fix History `key` prop — include `_id` (as string) in history API response
- [x] Add error state in History when `circuit-info` request fails
- [x] Add `conftest.py` and remove `sys.path.insert` from test file
- [ ] Add `TestClient`-based integration tests for auth + bio-local endpoints
- [x] Add Vite dev proxy to remove CORS friction in local dev
- [ ] Split `GroverPOC.jsx` into sub-components
- [x] Add `AbortController` to `HybridSearchPanel` fetch

### P3 — Low / Production hardening

- [ ] Docker Compose for local API + MongoDB orchestration
- [ ] `httpOnly` cookie strategy to remove JWT from `localStorage`
- [ ] Add security headers middleware (`X-Content-Type-Options`, `CSP`, `HSTS`)
- [ ] `useBackendStatus` health check: poll `/health` endpoint instead of `/docs`

---

## 11. Ratings Summary

| Dimension                                 | Score      | Notes                                                                    |
| ----------------------------------------- | ---------- | ------------------------------------------------------------------------ |
| **Architecture / Separation of concerns** | 8 / 10     | Clean router→service split; DB DI could be improved                      |
| **Security**                              | 5 / 10     | Fernet + rate limiting added; JWT secret default + credential redaction remain |
| **Backend code quality**                  | 7 / 10     | Readable, well-commented; duplicate operators in bio layer                |
| **Frontend code quality**                 | 7 / 10     | Polished dual-patient UI; step diagram mismatch; GroverPOC too large      |
| **Quantum implementation**                | 6 / 10     | Correct constrained Grover math; OOM risk at n=18; IBM/local inconsistency |
| **Biological framing**                    | 7 / 10     | Carrier vs. healthy control is now correct; marker specificity well handled |
| **Testing**                               | 4 / 10     | Unit tests exist; no bio-grover or integration tests                      |
| **DevOps / Deployability**                | 7 / 10     | CI/CD and Dockerfile in place; missing `.env.example`, no Docker Compose  |
| **Documentation**                         | 8 / 10     | Excellent inline docs; bio_grover.py docstrings are exemplary             |
| **Overall**                               | **6.5 / 10** | Impressive educational project; one P0 crash risk in quantum layer needs fixing |

---

_This review reflects the codebase as of March 4, 2026, including the BRCA1 dual-patient constrained Grover PoC additions._
