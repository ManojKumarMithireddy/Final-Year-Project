"""
pytest configuration — ensure backend package is importable
regardless of the directory pytest is invoked from.
"""
import sys
import os

# Insert the backend root so that `import auth`, `import services.grover`, etc. work
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
