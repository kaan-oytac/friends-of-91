"""Rigorous primality of every prime cited in a certificate.
Usage: python prove_primes.py cert.json [cert2.json ...]
verify.py uses a 12-base Miller-Rabin test, which is deterministic below 3.18e23 and a strong
probable-prime test above.  This script re-checks every integer that a certificate presents as a
prime (path primes, denominator factors, forbidden primes) with FLINT's fmpz_is_prime, which is a
proven primality test (BPSW followed by a Pocklington/APR-CL style proof), via python-flint.
Requires: pip install python-flint"""
import sys, json
from flint import fmpz
import flint

MR_LIMIT = 318665857834031151167461      # first 12 primes as bases are deterministic below this
for path in sys.argv[1:]:
    cert = json.load(open(path))
    cands = set()
    for br in cert['branches']:
        for nd in br['nodes']:
            for key in ('den_factors', 'den_partial'):
                cands.update(int(q) for q, _ in nd.get(key, []))
            cands.update(int(q) for q in nd.get('forbidden', []))
            cands.update(int(q) for q, _ in nd['path'])
    big = sorted(q for q in cands if q >= MR_LIMIT)
    bad = [q for q in cands if not fmpz(q).is_prime()]
    print(f"{path}: {len(cands)} distinct cited primes, {len(big)} above the deterministic Miller-Rabin range "
          f"(largest {len(str(max(cands)))} digits); python-flint {flint.__version__}, FLINT {flint.__FLINT_VERSION__}")
    if bad:
        print("  NOT PRIME:", bad); sys.exit(1)
    print("  all proved prime")
