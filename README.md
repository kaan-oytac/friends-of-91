# Friends of 91: certificates and code

Supplementary material for the paper *The Friends of 91 Live Beyond Googol Street* by Kaan Martin Oytac.

A positive integer m is a *friend* of n if m != n and sigma(m)/m = sigma(n)/n. This repository contains
the programs and machine-checkable certificates behind the paper's computational results:

* any friend of 91 exceeds 10^100;
* any friend of 91 has at least seven distinct prime factors;
* any friend of 14, 15, 46, 62, 94 exceeds 10^100, 10^99, 10^100, 10^101, 10^80 respectively.

All files are in the `repo/` folder. Run the commands below from inside it (`cd repo`).

## Checking a certificate

Requires Python 3.9 or later; the verifier uses only the standard library.

    python code/verify.py certificates/cert_91_1e100.json 100

The last line should read `Certified: any friend m != 91 satisfies m > 10^100.25`.
The other certificates are checked the same way:

    python code/verify.py certificates/hunt_cert_14_1e100.json 100
    python code/verify.py certificates/hunt_cert_15_1e100.json 99
    python code/verify.py certificates/hunt_cert_46_1e100.json 100
    python code/verify.py certificates/hunt_cert_62_1e100.json 101
    python code/verify.py certificates/hunt_cert_94_1e80.json 80

The verifier reports two bounds. The *leaf-product* bound is the minimum of 2^k * prod p^e over all
cap and unfactored leaves. The *strengthened* bound additionally multiplies each unfactored leaf by the
reduced denominator B of its residual, which must divide the remaining cofactor of m (property (P4) of
the paper); cap leaves are not strengthened. Both are exact integer comparisons. The verdict uses the
strengthened bound unless a third argument `plain` is given, in which case the leaf-product rule is used:

    python code/verify.py certificates/hunt_cert_46_1e100.json 63 plain

(For 46 the leaf-product bound is 10^63.998, so `... 64 plain` is correctly rejected; the strengthened
bound for the same certificate is 10^100.77.)

Every certificate must contain the root node of each alive 2-adic branch, and every other node must be
listed as a child (or cap leaf) by an internal parent node; a certificate reduced to a subset of the tree
is rejected. The adversarial tests in `code/test_verify.py` exercise this and ten other corruptions:

    python code/test_verify.py certificates/cert_91_1e100.json 100

## Contents

| file | purpose |
| --- | --- |
| `code/engine3.py` | search program that produces a size-bound certificate |
| `code/verify.py` | independent verifier for those certificates (shares no code with the search) |
| `code/test_verify.py` | adversarial tests: corrupted certificates that the verifier must reject |
| `code/prove_primes.py` | proves primality of every prime cited in a certificate with FLINT (needs `pip install python-flint`) |
| `code/hunt.py` | command-line wrapper around the search |
| `code/arith.py` | arithmetic helpers used by the search programs |
| `code/known_factors.txt` | external factorizations (YAFU, SymPy); each is re-checked by the verifier |
| `code/omega91.py` | fixed-number-of-primes search for 91 |
| `code/omega_check.py` | independent re-implementation of that search |
| `code/catalogue91.py` | numerical checks of the elementary structure lemmas for a friend of 91 |
| `code/ingest_log.py` | imports YAFU factor logs into `known_factors.txt` |
| `certificates/` | certificates for 91, 14, 15, 46, 62 and 94 |

The search programs need `sympy` (`pip install sympy`). The prime-count runs for 91 are reproduced with

    python code/omega91.py 6 30 3600
    python code/omega_check.py 6 30 3600

and, for seven primes, piece by piece as described in the paper, for example
`python code/omega91.py 7 30 72000 tail 11`.

## Changes since version 1

* `verify.py` now requires the root node of every alive branch and checks that every node is reachable
  from it through an internal parent that lists it (version 1 accepted a certificate reduced to a single
  leaf). It also reports the strengthened bound at unfactored leaves and uses exact integer comparisons.
  These defects and the strengthening were pointed out by an independent audit of version 1.
* `test_verify.py` added.
* `prove_primes.py` added: `python code/prove_primes.py certificates/cert_91_1e100.json` proves all 421 primes cited in
  the certificate of 91 (59 of them beyond the deterministic Miller-Rabin range) with FLINT's proven primality test.
* `omega91.py` and `omega_check.py`: `omega_check.py` no longer uses SymPy's `factorint(..., limit=3000)`, which could
  leave a composite number in SymPy's internal cache and crash a long run with "is not a prime factor"; it now does its
  own trial division. Both programs cap their factorization caches on long runs and append each verdict (or the full
  error) to `omega_log.txt`. The search logic and all reported state counts are unchanged.
* The certificates themselves are unchanged. Under the leaf-product rule the certificate for 46 gives
  10^63.998, not 10^64; under the strengthened rule it gives 10^100.77.

## AI assistance

The search programs, the verifier and the certificates were produced with the help of an AI coding assistant.
All results are reproducible from this repository.

## License

Code: MIT License (see `LICENSE`). Certificates and data: CC BY 4.0.
