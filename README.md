# Friends of 91: certificates and code

Supplementary material for the paper *The Friends of 91 Live Beyond Googol Street* by Kaan Martin Oytac.

A positive integer m is a *friend* of n if m != n and sigma(m)/m = sigma(n)/n. This repository contains
the programs and machine-checkable certificates behind the paper's computational results:

* any friend of 91 exceeds 10^100;
* any friend of 91 has at least seven distinct prime factors;
* any friend of 14, 15, 46, 62, 94 exceeds 10^70, 10^58, 10^64, 10^64, 10^55 respectively.

## Checking a certificate

Requires Python 3.9 or later; the verifier uses only the standard library.

    python code/verify.py certificates/cert_91_1e100.json 100

The last line should read `Certified: any friend m != 91 satisfies m > 10^100.25`.
The other certificates are checked the same way, for example

    python code/verify.py certificates/hunt_cert_14_1e100.json 70

## Contents

| file | purpose |
| --- | --- |
| `code/engine3.py` | search program that produces a size-bound certificate |
| `code/verify.py` | independent verifier for those certificates (shares no code with the search) |
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

## AI assistance

The search programs, the verifier and the certificates were produced with the help of an AI coding assistant.
All results are reproducible from this repository.

## License

Code: MIT License (see `LICENSE`). Certificates and data: CC BY 4.0.
