"""Independent verifier for engine3 certificates.
Usage: python verify.py cert_91_1e60.json 60 [plain]
Shares no code with the search engine: own primality test, own sigma, own enumeration.
It re-derives every claim from the raw arithmetic and reports the certified lower bound
on any friend m != n, or FAILS with the first unjustified step.

Bounds.  Every leaf lambda = ((p1,w1),...,(pt,wt)) of an alive branch k certifies
m >= 2^k * prod p_i^w_i for any friend m consistent with it ("leaf product").  At an
unfactored or free leaf the residual A/B of the leaf is known exactly, and (P4) gives
B | m', so the leaf also certifies m >= 2^k * prod p_i^w_i * B ("strengthened").  A cap
leaf stands for every exponent >= w_t, so its residual is not known and only the leaf
product is used.  By default the strengthened bound is reported and compared with the
claimed exponent X; the optional third argument "plain" uses the leaf product only.
All comparisons are exact integer comparisons; no floating point enters the verdict.

Structure.  Each alive branch must contain its root node (empty path), and every other
node must be the child of an internal node that lists its last exponent as "child" (or
"cap" for a cap leaf); unreachable or orphaned nodes are rejected.  Together with the
re-enumeration of exponents at every internal node this guarantees that the certificate
is the complete tree of Section 4, not a subset of it."""
import sys, json, math
from fractions import Fraction

def is_prime(n):
    if n < 2: return False
    for p in (2,3,5,7,11,13,17,19,23,29,31,37):
        if n % p == 0: return n == p
    d, s = n-1, 0
    while d % 2 == 0: d //= 2; s += 1
    for a in (2,3,5,7,11,13,17,19,23,29,31,37):   # deterministic for n < 3.18e23 (Sorenson-Webster); strong probable prime beyond
        x = pow(a, d, n)
        if x in (1, n-1): continue
        for _ in range(s-1):
            x = x*x % n
            if x == n-1: break
        else: return False
    return True

def sig(p, e): return (p**(e+1) - 1)//(p - 1)
def I(p, e): return Fraction(sig(p, e), p**e)
def v2(n):
    c = 0
    while n % 2 == 0: n //= 2; c += 1
    return c
def v2sig(p, e): return 0 if e % 2 == 0 else v2(p+1) + v2(e+1) - 1

def small_factor(n):
    f, q = {}, 2
    while q*q <= n:
        while n % q == 0: f[q] = f.get(q, 0)+1; n //= q
        q += 1
    if n > 1: f[n] = f.get(n, 0)+1
    return f

fails = []
def fail(msg): fails.append(msg)

def main():
    cert = json.load(open(sys.argv[1])); X = Fraction(sys.argv[2])
    plain = len(sys.argv) > 3 and sys.argv[3] == 'plain'
    n = cert['n']; cap = 10**cert['cap']
    tenX = None if X.denominator != 1 else 10**int(X)
    def above_claim(M):      # exact: is M > 10^X ?  (log10 comparison only when X is not an integer)
        return M > tenX if tenX is not None else math.log10(M) > X
    best = {'plain': None, 'strong': None}     # exact integer minima over all leaves
    def record(which, M):
        best[which] = M if best[which] is None else min(best[which], M)
    nf = small_factor(n)
    T = Fraction(1)
    for p, e in nf.items(): T *= I(p, e)
    if str(T) != cert['T']: fail(f"target mismatch {T} vs {cert['T']}")
    a, b = T.numerator, T.denominator
    # 2-adic split: m = 2^k * odd; verify which k survive
    alive_k = []
    k = v2(b)
    while True:
        r = Fraction(a * 2**k, b * (2**(k+1) - 1))
        if r < 1: break
        if r.denominator % 2 == 1: alive_k.append((k, r))
        k += 1
        if k > 60: fail("2-adic tail not finite"); break
    cert_k = {br['k'] for br in cert['branches']}
    for k, r in alive_k:
        if k not in cert_k: fail(f"alive 2-adic branch k={k} missing from certificate")
    for br in cert['branches']:
        k = br['k']
        if not br['complete']: fail(f"branch k={k} incomplete (timed out)")
        r0 = dict(alive_k).get(k)
        if r0 is None: fail(f"branch k={k} is not alive"); continue
        nodes = br['nodes']
        bypath = {tuple(map(tuple, nd['path'])): nd for nd in nodes}
        if len(bypath) != len(nodes): fail("duplicate node paths")
        if () not in bypath: fail(f"branch k={k}: root node missing"); continue
        for nd in nodes:
            path = tuple(map(tuple, nd['path']))
            if path:   # every non-root node must hang off an internal node that lists it
                parent = bypath.get(path[:-1])
                if parent is None or parent.get('kind') != 'internal':
                    fail(f"node {path}: parent missing or not internal"); continue
                p, e = path[-1]
                if parent.get('branch_prime') != p: fail(f"node {path}: parent branches on {parent.get('branch_prime')}, not {p}")
                want = 'cap' if nd['kind'] == 'cap' else 'child'
                lab = next((it[1] for it in parent['exponents'] if int(it[0]) == e), None)
                if lab != want: fail(f"node {path}: parent labels exponent {e} as {lab}, expected {want}")
            # recompute residual from the path
            r = r0
            for p, e in path: r /= I(p, e)
            P = 2**k
            for p, e in path: P *= p**e
            if nd['kind'] == 'cap':
                p, e = path[-1]
                if sig(p, e) <= cap: fail(f"cap node {path}: sigma not above cap")
                record('plain', P); record('strong', P)
                if not above_claim(P): fail(f"cap node {path}: leaf product 10^{math.log10(P):.4f} below claimed bound 10^{float(X)}")
                continue
            if str(r) != nd['r']: fail(f"node {path}: residual mismatch {r} vs {nd['r']}"); continue
            A, B = r.numerator, r.denominator
            used = {p for p, _ in path}
            kind = nd['kind']
            if kind == 'solution':
                m = 2**k
                for p, e in path: m *= p**e
                if r != 1: fail(f"solution node {path}: residual {r} != 1")
                if m != n: fail(f"FRIEND CLAIMED: {m}")
                continue
            if kind in ('unfactored', 'free'):
                record('plain', P); record('strong', P * B)      # B | m' by (P4), so m >= 2^k * prod p^e * B
                M = P if plain else P * B
                if not above_claim(M): fail(f"{kind} node {path}: certified 10^{math.log10(M):.4f} below claimed bound 10^{float(X)}")
                continue
            if kind == 'dead':
                why = nd['reason']
                if why == 'even denominator':
                    if B % 2: fail(f"{path}: denominator is odd");
                    continue
                if why == 'below 1':
                    if r >= 1: fail(f"{path}: residual not below 1")
                    continue
                fac = dict((int(q), int(e)) for q, e in nd.get('den_factors', nd.get('den_partial', [])))
                for q, e in fac.items():
                    if not is_prime(q) or B % q**e: fail(f"{path}: bad factor {q}^{e}")
                if why == 'used prime in denominator':
                    if not (set(fac) & used): fail(f"{path}: no used prime among den factors")
                    continue
                if why == 'forced primes exceed target':
                    lb = Fraction(1)
                    for q, e in fac.items(): lb *= I(q, e)
                    if lb <= r: fail(f"{path}: lower bound {float(lb):.6f} does not exceed {float(r):.6f}")
                    continue
                if why == 'no viable exponent':
                    pass  # falls through to the enumeration check below
                else:
                    fail(f"{path}: unknown reason {why}"); continue
            if kind not in ('internal', 'dead'): fail(f"{path}: unknown kind {kind}"); continue
            # full factorization of the denominator must be supplied and correct
            fac = dict((int(q), int(e)) for q, e in nd['den_factors'])
            prod = 1
            for q, e in fac.items():
                if not is_prime(q): fail(f"{path}: {q} not prime")
                prod *= q**e
            if prod != B: fail(f"{path}: den_factors do not multiply to denominator")
            if set(fac) & used: fail(f"{path}: used prime forced again, should be dead")
            lb = Fraction(1)
            for q, e in fac.items(): lb *= I(q, e)
            if lb > r: fail(f"{path}: should be dead by lower bound")
            p = nd['branch_prime']
            if p not in fac: fail(f"{path}: branch prime {p} not forced"); continue
            emin = fac[p]
            lb_others = lb / I(p, emin)
            jb = v2(A)
            F = nd.get('forbidden', [])
            for q in F:
                if not is_prime(q) or q in fac or q in used or A % q == 0 or lb * Fraction(q+1, q) <= r:
                    fail(f"{path}: {q} wrongly forbidden")
            # independent enumeration of exponents
            listed = {int(item[0]): item for item in nd['exponents']}
            e = emin
            children = 0
            while True:
                s = sig(p, e)
                item = listed.get(e)
                if item is None: fail(f"{path}: exponent {e} of {p} not accounted for"); break
                st = item[1]
                if s > cap:
                    if st != 'cap': fail(f"{path}: e={e} should be cap")
                    if tuple(list(path) + [(p, e)]) not in bypath: fail(f"{path}: cap node for e={e} missing")
                    children += 1; break
                if I(p, e) * lb_others > r:
                    if st != 'stop_size': fail(f"{path}: e={e} should be stop_size")
                    break
                if e % 2 == 1 and v2sig(p, e) > jb:
                    if st != 'skip_parity': fail(f"{path}: e={e} should be skip_parity")
                    e += 1; continue
                bad = next((q for q in F if s % q == 0), None)
                if bad is not None:
                    if st != 'skip_forbidden': fail(f"{path}: e={e} should be skip_forbidden")
                    e += 1; continue
                if st != 'child': fail(f"{path}: e={e} should be a child, certificate says {st}")
                if tuple(list(path) + [(p, e)]) not in bypath: fail(f"{path}: child node for e={e} missing")
                children += 1; e += 1
            if kind == 'dead' and children: fail(f"{path}: claimed no viable exponent but found {children}")
            if kind == 'internal' and not children: fail(f"{path}: internal node with no children")
    if fails:
        print("VERIFICATION FAILED:")
        for f in fails[:30]: print("  ", f)
        if len(fails) > 30: print("  ...", len(fails), "failures")
        sys.exit(1)
    print(f"VERIFIED: every node of the certificate for n={n} is justified.")
    print(f"Leaf-product bound:  any friend m != {n} satisfies m >= {best['plain']}  (10^{math.log10(best['plain']):.4f})")
    print(f"Strengthened bound:  any friend m != {n} satisfies m >= {best['strong']}  (10^{math.log10(best['strong']):.4f})")
    M = best['plain'] if plain else best['strong']
    print(f"Certified: any friend m != {n} satisfies m > 10^{math.log10(M):.2f}  (claimed 10^{float(X)}, rule: {'leaf product' if plain else 'strengthened'})")
    if not above_claim(M): print("WARNING: certified bound is below the claimed one"); sys.exit(1)

main()
