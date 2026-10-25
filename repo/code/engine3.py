"""engine3: certificate-producing search for friends of n (odd-part branches).
Every node is recorded with enough data for an independent verifier (verify.py).
Soundness fixes vs engine2: no pruning on used primes dividing sigma(p^e);
forbidden-prime pruning only when the forbidden prime does not divide the residual numerator."""
import time, json
from fractions import Fraction
from arith import factorize, sigma_pp, sigma, is_probable_prime, v2

class Timeout(Exception): pass

def v2sig(p, e): return 0 if e % 2 == 0 else v2(p+1) + v2(e+1) - 1
def I(p, e): return Fraction(sigma_pp(p, e), p**e)

def forbidden(r, Bf, used, upto=1000):
    lb = Fraction(1)
    for q, k in Bf.items(): lb *= I(q, k)
    A = r.numerator
    F = []
    q = 3
    while q < upto:
        if is_probable_prime(q) and q not in Bf and q not in used and A % q and lb * Fraction(q+1, q) > r: F.append(q)
        q += 2
    return F

def two_adic_branches(n):
    T = Fraction(sigma(n), n); a, b = T.numerator, T.denominator
    out, k = [], v2(b)
    while True:
        r = Fraction(a * 2**k, b * (2**(k+1) - 1))
        if r < 1: out.append(dict(k=k, r=str(r), status='dead: below 1')); break
        out.append(dict(k=k, r=str(r), status=('dead: even denominator' if r.denominator % 2 == 0 else 'alive')))
        k += 1
        if k > 40: out.append(dict(k=k, status='TAIL')); break
    return out

class Search:
    def __init__(self, n, cap, budget):
        self.n, self.cap, self.deadline = n, cap, time.time() + budget
        self.nodes = []
    def rec(self, path, r, k):
        if time.time() > self.deadline: raise Timeout
        node = dict(path=[list(x) for x in path], r=str(r))
        self.nodes.append(node)
        A, B = r.numerator, r.denominator
        if B % 2 == 0: node['kind'] = 'dead'; node['reason'] = 'even denominator'; return
        if r < 1: node['kind'] = 'dead'; node['reason'] = 'below 1'; return
        if r == 1:
            m = 2**k
            for p, e in path: m *= p**e
            node['kind'] = 'solution'; node['m'] = str(m); return
        used = {p for p, _ in path}
        small, rest, q = {}, B, 3
        while q < 100000 and q*q <= rest:
            while rest % q == 0: small[q] = small.get(q, 0)+1; rest //= q
            q += 2
        if any(p in used for p in small):
            node['kind'] = 'dead'; node['reason'] = 'used prime in denominator'; node['den_partial'] = sorted(small.items()); return
        lb = Fraction(1)
        for q, e in small.items(): lb *= I(q, e)
        if lb > r:
            node['kind'] = 'dead'; node['reason'] = 'forced primes exceed target'; node['den_partial'] = sorted(small.items()); return
        try: Bf = factorize(B)
        except RuntimeError:
            node['kind'] = 'unfactored'; node['cofactor'] = str(rest); return
        node['den_factors'] = sorted(Bf.items())
        if any(p in used for p in Bf): node['kind'] = 'dead'; node['reason'] = 'used prime in denominator'; return
        if B == 1: node['kind'] = 'free'; return
        lb = Fraction(1)
        for q, e in Bf.items(): lb *= I(q, e)
        if lb > r: node['kind'] = 'dead'; node['reason'] = 'forced primes exceed target'; return
        F = forbidden(r, Bf, used)
        jb = v2(A)
        best = None
        for p, emin in Bf.items():
            ex = self.exponents(p, emin, r, lb / I(p, emin), jb, F)
            cnt = len([x for x in ex if x[1] in ('child', 'cap')])
            if best is None or cnt < best[2]: best = (p, ex, cnt)
        p, ex, cnt = best
        node['branch_prime'] = p; node['forbidden'] = F; node['exponents'] = ex
        if cnt == 0:
            node['kind'] = 'dead'; node['reason'] = 'no viable exponent'; return
        node['kind'] = 'internal'
        for item in ex:
            e, status = item[0], item[1]
            if status == 'child': self.rec(path + [(p, e)], r / I(p, e), k)
            elif status == 'cap':
                self.nodes.append(dict(path=[list(x) for x in path] + [[p, e]], kind='cap', r=str(r), bound_exp=e))
    def exponents(self, p, emin, r, lb_others, jb, F):
        out, e = [], emin
        while True:
            s = sigma_pp(p, e)
            if s > self.cap: out.append((e, 'cap')); break
            if I(p, e) * lb_others > r: out.append((e, 'stop_size')); break
            if e % 2 == 1 and v2sig(p, e) > jb: out.append((e, 'skip_parity')); e += 1; continue
            bad = next((q for q in F if s % q == 0), None)
            if bad: out.append((e, 'skip_forbidden', bad)); e += 1; continue
            out.append((e, 'child')); e += 1
        return out

def run(n, cap_exp, budget):
    T = Fraction(sigma(n), n)
    brs = two_adic_branches(n)
    cert = dict(n=n, T=str(T), cap=cap_exp, two_adic=brs, branches=[])
    for b in brs:
        if b.get('status') != 'alive': continue
        S = Search(n, 10**cap_exp, budget)
        try: S.rec([], Fraction(b['r']), b['k']); done = True
        except Timeout: done = False
        cert['branches'].append(dict(k=b['k'], complete=done, nodes=S.nodes))
    return cert

if __name__ == '__main__':
    import sys
    n, cap_exp, budget = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3])
    cert = run(n, cap_exp, budget)
    json.dump(cert, open(f'cert_{n}_1e{cap_exp}.json', 'w'))
    for br in cert['branches']:
        kinds = {}
        for nd in br['nodes']: kinds[nd['kind']] = kinds.get(nd['kind'], 0) + 1
        print(f"n={n} k={br['k']} complete={br['complete']} nodes={kinds}")
