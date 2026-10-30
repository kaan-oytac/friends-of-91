"""omega91.py  --  prove that a friend of 91 has more than k distinct prime factors.

Nielsen / Thackeray style search adapted to target T = 16/13 (odd m, 13 | m, 3,5,7 do not divide m,
v2(sigma(m)) = 4 so at most four odd exponents).

State: on  = list of (p, a, exact)   exact exponent a, or exponent >= a when exact is False
       off = dict  q -> minimum exponent (q known to divide m, exponent not yet fixed)
       sup = dict  r -> amount of r supplied so far to sigma(m) by the sigma(p^a) of exact on primes
       j   = 2-adic budget used by exact odd exponents
       P   = every prime not in on/off and not forbidden that divides m exceeds P
Outcomes per leaf: dead, solution (a friend!), open (method could not close the branch).
A run with k, B returns counts; if open == 0 and no solution other than 91, omega(m) != k.

usage: python omega91.py k B_exponent [time_limit_seconds] [exponent_of_13 | tail]
"""
import sys, time, json
from fractions import Fraction
from arith import factorize, sigma_pp, is_probable_prime, v2
from sympy import primerange

T = Fraction(16, 13)
FORB = {3, 5, 7}
J = 4
def _val(n):
    f = {}
    for q, e in factorize(n).items():
        if q != 2: f[q] = e
    return f
def setup():
    global ANUM, DEN, J
    ANUM = _val(T.numerator)      # odd numerator primes absorb this much supply without dividing m
    DEN = _val(T.denominator)     # denominator primes: v_r(sigma(m)) = v_r(m) - DEN[r]
    J = v2(T.numerator)
setup()

def I(p, a): return Fraction(sigma_pp(p, a), p**a)
def Iinf(p): return Fraction(p, p - 1)
def v2sig(p, a): return 0 if a % 2 == 0 else v2(p + 1) + v2(a + 1) - 1

class Timeout(Exception): pass

from arith import _cache as _KNOWN, _rho_bounded
from sympy import perfect_power, primerange as _pr
_SMALL = list(_pr(3, 3000))
_QF = {}
def quick_factor(s):
    """Factor s as far as cheap methods allow. Returns (primes dict, list of leftover composites).
    Leftover composites have no prime factor below 3000 and are not prime."""
    if s in _QF: return _QF[s]
    if s in _KNOWN: r = (dict(_KNOWN[s]), []); _QF[s] = r; return r
    f, rest = {}, s
    while rest % 2 == 0: f[2] = f.get(2, 0) + 1; rest //= 2
    for q in _SMALL:
        if q*q > rest: break
        while rest % q == 0: f[q] = f.get(q, 0) + 1; rest //= q
    comps, stack = [], ([rest] if rest > 1 else [])
    while stack:
        x = stack.pop()
        if x in _KNOWN:
            for q, e in _KNOWN[x].items(): f[q] = f.get(q, 0) + e
            continue
        if is_probable_prime(x): f[x] = f.get(x, 0) + 1; continue
        d = _rho_bounded(x, 3000)
        if d: stack += [d, x // d]
        else: comps.append(x)
    if len(_QF) > 200000: _QF.clear()               # bounded memory on long runs
    r = (f, comps); _QF[s] = r; return r

def guaranteed(C):
    """number of distinct primes certainly dividing composite C (not prime, no small factors)"""
    return 1 if perfect_power(C) else 2

class Omega:
    def __init__(self, k, B, limit):
        self.k, self.B, self.deadline = k, B, time.time() + limit
        self.t0 = time.time()
        self.stats = dict(nodes=0, dead=0, open=0, solution=0, candidates=[], opens=[])
        self.dead_reasons = {}
        self.only2 = None    # restrict the second prime when 13 is in its tail
        self.only13 = None   # None = all; int = only that exact exponent of 13; 'tail' = only the tail
    def dead(self, why):
        self.stats['dead'] += 1; self.dead_reasons[why] = self.dead_reasons.get(why, 0) + 1
    def bounds(self, on, off):
        kmin = Fraction(1); kmax = Fraction(1)
        for p, a, ex in on:
            kmin *= I(p, a); kmax *= I(p, a) if ex else Iinf(p)
        for q, b in off.items():
            kmin *= I(q, b); kmax *= Iinf(q)
        return T / kmin, T / kmax   # (M, m_low): unknown part U lies in (m_low, M]
    def rec(self, on, off, sup, j, P, comps=()):
        if time.time() > self.deadline: raise Timeout
        self.stats['nodes'] += 1
        if self.stats['nodes'] % 20000 == 0:
            print(f"   ... {self.stats['nodes']} nodes, dead {self.stats['dead']}, open {self.stats['open']}, {time.time()-self.t0:.0f}s", flush=True)
        u = self.k - len(on) - len(off)
        if u < 0: self.dead('too many primes'); return
        if comps and sum(guaranteed(C) for C in comps) > u: self.dead('unfactored cofactors need more primes'); return
        M, mlow = self.bounds(on, off)
        if M < 1: self.dead('known part exceeds target'); return
        if u == 0 and not off:
            if all(ex for _, _, ex in on):
                if M == 1:
                    m = 1
                    for p, a, _ in on: m *= p**a
                    self.stats['solution'] += 1; self.stats['candidates'].append(('SOLUTION', m)); return
                self.dead('all exact, product not equal to target'); return
            if not (mlow < 1 <= M): self.dead('inexact primes cannot reach target'); return
            inex = [(p, a) for p, a, ex in on if not ex]
            if len(inex) == 1:
                p, a0 = inex[0]
                rest = Fraction(1)
                for q, a, ex in on:
                    if ex: rest *= I(q, a)
                need = T / rest
                a = a0
                while I(p, a) < need: a += 1
                if I(p, a) == need:
                    m = p**a
                    for q, b, ex in on:
                        if ex: m *= q**b
                    self.stats['solution'] += 1; self.stats['candidates'].append(('SOLUTION', m)); return
                self.dead('single inexact exponent has no solution'); return
            self.stats['open'] += 1; self.stats['opens'].append(('several inexact exponents', [(p, a, ex) for p, a, ex in on])); return
        if u == 0 and off:
            pass  # fall through: off primes still need exponents
        if u > 0 and M <= 1: self.dead('no room for further primes'); return
        if off:
            q = min(off); b = off[q]
            self.expand(q, b, on, {r: e for r, e in off.items() if r != q}, sup, j, P, comps)
            return
        # off empty, u >= 1: choose the smallest unknown prime p
        if mlow <= 1:
            self.stats['open'] += 1; self.stats['opens'].append(('no upper bound for next prime', [(p, a, ex) for p, a, ex in on], str(M), str(mlow))); return
        Bhigh = 1 + Fraction(u) / (mlow - 1)
        Blow = 1 / (M - 1) if M > 1 else 0
        lo = max(P, int(Blow) - 1, 10)
        used = {p for p, _, _ in on}
        second = (self.only2 is not None and len(on) == 1 and on[0][0] == 13 and not on[0][2])
        if second and len(on) == 1:
            print(f"   second-prime candidates in this tail: {[q for q in primerange(lo + 1, int(Bhigh) + 2) if q not in FORB and q not in used and Fraction(q) < Bhigh and Fraction(q + 1, q) <= M]}", flush=True)
        for p in primerange(lo + 1, int(Bhigh) + 2):
            if second and p != self.only2: continue
            if p in FORB or p in used or Fraction(p) >= Bhigh: continue
            if Fraction(p + 1, p) > M: continue
            self.expand(p, 1, on, {}, sup, j, p, comps)
    def strip(self, comps, primes, P, off):
        """Divide listed primes out of the pending composites, split them into pairwise coprime pieces,
        and turn any piece that is prime into a new off prime. Returns (comps, off) or None on contradiction.
        Pairwise coprimality makes the later count sum(guaranteed(C)) a sound lower bound on the number of
        distinct unlisted primes."""
        from math import gcd
        work = [C for C in comps]
        off = dict(off)
        listed = set(primes) | set(off)
        changed = True
        while changed:
            changed = False
            out = []
            for C in work:
                for q in listed:
                    while C % q == 0: C //= q
                if C == 1: continue
                if is_probable_prime(C):
                    if C in FORB: return None
                    if C not in off:
                        if C <= P: return None
                        off[C] = 1; listed.add(C); changed = True
                    continue
                out.append(C)
            # split pairwise until coprime
            i = 0
            while i < len(out):
                j = i + 1
                while j < len(out):
                    g = gcd(out[i], out[j])
                    if g > 1:
                        a, b = out[i], out[j]
                        pieces = {g, a // g, b // g}
                        pieces.discard(1)
                        out = out[:i] + out[i+1:j] + out[j+1:] + sorted(pieces)
                        changed = True
                        i = -1; break
                    j += 1
                i += 1
            work = out
        return tuple(sorted(set(work))), off
    def expand(self, p, b, on, off, sup, j, P, comps=()):
        """branch on the exact exponent a >= b of prime p (finite powers up to B, then the tail)."""
        a = b
        root = (p == 13 and not on and self.only13 is not None)
        while p**a <= self.B:
            if a % 2 == 1 and j + v2sig(p, a) > J: a += 1; continue
            if root and self.only13 != a: a += 1; continue
            s = sigma_pp(p, a)
            fs, newc = quick_factor(s)
            if newc: self.stats['unfactored'] = self.stats.get('unfactored', 0) + 1
            ok = True; new_off = dict(off); new_sup = dict(sup)
            for r, e in fs.items():
                if r == 2: continue
                if r in FORB: ok = False; break
                new_sup[r] = new_sup.get(r, 0) + e
                need = new_sup[r] - ANUM.get(r, 0)    # amount that m itself must carry
                if need <= 0: continue                 # absorbed by the numerator of T
                exact = next(((aa, ex) for q, aa, ex in on if q == r), None)
                if exact is not None:
                    aa, ex = exact
                    if ex and need > aa - DEN.get(r, 0): ok = False; break
                elif r in new_off:
                    new_off[r] = max(new_off[r], need)
                else:
                    if r <= P: ok = False; break      # a new prime below P is excluded by the ordering
                    new_off[r] = need
            if ok:
                listed = [q for q, _, _ in on] + [p] + list(new_off)
                st = self.strip(tuple(comps) + tuple(newc), listed, P, new_off)
                if st is None: self.dead('cofactor reduces to an excluded prime')
                else:
                    c2, o2 = st
                    self.rec(on + [(p, a, True)], o2, new_sup, j + v2sig(p, a), P, c2)
            a += 1
        # tail: exponent >= a, sigma unknown (relaxation)
        if root and self.only13 != 'tail': return
        st = self.strip(tuple(comps), [p], P, off)
        if st is None: self.dead('cofactor reduces to an excluded prime'); return
        c2, o2 = st
        self.rec(on + [(p, a, False)], o2, sup, j, P, c2)

def run(k, Bexp, limit, only13=None, only2=None):
    S = Omega(k, 10**Bexp, limit); S.only13 = only13; S.only2 = only2
    t0 = time.time(); done = True
    try: S.rec([], {13: 2}, {}, 0, 0)
    except Timeout: done = False
    return S, done, time.time() - t0

def _main():
    k = int(sys.argv[1]); Bexp = int(sys.argv[2]); limit = float(sys.argv[3]) if len(sys.argv) > 3 else 600
    only13 = None
    if len(sys.argv) > 4: only13 = 'tail' if sys.argv[4] == 'tail' else int(sys.argv[4])
    only2 = int(sys.argv[5]) if len(sys.argv) > 5 else None
    S, done, secs = run(k, Bexp, limit, only13, only2)
    st = S.stats
    sol = [c for c in st['candidates'] if c[0] == 'SOLUTION' and c[1] != 91]
    piece = '' if only13 is None else (' [piece: 13^%s%s]' % (only13 if only13 != 'tail' else 'tail', '' if only2 is None else ', next prime %d' % only2))
    verdict = ('omega != %d PROVED' % k + piece) if (done and st['open'] == 0 and not sol) else ('TIMEOUT' + piece if not done else 'INCONCLUSIVE' + piece)
    print(f"k={k} B=1e{Bexp} {verdict} in {secs:.1f}s: nodes={st['nodes']} dead={st['dead']} open={st['open']} solutions={st['solution']} (other than 91: {len(sol)})")
    for why, c in sorted(S.dead_reasons.items(), key=lambda x: -x[1]): print(f"   dead: {c:7d}  {why}")
    for o in st['opens'][:8]: print("   OPEN:", o)
    if sol: print("   SOLUTIONS:", sol)
    json.dump(dict(k=k, B=Bexp, verdict=verdict, stats={a: b for a, b in st.items() if a not in ('candidates', 'opens')}, opens=[str(o) for o in st['opens'][:200]], solutions=sol), open(f'omega91_k{k}_B{Bexp}' + ('' if only13 is None else f'_13e{only13}') + ('' if only2 is None else f'_p{only2}') + '.json', 'w'))
    return verdict

if __name__ == '__main__':
    import multiprocessing, traceback, datetime
    multiprocessing.set_start_method('spawn', force=True)
    try:
        v = _main()
        open('omega_log.txt', 'a').write(f"{datetime.datetime.now()}  {' '.join(sys.argv)}  ->  {v}\n")
    except BaseException:
        open('omega_log.txt', 'a').write(f"{datetime.datetime.now()}  {' '.join(sys.argv)}  ->  ERROR\n{traceback.format_exc()}\n")
        raise
