"""omega_check.py -- independent re-implementation of the fixed-omega search for friends of 91.

Written separately from omega91.py: iterative instead of recursive, its own bound formulas, its own
partial factorization (SymPy only), its own handling of pending composites. It shares no code with
omega91.py or arith.py. Verdicts must agree with omega91.py on every piece.

usage: python omega_check.py k B_exponent [time_limit] [exponent_of_13 | tail] [next_prime]
"""
import sys, time, json
from fractions import Fraction
from math import gcd
from sympy import isprime, factorint, primerange, perfect_power
from sympy.ntheory.factor_ import pollard_rho

TARGET = Fraction(16, 13)          # I(91)
EXCLUDED = (3, 5, 7)               # proved not to divide a friend of 91
TWO_ADIC = 4                       # v2(sigma(m)) = v2(16)
SPECIAL = 13                       # v13(sigma(m)) = v13(m) - 1

def sig(p, a): return (p**(a + 1) - 1) // (p - 1)
def abund(p, a): return Fraction(sig(p, a), p**a)
def ceiling(p): return Fraction(p, p - 1)
def two_adic_of_sigma(p, a):
    if a % 2 == 0: return 0
    s = sig(p, a); c = 0
    while s % 2 == 0: s //= 2; c += 1
    return c

_fcache = {}
_SMALL_PRIMES = list(primerange(2, 3000))
def partial_factor(n):
    """(dict prime->exp, list of composite leftovers with no factor < 3000)"""
    if n in _fcache: return _fcache[n]
    # own trial division (SymPy's factorint(limit=...) can poison its internal cache with composites)
    primes, comps, stack = {}, [], []
    rest = n
    for q in _SMALL_PRIMES:
        if q * q > rest: break
        while rest % q == 0: primes[q] = primes.get(q, 0) + 1; rest //= q
    if rest > 1: stack.append(rest)
    while stack:
        c = stack.pop()
        if isprime(c): primes[c] = primes.get(c, 0) + 1; continue
        d = pollard_rho(c, max_steps=4000)
        if d is None or d in (1, c): comps.append(c)
        else: stack += [d, c // d]
    if len(_fcache) > 200000: _fcache.clear()      # bounded memory on long runs
    res = (primes, comps); _fcache[n] = res; return res

def coprime_pieces(comps):
    """split a list of composites into pairwise coprime pieces (drop 1s)"""
    out = [c for c in comps if c > 1]
    again = True
    while again:
        again = False
        for i in range(len(out)):
            for j in range(i + 1, len(out)):
                g = gcd(out[i], out[j])
                if g > 1:
                    a, b = out[i], out[j]
                    out = [x for t, x in enumerate(out) if t not in (i, j)] + [x for x in {g, a // g, b // g} if x > 1]
                    again = True; break
            if again: break
    return out

def min_new_primes(comps):
    return sum(1 if perfect_power(c) else 2 for c in comps)

class Node:
    __slots__ = ('on', 'off', 'supply', 'twos', 'floor', 'pending')
    def __init__(self, on, off, supply, twos, floor, pending):
        self.on, self.off, self.supply, self.twos, self.floor, self.pending = on, off, supply, twos, floor, pending

def settle_pending(pending, on, off, floor):
    """divide listed primes out of pending composites; prime leftovers become off primes.
    Returns (pending, off) or None if a leftover prime is excluded or violates the ordering."""
    off = dict(off)
    listed = set(p for p, _, _ in on) | set(off)
    comps = list(pending)
    while True:
        new = []
        for c in comps:
            for q in listed:
                while c % q == 0: c //= q
            if c > 1: new.append(c)
        comps = coprime_pieces(new)
        found = False
        for c in list(comps):
            if isprime(c):
                comps.remove(c)
                if c in EXCLUDED or c <= floor: return None
                if c not in off: off[c] = 1; listed.add(c); found = True
        if not found: return tuple(comps), off

def check(k, Bexp, limit, piece13=None, piece2=None):
    B = 10**Bexp
    t0 = time.time()
    stats = {'nodes': 0, 'dead': 0, 'open': 0, 'solutions': [], 'opens': 0}
    stack = [Node((), {SPECIAL: 2}, {}, 0, 0, ())]
    while stack:
        if time.time() - t0 > limit: return stats, False, time.time() - t0
        nd = stack.pop(); stats['nodes'] += 1
        on, off, pending = nd.on, nd.off, nd.pending
        u = k - len(on) - len(off)
        if u < 0: stats['dead'] += 1; continue
        if pending and min_new_primes(pending) > u: stats['dead'] += 1; continue
        lo_known = Fraction(1); hi_known = Fraction(1)
        for p, a, exact in on:
            lo_known *= abund(p, a); hi_known *= (abund(p, a) if exact else ceiling(p))
        for q, b in off.items():
            lo_known *= abund(q, b); hi_known *= ceiling(q)
        U_max = TARGET / lo_known; U_min = TARGET / hi_known      # unknown part in (U_min, U_max]
        if U_max < 1: stats['dead'] += 1; continue
        if u == 0 and not off:
            inexact = [(p, a) for p, a, ex in on if not ex]
            if not inexact:
                if U_max == 1:
                    m = 1
                    for p, a, _ in on: m *= p**a
                    stats['solutions'].append(m)
                else: stats['dead'] += 1
                continue
            if not (U_min < 1 <= U_max): stats['dead'] += 1; continue
            if len(inexact) == 1:
                p, a0 = inexact[0]
                rest = Fraction(1)
                for q, a, ex in on:
                    if ex: rest *= abund(q, a)
                need = TARGET / rest; a = a0
                while abund(p, a) < need: a += 1
                if abund(p, a) == need:
                    m = p**a
                    for q, b, ex in on:
                        if ex: m *= q**b
                    stats['solutions'].append(m)
                else: stats['dead'] += 1
                continue
            stats['open'] += 1; continue
        if u > 0 and U_max <= 1: stats['dead'] += 1; continue
        # which prime to branch on
        if off:
            p = min(off); bmin = off[p]; rest_off = {q: b for q, b in off.items() if q != p}
            candidates = [(p, bmin, rest_off, nd.floor)]
        else:
            if U_min <= 1: stats['open'] += 1; continue
            hi = 1 + Fraction(u) / (U_min - 1)
            lo = Fraction(1) / (U_max - 1)
            start = max(nd.floor, int(lo) - 1, 10)
            used = {q for q, _, _ in on}
            candidates = []
            second_level = (piece2 is not None and len(on) == 1 and on[0][0] == SPECIAL and not on[0][2])
            for p in primerange(start + 1, int(hi) + 2):
                if p in EXCLUDED or p in used or Fraction(p) >= hi: continue
                if Fraction(p + 1, p) > U_max: continue
                if second_level and p != piece2: continue
                candidates.append((p, 1, {}, p))
        root_special = (not on and off == {SPECIAL: 2} and piece13 is not None)
        for p, bmin, rest_off, floor in candidates:
            a = bmin
            children = []
            while p**a <= B:
                if a % 2 == 1 and nd.twos + two_adic_of_sigma(p, a) > TWO_ADIC: a += 1; continue
                if root_special and piece13 != a: a += 1; continue
                primes, comps = partial_factor(sig(p, a))
                ok = True; new_off = dict(rest_off); new_sup = dict(nd.supply)
                for r, e in primes.items():
                    if r == 2: continue
                    if r in EXCLUDED: ok = False; break
                    new_sup[r] = new_sup.get(r, 0) + e
                    need = new_sup[r]
                    entry = next(((aa, ex) for q, aa, ex in on if q == r), None)
                    if entry is not None:
                        aa, ex = entry
                        limit_r = aa - (1 if r == SPECIAL else 0)
                        if ex and need > limit_r: ok = False; break
                    elif r in new_off: new_off[r] = max(new_off[r], need)
                    else:
                        if r <= floor: ok = False; break
                        new_off[r] = need
                if ok:
                    on2 = on + ((p, a, True),)
                    st = settle_pending(tuple(pending) + tuple(comps), on2, new_off, floor)
                    if st is None: stats['dead'] += 1
                    else:
                        pend2, off2 = st
                        children.append(Node(on2, off2, new_sup, nd.twos + two_adic_of_sigma(p, a), floor, pend2))
                a += 1
            if not (root_special and piece13 != 'tail'):
                on2 = on + ((p, a, False),)
                st = settle_pending(tuple(pending), on2, rest_off, floor)
                if st is None: stats['dead'] += 1
                else:
                    pend2, off2 = st
                    children.append(Node(on2, off2, nd.supply, nd.twos, floor, pend2))
            stack.extend(reversed(children))
    return stats, True, time.time() - t0

def _main():
    k = int(sys.argv[1]); Bexp = int(sys.argv[2]); limit = float(sys.argv[3]) if len(sys.argv) > 3 else 600
    piece13 = None
    if len(sys.argv) > 4: piece13 = 'tail' if sys.argv[4] == 'tail' else int(sys.argv[4])
    piece2 = int(sys.argv[5]) if len(sys.argv) > 5 else None
    stats, done, secs = check(k, Bexp, limit, piece13, piece2)
    others = [m for m in stats['solutions'] if m != 91]
    tag = '' if piece13 is None else f" [piece: 13^{piece13}" + ('' if piece2 is None else f", next prime {piece2}") + "]"
    verdict = (f"CHECK: omega != {k} PROVED" if done and stats['open'] == 0 and not others else ("CHECK: TIMEOUT" if not done else "CHECK: INCONCLUSIVE")) + tag
    print(f"{verdict} in {secs:.1f}s: nodes={stats['nodes']} dead={stats['dead']} open={stats['open']} solutions={stats['solutions']}")
    json.dump(dict(k=k, B=Bexp, verdict=verdict, stats=stats), open(f"check_k{k}_B{Bexp}" + ('' if piece13 is None else f"_13e{piece13}") + ('' if piece2 is None else f"_p{piece2}") + ".json", 'w'))
    return verdict

if __name__ == '__main__':
    import traceback, datetime
    try:
        v = _main()
        open('omega_log.txt', 'a').write(f"{datetime.datetime.now()}  {' '.join(sys.argv)}  ->  {v}\n")
    except BaseException:
        open('omega_log.txt', 'a').write(f"{datetime.datetime.now()}  {' '.join(sys.argv)}  ->  ERROR\n{traceback.format_exc()}\n")
        raise
