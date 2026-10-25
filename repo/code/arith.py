import random, subprocess
from math import gcd
from fractions import Fraction

def is_probable_prime(n):
    if n < 2: return False
    for p in (2,3,5,7,11,13,17,19,23,29,31,37):
        if n % p == 0: return n == p
    d, s = n-1, 0
    while d % 2 == 0: d //= 2; s += 1
    for a in (2,3,5,7,11,13,17,19,23,29,31,37):
        x = pow(a, d, n)
        if x in (1, n-1): continue
        for _ in range(s-1):
            x = x*x % n
            if x == n-1: break
        else: return False
    return True

def pollard_brent(n):
    if n % 2 == 0: return 2
    if n % 3 == 0: return 3
    while True:
        y, c, m = random.randrange(1,n), random.randrange(1,n), 128
        g = r = q = 1
        while g == 1:
            x = y
            for _ in range(r): y = (y*y+c) % n
            k = 0
            while k < r and g == 1:
                ys = y
                for _ in range(min(m, r-k)):
                    y = (y*y+c) % n; q = q*abs(x-y) % n
                g = gcd(q, n); k += m
            r *= 2
        if g == n:
            g = 1
            while g == 1:
                ys = (ys*ys+c) % n; g = gcd(abs(x-ys), n)
        if g != n: return g

_cache = {}

def _sympy_server(conn):
    import sympy
    while True:
        m = conn.recv()
        if m is None: break
        conn.send(dict(sympy.factorint(m)))

_W = {'proc': None, 'conn': None}
def _worker():
    import multiprocessing as mp
    if _W['proc'] is None or not _W['proc'].is_alive():
        a, b = mp.Pipe()
        p = mp.Process(target=_sympy_server, args=(b,), daemon=True)
        p.start(); b.close()
        _W['proc'], _W['conn'] = p, a
    return _W['conn']

def _sympy_factor_with_timeout(m, seconds):
    """Factor m with sympy in ONE long-lived child process (cheap on Windows); kill and respawn on timeout."""
    import sympy  # raises ImportError if missing
    conn = _worker()
    conn.send(m)
    if conn.poll(seconds):
        return conn.recv()
    _W['proc'].terminate(); _W['proc'].join(); _W['proc'] = None
    return None
try:
    import sympy as _SYMPY
except ImportError:
    _SYMPY = None

def _rho_bounded(n, maxiter):
    """Pollard-Brent with a hard iteration budget. Returns a nontrivial factor or None."""
    from math import gcd
    for c in (1, 3, 5):
        y, r, q, g, it = 2, 1, 1, 1, 0
        x = ys = y
        while g == 1 and it < maxiter:
            x = y
            for _ in range(r): y = (y*y + c) % n
            k = 0
            while k < r and g == 1:
                ys = y
                for _ in range(min(128, r - k)):
                    y = (y*y + c) % n; q = q * abs(x - y) % n
                g = gcd(q, n); k += 128
            it += r; r *= 2
        if g == n:
            g = 1
            while g == 1:
                ys = (ys*ys + c) % n; g = gcd(abs(x - ys), n)
        if 1 < g < n: return g
    return None

def _load_known():
    import os
    for d in (os.path.dirname(os.path.abspath(__file__)), os.getcwd()):
        p = os.path.join(d, "known_factors.txt")
        if os.path.exists(p):
            for line in open(p):
                xs = list(map(int, line.split()))
                if len(xs) >= 2:
                    f = {}
                    for q in xs[1:]: f[q] = f.get(q, 0) + 1
                    _cache[xs[0]] = f
_load_known()
def factorize(n, big_timeout=30):
    """return dict prime->exp; raises RuntimeError if a cofactor can't be factored in time"""
    if n in _cache: return dict(_cache[n])
    orig = n; f = {}
    for p in (2,3,5,7,11,13,17,19,23,29,31,37,41,43,47):
        while n % p == 0: f[p] = f.get(p,0)+1; n //= p
    stack = [n] if n > 1 else []
    while stack:
        m = stack.pop()
        if m == 1: continue
        if is_probable_prime(m): f[m] = f.get(m,0)+1; continue
        if m in _cache:
            for q, k in _cache[m].items(): f[q] = f.get(q,0)+k
            continue
        if m < 10**40 and _SYMPY is not None:   # in-process: fast and no time limit needed at this size
            for q, k in _SYMPY.factorint(m).items(): f[q] = f.get(q,0)+k
            continue
        if m >= 10**24:  # strip medium primes so the cofactor matches what external tools were given
            q, found = 49, False
            while q < 100000 and q*q <= m:
                if m % q == 0:
                    while m % q == 0: f[q] = f.get(q,0)+1; m //= q
                    found = True
                q += 2
            if found: stack.append(m); continue
        if m < 10**24:
            d = pollard_brent(m); stack += [d, m//d]; continue
        d = _rho_bounded(m, 60000)   # cheap in-process attempt; finds factors up to ~10^15 quickly
        if d:
            stack += [d, m//d]; continue
        try:
            fm = _sympy_factor_with_timeout(m, big_timeout)
            if fm is None: raise RuntimeError(f"unfactored {len(str(m))}-digit cofactor")
            for q, k in fm.items(): f[q] = f.get(q,0)+k
            continue
        except ImportError:
            pass
        import shutil
        if shutil.which("factor"):
            try:
                out = subprocess.run(["factor", str(m)], capture_output=True, text=True, timeout=big_timeout).stdout
                for t in out.split(":")[1].split(): f[int(t)] = f.get(int(t),0)+1
                continue
            except subprocess.TimeoutExpired:
                raise RuntimeError(f"unfactored {len(str(m))}-digit cofactor")
        raise RuntimeError(f"unfactored {len(str(m))}-digit cofactor (install sympy)")
    _cache[orig] = dict(f)
    return f

def sigma_pp(p, e): return (p**(e+1)-1)//(p-1)

def sigma(n):
    s = 1
    for p, e in factorize(n).items(): s *= sigma_pp(p, e)
    return s

def v2(n):
    c = 0
    while n % 2 == 0: n //= 2; c += 1
    return c
