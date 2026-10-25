"""Overnight friend hunt.   python hunt.py <n> <hours> [cap_exponent]
Runs the certificate engine (engine3) on n with a time budget. Any 'solution' node other than
n itself is a friend and is printed loudly. Appends a summary line to hunt_<n>.log.
The exponent cap (default 1e60) bounds how far each prime power is explored; the run may finish
early if the whole tree below the cap is exhausted, in which case raise the cap next time."""
import sys, time, json, math, multiprocessing
multiprocessing.set_start_method("spawn", force=True)
import engine3
from engine3 import run

if __name__ == '__main__':

    n = int(sys.argv[1]); hours = float(sys.argv[2]); cap = int(sys.argv[3]) if len(sys.argv) > 3 else 60
    t0 = time.time()
    cert = run(n, cap, hours * 3600)
    secs = time.time() - t0
    friends, kinds, lb = [], {}, None
    for br in cert['branches']:
        for nd in br['nodes']:
            kinds[nd['kind']] = kinds.get(nd['kind'], 0) + 1
            if nd['kind'] == 'solution' and int(nd['m']) != n: friends.append(nd['m'])
            if nd['kind'] in ('cap', 'unfactored', 'free'):
                lg = sum(e*math.log10(p) for p, e in nd['path']) + br['k']*math.log10(2)
                lb = lg if lb is None else min(lb, lg)
    complete = all(br['complete'] for br in cert['branches'])
    line = (f"n={n} cap=1e{cap} ran {secs/3600:.2f} h, tree {'exhausted' if complete else 'NOT exhausted (timed out)'}, "
            f"nodes {kinds}, " + (f"FRIEND FOUND: {friends}" if friends else "no friend found") +
            (f", every surviving leaf forces m > 10^{lb:.1f}" if lb is not None and complete else ""))
    print(line)
    open(f"hunt_{n}.log", "a").write(line + "\n")
    json.dump(cert, open(f"hunt_cert_{n}_1e{cap}.json", "w"))
    if friends:
        from fractions import Fraction
        print("\n*** FRIEND FOUND ***", friends)
        for m in friends: print(f"   check: sigma({m})/{m} = {Fraction(engine3.sigma(int(m)), int(m))}  target {cert['T']}")
