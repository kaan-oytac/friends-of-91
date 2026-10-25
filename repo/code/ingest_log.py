"""python ingest_log.py factor.log cert.json X  -> adds verified factorizations of needed cofactors to known_factors.txt"""
import re, math, json, sys
from arith import is_probable_prime
log = open(sys.argv[1]).read(); c = json.load(open(sys.argv[2])); X = float(sys.argv[3])
need = [nd['cofactor'] for nd in c['branches'][0]['nodes'] if nd['kind']=='unfactored' and sum(e*math.log10(p) for p,e in nd['path']) < X]
blocks = {}
for b in re.split(r'Starting factorization of ', log)[1:]:
    blocks.setdefault(int(b.split()[0]), []).append(b)
def primes_in(b): return [int(x) for x in re.findall(r'[pPcC]\d+ = (\d+)', b) if is_probable_prime(int(x))]
lines = []
for s in need:
    n = int(s)
    if n not in blocks: print("MISSING", len(s), "digits"); continue
    fs = []
    for b in blocks[n]: fs += primes_in(b)
    # composites YAFU split later appear as their own blocks; pull their primes too
    for b in blocks[n]:
        for cc in re.findall(r'[cC]\d+ = (\d+)', b):
            cc = int(cc)
            if not is_probable_prime(cc) and cc in blocks:
                for bb in blocks[cc]: fs += primes_in(bb)
    # assemble with multiplicity by repeated division
    rem, out = n, []
    for f in sorted(set(fs)):
        while rem % f == 0: rem //= f; out.append(f)
    if rem != 1 and is_probable_prime(rem): out.append(rem); rem = 1
    ok = rem == 1
    print(len(s), "digits:", "OK" if ok else f"PROBLEM remaining {len(str(rem))} digits", [len(str(f)) for f in out])
    if ok: lines.append(s + " " + " ".join(map(str, out)))
open('known_factors.txt','a').write("\n".join(lines) + ("\n" if lines else ""))
print(len(lines), "added of", len(need))
