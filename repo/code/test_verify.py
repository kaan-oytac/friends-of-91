"""Adversarial tests for verify.py: every corrupted certificate below must be rejected.
Usage: python test_verify.py [certificate.json] [X]    (default: cert_91_1e100.json 100)
Each mutant is written to a temporary file and verify.py is run on it as a subprocess;
the valid certificate must be accepted and every mutant must be rejected (non-zero exit)."""
import sys, json, copy, subprocess, tempfile, os

cert_path = sys.argv[1] if len(sys.argv) > 1 else 'cert_91_1e100.json'
X = sys.argv[2] if len(sys.argv) > 2 else '100'
here = os.path.dirname(os.path.abspath(__file__))
base = json.load(open(cert_path))

def run(cert, X=X):
    with tempfile.NamedTemporaryFile('w', suffix='.json', delete=False) as f:
        json.dump(cert, f); name = f.name
    try:
        r = subprocess.run([sys.executable, os.path.join(here, 'verify.py'), name, X], capture_output=True, text=True)
    finally:
        os.unlink(name)
    return r.returncode, (r.stdout + r.stderr).strip().splitlines()

def nodes(c): return c['branches'][0]['nodes']
def find(c, kind): return next(nd for nd in nodes(c) if nd['kind'] == kind)

mutants = {}
# 1. the whole tree replaced by one genuine unfactored leaf (the defect found in the audit of v1)
m = copy.deepcopy(base); m['branches'][0]['nodes'] = [copy.deepcopy(find(base, 'unfactored'))]; mutants['single leaf, no root'] = m
# 2. root removed, everything else kept
m = copy.deepcopy(base); m['branches'][0]['nodes'] = [nd for nd in nodes(m) if nd['path'] != []]; mutants['root removed'] = m
# 3. an unfactored leaf whose parent does not exist (orphan)
m = copy.deepcopy(base); nd = copy.deepcopy(find(base, 'unfactored')); nd['path'] = nd['path'] + [[1000003, 2]]; nodes(m).append(nd); mutants['orphan leaf'] = m
# 4. a child deleted (tree incomplete)
m = copy.deepcopy(base); ch = next(nd for nd in nodes(m) if nd['kind'] == 'internal' and nd['path']); nodes(m).remove(ch); mutants['child deleted'] = m
# 5. a child's exponent relabelled in its parent as skipped
m = copy.deepcopy(base); ch = next(nd for nd in nodes(m) if nd['kind'] == 'internal' and nd['path'])
par = next(nd for nd in nodes(m) if nd['path'] == ch['path'][:-1]); e = ch['path'][-1][1]
for it in par['exponents']:
    if int(it[0]) == e: it[1] = 'skip_parity'
mutants['exponent mislabelled'] = m
# 6. an internal node falsely declared dead
m = copy.deepcopy(base); nd = next(nd for nd in nodes(m) if nd['kind'] == 'internal'); nd['kind'] = 'dead'; nd['reason'] = 'below 1'; mutants['fake dead node'] = m
# 7. residual tampered at the root
m = copy.deepcopy(base); r = find(m, 'internal'); r['r'] = '17/13'; mutants['bad residual'] = m
# 8. a factor of a denominator replaced by a composite with the same product (p*q -> p*q written as one composite)
m = copy.deepcopy(base)
for nd in nodes(m):
    if nd.get('den_factors') and len(nd['den_factors']) >= 2:
        (p, a), (q, b) = nd['den_factors'][0], nd['den_factors'][1]
        if a == 1 and b == 1:
            nd['den_factors'] = [[p*q, 1]] + nd['den_factors'][2:]; break
mutants['composite in factorization'] = m
# 9. bound overclaimed by one order of magnitude
mutants['overclaimed bound'] = ('overclaim', copy.deepcopy(base))
# 10. a cap leaf removed
m = copy.deepcopy(base); nodes(m).remove(find(m, 'cap')); mutants['cap leaf deleted'] = m
# 11. unfactored leaf whose residual is edited to a smaller denominator (would weaken the strengthened bound)
m = copy.deepcopy(base); nd = find(m, 'unfactored'); nd['r'] = str(int(nd['r'].split('/')[0])) + '/3'; mutants['unfactored residual edited'] = m

ok = True
code, out = run(base)
print(('PASS' if code == 0 else 'FAIL'), 'valid certificate accepted:', out[-1] if out else '')
ok &= code == 0
for name, mut in mutants.items():
    if isinstance(mut, tuple):
        code, out = run(mut[1], str(int(float(X)) + 1))
    else:
        code, out = run(mut)
    rejected = code != 0
    print(('PASS' if rejected else 'FAIL'), f'{name}:', out[1].strip() if len(out) > 1 else (out[0] if out else ''))
    ok &= rejected
print('ALL TESTS PASSED' if ok else 'SOME TESTS FAILED')
sys.exit(0 if ok else 1)
