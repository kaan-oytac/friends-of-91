"""Machine-checked necessary conditions for a friend m of 91 (I = 16/13, i.e. 13*sigma(m) = 16*m)."""
from fractions import Fraction
from arith import *
T = Fraction(16, 13)
def I(p, e): return Fraction(sigma_pp(p, e), p**e)
print("L1  m odd: k=1 gives odd-part target", Fraction(16*2, 13*3), "< 1 ->", Fraction(32,39) < 1)
print("L2  13 | m. For odd e, 14 | sigma(13^e) so 7 | m and I(7)*I(13^e) >= 16/13 with equality only at e=1:",
      all(sigma_pp(13, e) % 14 == 0 and Fraction(8,7)*I(13, e) > T for e in range(3, 60, 2)), "; e=1 forces m=91 (residual 8/7, 7^1 only; 7^2 gives", I(7,2) > Fraction(8,7), ")")
print("L3  so m != 91 has 13^(2a) exactly, a>=1. Then 3,5,7 cannot divide m:")
for q in (3, 5, 7):
    print(f"      I({q})*I(13^2) = {float(Fraction(q+1,q)*I(13,2)):.4f} > 16/13 = {float(T):.4f}")
print("L4  v2(sigma(m)) = 4 exactly (16m = 13 sigma(m), m odd). Odd-exponent primes p^e contribute v2(p+1)+v2(e+1)-1 >= 1 each, so at most 4 of them, sum exactly 4.")
print("L5  3,5,7 do not divide sigma(m) either (they would have to divide m). Hence for every prime power p^e || m:")
print("      3: p=1 mod 3 needs 3 !| e+1 ; p=2 mod 3 needs e even")
print("      5: p=1 mod 5 needs 5 !| e+1 ; p=4 mod 5 needs e even ; p=2,3 mod 5 needs 4 !| e+1")
print("      7: p=1 mod 7 needs 7 !| e+1 ; p=6 mod 7 needs e even ; p=2,4 mod 7 needs 3 !| e+1 ; p=3,5 mod 7 needs 6 !| e+1")
# verify L5 numerically
ok = True
for p in [q for q in range(11, 2000) if is_probable_prime(q)]:
    for e in range(1, 25):
        s = sigma_pp(p, e)
        pred3 = (p % 3 == 1 and (e+1) % 3 == 0) or (p % 3 == 2 and e % 2 == 1)
        pred5 = (p % 5 == 1 and (e+1) % 5 == 0) or (p % 5 == 4 and e % 2 == 1) or (p % 5 in (2,3) and (e+1) % 4 == 0)
        pred7 = (p % 7 == 1 and (e+1) % 7 == 0) or (p % 7 == 6 and e % 2 == 1) or (p % 7 in (2,4) and (e+1) % 3 == 0) or (p % 7 in (3,5) and (e+1) % 6 == 0)
        if (s % 3 == 0) != pred3 or (s % 5 == 0) != pred5 or (s % 7 == 0) != pred7: ok = False
print("      numerically verified for p<2000, e<25:", ok)
print("L6  Consequence: an odd-exponent prime p must satisfy p = 1 mod 3, p != 4 mod 5, p != 6 mod 7.")
print("L7  13-adic: v13(sigma(m)) = 2a-1 >= 1, and sigma(13^2a) is coprime to 13, so other prime powers must supply 13^(2a-1).")
print("      p = 1 mod 13: v13(sigma(p^e)) = v13(e+1).  p = -1 mod 13 with e odd: v13 = v13(p+1)+v13((e+1)/2).  ord_13(p)=d>2: needs d | e+1.")
print("L8  smallest prime of m other than 13 is >= 11. Exponent of 11 (if present): e even, 5 !| e+1, 3 !| e+1 (11 = 2 mod 3, 1 mod 5, ord_7(11)=3):")
for e in range(1, 40):
    if e % 2 == 0 and (e+1) % 5 and (e+1) % 3:
        s = sigma_pp(11, e); f = factorize(s)
        print(f"      e={e}: sigma(11^e) = {s} = {f}", "-> dead (needs 3,5 or 7)" if any(q in f for q in (3,5,7)) else "")
print("L9  every prime q != 13 of m has v_q(sigma(m)) = v_q(m); every prime of sigma(m) other than 2 and 13 divides m.")
