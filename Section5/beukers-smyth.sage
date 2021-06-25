#
# This SageMath code verifies Proposition 5.7 using the Beukers-Smyth algorithm,
# as described in Remark 5.8 and Remark 5.40.
#

R.<x,y> = PolynomialRing(Rationals())
Q.<w> = PolynomialRing(Rationals())
S = [(-x,-y),(x^2,y^2),(-x^2,-y^2)]
sols = []
for n in [0] + [2..9]:
    f = (x^2*y+x*y^2+1)*(x+y+x^2*y^2) - n*x^2*y^2
    res = f.resultant(prod(f(i) for i in S))(0,w); assert(res != 0)
    res = res.cyclotomic_part()
    for (h1,null) in res.factor():
        K1.<a> = NumberField(h1); d = a.multiplicative_order()
        if d == Infinity: 
            continue
        QK.<t> = PolynomialRing(K1)
        for (h2,null) in f(a,t).factor():
            if h2.degree() == 0: continue
            K2.<b> = K1.extension(h2); e = b.multiplicative_order()
            if e == Infinity or e < d: continue
            if e%d == 0: c = e; u = b
            else:
                c = lcm(d,e)
                u = K2.unit_group().gen()
                u = K2(u)^(K2.number_of_roots_of_unity()//c)
            for i in range(c):
                if u^i == a: m1 = i
                if u^i == b: m2 = i
            m3 = -m1-m2
            l = sorted([sorted([i*m1%c, i*m2%c, i*m3%c]) \
                    for i in range(c) if gcd(i,c)==1])
            if not (c,l[0]) in sols: sols.append((c,l[0]))
ans = [(j[0]/i, j[1]/i, j[2]/i) for (i,j) in sorted(sols)]
print(ans)
