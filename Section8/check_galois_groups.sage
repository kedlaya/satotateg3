#
# This SageMath code confirms the assertions about Jacobians of curves, and their endomorphism fields,
# made in Section 8.2 and Section 8.3. The code uses Magma in order to compute some large Galois groups
# and some geometric automorphism groups of hyperelliptic curves.
#
# See also check_twisting_representations.sage for verification of the group presentations in Section 8.3.

def check_group_subgroup(u, id1, id2, D):
    # Given a polynomial over Q, check that its Galois group has GAP ID id1 and that it has a subgroup
    # with GAP ID id2 which is a quadratic field with discriminant D.
    L.<a> = u.splitting_field()
    G = L.galois_group()
    if (G.gap().IdGroup() != id1):
        return False
    l = G.gap().NormalSubgroups()
    for i in l:
        if i.IdGroup() == id2:
            L2 = G.subgroup(i.GeneratorsOfMagma()).fixed_field()[0]
            if L2.degree() == 2 and L2.disc() == D:
                return True
    return False
    
def ritzenthaler_romagny(f, g, h, c=1):
    ## Perform the Ritzenthaler-Romagny construction to compute the Prym of a 2:1 map from
    ## a genus-3 curve to a genus-1 curve.
    P = f.parent()
    f /= c
    g /= c
    A = Matrix([list(f)[::-1], list(g)[::-1], list(h)[::-1]])
    B = ~A.transpose()
    B.rescale_col(1, 2)
    (a1, b1, c1) = tuple(P(list(i)) for i in B.rows())
    return b1*(b1^2-a1*c1)

## Example 8.3 -- U(3)
G0 = magma.GU(3, 7)
l0 = [i for i in G0.Generators()]
l1 = [[[i[j,k] for k in range(1,4)] for j in range(1,4)] for i in l0]
F.<a> = GF(7^2)
l2 = [[[F(0) if i[j][k]==0 else a^(i[j][k].Log()) for k in range(3)] for j in range(3)] for i in l1]
l3 = [[[i[j][k].polynomial().list() for k in range(3)] for j in range(3)] for i in l2]
for i in l3:
    for j in range(3):
        for k in range(3):
            while len(i[j][k]) < 2:
                i[j][k].append(0)
t = [(a^i).matrix() for i in range(2)]
F0 = GF(7)
l4 = [block_matrix(F0, 3, 3, [i[j][k][0]*t[0]+i[j][k][1]*t[1] for j in range(3) for k in range(3)], subdivide=False) for i in l3]
G = GL(6, 7)
G1 = G.subgroup(l4)
assert G0.Order() == G1.order()
G2 = G1.gap().CentralizerInParent()
l5 = [i.matrix() for i in G2.GeneratorsOfMagma()]
G3 = magma.GL(6, 7).sub([magma(i) for i in l4+l5])
l = [H.get_magma_attribute('subgroup') for H in G3.MaximalSubgroups()]
P.<T> = F0[]
l2 = [[t[3].CharacteristicPolynomial().sage() for t in H.ConjugacyClasses()] for H in l]
l3 = [[t(T) for t in i] for i in l2]
Q.<x> = QQ[]
f = x^7 + 3*x^5 + 4*x^3 + 2*x
C = HyperellipticCurve(f)
primes = [3]
cp = [C.change_ring(GF(p)).frobenius_polynomial()(T) for p in primes]
assert not any(all(u in i for u in cp) for i in l3)
assert cp == [T^6 + 6*T^4 + 4*T^2 + 6]
print("Example 8.3 -- U(3) verified")

P.<x,y> = QQ[]
P1.<t> = QQ[]

## Example 8.4 -- SU(2) x USp(4)
C = Curve(y^2 - (4*x^4 - 7*x + 4))
J = Jacobian(C)
assert not J.has_cm()
print("Example 8.4 -- SU(2) x USp(4) verified")

## Example 8.4 -- N(U(1)) x USp(4)
C = Curve(y^2 - (x^4 - x^3 - 3*x^2 + x - 1))
J = Jacobian(C)
assert J.has_cm()
print("Example 8.4 -- N(U(1)) x USp(4) verified")

## Example 8.5 -- SU(2) x N(G_{3,3})
p = 23
u = t^3 - t^2 + t - 2
l = u.roots(GF(p))
l1 = [EllipticCurve([0,0,0,-1,i]).trace_of_frobenius() for (i, _) in l]
assert len(set(l1)) == 3
u = 2*t^7 + 4*t^6 - 7*t^4 + 4*t^3 - 4*t + 1
assert u == (t-1) * (2*t^3 + 1) * ((1+t)^3 - 2)
K.<a> = NumberField(t^3-2)
v = (-t+a^2)/(a*t+1)
assert v(v(t)) == t
assert u(v) * (a*t+1)^8 / 81 == u
C = Curve(y^2 - (x-a)*(x^2-x/a+1/a^2))
J = Jacobian(C)
assert not J.has_cm()
print("Example 8.5 -- SU(2) x N(G_{3,3}) verified")

## Example 8.6 -- N(U(1)) x N(G_{3,3})
C = Curve(y^2 - (x^4 + 2*x^3 + 4*x^2 + 4*x + 4))
J = Jacobian(C)
assert J.has_cm()
assert J.j_invariant() == 8000
u = t^5 + 2*t^4 + 4*t^3 + 4*t^2 + 4*t
K.<a> = QuadraticField(2)
u1 = u(a*(t+1)/(t-1))*(t-1)^6 / 16
assert u1 == (a+1)*t^6 + (a-1)*t^4 + (-a-1)*t^2 + (-a+1)
C = Curve(y^2 - ((a+1)*x^3 + (a-1)*x^2 + (-a-1)*x + (-a+1)))
J = Jacobian(C)
assert J.global_minimal_model() == EllipticCurve([0,-a+1,0,-a,0])
print("Example 8.6 -- N(U(1)) x N(G_{3,3}) verified")

## Example 8.7 -- SU(2) x F_{ac}
C = Curve(y^2 - (x^4 - 5*x^3 + 10*x^2 - 10*x + 5))
J = Jacobian(C)
assert not J.has_cm()
print("Example 8.7 -- SU(2) x F_{ac} verified")

## Example 8.8 -- SU(2) x F_{a,b}
C1 = Curve(2*y^2 - (x^3 + x^2 - 3*x + 1))
J1 = Jacobian(C1)
assert J1.has_cm()
assert J1.j_invariant() == 8000
C2 = Curve(2*y^2 - (x^3 - 3*x^2 + x + 1))
J2 = Jacobian(C2)
assert J2.has_cm()
assert J2.j_invariant() == 1728
print("Example 8.8 -- SU(2) x F_{a,b} verified")

## Example 8.9 -- N(U(1)) x F_{ac}
C = Curve(y^2 - (x^4 - 8*x^3 + 20*x^2 - 16*x + 2))
J = Jacobian(C)
assert J.has_cm()
assert J.j_invariant() == 8000
print("Example 8.9 -- N(U(1)) x F_{ac} verified")

## Example 8.9 -- N(U(1)) x N(U(1)) x N(U(1))
C1 = Curve(3*y^2 + (-6*x^2 + 6)*y + (2*x^4 - 12*x^2 + 6))
J1 = Jacobian(C1)
assert J1.has_cm()
j1 = J1.j_invariant()
C2 = Curve(2*y^2 + (-6*x^2 - 12)*y + (3*x^4 + 6*x^2 + 6))
J2 = Jacobian(C2)
assert J2.has_cm()
j2 = J2.j_invariant()
C3 = Curve(6*y^2 + (6*x^2 - 12)*y + (3*x^4 - 6*x^2 + 2))
J3 = Jacobian(C3)
assert J3.has_cm()
j3 = J3.j_invariant()
assert j1 != j2 and j1 != j2 and j2 != j3
print("Example 8.9 -- N(U(1)) x N(U(1)) x N(U(1)) verified")

## Example 8.10 -- SU(2) x J(E_4)
C = Curve(y^2 - (x^4 + x^2 + 2))
J = Jacobian(C)
assert not J.has_cm()
print("Example 8.10 -- SU(2) x J(E_4) verified")

## Example 8.10 -- N(U(1)) x J(E_4)
C = Curve(y^2 - (x^4 + 2*x^3 + 4*x - 4))
J = Jacobian(C)
assert J.has_cm()
assert J.j_invariant() == 1728
L.<a> = (t^4 - 6*t^2 + 3).splitting_field()
C = HyperellipticCurve(t^5 + 2*t^4 + 4*t^2 - 4*t)
G = magma(C).BaseChange(L).AutomorphismGroup()
assert list(G.IdentifyGroup()) == [8, 3]
print("Example 8.10 -- N(U(1)) x J(E_4) verified")

## Example 8.11 -- SU(2) x J(E_6)
C = Curve(y^2 - (x^4 + 8*x^3 + 18*x^2 + 16*x - 4))
J = Jacobian(C)
assert not J.has_cm()
u0 = t^3 - 2*t^2 + t - 16/27
L.<a> = NumberField(u0)
P2.<y> = L[]
C = HyperellipticCurve(y^5 + 8*y^4 + 18*y^3 + 16*y^2 - 4*y)
J = magma(C).Jacobian()
J1 = J.RichelotIsogenousSurfaces()[1]
C1 = J1.Curve()
C2 = HyperellipticCurve(a*(y^5 + 8*y^4 + 18*y^3 + 16*y^2 - 4*y))
assert magma.IsIsomorphic(magma(C2), C1) == magma(True)
primes = [5, 13, 17]
l = [C.change_ring(GF(p)).zeta_function().numerator() for p in primes]
assert all(i[1] != 0 for i in l)
assert all(any(legendre_symbol(-x, p) == -1 for p in primes) for x in divisors(6) if x != 1)
L1.<a1> = QuadraticField(-1)
assert l[0].change_ring(L1).is_irreducible()
print("Example 8.11 -- SU(2) x J(E_6) verified")

## Example 8.11 -- N(U(1)) x J(E_6)
C = Curve(2*y^2 + (4*x^2-6)*y + 4*x^4 + 6*x^3 + x + 3)
J = Jacobian(C)
assert J.has_cm()
assert J.j_invariant() == 1728
u = 4*t^4 + 6*t^3 + t + 3
L.<a> = NumberField(u//(4*(t+1)))
f = (t+1)*(t-a)
h = u.change_ring(L) // f
g = -4*t^2 + 6
v = ritzenthaler_romagny(f, g, h, 2)
C = HyperellipticCurve(v)
J = magma(C).Jacobian()
J1 = J.RichelotIsogenousSurfaces()[1]
C1 = J1.Curve().sage()
v1 = C1.hyperelliptic_polynomials()[0]
v1 = (56*v1/v1[6]).change_ring(QQ)
v1 = v1(t+2)
assert v1 == 56*t^6 + 192*t^5 + 132*t^4 - 256*t^3 - 390*t^2 - 48*t + 99
C1 = HyperellipticCurve(v1)
assert magma.IsIsomorphic(magma(C), magma(C1.change_ring(L))) == magma(False)
r = C.igusa_clebsch_invariants()[0] / C1.igusa_clebsch_invariants()[0]
assert r.is_square()
primes = [5, 11, 31, 53]
l = [C1.change_ring(GF(p)).zeta_function().numerator() for p in primes]
assert all(i[1] != 0 for i in l)
assert all(any(legendre_symbol(-x, p) == -1 for p in primes) for x in divisors(42) if x != 6)
L1.<a1> = QuadraticField(-6)
assert l[0].change_ring(L1).is_irreducible()
u1 = t^6 + 7*t^4 + 7*t^2 - 63
L1.<a1> = u1.splitting_field()
u0 = (t^2+1) * u1
L.<a> = u0.splitting_field()
u1 = t^12 + 10*t^10 + 177*t^8 + 244*t^6 + 172*t^4 + 96*t^2 + 36
L1.<a1> = u1.splitting_field()
assert len(u.roots(L1)) == u.degree() and len(u1.roots(L)) == u1.degree()
L2.<a2> = v1.splitting_field()
a3 = (u//(4*(t+1))).roots(L2)[0][0]
f = (t+1)*(t-a3)
h = u.change_ring(L2) // f
g = -4*t^2 + 6
v = ritzenthaler_romagny(f, g, h, 2)
assert len(v.roots(L2)) == v.degree()
print("Example 8.11 -- N(U(1)) x J(E_6) verified")

## Example 8.13 - SU(2) x J(D_6)
C = Curve(24*y^2 + (-8*x^2 + 24*x)*y + x^4 - 12*x - 9)
J = Jacobian(C)
assert not J.has_cm()
P1.<t> = QQ[]
u0 = t^6 - 3*t^2 - 4
L.<a> = NumberField(u0)
u1 = (t^4 - 12*t - 9).change_ring(L)
f = u1.factor()[0][0]
h = u1.factor()[1][0]
g = (8*t^2 - 24*t)
u = ritzenthaler_romagny(f, g, h, 24)
r1 = u.roots()[0][0]
r2 = u.roots()[1][0]
u = (u((r2*t-r1)/(t-1)) * (t-1)^6).numerator()
u = u(t*u[4]/u[5]/8)               
u /= u[5]
assert (u == t^5 + 8*t^4 + 16*t^2 - 4*t)
print("Example 8.13 -- SU(2) x J(D_6) verified")

P.<x> = QQ[]

## Example 8.14 -- SU(2) x J(O)
u0 = x^3 + 2*x + 1
L3.<alpha> = NumberField(u0)
a = -4-3*alpha^2
b = 2*alpha^2-3*alpha + 4
F1.<r1> = L3.extension(x^2 - b)
F1.<c> = F1.absolute_field()
u = F1.defining_polynomial()
assert check_group_subgroup(u, [48, 48], [24, 12], -4)
print("Example 8.14 -- SU(2) x J(O) verified")

## Example 8.16 -- M[S_4]
u0 = x^3 - x^2 + x - 2
L3.<a> = u0.splitting_field()
(a1, a2, a3) = tuple(i[0] for i in u0.roots(L3))
F.<b> = L3.extension(x^2 - a1/a2)
F1.<c> = F.absolute_field()
u = F1.defining_polynomial()
G = u.galois_group(algorithm='magma')
assert G.gap().IdGroup() == [24, 12]
print("Example 8.16 -- M[S_4] verified")

## Example 8.17 -- J_s(B(1,12))
u = x^8 - 12*x^6 + 12*x^4 - 52*x^2 - 3
assert check_group_subgroup(u, [16, 7], [8, 3], -3)
u *= x^3 - 2
assert check_group_subgroup(u, [48, 15], [24, 10], -3)
print("Example 8.17 -- J_s(B(1,12)) verified")

## Example 8.19 -- J(B(3,4;4))
u = x^8 - 7*x^4 + 28
assert check_group_subgroup(u, [16, 7], [8, 1], -4)
assert check_group_subgroup(u, [16, 7], [8, 3], -7)
u *= x^3 - x^2 + 2*x - 3
assert check_group_subgroup(u, [48, 15], [24, 1], -4)
print("Example 8.19 -- J(B(3,4;4)) verified")

## Example 8.20 -- J(B(3,4))
#u = cyclotomic_polynomial(24)(x)*(x^4-3)*(x^3+3*x+2)
u = x^12 - 3*x^8 + 24*x^4 - 48
assert check_group_subgroup(u, [48, 38], [24, 5], -4)
print("Example 8.20 -- J(B(3,4)) verified")

## Example 8.21 -- J_s(B(3,4))
u = x^8 - 8*x^4 + 25
assert check_group_subgroup(u, [16, 13], [8, 2], -4)
assert check_group_subgroup(u, [16, 13], [8, 4], -40)
u *= x^3 - x^2 + 5*x + 15
assert check_group_subgroup(u, [48, 41], [24, 5], -4)
print("Example 8.21 -- J_s(B(3,4)) verified")

## Example 8.23 -- J(B(O,1))
u = x^16 - 44*x^12 - 308*x^10 - 990*x^8 - 1936*x^6 - 2662*x^4 + 9196*x^2 + 20449
G = u.galois_group(algorithm='magma')
assert G.gap().IdGroup() == [96, 193]
L = magma(u).SplittingField()
u1 = P(seq(L.DefiningPolynomial().Coefficients()))
L1.<a1> = NumberField(u1)
u2 = x^6 + 4*x^4 - 22
assert len(u2.roots(L1)) == u2.degree()
print("Example 8.23 -- J(B(O,1)) verified")

## Example 8.24 -- J_s(B(T,3))
u = x^8 - 6*x^4 + 4*x^2 - 3
assert check_group_subgroup(u, [48, 29], [24, 3], -3)
L0.<a> = NumberField(u)
L1.<b> = L0.extension(x^3 - 2)
L2.<c> = L1.absolute_field()
u = L2.defining_polynomial()
G = u.galois_group(algorithm='magma')
assert G.gap().IdGroup() == [144, 125]
print("Example 8.24 -- J_s(B(T,3)) verified")

## Example 8.27 -- J(B(T,3))
u = P([1,0,0,0,-18,0,101,0,99,0,-1098,0,3043,0,-738,0,81])
assert check_group_subgroup(u, [48, 33], [24, 3], -3)
L0.<a> = NumberField(u)
L1.<b> = L0.extension(x^3 - 2)
L2.<c> = L1.absolute_field()
u = L2.defining_polynomial()
## Since this polynomial of degree 48, which is too big for the Gap classification of 
## transitive groups, we need special syntax to retrieve it as a permutation group.
G = magma(u).GaloisGroup()
l = list(G.Generators())
for i in range(len(l)):
    l[i] = l[i].CycleDecomposition()
    l[i] = tuple(tuple(j) for j in l[i])
G = PermutationGroup(l)
assert G.gap().IdGroup() == [144, 127]
print("Example 8.27 -- J(B(T,3)) verified")

## Example 8.28 -- J(B(O,2))
u0 = x^3 + 2*x + 1
L3.<alpha> = NumberField(u0)
a = -4-3*alpha^2
b = 2*alpha^2-3*alpha + 4
F1.<r1> = L3.extension(x^2 - b)
F2.<r2> = F1.extension(x^4 - (a - 2*r1))
F1.<c> = F2.absolute_field()
u = F1.defining_polynomial()
G = u.galois_group(algorithm='magma')
assert G.gap().IdGroup() == [192, 988]
print("Example 8.28 -- J(B(O,2)) verified")

## Example 8.30 -- J(D(4,4))
u0 = x^3 - x^2 + x - 2
L3.<a> = u0.splitting_field()
(a1, a2, a3) = tuple(i[0] for i in u0.roots(L3))
F.<b> = L3.extension(x^4 - a1/a2)
F1.<c> = F.absolute_field()
u = F1.defining_polynomial()
G = u.galois_group(algorithm='magma')
assert G.gap().IdGroup() == [192, 956]
print("Example 8.30 -- J(D(4,4)) verified")

## Example 8.32 -- J(D(6,6))
u0 = x^3 - x^2 + x - 2
L3.<a> = u0.splitting_field()
(a1, a2, a3) = tuple(i[0] for i in u0.roots(L3))
F.<b> = L3.extension(x^6 - a1/a2)
F1.<c> = F.absolute_field()
u = F1.defining_polynomial()
G = u.galois_group(algorithm='magma')
assert G.gap().IdGroup() == [432, 523]
print("Example 8.32 -- J(D(6,6)) verified")

## Example 8.34 -- J(E(216))
P0.<a,b,k0,k1> = PolynomialRing(QQ, 4)
P1.<x> = PolynomialRing(P0)
u = (a*x+b)^3 - (x^4 + x+1) + (x^2 + k1*x + k0)^2
I = P0.ideal(u.coefficients())
f = I.elimination_ideal([b,k0,k1]).gens()[0]
f = P(f.univariate_polynomial())
g = P(list(f)[::3])
assert g(x^3) == f
G = g.galois_group()
assert G.gap().IdGroup() == [432, 734]
K.<z> = NumberField(g)
assert not z.is_nth_power(3)
G = f.galois_group(algorithm='magma')
assert G.gap().IdGroup() == [1296, 2891]
print("Example 8.34 -- J(E(216)) verified")

print("All checks completed")

