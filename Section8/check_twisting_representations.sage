#
# This SageMath code confirms that the twisting constructions in Section 8.3
# give the specified groups and index-2 subgroups.
#
# See also check_galois_groups.sage for verification of the Galois group assertions in Section 8.3.

def check_twisting_rep(G, h=None, swap=None):
    # For G a matrix group defined over an imaginary quadratic field,
    # check that the action of complex conjugation followed by conjugation by h preserves G,
    # and compute the class of the resulting extension group.
    #
    # If swap is specified, we take h to be the 2x2 matrix [[0, swap], [1, 0]].

    if swap:
        h1 = Matrix([[0,swap],[1,0]]) 
    else:
        h1 = identity_matrix(G.degree())
    if h is not None:
        h1 *= h 
    assert all(~h1*g.matrix().conjugate()*h1 in G for g in G.gens())
    # Construct a permutation representation of the semidirect product.
    Glist = [g.matrix() for g in G]
    n = len(Glist)
    l = [[Glist.index((g*g1).matrix()) for g1 in G] for g in G.gens()]
    l.append([Glist.index(~h1*g.matrix().conjugate()*h1) for g in G])
    G1 = PermutationGroup(l, domain=list(range(n)))
    return (G.gap().IdGroup(), G1.gap().IdGroup())

K3.<z3> = CyclotomicField(3)
K4.<z4> = CyclotomicField(4)
K7.<z7> = CyclotomicField(7)
K8.<z8> = CyclotomicField(8)
s2 = z8^3 + z8
s7 = 1 + 2*(z7 + z7^2 + z7^4)

## Example 8.17 -- J_s(B(1, 12))
g1 = z3*identity_matrix(2)
g2 = Matrix([[-1,-2],[0,1]])
g3 = Matrix([[-1,0],[1,1]])
G = MatrixGroup([g1,g2,g3])
assert check_twisting_rep(G, swap=2) == ([24, 10], [48, 15])
print("Example 8.17 confirmed")

## Example 8.19 -- J(B(3, 4; 4))
g1 = z4*identity_matrix(2)
g2 = Matrix([[-z4-1,-3],[z4,z4+1]])
g3 = Matrix([[1,3],[-1,-2]])
G = MatrixGroup([g1,g2,g3])
assert check_twisting_rep(G, swap=3) == ([24, 1], [48, 15])
print("Example 8.18 confirmed")

## Example 8.20 -- J(B(3, 4))
g1 = z4*identity_matrix(2)
g2 = Matrix(K4, [[0,1],[1,0]])
g3 = Matrix(K4, [[0,-1],[1,-1]])
G = MatrixGroup([g1,g2,g3])
assert check_twisting_rep(G) == ([24, 5], [48, 38])
print("Example 8.20 confirmed")

## Example 8.21 -- J_s(B(3, 4))
g1 = z4*identity_matrix(2)
g2 = Matrix(K4, [[1,0],[-1,-1]])
g3 = Matrix(K4, [[1,3],[-1,-2]])
G = MatrixGroup([g1,g2,g3])
assert check_twisting_rep(G, swap=3) == ([24, 5], [48, 41])
print("Example 8.21 confirmed")

## Example 8.23 -- J(B(O, 1))
g1 = Matrix([[-1, 0], [1-s2, 1]])
g2 = Matrix([[-1, -1], [1, 0]])
g2 = Matrix([[-s2, 1], [1+s2, -1+s2]])
G = MatrixGroup([g1, g2])
assert check_twisting_rep(G, swap=1) == ([48, 29], [96, 193])
print("Example 8.22 confirmed")

## Example 8.24 -- J_s(B(T, 3))
g1 = z3*identity_matrix(2)
g2 = Matrix([[0,-1],[1,0]])
g3 = Matrix([[1+z3,-1],[0,-z3]])
G = MatrixGroup([g1,g2,g3])
assert check_twisting_rep(G, swap=1) == ([72, 25], [144, 125])
print("Example 8.23 confirmed")

## Example 8.27 -- J(B(T, 3))
g1 = z3*identity_matrix(2)
g2 = Matrix([[-1,-2],[1,1]])
g3 = Matrix([[1+z3,2*z3],[0,-z3]])
G = MatrixGroup([g1,g2,g3])
assert check_twisting_rep(G, swap=2) == ([72, 25], [144, 127])
print("Example 8.27 confirmed")

## Example 8.28 -- J(B(O, 2))
g1 = Matrix([[z4,0],[-z4,1]])
g2 = Matrix([[1,1],[0,z4]])
G = MatrixGroup([g1, g2])
assert check_twisting_rep(G, swap=1) == ([96, 67], [192, 988])
print("Example 8.28 confirmed")

## Example 8.31 -- J(E(168))
g1 = Matrix([[-1,0,(1-s7)/2],[0,-1,(-1-s7)/2],[0,0,1]])
g2 = Matrix([[2,(-1+s7)/2,(-3+s7)/2],[1,0,0],[(1+s7)/2,-1,-1]])
G = MatrixGroup([g1, g2])
h = Matrix([[0,1,0],[1,0,0],[0,0,-1]])
assert check_twisting_rep(G, h) == ([168, 42], [336, 208])
print("Example 8.31 confirmed")

## Example 8.34 -- J(E(216))
g1 = Matrix([[1,z3^2,0],[0,z3^2,0],[0,0,1]])
g2 = Matrix([[z3,0,0],[z3,1,0],[z3,0,1]])
g3 = Matrix([[1,0,z3^2],[0,1,0],[0,0,z3^2]])
G = MatrixGroup([g1, g2, g3])
h = Matrix([[1,0,0],[1,-1,0],[1,0,-1]])
assert check_twisting_rep(G, h) == ([648, 533], [1296, 2891])
print("Example 8.34 confirmed")

print("All checks completed")

