#
# This SageMath code builds the finite subgroups of SU(3).C2 described in Section 5.
# It is a dependency for some of the other code in this directory.
#

import re

N = 504; K.<z> = CyclotomicField(N)
(z3,z4,z6,z7,z8,z9) = (z^(N/3), z^(N/4), z^(N/6), z^(N/7), z^(N/8), z^(N/9))

def f(M):
     s = M.trace(); s *= conjugate(s) 
     return s.is_integer()
     
def verify_rationality(G):
    return all(f(g.matrix()) for g in G)

def embed(M):
     return matrix.block(2, 2, [M, 0, 0, M.conjugate()])

def enlarge_group(G, l):
     return MatrixGroup([g.matrix() for g in G.gens()] + l)

def strip_name(name):
     m = re.match(r"^J(_n|_s|)\((.*)\)$", name)
     if m: return(m.group(2))
     return name
     
small_list = []
A_data = {(1,1): [(1/3,1/3,1/3)],
          (1,2): [(1/3,5/6,5/6)],
          (1,3): [(1/9,4/9,4/9)],
          (1,4,1): [(1/6,5/12,5/12)],
          (1,4,2): [(1/3,1/12,7/12)],
          (1,6,1): [(1/9,17/18,17/18)],
          (1,6,2): [(1/18,2/9,13/18)],
          (1,7): [(1/21,16/21,4/21)],
          (1,8,1): [(1/6,1/24,19/24)],
          (1,8,2): [(1/12,5/24,17/24)],
          (1,12): [(1/9,7/36, 25/36)],
          (2,2): [(0,1/2,1/2), (1/6,1/6,2/3)],
          (2,4): [(1/2,0,1/2), (1/6,5/12,5/12)],
          (2,6): [(1/2,0,1/2), (1/9,17/18,17/18)],
          (3,1): [(1/3,0,2/3), (1/3,1/3,1/3)],
          (3,2): [(1/3,0,2/3), (2/3,1/6,1/6)],
          (3,3): [(0,1/3,2/3), (1/9,1/9,7/9)],
          (3,4): [(0,1/3,2/3), (1/6,5/12,5/12)],
          (3,6): [(0,1/3,2/3), (1/9,5/18,11/18)],
          (4,4): [(0,1/4,3/4), (1/12,1/12,5/6)],
          (6,2): [(1/6,0,5/6), (2/3,1/6,1/6)],
          (6,6): [(0,1/6,5/6), (1/18,1/18,8/9)]}
A_groups = {}
for i in sorted(A_data.keys()):
     G = MatrixGroup([diagonal_matrix([z^(N*j) for j in m]) for m in A_data[i]])
     if len(i) == 2:
          A_str = "A(%d, %d)" % i
     else:
          A_str = "A(%d, %d)_%d" % i
     A_groups[i] = G
     small_list.append((A_str, G))
     
T_mat = {i: Matrix([[-z^(-N/i),0,0],[0,0,z^(N/(2*i))],[0,z^(N/(2*i)),0]])
     for i in (1,2,4)}
B_keys = [((1,4,2), 1), ((1,8,1), 1), ((1,12), 1),
          ((2,4), 1), ((3,1), 1), ((3,2), 1), ((3,3), 1), ((3,4), 1),
          ((3,6), 1), ((4,4), 1), ((6,2), 1), ((6,6), 1),
          ((1,4,2), 2), ((1,12), 2), ((3,2), 2), ((3,6), 2),
          ((2,4), 4), ((3,4), 4)]
for (i,k) in B_keys:
     G = enlarge_group(A_groups[i], [T_mat[k]])
     if k==1 and len(i) == 2:
          B_str = "B(%d, %d)" % i
     elif k == 1:
          B_str = "B(%d, %d)_%d" % i
     elif len(i) == 2:
          B_str = "B(%d, %d; %d)" % (i[0], i[1], k)
     else:
          B_str = "B(%d, %d; %d)_%d" % (i[0], i[1], k, i[2])
     small_list.append((B_str, G))

BT_gens = [Matrix([[1,0,0],[0,0,1],[0,-1,0]]),
           Matrix([[1,0,0],[0,1/2+z4/2,1/2+z4/2],[0,-1/2+z4/2,1/2-z4/2]])]
for n in [1,2,3]:
    u = z^(N//(6*n))
    d = diagonal_matrix([u^(-2),u,u])
    G = MatrixGroup(BT_gens + [d])
    small_list.append(("B(T, %d)" %n, G))
    if n==3:
         G = MatrixGroup([BT_gens[0], BT_gens[1]*d])
         small_list.append(("B(T, 1; 1)", G))
for n in [1,2]:
    u = z^(N//(12*n))
    G = MatrixGroup(BT_gens + [diagonal_matrix([u^2,~u*z8,~u/z8])])
    small_list.append(("B(O, %d)" %n, G))

C_keys = [(1,7), (2,2), (3,1), (3,3), (4,4), (6,2), (6,6)]
S_mat = Matrix([[0,0,1],[1,0,0],[0,1,0]])
for i in C_keys:
     G = enlarge_group(A_groups[i], [S_mat])
     small_list.append(("C(%d, %d)" %i, G))

D_keys = [(2,2), (3,1), (3,3), (4,4), (6,2), (6,6)]
for i in D_keys:
     G = enlarge_group(A_groups[i], [S_mat, T_mat[1]])
     small_list.append(("D(%d, %d)" %i, G))

g1 = diagonal_matrix([1, z3, z3^2])
g2 = ~(z3-z3^2)*Matrix([[1,1,1],[1,z3,z3^2],[1,z3^2,z3]])
g3 = diagonal_matrix([z9^2, z9^2, z9^5])
small_list.append(("E(36)", MatrixGroup([g1, g2])))
small_list.append(("E(72)", MatrixGroup([g1, g2, g3*g2*~g3])))
small_list.append(("E(216)", MatrixGroup([g1, g2, g3])))

h = (z7+z7^2+z7^4-z7^3-z7^5-z7^6)/7 #equals sqrt(-7)
(a,b,c) = (z7^4-z7^3,z7^2-z7^5,z7-z7^6)
small_list.append(("E(168)", MatrixGroup([S_mat, diagonal_matrix([z7,z7^2,z7^4]),
                                          h*Matrix([[a,b,c],[b,c,a],[c,a,b]]),
                                          diagonal_matrix([z3,z3,z3])])))

verify = False

big_list1 = []
for (name, G) in small_list:
     if verify: assert verify_rationality(G)
     G2 = MatrixGroup([embed(-identity_matrix(3))] + [embed(g.matrix()) for g in G.gens()])
     big_list1.append((name, G2))

prime_data_A = ['A(1, 4)_2', 'A(1, 8)_1', 'A(1, 8)_2', 'A(1, 12)',
               'A(2, 2)', 'A(2, 4)', 'A(2, 6)', 'A(3, 1)', 'A(3, 2)',
               'A(3, 3)', 'A(3, 4)', 'A(3, 6)', 'A(4, 4)', 'A(6, 2)', 'A(6, 6)',
               'C(2, 2)', 'C(3, 3)', 'C(4, 4)', 'C(6, 2)', 'C(6, 6)']
prime_data_B = ['B(3, 2)', 'B(3, 4)', 'B(3, 6)',
                'B(3, 2; 2)', 'B(3, 6; 2)', 'B(3, 4; 4)']
prime_data_BT = ['B(1, 4)_2', 'B(1, 12)',
                 'B(2, 4)', 'B(1, 4; 2)_2', 'B(1, 12; 2)', 'B(2, 4; 4)',
                 'B(T, 1)', 'B(T, 2)', 'B(T, 3)', 'B(T, 1; 1)']

J = matrix.block(2,2,[0,identity_matrix(3),-identity_matrix(3),0])
big_list2 = []
big_list3 = []
for (name, G2) in big_list1:
     if name != 'B(T, 1; 1)':
          G3 = enlarge_group(G2, [J])
          Jname = "J(%s)" %name
          big_list2.append((Jname, G3))
     if name in prime_data_A:
          G3 = enlarge_group(G2, [J*embed(T_mat[1])])
     elif name in prime_data_B:
          G3 = enlarge_group(G2, [J*embed(diagonal_matrix([-1,1,-1]))])
     elif name in prime_data_BT:
          G3 = enlarge_group(G2, [J*embed(diagonal_matrix([z4,z4,-1]))])
     else:
          continue
     Jname = "J_s(%s)" %name
     big_list3.append((Jname, G3))

prime_data_A2a = ['A(1, 2)', 'A(1, 4)_1', 'A(1, 6)_1', 'A(3, 2)', 'A(3, 4)', 'A(3, 6)']
prime_data_A2b = ['A(1, 4)_2', 'A(1, 12)', 'A(2, 4)']
big_list4 = []
for (name, G2) in big_list1:
     if name in prime_data_A2a:
          G3 = enlarge_group(G2, [J*embed(Matrix([[1,0,0],[0,0,1],[0,-1,0]]))])
     elif name in prime_data_A2b:
          G3 = enlarge_group(G2, [J*embed(Matrix([[z4,0,0],[0,0,z4],[0,1,0]]))])
     elif name == 'E(36)':
          G3 = enlarge_group(G2, [J*embed(g3*g2*~g3)])
     else:
          continue
     Jname = "J_n(%s)" %name
     big_list4.append((Jname, G3))

big_list = big_list1 + big_list2 + big_list3 + big_list4
groupdict = dict(big_list)

print("Groups computed")
groups_computed = True


