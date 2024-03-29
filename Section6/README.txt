This folder contains 18 ".m" files (magma files) and the folder "Output", containing 137 ".txt" files (text files). The magma files perform the computations associated to Section 6 of 

"Sato-Tate groups of abelian threefolds", cited [FKS21] from now on.

The folder "Output" contains the data text files generated by the magma files. The data is relative to the 433 Sato-Tate groups classified in [FKS21].

With two exceptions, the text files can be generated by executing the magma file "execute.m". The exceptions are: the files "generators_*.txt", which give generators in TeX form for the groups of components of the Sato-Tate groups; the file "lmfdb_total.txt", described below. The execution of the code took around 6 hours on an intel i9-9960X with 128MB memory.

1) DESCRIPTION OF THE TEXT FILES (OUTPUT):

They are all magma readable files containing a list called "L" of diferent types of magma objects.

A) For each absolute type * in [A,B,...,M,N], as defined in [FKS21], the folder "Output" contains files

"names_*.txt"
"a1a2a3_*.txt"
"diagonals_*txt"
"comp_diagonals_*.txt"
"genenerators_*.txt"
"labels_*.txt"
"lattice_*.txt"
"moments_*.txt"
"zvectors_*.txt"

The file "names_*.txt" contains the list of the names (as appear on the LMFDB) of the Sato-Tate groups of type *. This sets an order of the groups which is again used for "a1a2a3_*.txt",...,"z_vectors.txt".

The file "a1a2a3_*.txt" contains the list of "a1a2a3-moments" of the Sato-Tate groups of type * (in the order prescribed by "names_*.txt"). For a group G in the classification, the "a1a2a3-moments" is a 36-tuple of integers containing the first 12 moments of a1, of a2, and of a3.

The file "diagonals_*.txt" contains the list of the 3-diagonals of character norms (as defined in [Sec. 6, FKS21]) of the Sato-Tate groups of type * (in the order prescribed by "names_*.txt"). The 3-diagonal of the group G is a 20-tuple of nonnegative integer numbers: it contains the expectation of the restriction to G of the 20 characters of USp(6) listed on [FKS21,Table 6] (according to that order).  

The file "comp_diagonals_*.txt" contains the list of distinct 3-diagonals of character norms that arise among the connected components of the Sato-Tate groups of type *. For a connected component, the 3-diagonal of character norms is a 20-tuple of integers, now possibly negative. 

The file "generators_*.txt" contains the list of "generators lists" of the Sato-Tate groups of type * (in the order prescribed by "names_*.tx"). For a group G in the classification, its "generators list" is a list of TeX strings describing matrix generators for the group of components of G.

The file "labels_*.txt" contains the list of the GAP identification labels of the groups of components of the Sato-Tate groups of type * (in the order prescribed by "names_*.txt"). 

The file "lattice_*.txt" contains the list of "subgroup lists" of the Sato-Tate groups of type * (in the order prescribed by "names_*.tx"). For a group G in the classification, its "subgroup list" is the list of names of the maximal subgroups of finite index of G (up to conjugacy and with multiplicities).

The file "moments_*.txt" contains the list of the 12-simplices of moments (as defined in [Sec. 6, FKS21]) of the Sato-Tate groups of type * (in the order prescribed by "names_*.txt"). The 12-simplex of a group G is a 56-tuple of nonnegative integer numbers: it contains the expectation M_{i,j,k}(G) of a1^i*a2^j*a3^k over G, for i+2*j+3*k <=12. The expectations are ordered by increasing weight i+2*j+3*k and, within the same weight, lexicographically in (i,j,k).
 
The file "zvectors_*.txt" contains a list of 23-tuples of rational numbers for each of the Sato-Tate groups of type * (in the order prescribed by "names_*.txt"). This tuple contains enough information to reconstruct the matrix of point densities (as defined in [Sec. 6, FKS21]). More precisely, this tuple is

(z_1,z_2^(-1),z_2^0,z_2^1,z_2^2,z_2^3,z_3,z_{12}^(-1),z_{12}^0,z_{12}^1,z_{12}^2,z_{12}^3,
z_{13},z_{23}^(-1),z_{23}^0,z_{23}^1,z_{23}^2,z_{23}^3,z_{123}^(-1),z_{123}^0,z_{123}^1,z_{123}^2,z_{123}^3) 


B) The files

"names_total.txt"
"a1a2a3_total.txt"
"diagonals_total.txt"
"generators_total.txt"
"labels_total.txt"
"lmfdb_total.txt"
"lattice_total.txt"
"moments_total.txt"
"zvectors_total.txt"

respectively contain lists (of length 433) of the names, a1a2a3-moments, 3-diagonals of character norms, GAP labels for the groups of components, lists of maximal subgroups, LMFDB labels, 12-simplices of moments, and z-vectors for all groups of the classification. The groups are ordered first lexicographically in their type, and within the type *, according to the order prescribed by "names_*.txt".  Note that LMFDB labels exist only for the 410 groups that actually arise as Sato-Tate groups of abelian varieties, for the other 23 the corresponding entry is set to "None".

The file "comp_diagonals_total.txt" collects all the 3-diagonals of character norms arising among the connected components of the groups in the classification (with multiplicity, if they arise for distinct connected components of the identity).

The file "orthogonality_relations_total.txt" contains the list of "orthogonality relations matrices" for each of the groups in the classification. For a group G in the classification, its "orthogonality relations matrix" contains the expectations over G of the products chi_i*chi_j, where chi_i and chi_j run through the first 15 characters listed on [FKS21, Table 6] (it is stored in the form of a 225-tuple of nonnegtive integers). The diagonal of the "orthogonality relations matrix" for G is thus the 2-diagonal of character norms of G. 


2) DESCRIPTION OF THE MAGMA FILES

A) We now describe the magma files. For each type *, the files "names_*.txt", "a1a2a3_*.txt", "comp_diagonals_*.txt", "diagonals_*.txt", "labels_*.txt", "lattice_*.txt", "moments_*.txt", and "zvectors_*.txt" are generated by the file "type_*.m". The file "type_*.m" first builds matrix generators for the groups of components of the Sato-Tate groups of type *. Then the function "generate_data" computes the data associated to the groups of type * (types B, E, F are treated slightly differently). The boolean parameters in "generate_data" can be set to generate only certain specific data (e.g. moments, diagonals, subgroups lattice, etc). At the top of some of the "type_*.m" files certain booleans can be activated to perform additional computations (for example, "type_L.m" allows to perform the verification of [Remark 4.1, FKS21]). 

The file "execute.m" calls "type_*.m" for every absolute type *. It suffices to load this file in Magma to generate all the text files in the folder (except the text files "generators_*.txt" that give generators in TeX form for the groups of components of the Sato-Tate groups).


B) There are the following additional magma files:

"common_functions.m"
"charactersUSp2g.m"
"Characters_in_terms_of_coef.m"

The file "common_functions.m" is called by "type_*.m" for every type *. It contains the common basic functions used to compute moments, labels, subgroup lattices, and z-vectors for all types. It asssumes defined functions "momentsp" and "label_group" which are specific to each type. It also contains functions "genus2_E" and "genus2_F" that constructs the Sato-Tate groups of dimension 2 of type E and F. 

The file "charactersUSp2g.m" contains the functions necessary to work with the characters of USp(2g). It assumes that a value of "g" has been given and that the field of fractions "FR" of a polynomial ring with (at least) g+1 variables has been defined. It is called by "type_A.m" (with g=3), "type_C.m" (with g=2) and "type_D.m" (with g=2).

The file "Characters_in_terms_of_coef.m" expresses the characters of USp(6) in terms of the coefficients of a generic characteristic polynomial. It calls the file "charactersUSp2g.m" and it is used to determine [Table 6, FKS21]. 


