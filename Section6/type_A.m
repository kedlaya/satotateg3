//This file builds the groups of type A.
//It then produces the files:
// names_A.txt
// moments_A.txt
// diagonals_A.txt
// labels_A.txt
// zvectors_A.txt
// a1a2a3_A.txt
// lattice_A.txt
// comp_diagonals_A.txt

K:= Rationals();

GK2:=GeneralLinearGroup(2,K);
MK2:=MatrixRing(K,2);
GK3:=GeneralLinearGroup(3,K);
MK3:=MatrixRing(K,3);
GK4:=GeneralLinearGroup(4,K);
MK4:=MatrixRing(K,4);
GK6:=GeneralLinearGroup(6,K);
MK6:=MatrixRing(K,6);

g:=3;
PR<u,v,w,T>:=PolynomialRing(K,4);
FR:=FieldOfFractions(PR);
AssignNames(~FR,["u","v","w","T"]);

U6:=Matrix(FR,6,6,[[u,0,0,0,0,0],[0,v,0,0,0,0],[0,0,w,0,0,0],[0,0,0,1/u,0,0],[0,0,0,0,1/v,0],[0,0,0,0,0,1/w]]);
I6:=GK6!IdentityMatrix(K,6);

//PR, FR need to be defined before calling "CharactersUSp2g.m"
load "CharactersUSp2g.m";


avervalue:=function(q)
	return mult_id(q);
end function;

label_group:=function(G)
	return IdentifyGroup(G);
end function;

//We can call "common_functions.m" only once "avervalue" and "label_group" have been defined.

load "common_functions.m";

groups_list:=[sub<GK6|I6>];

gn:=["USp(6)"];

//We generate "names_A.txt", "diagonals_A.txt", "moments_A.txt", "zvectors_A.txt", "labels_A.txt", "a1a2a3_A.txt", "comp_diagonals_A.txt"

generate_data("A", gn, groups_list, true,true,true,true, true,true);

//We next generate the subgroup lattice
load "diagonals_A.txt";
gdi:=L;
load "names_A.txt";
gn:=L;
load "labels_A.txt";
glab:=L;
	
glat:=subgroups_lattice(groups_list,gdi,gn,glab);
print_in_file(glat,"lattice_A.txt");



