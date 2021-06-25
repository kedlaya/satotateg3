//This file builds the groups of type C.
//It then produces the files:
// names_C.txt
// moments_C.txt
// diagonals_C.txt
// labels_C.txt
// zvectors_C.txt
// a1a2a3_C.txt
// lattice_C.txt
// comp_diagonals_C.txt

K:= Rationals();

GK2:=GeneralLinearGroup(2,K);
MK2:=MatrixRing(K,2);
GK3:=GeneralLinearGroup(3,K);
MK3:=MatrixRing(K,3);
GK4:=GeneralLinearGroup(4,K);
MK4:=MatrixRing(K,4);
GK6:=GeneralLinearGroup(6,K);
MK6:=MatrixRing(K,6);

g:=2;
PR<u,v,T,w>:=PolynomialRing(K,4);
FR:=FieldOfFractions(PR);
AssignNames(~FR,["u","v","T","w"]);

U6:=Matrix(FR,6,6,[[w,0,0,0,0,0],[0,1/w,0,0,0,0],[0,0,u,0,0,0],[0,0,0,1/u,0,0],[0,0,0,0,v,0],[0,0,0,0,0,1/v]]);
I6:=GK6!IdentityMatrix(K,6);

//Both g and PR,FR need to be defined before calling "CharactersUSp2g.m"
load "CharactersUSp2g.m";

coef_var:=function(p,var,d)
	dd:=Degree(Denominator(p),var);
	cd:=Coefficient(Denominator(p), var, dd);
	cn:=Coefficient(Numerator(p),var,d+dd);
	return cn/cd;
end function;

avervalue:=function(q)
	coefw:=coef_var(q,w,0)-coef_var(q,w,2);
	return mult_id(coefw);
end function;

label_group:=function(G)
	return IdentifyGroup(G);
end function;

//We can call "common_functions.m" only once "avervalue" and "label_group" have been defined.

load "common_functions.m";

groups_list:=[sub<GK6|I6>];

gn:=["SU(2)xUSp(4)"];

//We generate "names_C.txt", "moments_C.txt", "diagonals_C.txt", "zvectors_C.txt", "labels_C.txt", "a1a2a3_C.txt", "comp_diagonals_C.txt"

generate_data("C", gn, groups_list, true,true,true,true,true, true);

//We generate the subgroup lattice
load "diagonals_C.txt";
gdi:=L;
load "names_C.txt";
gn:=L;
load "labels_C.txt";
glab:=L;
	
glat:=subgroups_lattice(groups_list,gdi,gn,glab);
print_in_file(glat,"lattice_C.txt");


