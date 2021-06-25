//This file builds the groups of type D.
//It then produces the files:
// names_D.txt
// moments_D.txt
// diagonals_D.txt
// labels_D.txt
// zvectors_D.txt
// a1a2a3_D.txt
// lattice_D.txt
// comp_diagonals.txt

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
mI1:=GK6!DiagonalMatrix(K,[-1,-1,1,1,1,1]);

//Both g and PR,FR need to be defined before calling "CharactersUSp2g.m"
load "CharactersUSp2g.m";

coef_var:=function(p,var,d)
	dd:=Degree(Denominator(p),var);
	cd:=Coefficient(Denominator(p), var, dd);
	cn:=Coefficient(Numerator(p),var,d+dd);
	return cn/cd;
end function;

avervalue:=function(q)
	coefw:=coef_var(q,w,0);
	return mult_id(coefw);
end function;

label_group:=function(G)
	if mI1 in G then return IdentifyGroup(quo<G|mI1>); end if;
	return IdentifyGroup(G);
end function;

//We can call "common_functions.m" only once "avervalue" and "label_group" have been defined.

load "common_functions.m";

N6:=embed6(MK2!Matrix(K,2,2,[0,-1,1,0]),MK4!DiagonalMatrix(K,[1,1,1,1]));
groups_list:=[sub<GK6|I6>,sub<GK6|N6>];

gn:=["U(1)xUSp(4)", "N(U(1)xUSp(4))"];


//We finally generate "names_D.txt", "moments_D.txt", "diagonals_D.txt", "zvectors_D.txt", "labels_D.txt", "a1a2a3_D.txt", "comp_diagonals_D.txt"

generate_data("D", gn, groups_list, true,true,true,true,true,true);

//We generate the subgroup lattice
load "diagonals_D.txt";
gdi:=L;
load "names_D.txt";
gn:=L;
load "labels_D.txt";
glab:=L;
	
glat:=subgroups_lattice(groups_list,gdi,gn,glab);
print_in_file(glat,"lattice_D.txt");



