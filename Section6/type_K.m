//This file builds the groups of type K.
//It then produces the files:
// names_K.txt
// moments_K.txt
// diagonals_K.txt
// labels_K.txt
// zvectors_K.txt
// a1a2a3_K.txt
// lattice_K.txt
// comp_diagonals_K.txt

compute_lattice:=true;//If true, the code additionally computes the subgroup lattice

N := 24; K<zz> := CyclotomicField(N);

GK2:=GeneralLinearGroup(2,K);
MK2:=MatrixRing(K,2);
GK3:=GeneralLinearGroup(3,K);
MK3:=MatrixRing(K,3);
GK4:=GeneralLinearGroup(4,K);
MK4:=MatrixRing(K,4);
GK6:=GeneralLinearGroup(6,K);
MK6:=MatrixRing(K,6);

PR<u,v>:=PolynomialRing(K,2);
FR:=FieldOfFractions(PR);
AssignNames(~FR,["u","v"]);

U6:=Matrix(FR,6,6,[[v,0,0,0,0,0],[0,1/v,0,0,0,0],[0,0,u,0,0,0],[0,0,0,u,0,0],[0,0,0,0,1/u,0],[0,0,0,0,0,1/u]]);
mI4:=GK6!DiagonalMatrix(K,[1,1,-1,-1,-1,-1]);

label_group:=function(G)
	if (mI4 in G) then return IdentifyGroup(quo<G|mI4>); end if;
	return IdentifyGroup(G);
end function;

coef_var:=function(p,var,d)
	dd:=Degree(Denominator(p),var);
	cd:=Coefficient(Denominator(p), var, dd);
	cn:=Coefficient(Numerator(p),var,d+dd);
	return cn/cd;
end function;

avervalue:=function(q)
	coefv:=coef_var(q,v,0)-coef_var(q,v,2);
	return coef_var(coefv,u,0);
end function;

//We can call "common_functions.m" only once "avervalue" and "label_group" have been defined.

load "common_functions.m";

big_list:=[];//This will contain the names of the 32 groups of type K
big_dict:=AssociativeArray();//This will contain the 32 groups of type K
groups_list:=[];

small_list:=genus2_F(K)[1];
small_dict:=genus2_F(K)[2];

for i in small_list do
	L:=[embed6(MK2!DiagonalMatrix(K,[1,1]),MK4!g): g in small_dict[i]];
	big_dict[i]:=sub<GK6|L>;
	groups_list:=Append(groups_list,sub<GK6|L>);
	big_list:=Append(big_list,i);
end for;

//We set tex names

gn:=[];
for i in [1..#big_list] do
	gn:=Append(gn,"SU(2)x" cat big_list[i]);
end for;

//We finally generate "names_K.txt", "moments_K.txt", "diagonals_K.txt", "zvectors_K.txt", "labels_K.txt", "a1a2a3_K.txt", "comp_diagonals_K.txt"

generate_data("K", gn, groups_list, true,true,true,true,true,true);

//We next generate the subgroup lattice
load "diagonals_K.txt";
gdi:=L;
load "labels_K.txt";
glab:=L;
load "names_K.txt";
gn:=L;
if compute_lattice then
	glat:=subgroups_lattice(groups_list,gdi,gn,glab);
	print_in_file(glat,"lattice_K.txt");
end if;


