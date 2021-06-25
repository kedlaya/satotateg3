//This file builds the groups of type I.
//It then produces the files:
// names_I.txt
// moments_I.txt
// diagonals_I.txt
// labels_I.txt
// zvectors_I.txt
// a1a2a3_I.txt
// lattice_I.txt
// comp_diagonals_I.txt

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

U6:=Matrix(FR,6,6,[[v,0,0,0,0,0],[0,1/v,0,0,0,0],[0,0,u,0,0,0],[0,0,0,1/u,0,0],[0,0,0,0,1/u,0],[0,0,0,0,0,u]]);
mI4:=GK6!DiagonalMatrix(K,[1,1,-1,-1,-1,-1]);
mI2:=GK6!DiagonalMatrix(K,[-1,-1,1,1,1,1]);

label_group:=function(G)
	if (mI2 in G) and (mI4 in G) then return IdentifyGroup(quo<G|mI2,mI4>); end if;
	if (mI2 in G) then return IdentifyGroup(quo<G|mI2>); end if;
	if (mI4 in G) then return IdentifyGroup(quo<G|mI4>); end if;
	if (mI4*mI2 in G) then return IdentifyGroup(quo<G|mI4*mI2>); end if;
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
	return coef_var(coefv,u,0)-coef_var(coefv,u,2);;
end function;


//We can call "common_functions.m" only once "avervalue" and "label_group" have been defined.

load "common_functions.m";

big_list:=[];//This will contain the name of the 10 groups of type I
big_dict:=AssociativeArray();//This will contain the 10 groups of type I
groups_list:=[];

small_list:=genus2_E(K)[1];
small_dict:=genus2_E(K)[2];


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

//We finally generate "names_I.txt", "moments_I.txt", "diagonals_I.txt", "zvectors_I.txt", "labels_I.txt", "a1a2a3_I.txt", "comp_diagonals_I.txt"

generate_data("I", gn, groups_list, true,true,true,true,true,true);


//We next generate the subgroup lattice
load "diagonals_I.txt";
gdi:=L;
load "labels_I.txt";
glab:=L;
load "names_I.txt";
gn:=L;

glat:=subgroups_lattice(groups_list,gdi,gn,glab);
print_in_file(glat,"lattice_I.txt");


