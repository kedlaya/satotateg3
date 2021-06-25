//This file builds the groups of type H.
//It then produces the files:
// names_H.txt
// moments_H.txt
// diagonals_H.txt
// labels_H.txt
// zvectors_H.txt
// a1a2a3_H.txt
// lattice_H.txt
// comp_diagonals_H.txt

compute_lattice:=true;//If true, the code additionally computes the subgroup lattice

K:= Rationals();

GK2:=GeneralLinearGroup(2,K);
MK2:=MatrixRing(K,2);
GK3:=GeneralLinearGroup(3,K);
MK3:=MatrixRing(K,3);
GK4:=GeneralLinearGroup(4,K);
MK4:=MatrixRing(K,4);
GK6:=GeneralLinearGroup(6,K);
MK6:=MatrixRing(K,6);

PR<u,v,w>:=PolynomialRing(K,3);
FR:=FieldOfFractions(PR);
AssignNames(~FR,["u","v","w"]);

U6:=Matrix(FR,6,6,[[u,0,0,0,0,0],[0,1/u,0,0,0,0],[0,0,v,0,0,0],[0,0,0,1/v,0,0],[0,0,0,0,w,0],[0,0,0,0,0,1/w]]);
mI1:=GK6!DiagonalMatrix(K,[-1,-1,1,1,1,1]);
mI2:=GK6!DiagonalMatrix(K,[1,1,-1,-1,1,1]);
mI3:=GK6!DiagonalMatrix(K,[1,1,1,1,-1,-1]);

coef_var:=function(p,var,d)
	dd:=Degree(Denominator(p),var);
	cd:=Coefficient(Denominator(p), var, dd);
	cn:=Coefficient(Numerator(p),var,d+dd);
	return cn/cd;
end function;

label_group:=function(G)
	L:=[];
	for l in [mI1, mI2, mI3, mI1*mI2, mI1*mI3, mI2*mI3, mI1*mI2*mI3] do	
		if (l in G) then L:=Append(L,l); end if;
	end for;
	return IdentifyGroup(quo<G|L>);
end function;

avervalue:=function(q)
	coefw:=coef_var(q,w,0);
	coefv:=coef_var(coefw,v,0);
	return coef_var(coefv,u,0);
end function;

//We can call "common_functions.m" only once "avervalue" and "label_group" have been defined.

load "common_functions.m";

a:=GK6!Matrix(K,6,6,[[0,1,0,0,0,0],[-1,0,0,0,0,0],[0,0,1,0,0,0],[0,0,0,1,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1]]);

b:=GK6!Matrix(K,6,6,[[1,0,0,0,0,0],[0,1,0,0,0,0],[0,0,0,1,0,0],[0,0,-1,0,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1]]);

c:=GK6!Matrix(K,6,6,[[1,0,0,0,0,0],[0,1,0,0,0,0],[0,0,1,0,0,0],[0,0,0,1,0,0],[0,0,0,0,0,1],[0,0,0,0,-1,0]]);

t:=GK6!Matrix(K,6,6,[[0,0,1,0,0,0],[0,0,0,1,0,0],[1,0,0,0,0,0],[0,1,0,0,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1]]);

s:=GK6!Matrix(K,6,6,[[0,0,0,0,1,0],[0,0,0,0,0,1],[1,0,0,0,0,0],[0,1,0,0,0,0],[0,0,1,0,0,0],[0,0,0,1,0,0]]);

generators:=[[IdentityMatrix(K,6)],[a],[a*b],[a*b*c],[t],[c*t],[s],[a*t],[a*c*t],[a,b],[a,b*c],[a*b,b*c],[c,t],[a*b,t],[a*b,c*t],[a*b*c,t],[a*b*c,s],[s,t],[a*b*c*t,s],[c,a*t],[a,b,c],[a*b,c,t],[a,b,t],[a*b,b*c,t],[a,b,c*t],[a*b,b*c,c*t],[a*b*c,s,t],[a*b,b*c,s],[a,b,c,t],[a,b,c,s],[a*b,b*c,s,t],[a*b,b*c,a*t,s],[a,b,c,s,t]];


groups_list:=[];//This will contain the 33 groups of type H
for l in generators do
	groups_list:=Append(groups_list,sub<GK6|l>);
end for;

//We set tex names

gn:=["H","H_a", "H_{ab}", "H_{abc}", "H_t", "H_{ct}", "H_s", "H_{at}", "H_{act}", "H_{a,b}", "H_{a,bc}", "H_{ab,bc}", "H_{c,t}", "H_{ab,t}", "H_{ab,ct}", "H_{abc,t}", "H_{abc,s}", "H_{s,t}", "H_{abct,s}", "H_{c,at}", "H_{a,b,c}", "H_{ab,c,t}", "H_{a,b,t}", "H_{ab,bc,t}", "H_{a,b,ct}", "H_{ab,bc,ct}", "H_{abc,s,t}", "H_{ab,bc,s}", "H_{a,b,c,t}", "H_{a,b,c,s}", "H_{ab,bc,s,t}", "H_{ab,bc,at,s}", "H_{a,b,c,s,t}"];


//We generate "names_H.txt", "moments_H.txt", "diagonals_H.txt", "zvectors_H.txt", "labels_H.txt", "a1a2a3_H.txt", "comp_diagonals_H.txt"

generate_data("H", gn, groups_list, true,true,true,true, true, true);


//We next generate the subgroup lattice
load "diagonals_H.txt";
gdi:=L;
load "names_H.txt";
gn:=L;
load "labels_H.txt";
glab:=L;
if compute_lattice then
	
	glat:=subgroups_lattice(groups_list,gdi,gn,glab);
	print_in_file(glat,"lattice_H.txt");
end if;


