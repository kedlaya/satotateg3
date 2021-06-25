//This file builds the groups of type G.
//It then produces the files:
// names_G.txt
// moments_G.txt
// diagonals_G.txt
// labels_G.txt
// zvectors_G.txt
// a1a2a3_G.txt
// lattice_G.txt
// comp_diagonals_G.txt

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
I6:=GK6!DiagonalMatrix(K,[1,1,1,1,1,1]);

coef_var:=function(p,var,d)
	dd:=Degree(Denominator(p),var);
	cd:=Coefficient(Denominator(p), var, dd);
	cn:=Coefficient(Numerator(p),var,d+dd);
	return cn/cd;
end function;

avervalue:=function(q)
	coefw:=coef_var(q,w,0)-coef_var(q,w,2);
	coefv:=coef_var(coefw,v,0);
	return coef_var(coefv,u,0);
end function;

label_group:=function(G)
	L:=[];
	for l in [mI1, mI2, mI3, mI1*mI2, mI1*mI3, mI2*mI3, mI1*mI2*mI3] do	
		if (l in G) then L:=Append(L,l); end if;
	end for;
	return IdentifyGroup(quo<G|L>);
end function;

//We can call "common_functions.m" only once "avervalue" and "label_group" have been defined.

load "common_functions.m";

a:=GK6!Matrix(K,6,6,[[0,1,0,0,0,0],[-1,0,0,0,0,0],[0,0,1,0,0,0],[0,0,0,1,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1]]);

b:=GK6!Matrix(K,6,6,[[1,0,0,0,0,0],[0,1,0,0,0,0],[0,0,0,1,0,0],[0,0,-1,0,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1]]);

c:=GK6!Matrix(K,6,6,[[0,0,1,0,0,0],[0,0,0,1,0,0],[-1,0,0,0,0,0],[0,-1,0,0,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1]]);

generators:=[[I6],[a],[a*b],[a*c],[a,b],[c],[a*b,c],[a,b,c] ];

groups_list:=[];//This contains the 8 groups of type G
for l in generators do
	groups_list:=Append(groups_list, sub<GK6|l>);
end for;

//We set tex names

gn:=["FxSU(2)",
"F_axSU(2)",
"F_{ab}xSU(2)",
"F_{ac}xSU(2)",
"F_{a,b}xSU(2)",
"F_cxSU(2)",
"F_{ab,c}xSU(2)",
"F_{a,b,c}xSU(2)"
];

//We finally generate "names_G.txt", "moments_G.txt", "diagonals_G.txt", "zvectors_G.txt", "labels_G.txt", "a1a2a3_G.txt", "comp_diagonals_G.txt"

generate_data("G", gn, groups_list, true,true,true,true,true,true);

load "diagonals_G.txt";
gdi:=L;
load "names_G.txt";
gn:=L;
load "labels_G.txt";
glab:=L;
	
glat:=subgroups_lattice(groups_list,gdi,gn,glab);
print_in_file(glat,"lattice_G.txt");

