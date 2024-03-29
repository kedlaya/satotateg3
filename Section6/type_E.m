//This file builds the groups of type E.
//It then produces the files:
// names_E.txt
// moments_E.txt
// diagonals_E.txt
// labels_E.txt
// zvectors_E.txt
// a1a2a3_E.txt
// lattice_E.txt
// comp_diagonals_E.txt

K:= Rationals();

GK2:=GeneralLinearGroup(2,K);
MK2:=MatrixRing(K,2);
GK3:=GeneralLinearGroup(3,K);
MK3:=MatrixRing(K,3);
GK4:=GeneralLinearGroup(4,K);
MK4:=MatrixRing(K,4);
GK6:=GeneralLinearGroup(6,K);
MK6:=MatrixRing(K,6);

PR<u,v,w>:=PolynomialRing(Rationals(),3);
FR:=FieldOfFractions(PR);
AssignNames(~FR,["u","v","w"]);

U6:=Matrix(FR,6,6,[[u,0,0,0,0,0],[0,1/u,0,0,0,0],[0,0,v,0,0,0],[0,0,0,1/v,0,0],[0,0,0,0,w,0],[0,0,0,0,0,1/w]]);
mI1:=GK6!DiagonalMatrix(K,[-1,-1,1,1,1,1]);
mI2:=GK6!DiagonalMatrix(K,[1,1,-1,-1,1,1]);
mI3:=GK6!DiagonalMatrix(K,[1,1,1,1,-1,-1]);
I6:=GK6!IdentityMatrix(K,6);

coef_var:=function(p,var,d)
	dd:=Degree(Denominator(p),var);
	cd:=Coefficient(Denominator(p), var, dd);
	cn:=Coefficient(Numerator(p),var,d+dd);
	return cn/cd;
end function;

avervalue:=function(q)
	coefw:=coef_var(q,w,0)-coef_var(q,w,2);
	coefv:=coef_var(coefw,v,0)-coef_var(coefw,v,2);
	return coef_var(coefv,u,0)-coef_var(coefv,u,2);
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

t:=GK6!Matrix(K,6,6,[[1,0,0,0,0,0],[0,1,0,0,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1],[0,0,1,0,0,0],[0,0,0,1,0,0]]);

s:=GK6!Matrix(K,6,6,[[0,0,0,0,1,0],[0,0,0,0,0,1],[1,0,0,0,0,0],[0,1,0,0,0,0],[0,0,1,0,0,0],[0,0,0,1,0,0]]);

//We generate "names_E.txt", "moments_E.txt", "diagonals_E.txt", "zvectors_E.txt", "labels_E.txt". E and the remaining groups generated by [t],[s],[s,t] need to be treated sparately

p:=CharacteristicPolynomial(U6);
momI:=moments_simplexp(p);
diaI:=diagonalp(p);
zvectorI:=z_vector_p(p);
a1a2a3I:=a1a2a3p(p);

PR<u,v>:=PolynomialRing(K,2);
FR:=FieldOfFractions(PR);
AssignNames(~FR,["u","v"]);

PFR<T>:=PolynomialRing(FR);

avervalue:=function(q)
	coefu:=coef_var(q,u,0)-coef_var(q,u,2);
	return coef_var(coefu,v,0)-coef_var(coefu,v,2);
end function;

load "common_functions.m";

p:=(1-(u+1/u)*T+T^2)*(1-(v+1/v)*T^2+T^4);
momt:=moments_simplexp(p);
diat:=diagonalp(p);
zvectort:=z_vector_p(p);
a1a2a3t:=a1a2a3p(p);

p:=(1-(v+1/v)*T^3+T^6);
moms:=moments_simplexp(p);
dias:=diagonalp(p);
zvectors:=z_vector_p(p);
a1a2a3s:=a1a2a3p(p);

momentsI:=momI;
momentsC2:=(momI+momt)/2;
momentsA3:=(momI+2*moms)/3;
momentsS3:=(momI+2*moms+3*momt)/6;

diagonalI:=diaI;
diagonalC2:=(diaI+diat)/2;
diagonalA3:=(diaI+2*dias)/3;
diagonalS3:=(diaI+2*dias+3*diat)/6;

zvectorI:=zvectorI;
zvectorC2:=(zvectorI+zvectort)/2;
zvectorA3:=(zvectorI+2*zvectors)/3;
zvectorS3:=(zvectorI+2*zvectors+3*zvectort)/6;

a1a2a3I:=a1a2a3I;
a1a2a3C2:=(a1a2a3I+a1a2a3t)/2;
a1a2a3A3:=(a1a2a3I+2*a1a2a3s)/3;
a1a2a3S3:=(a1a2a3I+2*a1a2a3s+3*a1a2a3t)/6;

gn:=["E","E_t", "E_s", "E_{s,t}"];

gm:=[momentsI, momentsC2, momentsA3, momentsS3];
gm:=format_good(gm);

gdi:=[diagonalI, diagonalC2, diagonalA3, diagonalS3];
gdi:=format_good(gdi);

gz:=[zvectorI, zvectorC2, zvectorA3, zvectorS3];
gz:=format_good(gz);

glab:=[<1,1>, <2,1>, <3,1>, <6,1>];

ga1a2a3:=[a1a2a3I, a1a2a3C2, a1a2a3A3, a1a2a3S3];
ga1a2a3:=format_good(ga1a2a3);

gcomp:=[diaI,diat,dias];
gcomp:=format_good(gcomp);

print_in_file_str(gn,"names_E.txt");

print_in_file(gm,"moments_E.txt");

print_in_file(gdi,"diagonals_E.txt");

print_in_file(gz,"zvectors_E.txt");

print_in_file(glab,"labels_E.txt");

print_in_file(ga1a2a3,"a1a2a3_E.txt");

print_in_file(gcomp,"comp_diagonals_E.txt");

//We generate the subgroup lattice

glat:=[
[],
["\"E\""],
["\"E\""],
["\"E_s\"", "\"E_t\""]
];
print_in_file(glat,"lattice_E.txt");

