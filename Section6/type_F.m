//This file builds the groups of type F.
//It then produces the files:
// names_F.txt
// moments_F.txt
// diagonals_F.txt
// labels_F.txt
// zvectors_F.txt
// a1a2a3_F.txt
// lattice_F.txt
// comp_diagonals_F.txt

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
	coefv:=coef_var(coefw,v,0)-coef_var(coefw,v,2);
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

t:=GK6!Matrix(K,6,6,[[1,0,0,0,0,0],[0,1,0,0,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1],[0,0,1,0,0,0],[0,0,0,1,0,0]]);

//We generate "names_F.txt", "moments_F.txt", "diagonals_F.txt", "zvectors_F.txt", "labels_F.txt". The three groups involving the generators [t],[a*t],[a,t] need to be treated separately.

p:=CharacteristicPolynomial(U6);
momI:=moments_simplexp(p);
diaI:=diagonalp(p);
zvectorI:=z_vector_p(p);
a1a2a3I:=a1a2a3p(p);

p:=CharacteristicPolynomial(Parent(U6)!a*U6);
moma:=moments_simplexp(p);
diaa:=diagonalp(p);
zvectora:=z_vector_p(p);
a1a2a3a:=a1a2a3p(p);

PR<u,v>:=PolynomialRing(K,2);
FR:=FieldOfFractions(PR);
AssignNames(~FR,["u","v"]);

PFR<T>:=PolynomialRing(FR);

avervalue:=function(q)
	coefu:=coef_var(q,u,0);
	return coef_var(coefu,v,0)-coef_var(coefu,v,2);
end function;

load "common_functions.m";

p:=(1-(u+1/u)*T+T^2)*(1-(v+1/v)*T^2+T^4);
momt:=moments_simplexp(p);
diat:=diagonalp(p);
zvectort:=z_vector_p(p);
a1a2a3t:=a1a2a3p(p);

p:=(1+T^2)*(1-(v+1/v)*T^2+T^4);
momat:=moments_simplexp(p);
diaat:=diagonalp(p);
zvectorat:=z_vector_p(p);
a1a2a3at:=a1a2a3p(p);

momentsI:=momI;
momentsC2:=(momI+moma)/2;
momentsS2:=(momI+momt)/2;
momentsC2x_C2S2:=(momI+momat)/2;
momentsC2xS2:=(momI+moma+momt+momat)/4;

diagonalI:=diaI;
diagonalC2:=(diaI+diaa)/2;
diagonalS2:=(diaI+diat)/2;
diagonalC2x_C2S2:=(diaI+diaat)/2;
diagonalC2xS2:=(diaI+diaa+diat+diaat)/4;

zvectorI:=zvectorI;
zvectorC2:=(zvectorI+zvectora)/2;
zvectorS2:=(zvectorI+zvectort)/2;
zvectorC2x_C2S2:=(zvectorI+zvectorat)/2;
zvectorC2xS2:=(zvectorI+zvectora+zvectort+zvectorat)/4;

a1a2a3I:=a1a2a3I;
a1a2a3C2:=(a1a2a3I+a1a2a3a)/2;
a1a2a3S2:=(a1a2a3I+a1a2a3t)/2;
a1a2a3C2x_C2S2:=(a1a2a3I+a1a2a3at)/2;
a1a2a3C2xS2:=(a1a2a3I+a1a2a3a+a1a2a3t+a1a2a3at)/4;

gn:=["F", "F_a", "F_t" , "F_{at}", "F_{a,t}"];

gm:=[momentsI, momentsC2, momentsS2, momentsC2x_C2S2, momentsC2xS2];
gm:=format_good(gm);

gdi:=[diagonalI, diagonalC2, diagonalS2, diagonalC2x_C2S2, diagonalC2xS2];
gdi:=format_good(gdi);

gz:=[zvectorI, zvectorC2, zvectorS2, zvectorC2x_C2S2, zvectorC2xS2];
gz:=format_good(gz);

glab:=[<1,1>, <2,1>, <2,1>, <2,1>, <4,2>];

ga1a2a3:=[a1a2a3I, a1a2a3C2, a1a2a3S2, a1a2a3C2x_C2S2, a1a2a3C2xS2];
ga1a2a3:=format_good(ga1a2a3);

gcomp:=[diaI,diaa,diat,diaat];
gcomp:=format_good(gcomp);

print_in_file_str(gn,"names_F.txt");

print_in_file(gm,"moments_F.txt");

print_in_file(gdi,"diagonals_F.txt");

print_in_file(gz,"zvectors_F.txt");

print_in_file(glab,"labels_F.txt");

print_in_file(ga1a2a3,"a1a2a3_F.txt");

print_in_file(gcomp,"comp_diagonals_F.txt");

//We generate the subgroup lattice

glat:=[
[],
["\"F\""],
["\"F\""],
["\"F\""],
["\"F_a\"", "\"F_t\"", "\"F_{at}\""]
];

print_in_file(glat,"lattice_F.txt");

