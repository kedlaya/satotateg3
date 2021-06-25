//This file builds the groups of type B.
//It then produces the files:
// names_B.txt
// moments_B.txt
// diagonals_B.txt
// labels_B.txt
// zvectors_B.txt
// a1a2a3_B.txt
// lattice_B.txt
// comp_diagonals_B.txt

K:= Rationals();

GK2:=GeneralLinearGroup(2,K);
MK2:=MatrixRing(K,2);
GK3:=GeneralLinearGroup(3,K);
MK3:=MatrixRing(K,3);
GK4:=GeneralLinearGroup(4,K);
MK4:=MatrixRing(K,4);
GK6:=GeneralLinearGroup(6,K);
MK6:=MatrixRing(K,6);

//The general methodology will apply only to U(3). N(U(3)) will be treated separately.

PR<u,v,T,w>:=PolynomialRing(K,4);
FR:=FieldOfFractions(PR);
AssignNames(~FR,["u","v","T","w"]);

U6:=Matrix(FR,6,6,[[u*w,0,0,0,0,0],[0,v*w,0,0,0,0],[0,0,1/u*1/v*w,0,0,0],[0,0,0,1/u*1/w,0,0],[0,0,0,0,1/v*1/w,0],[0,0,0,0,0,u*v*1/w]]);
I6:=GK6!IdentityMatrix(K,6);
J6:=GK6!Matrix(K,6,6,[[0,0,0,1,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1],[-1,0,0,0,0,0],[0,-1,0,0,0,0],[0,0,-1,0,0,0]]);
mI6:=GK6!DiagonalMatrix(K,[-1,-1,-1,-1,-1,-1]);

//PR, FR need to be defined before calling "CharactersSU3.m"
load "CharactersSU3.m";

coef_var:=function(p,var,d)
	dd:=Degree(Denominator(p),var);
	cd:=Coefficient(Denominator(p), var, dd);
	cn:=Coefficient(Numerator(p),var,d+dd);
	return cn/cd;
end function;

avervalue:=function(q)
	coefw:=coef_var(q,w,0);
	return mult_idSU3(coefw);
end function;

label_group:=function(G)
	return IdentifyGroup(quo<G|mI6>);
end function;

//We can call "common_functions.m" only once "avervalue" and "label_group" have been defined.

load "common_functions.m";

//We generate "names_B.txt", "moments_B.txt", "diagonals_B.txt", "zvectors_B.txt", "labels_B.txt" for U(3) and N(U(3))

p:=CharacteristicPolynomial(U6);
momentsU3:=format_good_aux(moments_simplexp(p));
diagonalU3:=format_good_aux(diagonalp(p));
zvectorU3:=format_good_aux(z_vector_p(p));
a1a2a3U3:=format_good_aux(a1a2a3p(p));

//The data used to compute "momentsNU3" can be obtained by executing the file "computations_NU3.m"
momentsNU3:=<>;
for i in [1..#cc] do
	if not(cc[i][1] eq 0) or not(cc[i][3] eq 0) then
		momentsNU3:=Append(momentsNU3,momentsU3[i]/2);
	elif cc[i] eq [0,1,0] then momentsNU3:=Append(momentsNU3,1);
	elif cc[i] eq [0,2,0] then momentsNU3:=Append(momentsNU3,3);
	elif cc[i] eq [0,3,0] then momentsNU3:=Append(momentsNU3,11);
	elif cc[i] eq [0,4,0] then momentsNU3:=Append(momentsNU3,58);
	elif cc[i] eq [0,5,0] then momentsNU3:=Append(momentsNU3,382);
	elif cc[i] eq [0,6,0] then momentsNU3:=Append(momentsNU3,2922);
	end if;
end for;

diagonalJU3:=<
1,
0,
2-2*1+1,
0,
2,
0,
2-2*1+1,
9-2*4+2,
0,
9+4*2+1-4*4+2*2-4*1,
0,
9-4*4+4*2,
0,
0,
0,
0,
51-4*21+4*9,
0,
51+9*9+4*2-6*21+4*9-12*4,
0
>;

diagonalNU3:=<
1/2*(diagonalU3[1]+1),
1/2*(diagonalU3[2]),
1/2*(diagonalU3[3]+2-2*1+1),
1/2*(diagonalU3[4]),
1/2*(diagonalU3[5]+2),
1/2*(diagonalU3[6]),
1/2*(diagonalU3[7]+2-2*1+1),
1/2*(diagonalU3[8]+9-2*4+2),
1/2*(diagonalU3[9]),
1/2*(diagonalU3[10]+9+4*2+1-4*4+2*2-4*1),
1/2*(diagonalU3[11]),
1/2*(diagonalU3[12]+9-4*4+4*2),
1/2*(diagonalU3[13]),
1/2*(diagonalU3[14]),
1/2*(diagonalU3[15]),
1/2*(diagonalU3[16]),
1/2*(diagonalU3[17]+51-4*21+4*9),
1/2*(diagonalU3[18]),
1/2*(diagonalU3[19]+51+9*9+4*2-6*21+4*9-12*4),
1/2*(diagonalU3[20])
>;

//The only required information to compute the zvector for the component J(U(3)) is that a1 and a3 are 0 on J(U(3)), while a2 is non-constant on J(U(3)).

zvectorNU3:=<>;
for i in [1..23] do
	if i in [1,7,13] then
		zvectorNU3:=Append(zvectorNU3,1/2);
	else
		zvectorNU3:=Append(zvectorNU3,0);
	end if;
end for;
 
a1a2a3NU3:=<>;
for i in [1..36] do
	if i eq 13 then s:=1;
	elif i eq 14 then s:=3;
	elif i eq 15 then s:=11;
	elif i eq 16 then s:=58;
	elif i eq 17 then s:=382;
	elif i eq 18 then s:=2922;
	elif i in [19..24] then s:=-1; 
	else s:= a1a2a3U3[i]/2;
	end if;
	a1a2a3NU3:=Append(a1a2a3NU3,s);
end for;


gn:=["U(3)","N(U(3))"];

gm:=[momentsU3,momentsNU3];

gdi:=[diagonalU3,diagonalNU3];

gz:=[zvectorU3,zvectorNU3];

glab:=[<1,1>, <2,1>];

ga1a2a3:=<a1a2a3U3,a1a2a3NU3>;

gcomp:=[diagonalU3,diagonalJU3];

print_in_file_str(gn,"names_B.txt");

print_in_file(gm,"moments_B.txt");

print_in_file(gdi,"diagonals_B.txt");

print_in_file(gz,"zvectors_B.txt");

print_in_file(glab,"labels_B.txt");

print_in_file(ga1a2a3,"a1a2a3_B.txt");

print_in_file(gcomp, "comp_diagonals_B.txt");

//We next generate the subgroup lattice

glat:=[
[],
["\"U(3)\""]
];
print_in_file(glat,"lattice_B.txt");

