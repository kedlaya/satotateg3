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

g:=3;
PR:=PolynomialRing(K,4);

u:=PR.1;
v:=PR.2;
w:=PR.3;
T:=PR.4;

AssignNames(~PR,["u","v","w","T"]);
FR:=FieldOfFractions(PR);
AssignNames(~FR,["u","v","w","T"]);

U6:=Matrix(FR,6,6,[[u,0,0,0,0,0],[0,v,0,0,0,0],[0,0,w,0,0,0],[0,0,0,1/u,0,0],[0,0,0,0,1/v,0],[0,0,0,0,0,1/w]]);
I6:=GK6!IdentityMatrix(K,6);

load "CharactersUSp2g.m";

//label_group is only formally defined to be able to apply the general procedure

label_group:=function(G)
	return IdentifyGroup(G);
end function;

//We first consider U3

//The following function is based on Lemma 6.9.i)	
avervalue:=function(q)	
    we:=irred_dec(q)[1];
    mu:=irred_dec(q)[2];
    aver:=0;
    for i in [1..#we] do
        wei:=we[i];
        pes:=wei[1]+wei[2]+wei[3];
        if (wei[1] mod 2 eq 0) and (wei[2] mod 2 eq 0) and (wei[3] mod 2 eq 0) then
            aver:=aver+mu[i]; 
        end if;
    end for;
    return aver;
end function;

//We can call "common_functions.m" only once "avervalue" and "label_group" have been defined.

load "common_functions.m";	

p:=CharacteristicPolynomial(U6);
momentsU3:=format_good_aux(moments_simplexp(p));
diaU3:=diagonalp(p);
diagonalU3:=format_good_aux(diaU3);
zvectorU3:=format_good_aux(z_vector_p(p));
a1a2a3U3:=format_good_aux(a1a2a3p(p));

//We next consider NU3

//The following function is based on Lemma 6.9.ii)	
avervalue:=function(q)	
    we:=irred_dec(q)[1];
    mu:=irred_dec(q)[2];
    aver:=0;
    for i in [1..#we] do
        wei:=we[i];
        pes:=wei[1]+wei[2]+wei[3];
        if (pes mod 4 eq 0) and (wei[1] mod 2 eq 0) and (wei[2] mod 2 eq 0) and (wei[3] mod 2 eq 0) then
            aver:=aver+mu[i]; 
        end if;
    end for;
    return aver;
end function;

//We can call "common_functions.m" only once "avervalue" and "label_group" have been defined.

load "common_functions.m";	

p:=CharacteristicPolynomial(U6);
momentsNU3:=format_good_aux(moments_simplexp(p));
diaNU3:=diagonalp(p);
diagonalNU3:=format_good_aux(diaNU3);
zvectorNU3:=<>;
for i in [1..23] do
	if i in [1,7,13] then
		zvectorNU3:=Append(zvectorNU3,1/2);
	else
		zvectorNU3:=Append(zvectorNU3,0);
	end if;
end for;
a1a2a3NU3:=format_good_aux(a1a2a3p(p));


gn:=["U(3)","N(U(3))"];

gm:=[momentsU3,momentsNU3];

gdi:=[diagonalU3,diagonalNU3];

gz:=[zvectorU3,zvectorNU3];

glab:=[<1,1>, <2,1>];

ga1a2a3:=<a1a2a3U3,a1a2a3NU3>;

diaJU3:=2*diaNU3-diaU3;
diagonalJU3:=format_good_aux(diaJU3);

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
