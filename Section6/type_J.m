//This file builds the groups of type J.
//It then produces the files:
// names_J.txt
// moments_J.txt
// diagonals_ J.txt
// labels_J.txt
// zvectors_J.txt
// a1a2a3_J.txt
// lattice_J.txt
// comp_diagonals_J.txt

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
	coefv:=coef_var(q,v,0);
	return coef_var(coefv,u,0)-coef_var(coefv,u,2);
end function;

//We can call "common_functions.m" only once "avervalue" and "label_group" have been defined.

load "common_functions.m";

big_list:=[];//This will contain the names of the 31 groups of type J
big_dict:=AssociativeArray();//This will contain the 31 groups of type J

//Let G be a type E genus 2 group. There are 3 types of groups of type J. Those of the form:
//1) U(1)xG are denoted <su,G>.
//2) N(U(1))xG are denoted <nsu,G>.
//3) N(U(1))x_{C_2} G are denoted <H,G>, where H is the label of the kernel of the hom G-->C_2. 
//We build the 20 direct products

small_list:=genus2_E(K)[1];
small_dict:=genus2_E(K)[2];

N6:=embed6(MK2!Matrix(K,2,2,[0,1,-1,0]),MK4!DiagonalMatrix(K,[1,1,1,1]));

for i in small_list do
	name:=<"su",i>;
	L:=[embed6(MK2!DiagonalMatrix(K,[1,1]),MK4!g): g in small_dict[i]];
	big_dict[name]:=sub<GK6|L>;
	big_list:=Append(big_list,name);
end for;

for i in small_list do
	big_dict[<"nsu",i>]:=sub<GK6|big_dict[<"su",i>], N6>;
	big_list:=Append(big_list,<"nsu",i>);
end for;

//We build the 11 fiber products


fiber_products_list:=[];
fiber_products_dict:=AssociativeArray();

fiber_product:=function(G,H)
	elements_H:=[G!h : h in H];
	generators:=[];
	for g in G do
		if not (g in elements_H) then 
			gg:=embed6(MK2!Matrix(K,2,2,[0,1,-1,0]),MK4!g);
			generators:=Append(generators,gg); 
		end if;
	end for;
	return sub<GK6|generators>;
end function;

for i in [1,2,3,4,6] do
	nameH:="E_" cat IntegerToString(i);
	nameG:="J(E_" cat IntegerToString(i) cat ")";
	H:=small_dict[nameH];
	G:=small_dict[nameG];
	fiber_name:= <nameH,nameG>;
	big_dict[fiber_name]:=fiber_product(G,H);
	big_list:=Append(big_list, fiber_name);
end for;


for i in [1,2,3] do
	nameH:="E_" cat IntegerToString(i);
	nameG:="E_" cat IntegerToString(2*i);
	H:=small_dict[nameH];
	G:=small_dict[nameG];
	fiber_name:= <nameH, nameG>;
	big_dict[fiber_name]:=fiber_product(G,H);
	big_list:=Append(big_list, fiber_name);
end for;


for i in [1,2,3] do
	nameH:="J(E_" cat IntegerToString(i) cat ")";
	nameG:="J(E_" cat IntegerToString(2*i) cat ")";
	H:=small_dict[nameH];
	G:=small_dict[nameG];
	fiber_name:= <nameH, nameG>;
	big_dict[fiber_name]:=fiber_product(G,H);
	big_list:=Append(big_list, fiber_name);
end for;

groups_list:=[];
for i in big_list do
	groups_list:=Append(groups_list,big_dict[i]);
end for;


//We set tex names

gn:=[];
for i in [1..#big_list] do
	na:=big_list[i];
	if na[1] eq "su" then name:="J_1(" cat na[2] cat ")";
	elif na[1] eq "nsu" then name:="J_2(" cat na[2] cat ")";
	else name:="J(" cat na[2] cat "," cat na[1] cat ")";
	end if;
	
	gn:=Append(gn,name);
end for;

//We generate "names_J.txt", "moments_J.txt", "diagonals_J.txt", "zvectors_J.txt", "labels_J.txt", "a1a2a3_J.txt", "comp_diagonals_J.txt"

generate_data("J", gn, groups_list, true,true,true,true,true,true);


//We next generate the subgroup lattice
load "diagonals_J.txt";
gdi:=L;
load "labels_J.txt";
glab:=L;
load "names_J.txt";
gn:=L;
if compute_lattice then
	glat:=subgroups_lattice(groups_list,gdi,gn,glab);
	print_in_file(glat,"lattice_J.txt");
end if;




