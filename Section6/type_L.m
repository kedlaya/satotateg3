//This file builds the groups of type L.
//It then produces the files:
// names_L.txt
// moments_L.txt
// diagonals_L.txt
// labels_L.txt
// zvectors_L.txt
// a1a2a3_L.txt
// lattice_L.txt
// comp_diagonals_L.txt

verify_fiber_products:=false;//If true, the code additionally verifies that the 73 fiber products only yield 58 multisets of characteristic polynomials
verify_distinguished_by_moments:=false;//If true, the code addicitonally verifies that the fiber product of N(U(1)) and J(D_6) along J(C_6) and the fiber product of N(U(1)) and J(D_6) along D_6 have distinct 14-simplex of moments

compute_lattice:=true;//If true, the code additionally computes the subgroup lattice

N := 24; K := CyclotomicField(N);

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

U4:=Matrix(FR,4,4,[[u,0,0,0],[0,u,0,0],[0,0,1/u,0],[0,0,0,1/u]]);
U6:=Matrix(FR,6,6,[[v,0,0,0,0,0],[0,1/v,0,0,0,0],[0,0,u,0,0,0],[0,0,0,u,0,0],[0,0,0,0,1/u,0],[0,0,0,0,0,1/u]]);

coef_var:=function(p,var,d)
	dd:=Degree(Denominator(p),var);
	cd:=Coefficient(Denominator(p), var, dd);
	cn:=Coefficient(Numerator(p),var,d+dd);
	return cn/cd;
end function;

avervalue:=function(q)
	coefv:=coef_var(q,v,0);
	return coef_var(coefv,u,0);
end function;

mI4:=GK4!DiagonalMatrix(K,[-1,-1,-1,-1]);
mI46:=GK6!DiagonalMatrix(K,[1,1,-1,-1,-1,-1]);
mI26:=GK6!DiagonalMatrix(K,[-1,-1,1,1,1,1]);

label_group:=function(G)
	if (mI26 in G) and (mI46 in G) then return IdentifyGroup(quo<G|mI26,mI46>); end if;
	if (mI26 in G) then return IdentifyGroup(quo<G|mI26>); end if;
	if (mI46 in G) then return IdentifyGroup(quo<G|mI46>); end if;
	if (mI46*mI26 in G) then return IdentifyGroup(quo<G|mI46*mI26>); end if;
	return IdentifyGroup(G);
end function;

//We can call "common_functions.m" only once "avervalue" and "label_group" have been defined.

load "common_functions.m";

big_list:=[];//This list will contain the names of the 122 groups
big_dict:=AssociativeArray();//This contains the 6-dimensional groups

//Let G be a type F genus 2 group. There are 3 types of groups of type L. Those of the form:
//1) U(1)xG are denoted <u,G>.
//2) N(U(1))xG are denoted <nu,G>.
//3) N(U(1))x_{C_2} G are denoted <H,G>, where H is the label of the kernel of the hom G-->C_2. 
//We build the 64 direct products

small_list:=genus2_F(K)[1];
small_dict:=genus2_F(K)[2];

N6:=embed6(MK2!Matrix(K,2,2,[0,-1,1,0]),MK4!DiagonalMatrix(K,[1,1,1,1]));

for i in small_list do
	Name:=<"u",i>;
	L:=[embed6(MK2!DiagonalMatrix(K,[1,1]),MK4!g): g in small_dict[i]];
	big_dict[Name]:=sub<GK6|L>;
	big_list:=Append(big_list,Name);
end for;

for i in small_list do
	Name:=<"nu",i>; 
	iL:=<"u",i>;
	big_dict[Name]:=sub<GK6|big_dict[iL], N6>;
	big_list:=Append(big_list,Name);
end for;

//Next we build the 73 fiber products. These are stored in "fiber_products". In "big_list" and "big_dict" we only store the 58 giving rise to distinct multisets of characteristic polynomials

fiber_products:=[];

//Given a genus 2 group G, returns the multiset of characteristic polynomials of G (this uniquely identifies the group G)
char_polys_group:=function(G)
	return {* CharacteristicPolynomial(U4*(Parent(U4)!g)) : g in G *};
end function;

//We build a dictionary that stores the chracteristic polynomials of genus 2 groups
small_char_polys:=AssociativeArray();
for i in small_list do
	small_char_polys[i]:=char_polys_group(small_dict[i]);
end for;

//Given a genus 2 group G that contains mI4, finds the index 2 subgroups that contain mI4
index_2_subgroups:=function(G)
	L:=MaximalSubgroups(G);
	LL:=[];
	n:=Order(G);
	for l in L do
		H:=l`subgroup;
		if Order(H) eq n/2 and (mI4 in H) then LL:=Append(LL,H); end if;
	end for;
	return LL;
end function;

//Given a genus 2 group H, determines its name
identify_group:=function(H)
	char_polys_H:=char_polys_group(H);
	for j in Keys(small_char_polys) do
		if char_polys_H eq small_char_polys[j] then name:=j; end if;
	end for;
	return name;
end function;

for str in small_list do
	G:=small_dict[str];
	LL:=index_2_subgroups(G);
	for i in [1..#LL] do
		H:=LL[i];
		name:=identify_group(H);
		fiber_name:=<name,str>;
		elements_H:=[G!h : h in H];
		gens_fiber_prod:=[];
		for g in G do
			if not g in elements_H then
				g6:=embed6(MK2!Matrix(K,2,2,[0,-1,1,0]),MK4!g); 
				gens_fiber_prod:=Append(gens_fiber_prod,g6); 
			end if;
		end for;
		if not(fiber_name in big_list) then 
			big_dict[fiber_name]:=sub<GK6|gens_fiber_prod>; 
			big_list:=Append(big_list,fiber_name);
		end if;
		fiber_products:=Append(fiber_products,sub<GK6|gens_fiber_prod>);
	end for;
end for;

groups_list:=[];
for i in big_list do
	groups_list:=Append(groups_list,big_dict[i]);
end for;

//We set tex names

gn:=[];
for i in [1..#big_list] do
	na:=big_list[i];
	if na[1] eq "u" then name:="L_1(" cat na[2] cat ")";
	elif na[1] eq "nu" then name:="L_2(" cat na[2] cat ")";
	else name:="L(" cat na[2] cat "," cat na[1] cat ")";
	end if;
	gn:=Append(gn,name);
end for;

//We finally generate "names_L.txt", "moments_L.txt", "diagonals_L.txt", "zvectors_L.txt", "labels_L.txt", "a1a2a3_L.txt", "comp_diagonals_L.txt"

generate_data("L", gn, groups_list,true,true,true,true,true,true);

//We next generate the subgroup lattice
load "diagonals_L.txt";
gdi:=L;
load "labels_L.txt";
glab:=L;
load "names_L.txt";
gn:=L;
if compute_lattice then
	glat:=subgroups_lattice(groups_list,gdi,gn,glab);
	print_in_file(glat,"lattice_L.txt");
end if;

//Confirmation that the fiber products with respect to conjugate subgroups in USp(4) give the same moments
if verify_fiber_products eq true then

	//Given a genus 3 group G, returns the multiset of characteristic polynomials of G (this uniquely identifies the group G)
	char_polys_group6:=function(G)
		return {* CharacteristicPolynomial(U6*Parent(U6)!g) : g in G *};
	end function;

	L:=[];
	for G in fiber_products do
		cpols:=char_polys_group6(G);
		if not cpols in L then L:=Append(L,cpols); end if;
	end for;
	assert #L eq 58;
end if;

//Confirmation that the fiber product of N(U(1)) and J(D_6) along J(C_6) and the fiber product of N(U(1)) and J(D_6) along D_6 have distinct 14-simplex of moments
if verify_distinguished_by_moments then
	
	moment_ijk_G:=function(G,i,j,k)
		moment:=0;
		o:=0;
		for g in G do
			p:=CharacteristicPolynomial(Parent(U6)!g*U6);
			p2:=Coefficients(p);
			p3:=[p2[2],p2[3],p2[4]];
			moment:=moment+momentsp(p3,i,j,k);
			o:=o+1;
		end for;
		return moment/o;
	end function;

	cc:=hom_simplex(14);

	G1:=big_dict[<"J(C_6)","J(D_6)">];
	G2:=big_dict[<"D_6","J(D_6)">];
	for c in cc do
		if not (moment_ijk_G(G1,c[1],c[2],c[3]) eq moment_ijk_G(G2,c[1],c[2],c[3])) then 
			print c;
		end if;
	end for;
end if;


