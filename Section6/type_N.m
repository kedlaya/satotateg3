//This file builds the groups of type N.
//It then produces the files:
// names_N.txt
// moments_N.txt
// diagonals_N.txt
// labels_N.txt
// zvectors_N.txt
// a1a2a3_N.txt
// lattice_N.txt
// comp_diagonals_N.txt

compute_lattice:=true;//If true, the code additionally computes the subgroup lattice

N := 504; K:= CyclotomicField(N);
z:=K.1;

z3:=z^(Integers()!(N/3)); 
z4:=z^(Integers()!(N/4)); 
z7:=z^(Integers()!(N/7)); 
z8:=z^(Integers()!(N/8)); 
z9:=z^(Integers()!(N/9));

GK2:=GeneralLinearGroup(2,K);
MK2:=MatrixRing(K,2);
GK3:=GeneralLinearGroup(3,K);
MK3:=MatrixRing(K,3);
GK4:=GeneralLinearGroup(4,K);
MK4:=MatrixRing(K,4);
GK6:=GeneralLinearGroup(6,K);
MK6:=MatrixRing(K,6);

PR:=PolynomialRing(K,1);
u:=PR.1;
FR:=FieldOfFractions(PR);
AssignNames(~FR,["u"]);

//We have empirically verified that working over the Laurent series ring is slightly faster than working over the the field of fractions of K[u]. For types with more than one variable however, we also empirically verified that working over the field of fractions is much faster than working over a multivariate Laurent series ring.
//FR:=LaurentSeriesRing(K);
//u:=FR.1;

U6:=DiagonalMatrix(FR,[u,u,u,1/u,1/u,1/u]);

dd3:=GK6!DiagonalMatrix(K,[z^(Integers()!(N/3)),z^(Integers()!(N/3)),z^(Integers()!(N/3)),z^(Integers()!(-N/3)),z^(Integers()!(-N/3)),z^(Integers()!(-N/3))]);
dd2:=GK6!DiagonalMatrix(K,[-1,-1,-1,-1,-1,-1]);

label_group:=function(G)
	L:=[];
	for l in [dd2, dd3] do	
		if (l in G) then L:=Append(L,l); end if;
	end for;	
	Q :=quo<G|L>;	
	return IdentifyGroup(Q);
end function;

coef_var:=function(p,var,d)
	dd:=Degree(Denominator(p),var);
	cd:=Coefficient(Denominator(p), var, dd);
	cn:=Coefficient(Numerator(p),var,d+dd);
	return cn/cd;
end function;

avervalue:=function(q)
	//return Coefficient(q,0); In case one wanted to work of the the Laurent series ring
	return coef_var(q,u,0);
end function;

//We can call "common_functions.m" only once "avervalue" and "label_group" have been defined.

load "common_functions.m";

generators_c:=AssociativeArray();//Dictionary of the 3x3 generators of the groups in the centralizer
generators:=AssociativeArray();//Dictionary of the generators of the groups
groups_names_c:=[];//List of the names of the groups in the centralizer
groups_names:=[];//List of the names of the groups
  
verify_rationality:=function(G)
	for g in G do
		t:=Trace(g);
		if not t*ComplexConjugate(t) in IntegerRing() then return false; end if;
	end for;
	return true;
end function;

union_dicts:=function(L)
	union:=AssociativeArray();
	for d in L do
		for i in Keys(d) do
			union[i]:=d[i];
		end for;
	end for;
	return union;
end function;

A_str:=function(i)
	if #i eq 2 then Astr:="A(" cat IntegerToString(i[1])cat "," cat IntegerToString(i[2]) cat ")"; end if;
	if #i eq 3 then Astr:="A(" cat IntegerToString(i[1]) cat "," cat IntegerToString(i[2]) cat ")_" cat IntegerToString(i[3]); end if;
	return Astr;
end function;

B_str:=function(l)
	if l[2][1] eq 1 and #l[1] eq 2 then
		Bstr:="B(" cat IntegerToString(l[1][1])cat "," cat IntegerToString(l[1][2]) cat ")";
	elif l[2][1] eq 1 then
		Bstr:="B(" cat IntegerToString(l[1][1]) cat "," cat IntegerToString(l[1][2]) cat ")_" cat IntegerToString(l[1][3]);
	elif #l[1] eq 2 then
		Bstr:="B(" cat IntegerToString(l[1][1])cat "," cat IntegerToString(l[1][2]) cat ";" cat IntegerToString(l[2][1]) cat ")";
	else
		Bstr:="B(" cat IntegerToString(l[1][1]) cat "," cat IntegerToString(l[1][2]) cat ";" cat IntegerToString(l[2][1]) cat ")_" cat IntegerToString(l[1][3]);
	end if;
	return Bstr;
end function;

A_data:= [
[<1/3,1/3,1/3>],
[<1/3,5/6,5/6>],
[<1/9,4/9,4/9>],
[<1/6,5/12,5/12>],
[<1/3,1/12,7/12>],
[<1/9,17/18,17/18>],
[<1/18,2/9,13/18>],
[<1/21,16/21,4/21>],
[<1/6,1/24,19/24>],
[<1/12,5/24,17/24>],
[<1/9, 7/36, 25/36>],
[<0,1/2,1/2>, <1/6,1/6,2/3>],
[<1/2,0,1/2>, <1/6,5/12,5/12>],
[<1/2,0,1/2>, <1/9,17/18,17/18>],
[<1/3,0,2/3>, <1/3,1/3,1/3>],
[<1/3,0,2/3>, <2/3,1/6,1/6>],
[<0,1/3,2/3>, <1/9,1/9,7/9>],
[<0,1/3,2/3>, <1/6,5/12,5/12>],
[<0,1/3,2/3>, <1/9,5/18,11/18>],
[<0,1/4,3/4>, <1/12,1/12,5/6>],
[<1/6,0,5/6>, <2/3,1/6,1/6>],
[<0,1/6,5/6>, <1/18,1/18,8/9>]
];

A_keys:=<<1,1>, <1,2>, <1,3>, <1,4,1>, <1,4,2>, <1,6,1>, <1,6,2>, <1,7>, <1,8,1>, <1,8,2>, <1,12>, <2,2>, <2,4>, <2,6>, 
<3,1>, <3,2>, <3,3>, <3,4>, <3,6>, <4,4>, <6,2>, <6,6>>;
 
for i in [1..#A_keys] do
	Astr:=A_str(A_keys[i]);
	L:=[];
	for m in A_data[i] do
		D:=GK3!DiagonalMatrix(K,[z^(Integers()!(N*m[1])),z^(Integers()!(N*m[2])),z^(Integers()!(N*m[3]))]);
		L:=Append(L,D);
	end for;
	generators_c[Astr]:=L;
	groups_names_c:=Append(groups_names_c,Astr);
end for;


T_mat:=AssociativeArray();
for i in [1,2,4] do
	T_mat[i]:=Matrix(K,3,3,[-z^(-Integers()!(N/i)),0,0,0,0,z^(Integers()!(N/(2*i))),0,z^(Integers()!(N/(2*i))),0]);
end for;

B_keys:=<<<1,4,2>, <1>>, <<1,8,1>, <1>>, <<1,12>, <1>>,
          <<2,4>, <1>>, <<3,1>, <1>>, <<3,2>, <1>>, <<3,3>, <1>>, <<3,4>, <1>>,
          <<3,6>, <1>>, <<4,4>, <1>>, <<6,2>, <1>>, <<6,6>, <1>>,
          <<1,4,2>, <2>>, <<1,12>, <2>>, <<3,2>, <2>>, <<3,6>, <2>>,
          <<2,4>, <4>>, <<3,4>, <4>>>;


for l in B_keys do
	Astr:=A_str(l[1]);
	Bstr:=B_str(l);
	generators_c[Bstr]:=Append(generators_c[Astr], T_mat[l[2][1]]);
	groups_names_c:=Append(groups_names_c,Bstr);
end for;

BT_gens := [Matrix(K,3,3,[1,0,0,0,0,1,0,-1,0]),
           Matrix(K,3,3,[1,0,0,0,1/2+z4/2,1/2+z4/2,0,-1/2+z4/2,1/2-z4/2])];

for n in [1,2,3] do
	uu:= z^(Integers()!(N/(6*n)));
	d:= DiagonalMatrix(K,[uu^(-2),uu,uu]);
	BTstr:="B(T," cat IntegerToString(n) cat ")";	
	generators_c[BTstr]:= BT_gens cat [d];
	groups_names_c:=Append(groups_names_c,BTstr);
	if n eq 3 then
	        generators_c["B(T,1;1)"]:=[BT_gens[1], BT_gens[2]*d];
		groups_names_c:=Append(groups_names_c,"B(T,1;1)");	
	end if;	
end for;

for n in [1,2] do
	uu:= z^(Integers()!(N/(12*n)));
	d:=DiagonalMatrix(K,[uu^2,uu^(-1)*z8,uu^(-1)/z8]);
	BOstr:="B(O," cat IntegerToString(n) cat ")";	
	generators_c[BOstr]:= BT_gens cat [d];
	groups_names_c:=Append(groups_names_c,BOstr);
end for;


C_keys := <<1,7>, <2,2>, <3,1>, <3,3>, <4,4>, <6,2>, <6,6>>;
S_mat := Matrix(K,3,3,[0,0,1,1,0,0,0,1,0]);

for i in C_keys do
	Cstr:="C(" cat IntegerToString(i[1])cat "," cat IntegerToString(i[2]) cat ")";
	generators_c[Cstr]:=generators_c[A_str(i)] cat [S_mat];
	groups_names_c:=Append(groups_names_c,Cstr);
end for;

D_keys :=<<2,2>, <3,1>, <3,3>, <4,4>, <6,2>, <6,6>>;
for i in D_keys do
	Dstr:="D(" cat IntegerToString(i[1])cat "," cat IntegerToString(i[2]) cat ")";     
	generators_c[Dstr]:= generators_c[A_str(i)] cat [S_mat, T_mat[1]];
	groups_names_c:=Append(groups_names_c,Dstr);
end for;

g1 := DiagonalMatrix(K,[1, z3, z3^2]);
g2 := 1/(z3-z3^2)*Matrix(K,3,3,[1,1,1,1,z3,z3^2,1,z3^2,z3]);
g3 := DiagonalMatrix(K,[z9^2, z9^2, z9^5]);
h := (z7+z7^2+z7^4-z7^3-z7^5-z7^6)/7; //equals sqrt(-7)
a:=z7^4-z7^3;
b:=z7^2-z7^5;
c:=z7-z7^6;
g4:=DiagonalMatrix(K,[z7,z7^2,z7^4]);
g5:=h*Matrix(K,3,3,[a,b,c,b,c,a,c,a,b]);
g6:=DiagonalMatrix(K,[z3,z3,z3]);

generators_c["E(36)"]:=[g1, g2];
generators_c["E(72)"]:=[g1, g2, g3*g2*g3^(-1)];
generators_c["E(216)"]:=[g1, g2, g3];
generators_c["E(168)"]:=[S_mat, g4,g5,g6];

groups_names_c:=groups_names_c cat ["E(36)","E(72)","E(216)","E(168)"];

for i in groups_names_c do
	L:=generators_c[i];
	LL:=[];
	for g in L do
		LL:=Append(LL,GK6!embed(MK3!g));
	end for;
	generators[i]:=LL;
end for;

prime_data_A := ["A(1,4)_2", "A(1,8)_1", "A(1,8)_2", "A(1,12)",
               "A(2,2)", "A(2,4)", "A(2,6)", "A(3,1)", "A(3,2)",
               "A(3,3)", "A(3,4)", "A(3,6)", "A(4,4)", "A(6,2)", "A(6,6)",
               "C(2,2)", "C(3,3)", "C(4,4)", "C(6,2)", "C(6,6)"];

prime_data_B:= ["B(3,2)", "B(3,4)", "B(3,6)",
                "B(3,2;2)", "B(3,6;2)", "B(3,4;4)"];

prime_data_BT := ["B(1,4)_2", "B(1,12)",
                 "B(2,4)", "B(1,4;2)_2", "B(1,12;2)", "B(2,4;4)",
                 "B(T,1)", "B(T,2)", "B(T,3)", "B(T,1;1)"];

J:=GK6!BlockMatrix(2,2,[0,IdentityMatrix(K,3),-IdentityMatrix(K,3),0]);
groups_names:=groups_names_c;

for i in groups_names_c do
	if not i eq "B(T,1;1)" then
		Jname := "J(" cat i cat ")";
		generators[Jname]:=generators[i] cat [J];
		groups_names:=Append(groups_names,Jname);
	end if;
	if i in prime_data_A then
		Jname:= "J_s(" cat i cat ")";
		generators[Jname]:=generators[i] cat [J*(GK6!embed(MK3!T_mat[1]))];
		groups_names:=Append(groups_names,Jname);
	elif i in prime_data_B then
		Jname:= "J_s(" cat i cat ")";
		generators[Jname]:=generators[i] cat [J*(GK6!embed(MK3!DiagonalMatrix(K,[-1,1,-1])))];
		groups_names:=Append(groups_names,Jname);
	elif i in prime_data_BT then
		Jname:= "J_s(" cat i cat ")";
		generators[Jname]:=generators[i] cat [J*(GK6!embed(MK3!DiagonalMatrix(K,[z4,z4,-1])))];
		groups_names:=Append(groups_names,Jname);
	end if;
end for;                                         

prime_data_A2a:= ["A(1,2)", "A(1,4)_1", "A(1,6)_1", "A(3,2)", "A(3,4)", "A(3,6)"];
prime_data_A2b:= ["A(1,4)_2", "A(1,12)", "A(2,4)"];

for i in groups_names_c do
	if i in prime_data_A2a then
		Jname := "J_n(" cat i cat ")";
		generators[Jname]:=generators[i] cat [J*(GK6!embed(Matrix(K,3,3,[1,0,0,0,0,1,0,-1,0])))];
		groups_names:=Append(groups_names,Jname);
	elif i in prime_data_A2b then
		Jname:= "J_n(" cat i cat ")";
		generators[Jname]:=generators[i] cat [J*(GK6!embed(Matrix(K,3,3,[z4,0,0,0,0,z4,0,1,0])))];
		groups_names:=Append(groups_names,Jname);
	elif i eq "E(36)" then
		Jname:= "J_n(" cat i cat ")";
		generators[Jname]:=generators[i] cat [J*(GK6!embed(MK3!(g3*g2*g3^(-1))))];
		groups_names:=Append(groups_names,Jname);
	end if;
end for;

groups_list:=[];
for i in groups_names do
	G := sub<GK6|generators[i]>;
	groups_list:=Append(groups_list,G);
end for;

gn:=groups_names;
gg:=[];
for i in gn do
	gg:=Append(gg,generators[i]);
end for;

//We finally generate "names_N.txt", "moments_N.txt", "diagonals_N.txt", "zvectors_N.txt", "labels_N.txt", "a1a1a3_N.txt", "comp_diagonals_N.txt"

generate_data("N", gn, groups_list, true, true, true, true, true, true);

//We next generate the subgroup lattice
load "diagonals_N.txt";
gdi:=L;
load "labels_N.txt";
glab:=L;
load "names_N.txt";
gn:=L;
if compute_lattice then
	glat:=subgroups_lattice(groups_list,gdi,gn,glab);
	print_in_file(glat,"lattice_N.txt");
end if;



