//This file builds the groups of type M.
//It then produces the files:
// names_M.txt
// moments_M.txt
// diagonals_M.txt
// labels_M.txt
// zvectors_M.txt
// a1a2a3_M.txt
// lattice_M.txt
// comp_diagonals_M.txt

N := 24; K:= CyclotomicField(N);
zz:=K.1;

zet:=AssociativeArray();
for i in [1,2,3,4,6,8,12] do
	zet[i]:=zz^(Integers()!(N/i)); 
end for;

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

U6:=Matrix(FR,6,6,[[u,0,0,0,0,0],[0,u,0,0,0,0],[0,0,u,0,0,0],[0,0,0,1/u,0,0],[0,0,0,0,1/u,0],[0,0,0,0,0,1/u]]);
mI:=GK6!DiagonalMatrix(K,[-1,-1,-1,-1,-1,-1]);

label_group:=function(G)
	if (mI in G) then return IdentifyGroup(quo<G|mI>); end if;
	return IdentifyGroup(G);
end function;

coef_var:=function(p,var,d)
	dd:=Degree(Denominator(p),var);
	cd:=Coefficient(Denominator(p), var, dd);
	cn:=Coefficient(Numerator(p),var,d+dd);
	return cn/cd;
end function;

avervalue:=function(q)
	return coef_var(q,u,0)-coef_var(q,u,2);
end function;

//We can call "common_functions.m" only once "avervalue" and "label_group" have been defined.

load "common_functions.m";

big_list:=[];//This list will contain the names of the groups
big_dict:=AssociativeArray();//This contains the 6-dimensional groups

//We construct the 11 groups

c_mat:=AssociativeArray();
for i in [1,2,3,4,6] do
	c_mat[i]:=Matrix(K,2,2,[(zet[i]+ComplexConjugate(zet[i]))/2,(zet[i]-ComplexConjugate(zet[i]))/(2*zet[4]),(-zet[i]+ComplexConjugate(zet[i]))/(2*zet[4]),(zet[i]+ComplexConjugate(zet[i]))/2]);
end for;

for i in [1,2,3,4,6] do
	name:=<"c",i>;
	big_dict[name]:=sub<GK6|embed26(c_mat[i])>;
	big_list:=Append(big_list,name);
end for;

s_mat:=Matrix(K,2,2,[1,0,0,-1]);

for i in [2,3,4,6] do
	name:=<"d",i>;
	big_dict[name]:=sub<GK6|embed26(c_mat[i]),embed26(s_mat)>;
	big_list:=Append(big_list,name);
end for;

g1:=Matrix(K,3,3,[[1,0,0],[1,-1,1],[0,0,1]]);

g2:=Matrix(K,3,3,[[1,0,0],[1,0,-1],[0,1,-1]]);

g3:=Matrix(K,3,3,[[0,0,-1],[1,0,-1],[0,1,-1]]);

g4:=Matrix(K,3,3,[[-1,1,0],[0,1,0],[0,1,-1]]);

big_dict[<"a",4>]:=sub<GK6|embed(g2),embed(g4)>;

big_dict[<"s",4>]:=sub<GK6|embed(g1),embed(g3)>;

big_list:=Append(big_list,<"a",4>);
big_list:=Append(big_list,<"s",4>);

groups_list:=[];
for i in big_list do
	groups_list:=Append(groups_list,big_dict[i]);
end for;


//We set tex names

gn:=[];
for i in [1..#big_list] do
	na:=big_list[i];
	if na[1] eq "c" then str1:="C_";
	elif na[1] eq "d" then str1:="D_";
	elif na[1] eq "a" then str1:="A_";
	elif na[1] eq "s" then str1:="S_";
	end if;
	name:="M(" cat str1 cat IntegerToString(na[2]) cat ")";
	gn:=Append(gn,name);
end for;

//We finally generate "names_M.txt", "moments_M.txt", "diagonals_M.txt", "zvectors_M.txt", "labels_M.txt", "a1a2a3_M.txt", "comp_diagonals_M.txt"

generate_data("M", gn, groups_list, true,true,true,true, true,true);

//We next generate the subgroup lattice
load "diagonals_M.txt";
gdi:=L;
load "labels_M.txt";
glab:=L;
load "names_M.txt";
gn:=L;

glat:=subgroups_lattice(groups_list,gdi,gn,glab);
print_in_file(glat,"lattice_M.txt");


