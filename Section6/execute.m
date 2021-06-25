//This file generates all the data ".txt" files. If the files "labels_*.txt", "names_*.txt", "moments_*.", "diagonals_*.txt", "zvectors_*.txt", "dvectors_*.txt", "lattice_*.txt" for * in [A,...,N] have been generated, one can comment the text comprised between //11 ...//22, and this file will create the files "labels_total.txt", "names_total.txt", "moments_total.txt", "zvectors_total.txt", "dvectors_total.txt", "lattice_total.txt".


load "type_A.m";
load "type_B.m";
load "type_C.m";
load "type_D.m";
load "type_E.m";
load "type_F.m";
load "type_G.m";
load "type_H.m";
load "type_I.m";
load "type_J.m";
load "type_K.m";
load "type_L.m";
load "type_M.m";
load "type_N.m";




print_in_file:=procedure(g_list,file_name)

	SetOutputFile(  file_name: Overwrite:=true);
	print "L:=[";
	for i in [1..#g_list] do
		print g_list[i];
		if not i eq #g_list then print ","; end if;
	end for;
	print "];";
	UnsetOutputFile();

end procedure;

print_in_file_str:=procedure(g_list,file_name)//This one is used to print strings

	SetOutputFile(  file_name: Overwrite:=true);
	print "L:=[";
	for i in [1..#g_list] do
		print "\"" cat g_list[i] cat "\"";
		if not i eq #g_list then print ","; end if;
	end for;
	print "];";
	UnsetOutputFile();

end procedure;

print_in_file_seq_str:=procedure(g_list,file_name)//This one is used to print strings

	SetOutputFile(  file_name: Overwrite:=true);
	print "L:=[";
	for s in [1..#g_list] do
		seq:=g_list[s];		
		str:="[";
		for j in [1..#seq] do
			str:=str cat "\"" cat seq[j] cat "\"";
			if not j eq #seq then str:=str cat ","; end if;
		end for;
		str:=str cat "]";
		print str;
		if not s eq #g_list then print ","; end if;
	end for;
	print "];";
	UnsetOutputFile();

end procedure;

hom_simplex:=function(m)
	hs:=[];
	for a in [0..m] do	
	for i in [0..a] do
		for j in [0..a-i] do
			for k in [0..m-i-2*j] do
				if not(i+j+k eq 0) and i+2*j+3*k eq a and (i+k) mod 2 eq 0 then 
				hs:=Append(hs,[i,j,k]);
				end if;
			end for;
		end for;
	end for;
	end for;
	return hs;
end function;

names:=[];
load "names_A.txt";
names:=names cat L;
load "names_B.txt";
names:=names cat L;
load "names_C.txt";
names:=names cat L;
load "names_D.txt";
names:=names cat L;
load "names_E.txt";
names:=names cat L;
load "names_F.txt";
names:=names cat L;
load "names_G.txt";
names:=names cat L;
load "names_H.txt";
names:=names cat L;
load "names_I.txt";
names:=names cat L;
load "names_J.txt";
names:=names cat L;
load "names_K.txt";
names:=names cat L;
load "names_L.txt";
names:=names cat L;
load "names_M.txt";
names:=names cat L;
load "names_N.txt";
names:=names cat L;

print_in_file_str(names, "names_total.txt");

moments:=[];
load "moments_A.txt";
moments:=moments cat L;
load "moments_B.txt";
moments:=moments cat L;
load "moments_C.txt";
moments:=moments cat L;
load "moments_D.txt";
moments:=moments cat L;
load "moments_E.txt";
moments:=moments cat L;
load "moments_F.txt";
moments:=moments cat L;
load "moments_G.txt";
moments:=moments cat L;
load "moments_H.txt";
moments:=moments cat L;
load "moments_I.txt";
moments:=moments cat L;
load "moments_J.txt";
moments:=moments cat L;
load "moments_K.txt";
moments:=moments cat L;
load "moments_L.txt";
moments:=moments cat L;
load "moments_M.txt";
moments:=moments cat L;
load "moments_N.txt";
moments:=moments cat L;

print_in_file(moments, "moments_total.txt");

diagonals:=[];
load "diagonals_A.txt";
diagonals:=diagonals cat L;
load "diagonals_B.txt";
diagonals:=diagonals cat L;
load "diagonals_C.txt";
diagonals:=diagonals cat L;
load "diagonals_D.txt";
diagonals:=diagonals cat L;
load "diagonals_E.txt";
diagonals:=diagonals cat L;
load "diagonals_F.txt";
diagonals:=diagonals cat L;
load "diagonals_G.txt";
diagonals:=diagonals cat L;
load "diagonals_H.txt";
diagonals:=diagonals cat L;
load "diagonals_I.txt";
diagonals:=diagonals cat L;
load "diagonals_J.txt";
diagonals:=diagonals cat L;
load "diagonals_K.txt";
diagonals:=diagonals cat L;
load "diagonals_L.txt";
diagonals:=diagonals cat L;
load "diagonals_M.txt";
diagonals:=diagonals cat L;
load "diagonals_N.txt";
diagonals:=diagonals cat L;

print_in_file(diagonals, "diagonals_total.txt");

comp_diagonals:=[];
load "comp_diagonals_A.txt";
comp_diagonals:=comp_diagonals cat L;
load "comp_diagonals_B.txt";
comp_diagonals:=comp_diagonals cat L;
load "comp_diagonals_C.txt";
comp_diagonals:=comp_diagonals cat L;
load "comp_diagonals_D.txt";
comp_diagonals:=comp_diagonals cat L;
load "comp_diagonals_E.txt";
comp_diagonals:=comp_diagonals cat L;
load "comp_diagonals_F.txt";
comp_diagonals:=comp_diagonals cat L;
load "comp_diagonals_G.txt";
comp_diagonals:=comp_diagonals cat L;
load "comp_diagonals_H.txt";
comp_diagonals:=comp_diagonals cat L;
load "comp_diagonals_I.txt";
comp_diagonals:=comp_diagonals cat L;
load "comp_diagonals_J.txt";
comp_diagonals:=comp_diagonals cat L;
load "comp_diagonals_K.txt";
comp_diagonals:=comp_diagonals cat L;
load "comp_diagonals_L.txt";
comp_diagonals:=comp_diagonals cat L;
load "comp_diagonals_M.txt";
comp_diagonals:=comp_diagonals cat L;
load "comp_diagonals_N.txt";
comp_diagonals:=comp_diagonals cat L;

print_in_file(comp_diagonals, "comp_diagonals_total.txt");

zvectors:=[];
load "zvectors_A.txt";
zvectors:=zvectors cat L;
load "zvectors_B.txt";
zvectors:=zvectors cat L;
load "zvectors_C.txt";
zvectors:=zvectors cat L;
load "zvectors_D.txt";
zvectors:=zvectors cat L;
load "zvectors_E.txt";
zvectors:=zvectors cat L;
load "zvectors_F.txt";
zvectors:=zvectors cat L;
load "zvectors_G.txt";
zvectors:=zvectors cat L;
load "zvectors_H.txt";
zvectors:=zvectors cat L;
load "zvectors_I.txt";
zvectors:=zvectors cat L;
load "zvectors_J.txt";
zvectors:=zvectors cat L;
load "zvectors_K.txt";
zvectors:=zvectors cat L;
load "zvectors_L.txt";
zvectors:=zvectors cat L;
load "zvectors_M.txt";
zvectors:=zvectors cat L;
load "zvectors_N.txt";
zvectors:=zvectors cat L;

print_in_file(zvectors, "zvectors_total.txt");

labels:=[];
load "labels_A.txt";
labels:=labels cat L;
load "labels_B.txt";
labels:=labels cat L;
load "labels_C.txt";
labels:=labels cat L;
load "labels_D.txt";
labels:=labels cat L;
load "labels_E.txt";
labels:=labels cat L;
load "labels_F.txt";
labels:=labels cat L;
load "labels_G.txt";
labels:=labels cat L;
load "labels_H.txt";
labels:=labels cat L;
load "labels_I.txt";
labels:=labels cat L;
load "labels_J.txt";
labels:=labels cat L;
load "labels_K.txt";
labels:=labels cat L;
load "labels_L.txt";
labels:=labels cat L;
load "labels_M.txt";
labels:=labels cat L;
load "labels_N.txt";
labels:=labels cat L;

print_in_file(labels, "labels_total.txt");


a1a2a3:=[];
load "a1a2a3_A.txt";
a1a2a3:=a1a2a3 cat L;
load "a1a2a3_B.txt";
a1a2a3:=a1a2a3 cat L;
load "a1a2a3_C.txt";
a1a2a3:=a1a2a3 cat L;
load "a1a2a3_D.txt";
a1a2a3:=a1a2a3 cat L;
load "a1a2a3_E.txt";
a1a2a3:=a1a2a3 cat L;
load "a1a2a3_F.txt";
a1a2a3:=a1a2a3 cat L;
load "a1a2a3_G.txt";
a1a2a3:=a1a2a3 cat L;
load "a1a2a3_H.txt";
a1a2a3:=a1a2a3 cat L;
load "a1a2a3_I.txt";
a1a2a3:=a1a2a3 cat L;
load "a1a2a3_J.txt";
a1a2a3:=a1a2a3 cat L;
load "a1a2a3_K.txt";
a1a2a3:=a1a2a3 cat L;
load "a1a2a3_L.txt";
a1a2a3:=a1a2a3 cat L;
load "a1a2a3_M.txt";
a1a2a3:=a1a2a3 cat L;
load "a1a2a3_N.txt";
a1a2a3:=a1a2a3 cat L;

print_in_file(a1a2a3, "a1a2a3_total.txt");

lattice:=[];
load "lattice_A.txt";
lattice:=lattice cat L;
load "lattice_B.txt";
lattice:=lattice cat L;
load "lattice_C.txt";
lattice:=lattice cat L;
load "lattice_D.txt";
lattice:=lattice cat L;
load "lattice_E.txt";
lattice:=lattice cat L;
load "lattice_F.txt";
lattice:=lattice cat L;
load "lattice_G.txt";
lattice:=lattice cat L;
load "lattice_H.txt";
lattice:=lattice cat L;
load "lattice_I.txt";
lattice:=lattice cat L;
load "lattice_J.txt";
lattice:=lattice cat L;
load "lattice_K.txt";
lattice:=lattice cat L;
load "lattice_L.txt";
lattice:=lattice cat L;
load "lattice_M.txt";
lattice:=lattice cat L;
load "lattice_N.txt";
lattice:=lattice cat L;

print_in_file_seq_str(lattice, "lattice_total.txt");



//We next find the groups that share the same 12 simplex
 
conflated_groups:=[];
seen:=[];
for i in [1..#moments] do
	for j in [1..i-1] do
	if moments[i] eq moments[j] then
		conflated_groups:=Append(conflated_groups,names[j]);
		conflated_groups:=Append(conflated_groups,names[i]);
		print j,i; 
	end if;
	end for;
end for;
 
print "The groups conflated by the 12-homogeneous simplex of moments are:";
print conflated_groups;

//We next compute the orthogonality relations for all groups using the 12-simplex of moments

PA:=PolynomialRing(Rationals(),3);

a1:=PA.1;
a2:=PA.2;
a3:=PA.3;

//We need to change the sign of a1 and a3 in the expressions obtained from "Characters_in_terms_of_coef.m" since the current convention is that the characteristic polynomial is 1+a1*T+a2*T^2+a3*T^3+a2*T^4+a1*T^5+T^6
chi:=AssociativeArray();
chi[1]:=PA!1;
chi[2]:=-a1;
chi[3]:=a2-1;
chi[4]:=-a3+a1;
chi[5]:=a1^2-a2;
chi[6]:=-a1*a2+a1+a3;
chi[7]:=a1*a3-a1^2-a2+1;
chi[8]:=a2^2-a1*a3-a2;
chi[9]:=-a2*a3+2*a1*a2-a1;
chi[10]:=a3^2-a2^2-a1*a3+2*a2-1;
chi[11]:=-a1^3 + 2*a1*a2 - a3;
chi[12]:=a1^2*a2 - a1^2 - a1*a3 - a2^2 + 2*a2;
chi[13]:=a1^3 - a1^2*a3 - 2*a1 + a2*a3;
chi[14]:=a1^2*a3 - a1*a2^2 + a2*a3 - a3;
chi[15]:=-2*a1^2*a2 + 2*a1^2 + a1*a2*a3 + a1*a3 - a3^2;

charN:=#Keys(chi);
cc:=hom_simplex(12);

ortrels_group:=function(g)
	ortrellist:=<>;
	
	for i in [1..charN] do
	for j in [1..charN] do
		ortrel:=Coefficient(Coefficient(Coefficient(chi[i]*chi[j], a1, 0), a2, 0), a3, 0);
		
		for c in [1..#cc] do
			e1:=cc[c][1];
			e2:=cc[c][2];
			e3:=cc[c][3];
			coef:=Coefficient(Coefficient(Coefficient(chi[i]*chi[j], a1, e1), a2, e2), a3, e3);
			ortrel:=ortrel+coef*moments[g][c];
		end for;
		ortrellist:=Append(ortrellist,ortrel);
	end for;
	end for;
	
	return ortrellist;
end function;

ortrels:=[];
for g in [1..#moments] do
	ortrels:=Append(ortrels,ortrels_group(g));
end for;



print_in_file(ortrels, "orthogonality_relations_total.txt");







