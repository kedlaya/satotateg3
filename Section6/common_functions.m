//This file contains functions used by several of the files type_*.m

//We start by defining some auxilliary functions

//Embedding that sends M to [[M,0],[0,Mbar]].

conj_matrix:=function(M)
	r:=Ncols(M);
	return Matrix(K, r, r, [ComplexConjugate(M[i][j]): i, j in [1..r]]);
end function;

embed:=function(M)
	return BlockMatrix(2,2,[M,0,0,conj_matrix(M)]);
end function;


//Embedding used in type M
embed36:=function(M)
	return BlockMatrix(3,3, [ MK2!DiagonalMatrix(K,[ M[i][j], M[i][j]]) : i,j in [1,2,3]]);
end function;

embed23_aux:=function(M,i,j)
	if i in [2,3] and j in [2,3] then
		return M[i-1][j-1];
	end if;
	if i eq 1 and j eq 1 then return 1; end if;
	return 0;
end function;

embed26:=function(M)
	mat:=Matrix(K,3,3, [ embed23_aux(M,i,j) : i,j in [1,2,3]] );
	return GK6!embed(mat);
end function;




//Embedding that sends a pair of matrices (M2,M4) to [[M2,0],[0,M4]].
 
aux_embed6:=function(i,j,M2,M4)
	if i in [1,2] and j in [1,2] then return M2[i][j]; end if;
	if i in [3,4,5,6] and j in [3,4,5,6] then return M4[i-2][j-2]; end if;
	return 0; 
end function;

embed6:=function(M2,M4)
	return GK6!Matrix(K, 6, 6, [aux_embed6(i,j,M2,M4): i, j in [1..6]]);
end function;

//We construct the homogeneous m-simplex for a given value m
m:=12;

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

cc:=hom_simplex(m);
zero_vector:=Zero(VectorSpace(Rationals(),#cc));

//Assumes the function "avervalue" has been defined

momentsp:=function(p,i,j,k)
	moment:=0;
	if (i+k) mod 2 eq 0 then
		q:=p[1]^i*p[2]^j*p[3]^k;
		moment:=avervalue(q);
	end if;	
	return moment;
end function;
 
moments_simplexp:=function(p)
	p2:=Coefficients(p);
	p3:=[p2[2],p2[3],p2[4]];
	simplex:=[];
	for c in cc do 
		simplex:=Append(simplex,momentsp(p3,c[1],c[2],c[3])); 
	end for;
	return Vector(Rationals(),simplex);
end function;


//Given a list of groups gl (of a given type *), returns the list of moments for each of the groups in gl in the same order.

groups_moments:=function(gl)
	char_polys:=[];
	char_poly_moments:=AssociativeArray();
	gm:=[];
	for i in [1..#gl] do
		G:=gl[i];
		moment:=zero_vector;
		o:=0;
		for g in G do
			p:=CharacteristicPolynomial(Parent(U6)!g*U6);
			if not p in char_polys then
				char_polys:=Append(char_polys,p);	
				char_poly_moments[p]:=moments_simplexp(p);
			end if;
			moment:=moment+char_poly_moments[p];
			o:=o+1;
		end for;
		gm:=Append(gm,moment/o);
	end for;
	return gm;
end function;

//Given a list of groups gl (of a given type *), returns the lists of moments of a1,a2,a3 for each of the groups in gl in the same order.

a1a2a3p:=function(p)
	p2:=Coefficients(p);
	p3:=[p2[2],p2[3],p2[4]];
	seq1:=[];
	seq2:=[];
	seq3:=[];
	for i in [1..12] do 
		seq1:=Append(seq1,momentsp(p3,i,0,0));
		seq2:=Append(seq2,momentsp(p3,0,i,0));
		seq3:=Append(seq3,momentsp(p3,0,0,i)); 
	end for;
	seq:=seq1 cat seq2 cat seq3;
	return Vector(Rationals(),seq);
end function;

zero_vector:=Zero(VectorSpace(Rationals(),36));
groups_a1a2a3:=function(gl)
	char_polys:=[];
	char_poly_a1a2a3:=AssociativeArray();
	ga1a2a3:=[];
	for i in [1..#gl] do
		G:=gl[i];
		seq:=zero_vector;
		o:=0;
		for g in G do
			p:=CharacteristicPolynomial(Parent(U6)!g*U6);
			if not p in char_polys then
				char_polys:=Append(char_polys,p);	
				char_poly_a1a2a3[p]:=a1a2a3p(p);
			end if;
			seq:=seq+char_poly_a1a2a3[p];
			o:=o+1;
		end for;
		ga1a2a3:=Append(ga1a2a3,seq/o);
	end for;
	return ga1a2a3;
end function;


//We next define the functions to compute the squares of the restrictions of the irreducible characters of USp(6) to each of the groups. We assume that "avervalue" has been defined.

diagonalp:=function(p)
	p2:=Coefficients(p);
	a1:=p2[2];
	a2:=p2[3];
	a3:=p2[4];
	chi:=[
	1,
	-a1,
	a2-1,
	-a3+a1,
	a1^2-a2,
	
	-a1*a2+a1+a3,
	a1*a3-a1^2-a2+1,
	a2^2-a1*a3-a2,
	-a2*a3+2*a1*a2-a1,
	a3^2-a2^2-a1*a3+2*a2-1,

	-a1^3 + 2*a1*a2 - a3,
	a1^2*a2 - a1^2 - a1*a3 - a2^2 + 2*a2,
	a1^3 - a1^2*a3 - 2*a1 + a2*a3,
	a1^2*a3 - a1*a2^2 + a2*a3 - a3,
	-2*a1^2*a2 + 2*a1^2 + a1*a2*a3 + a1*a3 - a3^2,

	-a1^2*a3 - a1*a2^2 + 4*a1*a2 + a1*a3^2 - 2*a1 - a2*a3, 
	a1^2*a2 - 2*a1*a2*a3 + a1*a3 + a2^3 - 2*a2^2 + a3^2,
	-a1^3 - a1^2*a3 + 2*a1*a2^2 + a1*a3^2 - a2^2*a3 - a2*a3 + a3,
	2*a1^2*a2 - a1^2 - 2*a1*a2*a3 - a1*a3 - a2^3 + 3*a2^2 + a2*a3^2 - 2*a2,
	a1^2*a3 - 3*a1*a2^2 + 2*a1*a2 + a1*a3^2 + 2*a2^2*a3 - 2*a2*a3 - a3^3 + a3
	];
	dia:=[1];
	for i in [2..20] do//20 is the length of chi
		dia:=Append(dia,avervalue(chi[i]^2));
	end for;
	return Vector(Rationals(),dia);
end function;

groups_diagonals:=function(gl)
	zero_vector:=Zero(VectorSpace(Rationals(),20));//20 is the length of chi
	char_polys:=[];
	char_poly_diagonals:=AssociativeArray();
	gdi:=[];
	for i in [1..#gl] do
		G:=gl[i];
		diagonal:=zero_vector;
		o:=0;
		for g in G do
			p:=CharacteristicPolynomial(Parent(U6)!g*U6);
			if not p in char_polys then
				char_polys:=Append(char_polys,p);	
				char_poly_diagonals[p]:=diagonalp(p);
			end if;
			diagonal:=diagonal+char_poly_diagonals[p];
			o:=o+1;
		end for;
		gdi:=Append(gdi,diagonal/o);
	end for;
	return gdi;
end function;

components_diagonals:=function(gl)
	char_polys:=[];
	components_diagonals:=[];
	for i in [1..#gl] do
		G:=gl[i];
		for g in G do
			p:=CharacteristicPolynomial(Parent(U6)!g*U6);
			if not p in char_polys then
				char_polys:=Append(char_polys,p);	
				diagonal:=diagonalp(p);
				if not diagonal in components_diagonals then
    				components_diagonals:=Append(components_diagonals,diagonal);
                end if;
			end if;
		end for;
	end for;
	return components_diagonals;
end function;

//We next define the functions to compute the zvectors
E:=VectorSpace(Rationals(),23);

z_vector_p:=function(p)
	p:=Coefficients(p);
	p:=[p[2],p[3],p[4]];
	zvp:=Zero(E);	
	
	if p[2] eq -1 then zvp[2]:=1; 
		elif p[2] eq  0 then zvp[3]:=1; 	
		elif p[2] eq  1 then zvp[4]:=1; 
		elif p[2] eq  2 then zvp[5]:=1; 	
		elif p[2] eq  3 then zvp[6]:=1;
	end if;

	if p[1] eq 0 then
		zvp[1]:=1; 
		if   p[2] eq -1 then zvp[8]:=1;  
		elif p[2] eq  0 then zvp[9]:=1;  	
		elif p[2] eq  1 then zvp[10]:=1; 
		elif p[2] eq  2 then zvp[11]:=1; 	
		elif p[2] eq  3 then zvp[12]:=1;
		end if; 
		if p[3] eq 0 then 
			zvp[13]:=1; 
			if   p[2] eq -1 then zvp[19]:=1; 
			elif p[2] eq  0 then zvp[20]:=1; 	
			elif p[2] eq  1 then zvp[21]:=1; 
			elif p[2] eq  2 then zvp[22]:=1; 
			elif p[2] eq  3 then zvp[23]:=1; 
			end if;	
		end if;
	end if;
	
	if p[3] eq 0 then
		zvp[7]:=1;
		if   p[2] eq -1 then zvp[14]:=1;
		elif p[2] eq  0 then zvp[15]:=1;	
		elif p[2] eq  1 then zvp[16]:=1; 
		elif p[2] eq  2 then zvp[17]:=1; 	
		elif p[2] eq  3 then zvp[18]:=1;
		end if; 	
	end if;
	
	assert not(p[2] in [-15..-2]);
	assert not(p[2] in [4..-15]);
	return zvp;

end function;

//Given a list of groups gl (of a given type *), returns the list of zvectors for each of the groups in gl in the same order.

groups_zvectors:=function(gl)
	char_polys:=[];
	char_poly_zvectors:=AssociativeArray();
	gz:=[];

	for i in [1..#gl] do
		G:=gl[i];
		zv:=Zero(E);
		o:=0;
		for g in G do
			p:=CharacteristicPolynomial(Parent(U6)!g*U6);
			if not p in char_polys then
				char_polys:=Append(char_polys,p);	
				char_poly_zvectors[p]:=z_vector_p(p);
			end if;
			zv:=zv+char_poly_zvectors[p];		
			o:=o+1;
		end for;
		gz:=Append(gz,zv/o);
	end for;
	return gz;
end function;


//Given a list of groups gl (of a given type *), returns the list of GAP labels of their groups of components for each of the groups in gl in the same order.

groups_labels:=function(gl)

	glab:=[];
	for i in [1..#gl] do
		lab:=label_group(gl[i]);
		glab:=Append(glab,lab);
	end for;
	return glab;
end function;

//Given K<zz>=cyclotomic field of 24 roots of unity, constructs the ST of genus 2 with identity component U(1). All presentations contain -1. Returns: a list "small_list" containing all names and a dictionary "small_dict" containing all groups.

genus2_F:=function(K)

	small_list:=[];//This will contain the names of the 32 genus 2 groups of type F
	small_dict:=AssociativeArray();//This will contain the 32 groups
	small_generators:=AssociativeArray();//This will contain the generators of the 32 groups	

	zz:=K.1;
	z:=AssociativeArray();
	for i in [1,2,3,4,6,8,12] do
		z[i]:=zz^(Integers()!(24/i)); 
	end for;
	
	zeta_mat:=AssociativeArray();

	for i in [1,2,3,4,6] do
		Cstr:="C_" cat IntegerToString(i);
		zeta_mat[2*i]:=GK4!embed(Matrix(K,2,2,[z[2*i],0,0,ComplexConjugate(z[2*i])]));
		small_dict[Cstr]:=sub<GK4|zeta_mat[2*i]>;
		small_generators[Cstr]:=[zeta_mat[2*i]];
		small_list:=Append(small_list, Cstr);
	end for;

	j_mat:=GK4!embed(Matrix(K,2,2,[0,1,-1,0]));

	for i in [2,3,4,6] do
		Dstr:="D_" cat IntegerToString(i);
		small_dict[Dstr]:=sub<GK4|j_mat,zeta_mat[2*i]>;
		small_generators[Dstr]:=[j_mat,zeta_mat[2*i]];
		small_list:=Append(small_list, Dstr);
	end for;

	t_mat:=GK4!embed(Matrix(K,2,2,[(1+z[4])/2,(1+z[4])/2,(-1+z[4])/2,(1-z[4])/2]));

	small_dict["T"]:=sub<GK4|j_mat,t_mat>;
	small_generators["T"]:=[j_mat,t_mat];
	small_dict["O"]:=sub<GK4|j_mat,t_mat,zeta_mat[8]>;
	small_generators["O"]:=[j_mat,t_mat,zeta_mat[8]];


	J:=BlockMatrix(2,2,[0,Matrix(K,2,2,[0,1,-1,0]),Matrix(K,2,2,[0,-1,1,0]),0]);


	small_list:=Append(small_list, "T");
	small_list:=Append(small_list, "O");

	for i in small_list do
		Jname:="J(" cat i cat ")";
		small_dict[Jname]:=sub<GK4| small_dict[i], J>;
		small_generators[Jname]:=small_generators[i] cat [J];
		small_list:=Append(small_list,Jname);
	end for;

	for i in [2,4,6] do
		C1str:= "C_{" cat IntegerToString(i) cat ",1}";
		small_dict[C1str]:=sub<GK4| J*zeta_mat[2*i]>;
		small_generators[C1str]:=[J*zeta_mat[2*i]];
		small_list:=Append(small_list,C1str);
	end for;
	
	for i in [2,4,6] do
		D1str:="D_{" cat IntegerToString(i) cat ",1}";
		small_dict[D1str]:=sub<GK4| J*zeta_mat[2*i],j_mat>;
		small_generators[D1str]:=[ J*zeta_mat[2*i],j_mat];
		small_list:=Append(small_list,D1str);
	end for;
	
	for i in [3,4,6] do
		D2str:="D_{" cat IntegerToString(i) cat ",2}";
		small_dict[D2str]:=sub<GK4| zeta_mat[2*i],J*j_mat>;
		small_generators[D2str]:=[zeta_mat[2*i],J*j_mat];
		small_list:=Append(small_list,D2str);
	end for;

	small_dict["O_1"]:=sub<GK4|j_mat,t_mat,J*zeta_mat[8]>;
	small_generators["O_1"]:=[j_mat,t_mat,J*zeta_mat[8]];
	small_list:=Append(small_list, "O_1");

	return <small_list, small_dict,small_generators>;
end function;



//Given K<zz>=cyclotomic field of 24 roots of unity, constructs the ST of genus 2 with identity component SU(2). All presentations contain -1. Returns: a list "small_list" containing all names and a dictionary "small_dict" containing all groups.

genus2_E:=function(K)
	
	small_list:=[];//This will contain the names of the 10 genus 2 groups of type E
	small_dict:=AssociativeArray();//This will contain the 10 groups
	
	zz:=K.1;
	z:=AssociativeArray();
	for i in [1,2,3,4,6,8,12] do
		z[i]:=zz^(Integers()!(24/i)); 
	end for;

	exp_mat:=AssociativeArray();
	
	for i in [1,2,3,4,6] do
		Estr:="E_" cat IntegerToString(i);
		exp_mat[2*i]:=GK4!embed(MK2!Matrix(K,2,2,[z[2*i],0,0,z[2*i]]));
		small_dict[Estr]:=sub<GK4|exp_mat[2*i]>;
		small_list:=Append(small_list, Estr);
	end for;

	J:=BlockMatrix(2,2,[0,MK2!Matrix(K,2,2,[0,1,-1,0]),MK2!Matrix(K,2,2,[0,-1,1,0]),0]);

	for i in small_list do
		Jname:="J(" cat i cat ")";
		small_dict[Jname]:=sub<GK4| small_dict[i], J>;
		small_list:=Append(small_list,Jname);
	end for;

	return <small_list, small_dict>;
end function;


//The following are the printing functions that we use.

print_in_file:=procedure(g_list,file_name)

	SetOutputFile( file_name: Overwrite:=true);
	print "L:=[";
	for i in [1..#g_list] do
		print g_list[i];
		if not i eq #g_list then print ","; end if;
	end for;
	print "];";
	UnsetOutputFile();

end procedure;

print_in_file_str:=procedure(g_list,file_name)//This one is used to print strings

	SetOutputFile( file_name: Overwrite:=true);
	print "L:=[";
	for i in [1..#g_list] do
		print "\"" cat g_list[i] cat "\"";
		if not i eq #g_list then print ","; end if;
	end for;
	print "];";
	UnsetOutputFile();

end procedure;


//Auxilliary function for printing
format_good_aux:=function(m)
	mm:=<>;
	for r in [1..Ncols(m)] do
		mm:=Append(mm,m[r]);
	end for;
	return mm;
end function;

format_good:=function(L)
	LL:=[];
	for l in L do
		LL:=Append(LL,format_good_aux(l));
	end for;
	return LL;
end function;

//Given the type "letter", a list of group names "gn", and a list of groups "groups_list", computes and prints moments, diagonals, zvectors and GAP labels according to booleans bm, bdi, bz, blab, and ba1a2a3.

generate_data:=procedure(letter, gn, groups_list, bm, bdi, bz, blab, ba1a2a3, bcomp)

	file_name:="names_" cat letter cat ".txt";
	print_in_file_str(gn,file_name);

	if bm then
		file_name:="moments_" cat letter cat ".txt";
		gm:=groups_moments(groups_list);
		gm:=format_good(gm);
		print_in_file(gm,file_name);
	end if;

	if bdi then
		file_name:="diagonals_" cat letter cat ".txt";
		gdi:=groups_diagonals(groups_list);
		gdi:=format_good(gdi);
		print_in_file(gdi,file_name);
	end if;

	if bz then
		file_name:="zvectors_" cat letter cat ".txt";
		gz:=groups_zvectors(groups_list);
		gz:=format_good(gz);
		print_in_file(gz,file_name);
	end if;
	
	if blab then
		file_name:="labels_" cat letter cat ".txt";
		glab:=groups_labels(groups_list);
		print_in_file(glab,file_name);
	end if;

	if ba1a2a3 then
		file_name:="a1a2a3_" cat letter cat ".txt";
		ga1a2a3:=groups_a1a2a3(groups_list);
		ga1a2a3:=format_good(ga1a2a3);
		print_in_file(ga1a2a3,file_name);
	end if;

    if bcomp then
        file_name:="comp_diagonals_" cat letter cat ".txt";
        gcomp:=components_diagonals(groups_list);
        gcomp:=format_good(gcomp);
        print_in_file(gcomp,file_name);
    end if;

end procedure;



//The next function identifies one group in the classification using the irreducible characters values and the label in the two cases of ambiguity.
identify_subgroup:=function(H,gdi,gn,glab)
	diagH:=groups_diagonals([H])[1];
	str:="Missing diagonal";
	for i in [1..#gdi] do
		if gdi[i] eq format_good_aux(diagH) then
			if not(gn[i] eq "J(C(3,3))") and not(gn[i] eq "J_s(C(3,3))") then	
				str:=gn[i];
			else
				if label_group(H) eq <54, 5> then str:="J(C(3,3))"; end if;
				if label_group(H) eq <54, 8> then str:="J_s(C(3,3))"; end if;
			end if; 
		end if;
	end for;
	return str;
end function;
	
//The next function computes the maximal subgroups of each of the groups in gl.
subgroups_lattice:=function(gl,gdi,gn,glab)
	glat:=[];
	for i in [1..#gl] do
		G:=gl[i];
		L:=MaximalSubgroups(G);	
		max_subgroups:=[];
		for l in L do
			H:=l`subgroup;
			if not(label_group(H) eq glab[i]) then//Rules out the possibility that H and G are the same ST group			
				max_subgroups:=Append(max_subgroups, "\"" cat identify_subgroup(H,gdi,gn,glab) cat "\"");	
			end if;
		end for;
		glat:=Append(glat,max_subgroups);	
	end for;
	return glat;
end function;



