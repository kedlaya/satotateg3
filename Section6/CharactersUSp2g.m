//In this file, we compute the irreducible characters of USp(2*g).
//This file assumes that a value for g has been set and that the fraction field FR of a polynomial ring PR with (at least) g+1 variables has been defined.


//First we compute the coefficients of the characteristic polynomial of an element in the maximal torus.

coefs_char_poly:=function()
	char_poly:=1;
	for i in [1..g] do
		char_poly:=char_poly*(PR.(g+1)-FR!PR.i)*(PR.(g+1)-(1 / FR!PR.i));
	end for;

	L:=Coefficients(Numerator(char_poly),PR.(g+1));
	cc:=[];
	for i in [1..g] do
		cc:=Append(cc,L[2*g+1-i]/Denominator(char_poly));
	end for;
	return cc;
end function;

//Given nonnegative integers l,r, returns the decreasing sequences of nonnegative integers of length l starting with an integer at most r.
seqlr:=function(l,r)
	if l eq 0 then R:=[[]];	
	else
		if r eq 0 then 
			zz:=[];
			for i in [1..l] do
				zz:=Append(zz,0);
			end for;
			R:=[zz]; 
		else 
			R:=[];
			for i in [0..r] do
				Ri:=$$(l-1,r-i);
				for req in Ri do
					R:=Append(R, [r-i] cat req);
				end for; 
			end for;
		end if;
	end if;
	return R;
end function;

//Given a characteristic function p of USp(2*g), finds a simplex containing its dominant weights.

weights:=function(p);
	d:=Degree(Numerator(p),PR.1) - Degree(Denominator(p),PR.1);
	return seqlr(g,d);
end function;


//Given a characteristic function pol of USp(2*g) and a g-tuple of decreasing nonnegative integers ww (a dominant weight), finds the multitplicity of ww in pol. 
 
coef_of:=function(pol,ww)
	coef:=Numerator(pol);
	for i in [1..g] do
		a:=ww[i];
		b:=Degree(Denominator(pol),PR.i);
		coef:=Coefficient(coef,PR.i,a+b);
	end for;
	return coef;
end function;

//Given a characteristic function p of USp(2*g), finds the multiplicities of all of its dominant weights.
multiplicities:=function(p)
	we:=weights(p);
	mu:=[];
	for ww in we do
		mu:=Append(mu,coef_of(p,ww));
	end for;
	return mu;
end function;

//We next construct the Weyl group as a group of matrices.
GZ:=GeneralLinearGroup(g,Integers());

entA:=function(i,j,l)
	if i eq l then
		if j eq (l+1) then e:=1;
		else e:=0; end if;
	elif i eq (l+1) then
		if j eq l then e:=1;
		else e:=0; end if;
	else
		if i eq j then e:=1;
		else e:=0; end if;	
	end if;
		
	return e;
end function;

entB:=function(i,j,l)
	if i eq j then
		if i eq l then e:=-1; 
		else e:=1; end if;
	else
		e:=0;
	end if;
	return e;
end function;

gene:=[];
for l in [1..g-1] do
	A:=Matrix(Rationals(),g,g,[<i,j,entA(i,j,l)>:i,j in [1..g]]);
	gene:=Append(gene,A);
end for;
for l in [1..g] do
	B:=Matrix(Rationals(),g,g,[<i,j,entB(i,j,l)>:i,j in [1..g]]);
	gene:=Append(gene,B);
end for;

Weyl:=sub<GZ|gene>;

//The next function will be used to obtain the numerator and denominator in Weyl's character formula.
DW:=function(ww)
	V:=Matrix(Integers(),g,1,[<i,j,ww[i]>: i in [1..g], j in [1]]);   
	D:=0;
	for gg in Weyl do
		exps:=gg*V;
		prod:=1;
		for i in [1..g] do
			prod:=prod*(FR!PR.i)^(exps[i][1]);
		end for;
        	D:=D+Determinant(gg)*prod;
	end for;
	return D;
end function;

//Given a dominant weight ww, returns the character of the representation of USp(2g) with highest weight ww.
character:=function(ww)
	wwa:=[]; wwg:=[];
	for i in [1..g] do wwa[i]:=ww[i]+g+1-i; wwg[i]:=g+1-i; end for;
	return DW(wwa)/DW(wwg);
end function;

//Given a characteristic function p of USp(2*g), returns the multiplicity of the identity in p.
mult_id:=function(p)
	if p eq 0 then return 0; end if;
	we:=weights(p);	
	mu:=multiplicities(p);
	val:=p;
	for i in [1..#we-1] do
		if not (coef_of(val,we[i]) eq 0) then
			ch:=character(we[i]);
			k:=coef_of(val,we[i])/coef_of(ch,we[i]);
			val:=val-k*ch;	
		end if;
	end for;
	return val;
end function;




