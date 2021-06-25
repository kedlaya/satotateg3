//This file computes the inhomogeneous m-simplex of moments of N(U(3)) via the dual Cauchy formula -DCF- for N(U(3)) in USp(6). In order to compute the right-hand side of DCF for a given value of m, the first task will consist on determining the multiplicity of the identity in a character of USp(6) when restricted to U(3). We provide a computational method that confirms Lemma 6.10.i)

m:=5;
g:=3;
PR:=PolynomialRing(Rationals(),4);

u:=PR.1;
v:=PR.2;
w:=PR.3;
T:=PR.4;

AssignNames(~PR,["u","v","w","T"]);
FR:=FieldOfFractions(PR);
AssignNames(~FR,["u","v","w","T"]);

load "CharactersUSp2g.m";
load "CharactersSU3.m";


//We next restrict a central function of USp(6) to U(3). The term u^a*v^b*w^c is mapped to u^(a-c)*v^(b-c)*w^(a+b+c)
restU3:=function(p)
	M:=Monomials(Numerator(p));
	C:=Coefficients(Numerator(p));
	p1:=0;
	for i in [1..#M] do
		m:=M[i];
		co:=C[i];
		a:=Degree(m,u);
		b:=Degree(m,v);
		c:=Degree(m,w);
		p1:=p1+co*(FR!u)^(a-c)*(FR!v)^(b-c)*(FR!w)^(a+b+c);
	end for; 
	a:=Degree(Denominator(p),u);
	b:=Degree(Denominator(p),v);
	c:=Degree(Denominator(p),w);
	p2:=(FR!u)^(a-c)*(FR!v)^(b-c)*(FR!w)^(a+b+c);  	
	return p1/p2;
end function;

//Given a characteristic function of U(3), takes its constant term wrt the variable w.
restSU3:=function(p)
	dw:=Degree(Denominator(p),w);
	dv:=Degree(Denominator(p),v);
	du:=Degree(Denominator(p),u);
	return Coefficient(Numerator(p),w,dw)/(u^du*v^dv);
end function;

//We next compute the right-hand side of the dual Cauchy formula for N(U(3)) in USp(6) (for the value m). Given a dominant weight ww of USp(2*m) that is a subpartition of the rectangular partition (3^m)=(3,...(m)...,3), we want to compute the value "mHNU3(ww)". Let us describe mHNU3(ww). Set wwtilde=(m-wwp[3],m-wwp[2],m-wwp[1]), where wwp is the transpose partition of ww. Then mHNU3(ww) is defined as the multiplicity of the identity in the restriction to N(U(3)) of the representation of USp(6) with highest weight wwtilde. Along the way, we also compute the dual Cauchy formula for U(3), by computing the corresponding quantity mHU3. 

//First compute the transpose of a subpartition ww of (3,...,(m)...,3).
ultimge:=function(ww,j)
	i:=1;
	while i le #ww and ww[i] ge j do
		i:=i+1;
	end while;
	return i-1;
end function;

transposa:=function(ww)
	trans:=[];
	for j in [1..3] do
		trans:=Append(trans, ultimge(ww,j));
	end for;
	return trans;
end function;

tilde:=function(ww)
	wwp:=transposa(ww);
	return [m-wwp[3],m-wwp[2],m-wwp[1]];
end function;

//ww is assumed to be of length m (completed with zeroes if necessary)
mHU3:=function(ww)
	cha:=character(tilde(ww));
	h:=restSU3(restU3(cha));
	return Integers()!mult_idSU3(h);
end function;

//We next compute the right-hand side of the dual Cauchy formula. First we compute the index set for the sum
Iset:=seqlr(m,3);
//We next compute the mHvector relative to U(3).
mHvectorU3:=[];
for ww in Iset do
	mHvectorU3:=Append(mHvectorU3,mHU3(ww));
end for;
//Note that the values of mH are only 0 or 1, as predicted by Lemma 6.10.i).

//The following function uses Lemma 6.10.ii).
mHNU3:=function(ww)
	tww:=tilde(ww);
	pes:=tww[1]+tww[2]+tww[3];
	mult:=0;
	if (pes mod 4 eq 0) and (tww[1] mod 2 eq 0) and (tww[2] mod 2 eq 0) and (tww[3] mod 2 eq 0) then mult:=1; end if;
	return mult;
end function;

//The mHvector relative to N(U(3))
mHvectorNU3:=[];
for ww in Iset do
	mHvectorNU3:=Append(mHvectorNU3,mHNU3(ww));
end for;


mHvector:=mHvectorNU3;

//From now on "character(.)" returns characters of USp(2*m), not of USp(6).
PR:=PolynomialRing(Rationals(),m+1);
FR:=FieldOfFractions(PR);
g:=m;

load "CharactersUSp2g.m";

S:=0;
for i in [1..#Iset] do
	if not(mHvector[i] eq 0) then
		newterm:=mHvector[i]*character(Iset[i]);
		S:=S+newterm;
	end if;
end for;

prod:=1;
for i in [1..m] do
	prod:=prod*(PR.i)^3;
end for;

RHSDCF:=prod*S;

coef_monomial:=function(p,e)
	val:=Numerator(p);
	for i in [1..m] do
		dden:=Degree(Denominator(p),PR.i);
		val:=Coefficient(val,PR.i,e[i]+dden);
	end for;
	return val;
end function;

//Builds a list of (lexicographically ordered) l-tuples whose entries add up to at most m
msimplex:=function(m,l)
	if l eq 0 then R:=[[]]; 
	else
		R:=[];
		for i in [0..m] do
			Ri:=$$(m-i,l-1);
			for req in Ri do
				R:=Append(R, [i] cat req);
			end for;
		end for;
	end if;
	return R;
end function;

//An easy calculation shows that the moment of a_1^i*a_2^j*a_3^k is given by the coefficient of (T_1*...*T_i)*(T_{i+1}*...*T_{i+j})^2*(T_{i+j+1}*...*T_{i+j+k})^3 of the RHSDCF. The next function computes [0...(m-i-j-k)...0,1,...(i)...,1,2...(j)...2,3...(k)...3] from (i,j,k) 
Texponent:=function(ww)
	i:=ww[1];
	j:=ww[2];
	k:=ww[3];
	R:=[];
	for l in [1..(m-i-j-k)] do
		R:=Append(R,0);
	end for;
	for l in [1..i] do
		R:=Append(R,1);
	end for;
	for l in [1..j] do
		R:=Append(R,2);
	end for;
	for l in [1..k] do
		R:=Append(R,3);
	end for;
	return R;
end function;

inhomsimplexNU3:=<>;
for ee in msimplex(m,3) do
	mom:=coef_monomial(RHSDCF,Texponent(ee));
	inhomsimplexNU3:=Append(inhomsimplexNU3,mom);
end for;

inhomsimplexNU3;


