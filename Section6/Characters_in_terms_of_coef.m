//In this file, we express the irreducible characters of USp(6) in terms of the characteristic polynomial coefficient characters

g:=3;
K:=Rationals();
PR<u,v,w,T,a1,a2,a3>:=PolynomialRing(K,2*g+1);
FR:=FieldOfFractions(PR);
AssignNames(~FR,["u","v","w","T","a1","a2","a3"]);

load "CharactersUSp2g.m";

//Given a g-tuple ww, returns the character ch=a1^(ww[1]-ww[2])*...*a(g-1)^(ww[g-1]-w[g])*ag^(ww[g]). The theory of representations of USp(2g) ensures that the representation with highest weight ww is a subrepresenation with multiplicity 1 of the representation afforded by ch
power_ais_char:=function(ww)
	cc:=coefs_char_poly();
	ch:=1;
	for i in [1..(g-1)] do
		ch:=ch*((-1)^i*cc[i])^(ww[i]-ww[i+1]);
	end for;
	return ch*((-1)^g*cc[g])^ww[g];
end function;

//Given a g-tuple ww, returns the formal expression a1^(ww[1])*...*ag^(ww[g])
power_ais:=function(ww)
	prod:=1;
	for i in [1..(g-1)] do
		prod:=prod*(PR.(g+1+i))^(ww[i]-ww[i+1]);
	end for;
	return prod*PR.(2*g+1)^ww[g];
end function;

//Given a characteristic function p of USp(2*g), expresses it in terms of the characteristic polynomial coefficient characters

express_in_ais:=function(p)
	if p eq 0 then ais:=0; 
	else
		we:=weights(p);	
		mu:=multiplicities(p);
		val:=p;
		ais:=0;
		for i in [1..#we] do
			if not (coef_of(val,we[i]) eq 0) then
				ch:=power_ais_char(we[i]);
				term:=power_ais(we[i]);
				k:=coef_of(val,we[i])/coef_of(ch,we[i]);
				k:=Integers()!k;
				val:=val-k*ch;
				ais:=ais+k*term;	
			end if;
		end for;
	end if;
	return ais;
end function;


for ww in seqlr(g,3) do
	print ww;
	express_in_ais(character(ww));
end for;

