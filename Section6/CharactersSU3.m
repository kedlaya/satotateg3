//In this file, we compute the irreducible characters of SU(3).
//This file assumes that the fraction field FR of a polynomial ring PR with variables u,v has been defined.


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

//Given a characteristic function p of SU(3), finds a simplex containing its dominant weights.
weightsSU3:=function(p);
	we:=[];
	d:=Degree(Numerator(p),u) - Degree(Denominator(p),u);
	return seqlr(2,d);
end function;


//Given a characteristic function pol of SU(3) and a 2-tuple of decreasing nonnegative integers ww (a dominant weight), finds the multiplicity of ww in pol 
coef_ofSU3:=function(pol,ww)
	a:=ww[1];
	b:=ww[2];
	dqu:=Degree(Denominator(pol),u);
	dqv:=Degree(Denominator(pol),v);
	c1:=Coefficient(Numerator(pol),u,a+dqu);
	c2:=Coefficient(c1,v,b+dqv);
	return c2;    
end function;

//Given a characteristic function p of SU(3), finds the multiplicities of all of its dominant weights.
multiplicitiesSU3:=function(p)
	we:=weightsSU3(p);
	mu:=[];
	for ww in we do
		mu:=Append(mu,coef_ofSU3(p,ww));
	end for;
	return mu;
end function;


//Returns numerator and denominator in Weyl's character formula for SU(3)
DWSU3:=function(ww)
	a:=ww[1];
	b:=ww[2];
	D:=(FR!u)^a*(FR!v)^b+(FR!u)^(b-a)*(FR!v)^(-a)+(FR!u)^(-b)*(FR!v)^(a-b)-(FR!u)^(-a)*(FR!v)^(b-a)-(FR!u)^(a-b)*(FR!v)^(-b)-(FR!u)^b*(FR!v)^a;
	return D;
end function;

//Given a dominant weight ww, returns the character of the representation of SU(3) with highest weight ww
characterSU3:=function(ww)
	a:=ww[1];
	b:=ww[2];
	return DWSU3([a+2,b+1])/DWSU3([2,1]) ;
end function;

//Given a characteristic function p of SU(3), returns the multiplicity of the identity in p.
mult_idSU3:=function(p)
	if p eq 0 then return 0; end if;
	we:=weightsSU3(p);	
	mu:=multiplicitiesSU3(p);
	val:=p;
	for i in [1..#we-1] do
		if not (coef_ofSU3(val,we[i]) eq 0) then
			ch:=characterSU3(we[i]);
			k:=coef_ofSU3(val,we[i])/coef_ofSU3(ch,we[i]);
			val:=val -k*ch;	 
		end if;
	end for;
	return val;
end function;

