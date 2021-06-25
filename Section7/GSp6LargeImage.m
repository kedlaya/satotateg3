Sp6primes := [2,3,5];
Sp4primes := [2,3,5,7,11,13];
Sp2primes := [5,7,11,13,17,19,23,29,31];  // Can't use 2 or 3 for Sp(2), surjective mod 2,3 image is not sufficient

M6 := [[H`subgroup: H in MaximalSubgroups(Sp(6,p))]:p in Sp6primes];
M4 := [[H`subgroup: H in MaximalSubgroups(Sp(4,p))]:p in Sp4primes];
M2 := [[H`subgroup: H in MaximalSubgroups(Sp(2,p))]:p in Sp2primes];

X6 := [*[{CharacteristicPolynomial(c[3]):c in ConjugacyClasses(H)}:H in M6[i]]: i in [1..#Sp6primes]*];
X4 := [*[{CharacteristicPolynomial(c[3]):c in ConjugacyClasses(H)}:H in M4[i]]: i in [1..#Sp4primes]*];
X2 := [*[{CharacteristicPolynomial(c[3]):c in ConjugacyClasses(H)}:H in M2[i]]: i in [1..#Sp2primes]*];

ZZ := Integers();

intrinsic TestSp6Image (D::RngIntElt,C::Crv:B:=200) -> BoolElt
{ Given an integer D divisible by all primes of bad reduction and a string C specifying a genus 3 curve, returns true if good Euler factors over p <= B prove large Galois image mod any of 2,3,5. }
    S6 := X6;
    for p in PrimesInInterval(1,B) do
        if IsDivisibleBy(ZZ!D,p) then continue; end if;
        try
            f := LPolynomial(ChangeRing(C,GF(p)));
        catch e
            continue;
        end try;
        if Degree(f) ne 6 then continue; end if;
        for i in [1..#Sp6primes] do
            if p mod Sp6primes[i] eq 1 then
                S6[i] := [s:s in S6[i]|ChangeRing(f,GF(Sp6primes[i])) in s];
                if #S6[i] eq 0 then return true; end if;
            end if;
        end for;
    end for;
//    print [#S6[i]:i in [1..#Sp6primes]];
    return false;
end intrinsic;

intrinsic TestSp4ImageGenus3 (D::RngIntElt,C::Crv:B:=200) -> BoolElt, MonStgElt
{ Given an integer D divisible by all primes of bad reduction and a genus 3 curve C, returns true if good Euler factors over p <= B prove Jac(C) is product of generic abelian surface and an elliptic curve (over Q) in the Cremona database. }
    S4 := X4;
    L2 :=[];
    for p in PrimesInInterval(1,B) do
        if IsDivisibleBy(ZZ!D,p) then continue; end if;
        f := LPolynomial(ChangeRing(C,GF(p)));
        if Degree(f) ne 6 then continue; end if;
        A := Factorization(f);
        if Degree(A[#A][1]) gt 4 or (Degree(A[1][1]) eq 3 and Degree(A[2][1]) eq 3) then return false,_,_; end if;
        if Degree(A[#A][1]) eq 4 then
            g := A[#A][1];
            Append(~L2,<p,ExactQuotient(f,g)>);
            for i in [1..#Sp4primes] do
                if p mod Sp4primes[i] eq 1 then
                    S4[i] := [s:s in S4[i]|ChangeRing(g,GF(Sp4primes[i])) in s];
                end if;
            end for;            
        end if;
    end for;
    if &and[#S4[i] gt 0: i in [1..#Sp4primes]] then return false,_,_; end if;
    DB := CremonaDatabase();
    E := &cat [[e: e in EllipticCurves(DB,N)| n eq 1 where _,_,n:=CremonaReferenceData(DB,e)] : N in Divisors(ZZ!D) | N ge 11 and N lt 400000];
    for r in L2 do
        E := [e : e in E | EulerFactor(e,r[1]) eq r[2]];
        if #E eq 0 then printf "Warning: could not find elliptic curve factor in Cremona database for curve C=%o with discriminant D=%o",C,D; return false,_,_; end if;
    end for;
    if #E gt 1 then printf "Warning: %o possible elliptic isogeny factors, need to increase prime bound B=%o\n", [CremonaReference(e):e in E],B; return false,_; end if;
    return true, CremonaReference(E[1]);
end intrinsic;

intrinsic TestSp4Image (D::RngIntElt,C::Crv:B:=200,SurjectivePrimes:=false,lpdata:=[]) -> BoolElt, MonStgElt
{ Given an integer D divisible by all primes of bad reduction and a genus 2 curve C, returns true if good Euler factors over p <= B prove Jac(C) is generic abelian surface. }
    S4 := X4;
    lp := AssociativeArray();
    for r in lpdata do lp[r[1]]:=r[2]; end for;
    for p in PrimesInInterval(1,B) do
        if IsDivisibleBy(ZZ!D,p) then continue; end if;
        f := IsDefined(lp,p) select lp[p] else LPolynomial(ChangeRing(C,GF(p)));
        if Degree(f) ne 4 then continue; end if;
        for i in [1..#Sp4primes] do
            if p mod Sp4primes[i] eq 1 then
                S4[i] := [s:s in S4[i]|ChangeRing(f,GF(Sp4primes[i])) in s];
                if not SurjectivePrimes and #S4[i] eq 0 then return true; end if;
            end if;
        end for;            
    end for;
    return &or[#S4[i]eq 0:i in [1..#Sp4primes]],[Sp4primes[i]:i in [1..#Sp4primes]|#S4[i] eq 0];
end intrinsic;

intrinsic TestEllipticFactor (D::RngIntElt, C::Crv:B:=200) -> BoolElt, MonStgElt
{ Given an integer D divisible by all primes of bad reduction and a string C specifying a genus 3 curve, returns true if there is a unique elliptic curve in the Cremona database whose L-poly divides that of C for all good p <= B. }
    L:=[];
    for p in PrimesInInterval(1,B) do
        if IsDivisibleBy(ZZ!D,p) then continue; end if;
        f := LPolynomial(ChangeRing(C,GF(p)));
        if Degree(f) ne 6 then continue; end if;
        A := Factorization(f);
        if Degree(A[#A][1]) gt 4 or (Degree(A[1][1]) eq 3 and Degree(A[2][1]) eq 3) then return false,_,_; end if;
        Append(~L,<p,f>);
    end for;
    DB := CremonaDatabase();
    E := &cat [[e: e in EllipticCurves(DB,N)| n eq 1 where _,_,n:=CremonaReferenceData(DB,e)] : N in Divisors(ZZ!D) | N ge 11 and N lt 400000];
    for r in L do
        E := [e : e in E | IsDivisibleBy(r[2],EulerFactor(e,r[1]))];
        if #E eq 0 then return false,_,_; end if;
    end for;
    if #E gt 1 then printf "Warning: %o possible elliptic isogeny factors, need to increase prime bound B=%o\n", [CremonaReference(e):e in E],B; return false,_; end if;
    return true, CremonaReference(E[1]);
end intrinsic;
    