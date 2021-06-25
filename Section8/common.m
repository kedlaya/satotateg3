intrinsic gg(f::RngUPolElt) -> Tup
{ The GAP id of the Galois group of f. }
    return IdentifyGroup(GaloisGroup(f));
end intrinsic;

intrinsic gg(K::FldNum) -> Tup
{ The GAP id of the Galois group of f. }
    return gg(DefiningPolynomial(K));
end intrinsic;

intrinsic nfdisc(f::RngUPolElt) -> RngIntElt
{ Given an irreducible polynomial with rational coefficients, returns the integer discriminant of the number field it defines. }
    return Integers()!Discriminant(RingOfIntegers(NumberField(PolynomialRing(Rationals())!f)));
end intrinsic;

intrinsic nfdisc(f::SeqEnum) -> RngIntElt
{ Given an irreducible polynomial with rational coefficients, returns the integer discriminant of the number field it defines. }
    return Integers()!Discriminant(RingOfIntegers(NumberField(PolynomialRing(Rationals())!f)));
end intrinsic;

intrinsic Elements(A::GrpAuto) -> SetEnum
{ Returns a list of the elements in the specified automorphism group. }
    n := #A;
    require n ne 0: "Automorphism group must be finite.";
    gens := Generators(A);
    S := gens;
    while true do
        S := {a*b: a in S, b in gens};
        if #S eq n then return S; end if;
    end while;
end intrinsic;

intrinsic ChangeRing(f::RngUPolElt, pi::Map) -> RngUPolElt
{ Given f = sum a_i*x^i returns sum pi(a_i)*x^i }
    return PolynomialRing(Codomain(pi))![pi(c):c in Coefficients(f)];
end intrinsic;

intrinsic NumberField (f::SeqEnum) -> FldNum
{ The number field defined by the polynoimal with the specified coefficients. }
    R<x> := PolynomialRing(Rationals());
    return NumberField(R!f);
end intrinsic;

// Given E:y^2=x^3+D and a prime p = 1 mod 3, returns the value of the Hecke character psi at the prime ideal (pi) above p with pi congruent to (2 mod 3)
// See [Sil94, II, Ex. 10.6]
intrinsic SexticPsi (D::RngIntElt,p::RngIntElt:z6:=0) -> FldNumElt
{ Returns the value of the Hecke character of y^2=x^3+D evaluated at the prime ideal (pi) of Z[zeta3] above p, with pi congruent to 2 modulo (3). }
    require D ne 0: "D must be nonzero";
    // compute a generator pi of prime above p with pi = 2 mod 3 (we make an arbitrary choice here by choosing a solution to the norm equation)
    require p mod 3 eq 1: "p must split in Z[zeta3])";
    _,a,b := NormEquation(3,4*p);
    if z6 eq 0 then z6 := NumberField([1,-1,1]).1; end if;
    c := (a-b) div 2;
    while c mod 3 ne 2 or b mod 3 ne 0 do x := -b; b +:= c; c:= x; end while;
    pi := c + b*z6;
    Fp := GF(p);
    vp := (Fp!D*4)^((p-1) div 6);
    a := [Fp!c:c in Eltseq(pi)];
    z6p := -a[1]/a[2];
    up := z6p; u := z6;
    while vp*up ne 1 do up *:= z6p; u *:= z6; end while;
    return -u*pi;
end intrinsic;

// Given E:y^2=x^3-Dx and a prime p = 1 mod 4, returns the value of the Hecke character psi at the prime ideal (pi) above p with pi congruent to 1 mod (2+2i)
// See [Sil94, II, Ex. 2.34]
intrinsic QuarticPsi (D::RngIntElt,p::RngIntElt:z4:=0) -> FldNumElt
{ Returns the value of the Hecke character of y^2=x^3-Dx evaluated at the prime ideal (pi) of Z[i] above p with pi congruent to 1 modulo (2+2i). }
    // compute a generator pi of prime above p with pi = 1 mod (2) (we make an arbitrary choice of solution to the norm equation)
    require p mod 4 eq 1: "p must split in Z[i]";
    _,a,b := NormEquation(1,p);
    if z4 eq 0 then R<x>:=PolynomialRing(Rationals()); z4 := NumberField(x^2+1).1; end if;
    while ((a-b) mod 4 ne 1 or b mod 2 ne 0) do x := a; a:= -b; b:= x; end while;
    pi := a + b*z4;
    Fp := GF(p);
    vp := (Fp!-D*4)^((p-1) div 4);
    a := [Fp!c:c in Eltseq(pi)];
    z4p := -a[1]/a[2];
    up := z4p; u := z4;
    while vp*up ne 1 do up *:= z4p; u *:= z4; end while;
    return u*pi;
end intrinsic;

intrinsic FrobeniusOrder(f::RngUPolElt, p::RngIntElt) -> RngIntElt
{ Returns the order of Frob_p in the Galois group of the splitting field of the polynomial f. }
    fp := ChangeRing(f,GF(p));
    if GCD(fp,Derivative(fp)) ne 1 then
        return Order(FrobeniusElement(NumberField(f),p));
    else
        return LCM([a[1]:a in DistinctDegreeFactorization(fp)]);
    end if;
end intrinsic;

intrinsic FrobeniusOrder(f::SeqEnum, p::RngIntElt) -> RngIntElt
{ Returns the order of Frob_p in the Galois group of the splitting field of the polynomial f. }
    R := PolynomialRing(Rationals());
    return FrobeniusOrder(R!f,p);
end intrinsic;

intrinsic FrobeniusOrder(K::FldNum, p::RngIntElt) -> RngIntElt
{ Returns the order of Frob_p in the Galois group of the splitting field of the polynomial f. }
    return FrobeniusOrder(DefiningPolynomial(K),p);
end intrinsic;

intrinsic FrobeniusElement(ctx::Tup, p::Any) -> RngIntElt
{ Given a tuple computed by FrobeniusContext and a prime number or prime ideal p unramified in the splitting field of f=ctx[1], returns the Frobenius element Frob_p as an index into the list of conjugacy classes ctx[4] for the permutation group ctx[3] corresponding to the splitting field of f=ctx[1]. }
    f := ctx[1];  S := ctx[2];
    if Degree(f) eq 1 then return 1; end if;
    if Degree(f) eq 2 and BaseRing(f) eq Integers() then n := KroneckerSymbol(Discriminant(f),p); return n eq -1 select 2 else n; end if;
    if Type(p) ne RngOrdIdl and Type(p) ne RngInt then K:=BaseRing(f); p := ideal<K|K!p>; end if;
    require IsPrime(p): "p must be a prime in the specified number field.";
    Fp, rp := ResidueClassField(p);
    RFp := PolynomialRing(Fp);
    fp := ChangeRing(f,rp);
    e := [i:i in [1..#ctx[4]]];
/*    try
        cc := Sort(Eltseq({* a[1]^^(Degree(a[2]) div a[1]):a in DistinctDegreeFactorization(fp)*}));
        e := [i:i in [1..#ctx[4]] | cc eq ctx[3][i][4]];
    catch err;
    end try;
    if #e eq 1 then return e[1],1; end if;*/
    Fpq := quo<RFp|fp>;
    xp := RFp!(Fpq.1^#Fp);
    tp := [rp(t):t in ctx[5]];
    for s in S do
        cp := Coefficients((xp*ChangeRing(s[1],rp)) mod fp);
        ap := &+[cp[i]*tp[i]:i in [1..#cp]];
        e := [i:i in e|Evaluate(ChangeRing(s[2][i],rp),ap) eq 0];
        if #e eq 1 then return e[1], s[1]; end if;
    end for;
    return 0, 0;
end intrinsic;

function _frobenius_test(ctx,P)
    S := {}; T := {};
    for p in P do
        c,h := FrobeniusElement(ctx,p);
        if c eq 0 then Include(~T,p); else Include(~S,h); end if;
    end for;
    return S,T;
end function;


intrinsic FrobeniusContext(f::RngUPolElt:Polys:=[],BadPrimes:=[],Verbose:=false) -> Tup
{ Given a polynomial f with coefficients in a number field precomputes data needed to compute Frobenius elements in the splitting field of f as a relative extension of its base field using the Dokchitser-Dokchitser algorithm. Returns a tuple <f,L,G,C,D,T>, where L is a list of pairs <h(x),[F_c]>, where F_c(X)=prod_s(X-sum_r(h(r)s(r)) with c ranging over the set C of conjugacy classes of G=Gal(f), s ranging over elements of C, and r ranging over the roots of f, G is the Galois group G of f (as a permutation group on deg(f) letters), C is a list of conjugacy classes of G (note the order may vary with each call to FrobeniusContext!), D is the norm of the discriminant of f (times any specified BadPrimes), and T is a list of traces of nth powers of the companion matrix of f from 0 to 2*deg(f)-2.}
    require IsMonic(f) and &and[IsIntegral(c):c in Coefficients(f)]: "Input most be a monic and integral";
    if Type(BadPrimes) eq RngIntElt then BadPrimes := [BadPrimes]; end if;
    K := BaseRing(f);
    if IsField(K) then K := RingOfIntegers(K); f:= ChangeRing(f,K); end if;
    T := [Trace(CompanionMatrix(f)^i):i in [0..2*Degree(f)-2]];
    if Verbose then printf "Computing splitting field with Galois action..."; t := Cputime(); end if;
    L,r,G := GaloisSplittingField(f);
    if Verbose then printf "%.3os\n", Cputime()-t; printf "Computing conjugacy classes..."; t := Cputime(); end if;
    C := ConjugacyClasses(G); C := [<c[1],c[2],c[3],Sort(CycleStructure(c[3]))>:c in C];
    if Degree(f) eq 1 then return <f,[],G,C,T>; end if;
    if Type(K) eq RngInt and Degree(f) eq 2 then return <f,[],G,C,T>; end if;
    RL<y> := PolynomialRing(L);  RK<X> := PolynomialRing(K);  RZ<x> := PolynomialRing(Integers());
    S := [];
    F := FieldOfFractions(K);
    if Verbose then printf "%.3os\n", Cputime()-t; printf "Generating list of test primes..."; t := Cputime(); end if;
    if F eq Rationals() then
        D := Integers()!Discriminant(f);
        for p in BadPrimes do D *:= p; end for;
        P := [p : p in PrimesInInterval(1,2^16) | not IsDivisibleBy(D,p)];
    else
        F := NumberField(F);
        D := Integers()!Norm(Discriminant(f));
        for p in BadPrimes do D *:= p; end for;
        P := &cat[[pp : pp in PrimeIdealsOverPrime(F,p)|Norm(pp) eq p]:p in PrimesInInterval(1,2^16) | not IsDivisibleBy(D,p)];
    end if;
    if Verbose then printf "%.3os\n", Cputime()-t; end if;
    if #Polys gt 0 then
        for g in Polys do
            h := RZ!g;
            if Verbose then printf "Checking h=%o...", h; t := Cputime(); end if;
            F := [RK!&*[y-&+[Evaluate(h,r[j])*r[j^g]:j in [1..#r]]:g in C[i][3]^G]:i in [1..#C]];
            if Verbose then printf "%.3os\n", Cputime()-t; end if;
            if #F eq #Set(F) then
                Append(~S,<h,F>);
                if Verbose then printf "Added h=%o to S.\n", h; end if;
            end if;
        end for;
        if Verbose then printf "Testing S of cardinality %o...", #S; t := Cputime(); end if;
        H,B := _frobenius_test(<f,S,G,C,T>,P);
        if Verbose then printf "%.3os, needed %o polys, found %o bad primes\n", Cputime()-t, #H, #B; end if;
        if Verbose then print Sort([n:n in {Norm(p):p in B}]); end if;
        assert #B eq 0;
        return <f,[r:r in S|r[1] in H],G,C,T>;
    else
        B := P; b:= #B; bt := Cputime();
        maxd := Degree(f)-1;
        E := {1..maxd};
        for k:=1 to maxd do
            if k le 2 then
                Ek := Sort([e:e in Subsets(E,k)],func<a,b|&+[2^i:i in a]-&+[2^i:i in b]>);
                Ek := [&+[x^i:i in e]:e in Ek];
            else
                // if no monomials or bionomials work, just pick random polys, these will only be used for a handful of primes
                Ek := [&+[Random([-maxd..maxd])*x^i: i in [0..maxd]]: j in [1..10]];
            end if;
            for h in Ek do
                if h eq x then continue; end if;
                if &+[Evaluate(h,a)*a:a in r] eq 0 then continue; end if; // quick abort
                if Verbose then printf "Checking h=%o...", h; t := Cputime(); end if;
                F := [RK!&*[y-&+[Evaluate(h,r[j])*r[j^g]:j in [1..#r]]:g in C[i][3]^G]:i in [1..#C]];
                if Verbose then printf "%.3os\n", Cputime()-t; end if;
                if #F eq #Set(F) then
                    Append(~S,<h,F>);
                    if Verbose then printf "Added h=%o to S.\n", h; end if;
                    if Verbose then printf "Testing S of cardinality %o...", #S; t := Cputime(); end if;
                    H,B := _frobenius_test(<f,S,G,C,T>,B);
                    if #B eq 0 then if Verbose then print "done!"; end if; return <f,S,G,C,T>; end if;
                    if Verbose then printf "%.3os, needed %o polys, found %o bad primes\n", Cputime()-t, #H, #{Norm(p):p in B}; end if;
                    if Verbose then print Sort([n:n in {Norm(p):p in B}]); end if;
                    if #B eq b then
                        Prune(~S);
                        if Cputime()-bt gt 60 then
                            printf "Problematic primes: %o.  You can exclude these using the BadPrimes parameter.\n", Sort([n:n in {Norm(p):p in B}]);
                            bt := Cputime();
                        end if;
                    else
                        b := #B;
                        bt := Cputime();
                    end if;
                end if;
            end for;
        end for;
        printf "Failed for primes: %o.  You can exclude these using the BadPrimes parameter.", Sort([n:n in {Norm(p):p in B}]);
        return <f,S,G,C,T>;
    end if;
    printf "Failed for primes: %o.  You can exclude these using the BadPrimes parameter.", Sort([n:n in {Norm(p):p in B}]);
    assert false;
end intrinsic;

intrinsic FrobeniusContextToString(ctx::Tup) -> MonStgElt
{ Converts the specified Frobneius context to a string that can be saved to a file. }
    function sprint(s) return Join(Split(Join(Split(Sprintf("%o",s)," "),""),"\n"),""); end function;
    return sprint(<<Eltseq(ctx[1]),Eltseq(DefiningPolynomial(FieldOfFractions(BaseRing(ctx[1]))))>,
        [<Eltseq(r[1]),[Eltseq(f):f in r[2]]>:r in ctx[2]],
        [Eltseq(g):g in Generators(ctx[3])],
        [<r[1],r[2],Eltseq(r[3]),r[4]> : r in ctx[4]],
        ctx[5]>);
end intrinsic;

intrinsic FrobeniusContextFromString(s::MonStgElt:K:=0) -> Tup
{ Reconstructs a Frobenius context from a string returned by FrobeniusContextTString. }
    raw := eval(s);
    assert #raw in [5,6];
    assert #raw[1] eq 2;
    RZ<x> := PolynomialRing(Integers());
    if Type(K) eq FldNum or Type(K) eq FldRat or Type(K) eq RngOrd or Type(K) eq RngInt then
        assert raw[1][2] eq Eltseq(DefiningPolynomial(K));
        if IsField(K) then K := RingOfIntegers(K); end if;
    else
        K := RingOfIntegers(NumberField(RZ!raw[1][2]));
    end if;
    RK<X> := PolynomialRing(K);
    f := RK!raw[1][1];
    S := [<RZ!r[1],[RK!g:g in r[2]]>:r in raw[2]];
    G := SymmetricGroup(Degree(f));
    G := sub<G|raw[3]>;
    C := [<r[1],r[2],G!r[3],r[4]>:r in raw[4]];
    return <f,S,G,C,raw[5]>;
end intrinsic;