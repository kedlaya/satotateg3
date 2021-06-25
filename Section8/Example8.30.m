Attach("common.m");
prefix := "ResE4a";
M<z4> := NumberField([1,0,1]);  OM := RingOfIntegers(M); RM := PolynomialRing(M);
H := sub<GL(3,M)|[0,1,0,0,0,1,1,0,0],[0,1,0,1,0,0,0,0,1],[z4,0,0,0,1,0,0,0,1]>;
assert IdentifyGroup(H) eq <384,5557>;
assert IdentifyGroup(WreathProduct(CyclicGroup(4),SymmetricGroup(3))) eq <384,5557>;

E:=EllipticCurve([-1,0]);
assert &and[Trace(QuarticPsi(1,p)) eq TraceOfFrobenius(E,p):p in PrimesInInterval(1,10000)|p mod 4 eq 1];

RZ<x> := PolynomialRing(Integers());
f := x^3-x^2+x-2;
assert gg(f) eq <6,1>;
D := Integers()!Discriminant(f);
K<a>:=NumberField(f);
E:=EllipticCurve([-a,1]);
P := PrimeDivisors(Integers()!Norm(Discriminant(E)));

h := x^12 + x^8 + 2*x^4 + 4; // endomorphism field L of E_{-4}^3 and Res E_{-4,alpha}
assert gg(h) eq <192,956>;

g := Evaluate(f,x^4);
G := GaloisGroup(g);
assert IdentifyGroup(G) eq <768, 1088009>;
P := PrimeDivisors(5*nfdisc(g));

filename := prefix cat "FrobeniusContext.txt";
if not OpenTest(filename, "r") then
  print "This should take 15-20 minutes...";
  time ctx := FrobeniusContext(ChangeRing(g,OM):Verbose); // f is irreducible over OM=Q(i)
  fp := Open(filename,"w"); Puts(fp,FrobeniusContextToString(ctx)); Flush(fp); delete fp;
else
  ctx := FrobeniusContextFromString(Read(filename):K:=OM);
end if;

// We only want to consider reps for which rho_M and rho_M^c' are O_M-equivalent, per Remark 7.10 (with V the identity)
V := GL(3,M)![1,0,0,0,1,0,0,0,1];
// Find an element s of G that induces the automorphism of H corresponding to Gal(kM/k) (this will be an order 2 element not in H)
GH := [K`subgroup: K in Subgroups(G:IndexEqual:=2) | IsIsomorphic(K`subgroup,H)]; assert #GH eq 1; GH := GH[1];
C := ConjugacyClasses(G); s := [c[3]:c in C|CycleStructure(c[3]) eq [<2,3>,<1,6>] and #(Conjugates(G,c[3]) meet Set(GH)) eq 0]; assert #s eq 1; s := s[1];
// We can index the reps rho_M with image H by automorphism of H, and we only need to distinguish them up to charpolys
A := Elements(AutomorphismGroup(H));
_,psi := IsIsomorphic(ctx[3],GH);  _,phi := IsIsomorphic(GH,H);
rho := func<a,g|a(phi(psi(g)))>;
rhoc := func<a,g|V^-1*cc(a(phi(s*psi(g)*s)))*V> where cc := func<m|GL(3,M)![ComplexConjugate(x):x in Eltseq(m)]>;
AA := [a:a in A| &and[rho(a,g) eq rhoc(a,g):g in Generators(ctx[3])]];
S1 := Sort([s:s in {[Eltseq(CharacteristicPolynomial(rho(a,c[3]))):c in ctx[4]]:a in AA}]);
assert #S1 eq 8 and S1[1] eq cc(S1[4]) and S1[2] eq cc(S1[3]) and S1[5] eq cc(S1[8]) and S1[6] eq cc(S1[7]) where cc := func<s|[[ComplexConjugate(c):c in r]:r in s]>;

S1 := [[Reverse(RM!h):h in s]:s in S1];
S2 := [3,-1,0,1,0,2];

function EulerTwist(p,i)
  if p mod 4 eq 1 then return RZ!Norm(Evaluate(S1[i][FrobeniusElement(ctx,psi)],psi*RM.1)) where psi:=QuarticPsi(1,p:z4:=z4);
  else return RZ![1,0,c*p,0,c*p^2,0,p^3] where c:= S2[FrobeniusOrder(g,p) div 2]; end if;
end function;

L529t := [<Coefficients(EulerTwist(5,i))[2..4],Coefficients(EulerTwist(29,i))[2..4]> : i in [1..8]];
assert Set(L529t) eq {<[-2,13,-16],[0,0,284]>,<[-2,-3,16],[0,0,-284]>,<[2,-3,-16],[0,0,284]>,<[2,-3,-16],[0,0,-284]>,
                     <[2,13,16],[0,0,284]>,<[-2,13,-16],[0,0,-284]>,<[2,13,16],[0,0,-284]>,<[-2,-3,16],[0,0,284]>};

R<T> := PolynomialRing(Integers());
function EulerRes(p:minpoly:=f)
    if p eq 83 then return (83*T^2+1)^3; end if;
    Fp:=GF(p);
    Rp<z>:=PolynomialRing(Fp);
    fp:=Rp!minpoly;
    S:=quo<Rp|fp>; t := S.1;
    if p mod 4 eq 1 then
        _,u,v := NormEquation(1,p);
        if IsEven(u) then s:=u; u:=v; v:=s; end if; // make u odd and v even
        if (u-1-v) mod 4 eq 2 then u := -u; end if; // make u-1 = v mod 4
        h := t^((p-1) div 4);
        A := [[u, 1, 1], [-u, -1, 1], [-v, u, v], [v, -u, v]];
        L := &*[R| (p*T^2 - 2*a[1]*T + 1)^Degree(GCD(fp,a[3]*Evaluate(h,z)-a[2])) : a in A];
        if Degree(L) eq 6 then return L; end if;
        if Degree(L) eq 0 then
            h := t^((p^3-1) div 4);  assert Degree(h) eq 0;
            s := u^3 - 3*u*v^2;  v := 3*u^2*v - v^3; u := s;
            A := [[u, 1, 1], [-u, -1, 1], [-v, u, v], [v, -u, v]];
            a := [a:a in A|a[3]*h eq a[2]][1];
            return p^3*T^6 - 2*a[1]*T^3 + 1;
        else
            assert Degree(L) eq 2;
            fp := ExactQuotient(fp,GCD(fp,Evaluate(h^4,z)-1));  assert Degree(fp) eq 2;
            S:=quo<Rp|fp>;  t := S.1;
            s := u^2-v^2; v := 2*u*v; u := s;
            A := [[u, 1, 1], [-u, -1, 1], [-v, u, v], [v, -u, v]];
            h := t^((p^2-1) div 4);  assert Degree(h) eq 0;
            a := [a:a in A|a[3]*h eq a[2]][1];
            return L*(p^2*T^4 - 2*a[1]*T^2 + 1);
        end if;
    else
        gp := GCD(fp,Evaluate(t^(p-1)-1,z));
        if Degree(gp) eq 0 then return p^3*T^6+1; end if;
        if Degree(gp) eq 3 then return (p*T^2+1)^3; end if;
        fp := ExactQuotient(fp,gp);  assert Degree(fp) eq 2;
        S:=quo<Rp|fp>;  t := S.1;
        h := t^((p^2-1) div 4);
        if h ne 1 and -h ne 1 then return (p*T^2+1)*(p^2*T^4+1); end if;
        return (p*T^2+1) * (h eq 1 select (p*T^2+1)^2 else (p*T^2-1)^2);
    end if;
end function;

twists := [MinimalPolynomial(c):c in [a,-a,a^3,-a^3,83^2*a,-83^2*a,83^2*a^3,-83^2*a^3]];
L529r := [<Coefficients(EulerRes(5:minpoly:=twists[i]))[2..4],Coefficients(EulerRes(29:minpoly:=twists[i]))[2..4]> : i in [1..8]];
assert Set(L529r) eq Set(L529t);
indexes := [Index(L529t,L529r[i]):i in [1..8]];

for p in PrimesInInterval(1,2^12) do if not p in P then assert &and[EulerRes(p:minpoly:=twists[i]) eq EulerTwist(p,indexes[i]): i in [1..8]]; end if; end for;
//print "Verified 8 twists match 8 restrictions.";

// Jac(C) does not match any of the twists/restrictions and we don't use it, it's just here for reference.
PP<u,v,w>:=ProjectiveSpace(Rationals(),2);
C:=Curve(PP,44*u^4+120*u^3*v+36*u^3*w+60*u^2*v*w+9*u^2*w^2-200*u*v^3+u*w^3-150*v^4-15*v^2*w^2);

if assigned jobs and assigned jobid and assigned n then
  jobs:=StringToInteger(jobs);  jobid:=StringToInteger(jobid); N:=2^StringToInteger(n);
  min := Floor(jobid*N/jobs)+1;  max := Floor((jobid+1)*N/jobs);
  t := Cputime();
  filename := Sprintf("lpdata_%o_%o.txt",prefix,jobid);
  A := [[p] cat Coefficients(EulerRes(p))[2..4]:p in PrimesInInterval(Floor(jobid*N/jobs)+1,Floor((jobid+1)*N/jobs)) | not p in P];
  fp := Open(filename,"w"); for a in A do Puts(fp,Sprintf("%o,%o,%o,%o",a[1],a[2],a[3],a[4])); end for; Flush(fp);
  printf "Wrote %o records to %o for p in [%o,%o] using %o CPU secs\n", #A, filename, min, max, Cputime()-t;
  exit;
end if;
