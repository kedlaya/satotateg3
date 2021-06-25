Attach("common.m");
prefix := "A4a";
M<z4> := NumberField([1,0,1]);  OM := RingOfIntegers(M); RM := PolynomialRing(M);
H := sub<GL(2,M)|[z4,0,0,z4],[0,1,1,0],[0,-1,1,-1]>;
assert IdentifyGroup(H) eq <24,5>;

E:=EllipticCurve([-2,0]);
assert &and[Trace(QuarticPsi(2,p:z4:=z4)) eq TraceOfFrobenius(E,p):p in PrimesInInterval(1,10000)|p mod 4 eq 1];

RZ<x> := PolynomialRing(Integers());
h := x^6 - 3*x^2 - 6;               // Endomorphism field K of A_{-4a} with ST group D_{6,1}
assert gg(h) eq <12,4>;
f := x^12 - 3*x^8 + 24*x^4 - 48;    // Endomorphism field L of E_{-4,2} x A_{-4a} with ST group J(B(3,4))
G := GaloisGroup(f);
assert IdentifyGroup(G) eq <48, 38>;
P := Set(PrimeDivisors(nfdisc(f)));

filename := prefix cat "FrobeniusContext.txt";
if not OpenTest(filename, "r") then
  print "This should take 10-15 seconds...";
  time ctx := FrobeniusContext(ChangeRing(f,OM):Verbose); // f is irreducible over M
  fp := Open(filename,"w"); Puts(fp,FrobeniusContextToString(ctx)); Flush(fp); delete fp;
else
  ctx := FrobeniusContextFromString(Read(filename):K:=OM);
end if;

// We only want to consider reps for which rho_M and rho_M^c' are O_M-equivalent, per Remark 7.10 (with V the identity)
V := GL(2,M)![1,0,0,1];
// Find an element s of G that induces the automorphism of H corresponding to Gal(kM/k) (this will be an order 2 element not in H)
GH := [K`subgroup: K in Subgroups(G:IndexEqual:=2) | IsIsomorphic(K`subgroup,H)]; assert #GH eq 1; GH := GH[1];
C := ConjugacyClasses(G); s := [c[3]:c in C|CycleStructure(c[3]) eq [ <2, 3>, <1, 6> ] and #(Conjugates(G,c[3]) meet Set(GH)) eq 0]; assert #s eq 1; s := s[1];
// We can index the reps rho_M with image H by automorphism of H, and we only need to distinguish them up to charpolys
A := Elements(AutomorphismGroup(H));
_,psi := IsIsomorphic(ctx[3],GH);  _,phi := IsIsomorphic(GH,H);
rho := func<a,g|a(phi(psi(g)))>;
rhoc := func<a,g|V^-1*cc(a(phi(s*psi(g)*s)))*V> where cc := func<m|GL(2,M)![ComplexConjugate(x):x in Eltseq(m)]>;
AA := [a:a in A| &and[rho(a,g) eq rhoc(a,g):g in Generators(ctx[3])]];
S1 := Sort([s:s in {[Eltseq(CharacteristicPolynomial(rho(a,c[3]))):c in ctx[4]]:a in AA}]);
assert #S1 eq 2 and S1[1] eq cc(S1[2]) where cc := func<s|[[ComplexConjugate(c):c in r]:r in s]>;

S1 := [[Reverse(RM!h):h in s]:s in S1];
S2 := [2,-2,-1,0,0,1];

function Euler(p,i)
  if p mod 4 eq 1 then return RZ!Norm(Evaluate(S1[i][FrobeniusElement(ctx,psi)],psi*RM.1)) where psi:=QuarticPsi(2,p:z4:=z4);
  else return RZ![1,0,S2[FrobeniusOrder(f,p) div 2]*p,0,p^2]; end if;
end function;

L29 := [Coefficients(Euler(29,i))[2..3] : i in [1..#S1]];
assert Set(L29) eq { [20,158], [-20,158] };
i := Index(L29,[20,158]);

E := EllipticCurve([-2,0]);
C := HyperellipticCurve (x^8 - x^7 - 7*x^6 - 7*x^5 - 35*x^4 - 35*x^3 - 21*x^2 + 3*x + 6);
for p in PrimesInInterval(1,200) do
    if p in P then continue; end if;
    Lp := Euler(p,i);
    if Lp*LPolynomial(ChangeRing(E,GF(p))) ne LPolynomial(ChangeRing(C,GF(p))) then
        printf "Failed at p=%o: %o != %o\n", p, Lp, LPolynomial(ChangeRing(C,GF(p))) / LPolynomial(ChangeRing(E,GF(p)));
        assert false;
    end if;
end for;

if assigned jobs and assigned jobid and assigned n then
  jobs:=StringToInteger(jobs);  jobid:=StringToInteger(jobid); N:=2^StringToInteger(n);
  min := Floor(jobid*N/jobs)+1;  max := Floor((jobid+1)*N/jobs);
  t := Cputime();
  filename := Sprintf("lpdata_%o_%o.txt",prefix,jobid);
  A := [[p] cat Coefficients(Euler(p,i))[2..3]:p in PrimesInInterval(min,max) | not p in P];
  fp := Open(filename,"w"); for a in A do Puts(fp,Sprintf("%o,%o,%o",a[1],a[2],a[3])); end for; Flush(fp);
  printf "Wrote %o records to %o for p in [%o,%o] using %o CPU secs\n", #A, filename, min, max, Cputime()-t;
  exit;
end if;
