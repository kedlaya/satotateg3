Attach("common.m");
prefix := "A4b";
M<z4> := NumberField([1,0,1]);  OM := RingOfIntegers(M); RM := PolynomialRing(M);
H := sub<GL(2,M)|[z4,0,-z4,1],[1,1,0,z4]>;
assert IdentifyGroup(H) eq <96,67>;

E:=EllipticCurve([-1,0]);
assert &and[Trace(QuarticPsi(1,p)) eq TraceOfFrobenius(E,p):p in PrimesInInterval(1,10000)|p mod 4 eq 1]; // sanity check

RZ<x> := PolynomialRing(Integers());
g := x^6 - 4*x^4 - 59; 
assert gg(g) eq <48,48>;                // minpoly(sqrt(b)), Endomorphism field K of A_{-4b} with ST group J(O)
f := x^24 - 40*x^16 + 430*x^12 - 240*x^8 + 344*x^4 - 59;
G := GaloisGroup(f);
assert IdentifyGroup(G) eq <192,988>;   // Endomorphism field L of E_{-4} x A_{-4b} with ST group J(B(O,2))
P := Set(PrimeDivisors(5*nfdisc(f)));   // 5 is problematic for FrobeniusContext, so we skip it

filename := prefix cat "FrobeniusContext.txt";
if not OpenTest(filename, "r") then
  print "This should take 30-60 seconds...";
  time ctx := FrobeniusContext(ChangeRing(f,OM):BadPrimes:=[5],Verbose); // f is irreducible over OM=Q(i)
  fp := Open(filename,"w"); Puts(fp,FrobeniusContextToString(ctx)); Flush(fp); delete fp;
else
  ctx := FrobeniusContextFromString(Read(filename):K:=OM);
end if;

// We only want to consider reps for which rho_M and rho_M^c' are O_M-equivalent, per Remark 7.10 (with V the identity)
V := GL(2,M)![0,1,1,0];
// Find an element s of G that induces the automorphism of H corresponding to Gal(kM/k) (this will be an order 2 element not in H)
GH := [K`subgroup: K in Subgroups(G:IndexEqual:=2) | IsIsomorphic(K`subgroup,H) and IsIsomorphic(NumberField(GaloisSubgroup(f,K`subgroup)),M)]; assert #GH eq 1; GH := GH[1];
C := ConjugacyClasses(G); s := [c[3]:c in C|CycleStructure(c[3]) eq [<2,11>,<1,2>] and #(Conjugates(G,c[3]) meet Set(GH)) eq 0]; assert #s eq 1; s := s[1];
// We can index the reps rho_M with image H by automorphism of H, and we only need to distinguish them up to charpolys
A := Elements(AutomorphismGroup(H));
_,psi := IsIsomorphic(ctx[3],GH);  _,phi := IsIsomorphic(GH,H);
rho := func<a,g|a(phi(psi(g)))>;
rhoc := func<a,g|V^-1*cc(a(phi(s*psi(g)*s)))*V> where cc := func<m|GL(2,M)![ComplexConjugate(x):x in Eltseq(m)]>;
AA := [a:a in A| &and[rho(a,g) eq rhoc(a,g):g in Generators(ctx[3])]];
S1 := Sort([s:s in {[Eltseq(CharacteristicPolynomial(rho(a,c[3]))):c in ctx[4]]:a in AA}]);
assert #S1 eq 4 and S1[1] eq cc(S1[4])and S1[2] eq cc(S1[3]) where cc := func<s|[[ComplexConjugate(c):c in r]:r in s]>;

S1 := [[Reverse(RM!h):h in s]:s in S1];
S2 := [2,-2,-1,0,0,1];

function Euler(p,i)
  if p mod 4 eq 1 then return RZ!Norm(Evaluate(S1[i][FrobeniusElement(ctx,psi)],psi*RM.1)) where psi:=QuarticPsi(1,p:z4:=z4);
  else return RZ![1,0,S2[FrobeniusOrder(f,p) div 2]*p,0,p^2]; end if;
end function;

// pick the representation matching Jac(C)/E
L61 := [Coefficients(Euler(61,i))[2..3] : i in [1..#S1]];
assert Set(L61) eq { [2,2], [-2,2], [22,242], [-22,242] };
i := Index(L61,[2,2]);

E := EllipticCurve("236a1");
P2<U,V,W> := ProjectiveSpace(Rationals(),2); C:=Curve(P2,-V^4+U^4+2*U^2*W^2+U*W^3);
for p in PrimesInInterval(1,100) do
  if p in P then continue; end if;
  Lp := Euler(p,i);
  if Lp*LPolynomial(ChangeRing(E,GF(p))) ne LPolynomial(ChangeRing(C,GF(p))) then
    printf "Failed at p=%o: %o != %o\n", p, Lp, LPolynomial(ChangeRing(C,GF(p))) / LPolynomial(ChangeRing(E,GF(p))); assert false;
  end if;
end for;

if assigned jobs and assigned jobid and assigned n then
  jobs:=StringToInteger(jobs);  jobid:=StringToInteger(jobid); N:=2^StringToInteger(n);
  min := Floor(jobid*N/jobs)+1;  max := Floor((jobid+1)*N/jobs);
  t := Cputime();
  filename := Sprintf("lpdata_%o_%o.txt",prefix,jobid);
  A := [[p] cat Coefficients(Euler(p,i))[2..3]:p in PrimesInInterval(Floor(jobid*N/jobs)+1,Floor((jobid+1)*N/jobs)) | not p in P];
  fp := Open(filename,"w"); for a in A do Puts(fp,Sprintf("%o,%o,%o",a[1],a[2],a[3])); end for; Flush(fp);
  printf "Wrote %o records to %o for p in [%o,%o] using %o CPU secs\n", #A, filename, min, max, Cputime()-t;
  exit;
end if;
