Attach("common.m");
prefix := "A32a";
M<z6> := NumberField([1,-1,1]);  OM := RingOfIntegers(M); RM := PolynomialRing(M);
H := sub<GL(2,M)|[z6-1,0,0,z6-1],[-1,-2,0,1],[-1,0,1,1]>;
assert IdentifyGroup(H) eq <24,10>;

E:=EllipticCurve([0,-1]);
assert &and[Trace(SexticPsi(-1,p)) eq TraceOfFrobenius(E,p):p in PrimesInInterval(1,10000)|p mod 3 eq 1];

RZ<x> := PolynomialRing(Integers());
f1 := x^8 - 12*x^6 + 12*x^4 - 52*x^2 - 3;
assert gg(f1) eq <16,7>;  // D8
f2 := x^3-2;
assert gg(f2) eq <6,1>;   // S3
h := x^4-4*x^2-3;
assert gg(h) eq <8,3>;                                 // D4 endomorphism field K of A_{-3[2]a} with ST group D_{4,1}
f := x^24 + 488*x^18 + 63384*x^12 - 5872*x^6 + 1372;
G := GaloisGroup(f);
assert IdentifyGroup(G) eq <48,15>;                    // endomorphism field L of E_{-3} x A_{-3[2]a} with ST group J_s(B(1,12))
P := Set(PrimeDivisors(nfdisc(f1)) cat PrimeDivisors(nfdisc(f2)));

filename := prefix cat "FrobeniusContext.txt";
if not OpenTest(filename, "r") then
  print "This should take about 10 seconds...";
  time ctx := FrobeniusContext(Factorization(ChangeRing(f,OM))[2][1]:Verbose); // f splits over OM=Q(z6)
  fp := Open(filename,"w"); Puts(fp,FrobeniusContextToString(ctx)); Flush(fp); delete fp;
else
  ctx := FrobeniusContextFromString(Read(filename):K:=OM);
end if;

// We only want to consider reps for which rho_M and rho_M^c' are O_M-equivalent, per Remark 7.12, which we are applying with m=2
cp := GL(2,M)![0,2,1,0];
// Find an element s of G that induces the automorphism of H corresponding to Gal(kM/k) (this will be an order 2 element not in H)
GH := [K`subgroup: K in Subgroups(G:IndexEqual:=2) | IsIsomorphic(K`subgroup,H)]; assert #GH eq 1; GH := GH[1];
C := ConjugacyClasses(G); s := [c[3]:c in C|Order(c[3]) eq 2 and #(Conjugates(G,c[3]) meet Set(GH)) eq 0]; assert #s eq 1; s := s[1];
// We can index the reps rho_M with image H by automorphism of H, and we only need to distinguish them up to charpolys
A := Elements(AutomorphismGroup(H));
_,psi := IsIsomorphic(ctx[3],GH);  _,phi := IsIsomorphic(GH,H);
rho := func<a,g|a(phi(psi(g)))>;
rhocp := func<a,g|cp^-1*cc(a(phi(s*psi(g)*s)))*cp> where cc := func<m|GL(2,M)![ComplexConjugate(x):x in Eltseq(m)]>;
AA := [a:a in A| &and[rho(a,g) eq rhocp(a,g):g in Generators(ctx[3])]];
S1 := Sort([s:s in {[Eltseq(CharacteristicPolynomial(rho(a,c[3]))):c in ctx[4]]:a in AA}]);
assert #S1 eq 2 and S1[1] eq cc(S1[2]) where cc := func<s|[[ComplexConjugate(c):c in r]:r in s]>;

S1 := [[Reverse(RM!h):h in s]:s in S1];
S2 := [2,-2,-1,0,0,1];

function Euler(p,n)
  if p mod 3 eq 1 then return RZ!Norm(Evaluate(S1[n][FrobeniusElement(ctx,psi)],psi*RM.1)) where psi:=SexticPsi(-1,p:z6:=z6);
  else return RZ![1,0,S2[FrobeniusOrder(f,p) div 2]*p,0,p^2]; end if;
end function;

// pick a particular representation to make the results consistent
L13 := [Coefficients(Euler(13,i))[2..3] : i in [1..#S1]];
assert Set(L13) eq { [0,1], [0,-23] };
i := Index(L13,[0,-23]);

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
