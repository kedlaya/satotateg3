Attach("common.m");
prefix := "A32b";
M<z6> := NumberField([1,-1,1]);  OM := RingOfIntegers(M); RM := PolynomialRing(M);
H := sub<GL(2,M)|[z6-1,0,0,z6-1],[-1,-2,1,1],[z6,2*z6-2,0,1-z6]>;
assert IdentifyGroup(H) eq <72,25>;

E:=EllipticCurve([0,-1]);
assert &and[Trace(SexticPsi(-1,p)) eq TraceOfFrobenius(E,p):p in PrimesInInterval(1,10000)|p mod 3 eq 1];

RZ<x> := PolynomialRing(Integers());
f1 := x^16 - 18*x^12 + 101*x^10 + 99*x^8 - 1098*x^6 + 3043*x^4 - 738*x^2 + 81;
assert gg(f1) eq <48,33>;
f2 := x^3 - 3;
assert gg(f2) eq <6,1>;
g := x^6 - 3*x^5 + 10*x^4 - 15*x^3 - 41*x^2 + 48*x + 27;
assert gg(g) eq <24,13>;                // Endomorphism field K of A_{-3[2]b} with ST group J(T)
f := x^48 - 909*x^42 + 353214*x^36 - 71099181*x^30 + 8065606914*x^24 - 508642563096*x^18 + 15439889573460*x^12 - 267151139172*x^6 + 3486784401;
G := GaloisGroup(f);
assert IdentifyGroup(G) eq <144, 127>; // Endomorphism field L of E_{-3} x A_{-3[2]b} with ST group J(B(T,3))
P := Set(PrimeDivisors(nfdisc(f1)*nfdisc(f2)));

filename := prefix cat "FrobeniusContext.txt";
if not OpenTest(filename, "r") then
  print "This should take about 60 seconds...";
  time ctx := FrobeniusContext(Factorization(ChangeRing(f,OM))[2][1]:Verbose); // f splits over M
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
assert #S1 eq 2;
assert S1[1] eq cc(S1[2]) where cc := func<s|[[ComplexConjugate(c):c in r]:r in s]>;

S1 := [[Reverse(RM!h):h in s]:s in S1];
S2 := [2,-2,-1,0,0,1];

RZ<T> := PolynomialRing(Integers());
function Euler(p,n)
  if p mod 3 eq 1 then return RZ!Norm(Evaluate(S1[n][FrobeniusElement(ctx,psi)],psi*RM.1)) where psi:=SexticPsi(-1,p:z6:=z6);
  else return RZ![1,0,S2[FrobeniusOrder(f1,p) div 2]*p,0,p^2]; end if; // at inert pirmes the order is determined by the splitting of f1
end function;

L7 := [Coefficients(Euler(7,i))[2..3] : i in [1..#S1]];
assert Set(L7) eq { [-1,-6], [5,18] };
i := Index(L7,[-1,-6]);

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
