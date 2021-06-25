prefix := "A8";
M<z> := NumberField([2,0,1]);  OM := RingOfIntegers(M); RM := PolynomialRing(M);
H := sub<GL(2,M)|[-1,0,1-z,1],[-1,-1,1,0]>;
assert IdentifyGroup(H) eq <48,29>;

E:=EllipticCurve(x^3 + x^2 - 3*x + 1);
Psi := function(p)
  if p mod 8 gt 4 then return 0; end if;
  _,a,b := NormEquation(2,p);
  if a mod 4 ne 1 then a := -a; end if; // we can always make a = 1 mod 4 (since a most be odd)
  if b mod 4 eq 3 then b := -b; end if; // we can always make b != 3 mod 4
  e := IsSquare(GF(p)!2*a/b) select 1 else -1;
  if b mod 4 eq 2 then e := -e; end if;
  return e*(a+b*z);
end function;

assert &and[Trace(Psi(p)) eq TraceOfFrobenius(E,p) : p in PrimesInInterval(3,10000)];

RZ<x> := PolynomialRing(Integers());
h := x^6 + 4*x^4 - 22;      // Endomorphism field K of A_{-8} with ST group J(O)
assert gg(h) eq <48,48>;
f := x^16 - 44*x^12 - 308*x^10 - 990*x^8 - 1936*x^6 - 2662*x^4 + 9196*x^2 + 20449;
G := GaloisGroup(f);
assert IdentifyGroup(G) eq <96,193>;   // Endomorphism field L of E_{-8} x A_{-8} with ST group J(B(O,1))
P := Set(PrimeDivisors(nfdisc(f)));

filename := prefix cat "FrobeniusContext.txt";
if not OpenTest(filename, "r") then
  print "This should take about 20 seconds...";
  time ctx := FrobeniusContext(ChangeRing(f,OM):Verbose); // f is irreducible over M
  fp := Open(filename,"w"); Puts(fp,FrobeniusContextToString(ctx)); Flush(fp); delete fp;
else
  ctx := FrobeniusContextFromString(Read(filename):K:=OM);
end if;

// We only want to consider reps for which rho_M and rho_M^c' are O_M-equivalent, per Remark 7.10 (with V the identity)
V := GL(2,M)![0,1,1,0];
// Find an element s of G that induces the automorphism of H corresponding to Gal(kM/k) (this will be an order 2 element not in H)
GH := [K`subgroup: K in Subgroups(G:IndexEqual:=2) | IsIsomorphic(K`subgroup,H) and IsIsomorphic(NumberField(GaloisSubgroup(f,K`subgroup)),M)]; assert #GH eq 1; GH := GH[1];
C := ConjugacyClasses(G); s := [c[3]:c in C|CycleStructure(c[3]) eq [<2,6>,<1,4>] and #(Conjugates(G,c[3]) meet Set(GH)) eq 0]; assert #s eq 1; s := s[1];
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
  if p mod 8 in [1,3] then return RZ!Norm(Evaluate(S1[i][FrobeniusElement(ctx,psi)],psi*RM.1)) where psi:=Psi(p);
  else return RZ![1,0,S2[FrobeniusOrder(f,p) div 2]*p,0,p^2]; end if;
end function;

L17 := [Coefficients(Euler(17,i))[2..3] : i in [1..#S1]];
assert Set(L17) eq { [8,32], [-8,32] };
i := Index(L17,[8,32]);

C:=HyperellipticCurve(x^6-5*x^4+10*x^3-5*x^2+2*x-1);
for p in PrimesInInterval(3,200) do if not p in P then assert LPolynomial(ChangeRing(C,GF(p))) eq Euler(p,i); end if; end for;

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
