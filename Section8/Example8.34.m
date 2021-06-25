Attach("common.m");
prefix := "TwistE3";
M<z6> := NumberField([1,-1,1]);  OM := RingOfIntegers(M); RM := PolynomialRing(M);
H := sub<GL(3,M)|[1,z3^2,0,0,z3^2,0,0,0,1],[z3,0,0,z3,1,0,z3,0,1],[1,0,z3^2,0,1,0,0,0,z3^2]> where z3:=z6-1;
assert IdentifyGroup(H) eq <648,533>;

E:=EllipticCurve([0,-1]);
assert &and[Trace(SexticPsi(-1,p)) eq TraceOfFrobenius(E,p):p in PrimesInInterval(1,10000)|p mod 3 eq 1];

RZ<x> := PolynomialRing(Rationals());
g := x^9 + 168*x^6 + 1080*x^5 - 636*x^3 + 864*x^2 - 432*x + 8;  // Endomorphism field of E_{-3}^3 with ST group J(E(216))
assert gg(g) eq <432,734>;
f := 512*Evaluate(g,x^3/2);
G := GaloisGroup(f);
assert IdentifyGroup(G) eq <1296,2891>;
P := PrimeDivisors(7*nfdisc(f));

filename := prefix cat "FrobeniusContext.txt";
if not OpenTest(filename, "r") then
  print "This should take 8-12 hours...";
  time ctx := FrobeniusContext(ChangeRing(f,OM):Verbose,BadPrimes:=[7]); // f is irreducible over Q(zeta_3)
  fp := Open(filename,"w"); Puts(fp,FrobeniusContextToString(ctx)); Flush(fp); delete fp;
else
  ctx := FrobeniusContextFromString(Read(filename):K:=OM);
end if;

// We only want to consider reps for which rho_M and rho_M^c' are O_M-equivalent, per Remark 7.10 (with V the identity)
V := GL(3,M)![1,0,0,1,-1,0,1,0,-1];
// Find an element s of G that induces the automorphism of H corresponding to Gal(kM/k) (this will be an order 2 element not in H)
GH := [K`subgroup: K in Subgroups(G:IndexEqual:=2) | IsIsomorphic(K`subgroup,H)]; assert #GH eq 1; GH := GH[1];
C := ConjugacyClasses(G); s := [c[3]:c in C|Order(c[3]) eq 2 and #(Conjugates(G,c[3]) meet Set(GH)) eq 0]; assert #s eq 1; s := s[1];
// We can index the reps rho_M with image H by automorphism of H, and we only need to distinguish them up to charpolys
printf "Computing large automorphism group, this may take 20-30 seconds..."; t:=Cputime();
A := Elements(AutomorphismGroup(H));
printf "computed %o automorphisms in %.3o secs.\n", #A, Cputime()-t;
_,psi := IsIsomorphic(ctx[3],GH);  _,phi := IsIsomorphic(GH,H);
rho := func<a,g|a(phi(psi(g)))>;
rhoc := func<a,g|V^-1*cc(a(phi(s*psi(g)*s)))*V> where cc := func<m|GL(3,M)![ComplexConjugate(x):x in Eltseq(m)]>;
AA := [a:a in A| &and[rho(a,g) eq rhoc(a,g):g in Generators(ctx[3])]];
S1 := Sort([s:s in {[Eltseq(CharacteristicPolynomial(rho(a,c[3]))):c in ctx[4]]:a in AA}]);
assert #S1 eq 6 and S1[1] eq cc(S1[5]) and S1[2] eq cc(S1[6]) and S1[3] eq cc(S1[4]) where cc := func<s|[[ComplexConjugate(c):c in r]:r in s]>;

S1 := [[Reverse(RM!h):h in s]:s in S1];
S2 := [3,-1,0,1,0,2];

function Euler(p,i)
  if p mod 6 eq 1 then return RZ!Norm(Evaluate(S1[i][FrobeniusElement(ctx,psi)],psi*RM.1)) where psi:=SexticPsi(-1,p:z6:=z6);
  else return RZ![1,0,c*p,0,c*p^2,0,p^3] where c:= S2[FrobeniusOrder(g,p) div 2]; end if;
end function;

L19 := [Coefficients(Euler(19,i))[2..4] : i in [1..#S1]];
assert Set(L19) eq { [1,-7,-26], [1,8,-11], [7,-7,-182], [7,56,259], [-8,8,88], [-8,56,-296] };
i := Index(L19,[1,-7,-26]);

if assigned jobs and assigned jobid and assigned n then
  jobs:=StringToInteger(jobs);  jobid:=StringToInteger(jobid); N:=2^StringToInteger(n);
  min := Floor(jobid*N/jobs)+1;  max := Floor((jobid+1)*N/jobs);
  t := Cputime();
  filename := Sprintf("lpdata_%o_%o.txt",prefix,jobid);
  A := [[p] cat Coefficients(Euler(p,i))[2..4]:p in PrimesInInterval(Floor(jobid*N/jobs)+1,Floor((jobid+1)*N/jobs)) | not p in P];
  fp := Open(filename,"w"); for a in A do Puts(fp,Sprintf("%o,%o,%o,%o",a[1],a[2],a[3],a[4])); end for; Flush(fp);
  printf "Wrote %o records to %o for p in [%o,%o] using %o CPU secs\n", #A, filename, min, max, Cputime()-t;
  exit;
end if;
