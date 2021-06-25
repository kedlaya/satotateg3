Attach("common.m");
prefix := "TwistE7";
M<z> := NumberField([7,0,1]);  OM := RingOfIntegers(M); RM := PolynomialRing(M);
H := sub<GL(3,M)|[-1,0,(1-z)/2,0,-1,(-1-z)/2,0,0,1],[2,(-1+z)/2,(-3+z)/2,1,0,0,(1+z)/2,-1,-1]>;
assert IdentifyGroup(H) eq <168,42>;

E := EllipticCurve([1,-1,0,-2,-1]); // 49a1

Psi := function(p)
    if p eq 2 then return 1; end if;
    if p mod 7 in [3,5,6] then return 0; end if;
    _,a,b := NormEquation(7,p);
    if (a-b) mod 4 ne 1 then if a mod 2 eq 0 then b := -b; else a:= -a; end if; end if;
    assert (a-b) mod 4 eq 1;
    e := IsSquare(GF(p)!a/b) select 1 else -1; // 6*B*gamma_3 = -(3402)^2*sqrt(-7) = a/b mod p, since sqrt(-7)=-a/b mod p
    if IsEven(a) then e := -e; end if;
    return e*(a+b*z);
end function;

assert &and[Trace(Psi(p)) eq TraceOfFrobenius(E,p):p in PrimesInInterval(1,10000)];

RZ<x> := PolynomialRing(Integers());
f := x^8 + 4*x^7 + 21*x^4 + 18*x + 9;   // Endomorphism field K of E_{-7}^3 with ST group J(E(168))
G := GaloisGroup(f);
assert IdentifyGroup(G) eq <336,208>;
P := PrimeDivisors(nfdisc(f));

filename := prefix cat "FrobeniusContext.txt";
if not OpenTest(filename, "r") then
  print "This should take about five minutes...";
  time ctx := FrobeniusContext(ChangeRing(f,OM):Verbose); // f is irreducible over OM=Q(i)
  fp := Open(filename,"w"); Puts(fp,FrobeniusContextToString(ctx)); Flush(fp); delete fp;
else
  ctx := FrobeniusContextFromString(Read(filename):K:=OM);
end if;

// We only want to consider reps for which rho_M and rho_M^c' are O_M-equivalent, per Remark 7.10 (with V the identity)
V := GL(3,M)![0,1,0,1,0,0,0,0,-1];
// Find an element s of G that induces the automorphism of H corresponding to Gal(kM/k) (this will be an order 2 element not in H)
GH := [K`subgroup: K in Subgroups(G:IndexEqual:=2) | IsIsomorphic(K`subgroup,H)]; assert #GH eq 1; GH := GH[1];
C := ConjugacyClasses(G); s := [c[3]:c in C|Order(c[3]) eq 2 and #(Conjugates(G,c[3]) meet Set(GH)) eq 0]; assert #s eq 1; s := s[1];
// We can index the reps rho_M with image H by automorphism of H, and we only need to distinguish them up to charpolys
A := Elements(AutomorphismGroup(H));
_,psi := IsIsomorphic(ctx[3],GH);  _,phi := IsIsomorphic(GH,H);
rho := func<a,g|a(phi(psi(g)))>;
rhoc := func<a,g|V^-1*cc(a(phi(s*psi(g)*s)))*V> where cc := func<m|GL(3,M)![ComplexConjugate(x):x in Eltseq(m)]>;
AA := [a:a in A| &and[rho(a,g) eq rhoc(a,g):g in Generators(ctx[3])]];
S1 := Sort([s:s in {[Eltseq(CharacteristicPolynomial(rho(a,c[3]))):c in ctx[4]]:a in AA}]);
assert #S1 eq 2 and S1[1] eq cc(S1[2]) where cc := func<s|[[ComplexConjugate(c):c in r]:r in s]>;

S1 := [[Reverse(RM!h):h in s]:s in S1];
S2 := [3,-1,0,1,0,2];

function Euler(p,i)
  if p mod 7 in [1,2,4] then return RZ!Norm(Evaluate(S1[i][FrobeniusElement(ctx,psi)],psi*RM.1)) where psi:=Psi(p);
  else return RZ![1,0,c*p,0,c*p^2,0,p^3] where c:= S2[FrobeniusOrder(f,p) div 2]; end if;
end function;

// pick the representation matching Jac(C)
L11 := [Coefficients(Euler(11,i))[2..4] : i in [1..#S1]];
assert Set(L11) eq { [-5,-3,57], [9,53,211] };
i := Index(L11,[9,53,211]);

R<u,v,w>:=ProjectiveSpace(Rationals(),2);
C:=Curve(R,2*u^3*v-2*u^3*w-3*u^2*w^2-2*u*v^3-2*u*w^3-4*v^3*w+3*v^2*w^2-v*w^3);
for p in PrimesInInterval(1,100) do if not p in P then assert Euler(p,i) eq LPolynomial(ChangeRing(C,GF(p))); end if; end for;

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
