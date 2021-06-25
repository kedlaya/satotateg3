Attach("common.m");
prefix := "ResE3a";
M<z6> := NumberField([1,-1,1]);  OM := RingOfIntegers(M); RM := PolynomialRing(M);
H := sub<GL(3,M)|[0,1,0,0,0,1,1,0,0],[0,1,0,1,0,0,0,0,1],[z6,0,0,0,1,0,0,0,1]>;
assert IdentifyGroup(H) eq <1296,1827>;
assert IdentifyGroup(WreathProduct(CyclicGroup(6),SymmetricGroup(3))) eq <1296,1827>;

E:=EllipticCurve([0,-1]);
assert &and[Trace(SexticPsi(-1,p)) eq TraceOfFrobenius(E,p):p in PrimesInInterval(1,10000)|p mod 3 eq 1];

RZ<x> := PolynomialRing(Rationals());
f := x^3-x^2+x-2;
K<a>:=NumberField(f);
E:=EllipticCurve([0,-a]);
P := PrimeDivisors(Integers()!Norm(Discriminant(E)));

g := Evaluate(f,x^6);
G := GaloisGroup(g);
assert TransitiveGroupIdentification(G) eq 396;                     // 18T396
h := x^36 + x^30 + 27*x^24 - 42*x^18 + 28*x^12 - 384*x^6 + 4096;    // Endomorphism field K of Res E_{-3,alpha} with ST group J(D(6,6))
assert gg(h) eq <432,523>;

filename := prefix cat "FrobeniusContext.txt";
if not OpenTest(filename, "r") then
  print "This should take 24-48 hours...";
  time ctx := FrobeniusContext(ChangeRing(g,OM):Verbose); // g is irreducible over Q(zeta_3)
  fp := Open(filename,"w"); Puts(fp,FrobeniusContextToString(ctx)); Flush(fp); delete fp;
else
  ctx := FrobeniusContextFromString(Read(filename):K:=OM);
end if;

// We only want to consider reps for which rho_M and rho_M^c' are O_M-equivalent, per Remark 7.10 (with V the identity)
V := GL(3,M)![1,0,0,0,1,0,0,0,1];
// Find an element s of G that induces the automorphism of H corresponding to Gal(kM/k) (this will be an order 2 element not in H)
GH := [K`subgroup: K in Subgroups(G:IndexEqual:=2) | IsIsomorphic(K`subgroup,H)]; assert #GH eq 1; GH := GH[1];
C := ConjugacyClasses(G); s := [c[3]:c in C|CycleStructure(c[3]) eq [ <2, 6>, <1, 6> ] and #(Conjugates(G,c[3]) meet Set(GH)) eq 0]; assert #s eq 1; s := s[1];
// We can index the reps rho_M with image H by automorphism of H, and we only need to distinguish them up to charpolys
printf "Computing large automorphism group, this will take a minute or so..."; t:=Cputime();
A := Elements(AutomorphismGroup(H));
printf "computed %o automorphisms in %.3o secs.\n", #A, Cputime()-t;
_,psi := IsIsomorphic(ctx[3],GH);  _,phi := IsIsomorphic(GH,H);
rho := func<a,g|a(phi(psi(g)))>;
rhoc := func<a,g|V^-1*cc(a(phi(s*psi(g)*s)))*V> where cc := func<m|GL(3,M)![ComplexConjugate(x):x in Eltseq(m)]>;
AA := [a:a in A| &and[rho(a,g) eq rhoc(a,g):g in Generators(ctx[3])]];
S1 := Sort([s:s in {[Eltseq(CharacteristicPolynomial(rho(a,c[3]))):c in ctx[4]]:a in AA}]);
assert #S1 eq 12 and S1[1] eq cc(S1[5]) and S1[2] eq cc(S1[6]) and S1[3] eq cc(S1[4]) and S1[7] eq cc(S1[11]) and S1[8] eq cc(S1[12]) and S1[9] eq cc(S1[10]) where cc := func<s|[[ComplexConjugate(c):c in r]:r in s]>;

S1 := [[Reverse(RM!h):h in s]:s in S1];
S2 := [3,-1,0,1,0,2];

function EulerTwist(p,i)
  if p mod 6 eq 1 then return RZ!Norm(Evaluate(S1[i][FrobeniusElement(ctx,psi)],psi*RM.1)) where psi:=SexticPsi(-1,p:z6:=z6);
  else return RZ![1,0,c*p,0,c*p^2,0,p^3] where c:= S2[FrobeniusOrder(g,p) div 2]; end if;
end function;

L13t := [ Coefficients(EulerTwist(13,i))[2..4] : i in [1..12]];
assert Set(L13t) eq { [2,-10,-46],[2,14,2],[5,-10,-115],[5,35,110],[7,14,7],[7,35,154],
                      [-2,-10,46],[-2,14,-2],[-5,-10,115],[-5,35,-110],[-7,14,-7],[-7,35,-154] };

R<T> := PolynomialRing(Integers());
function EulerRes(p:minpoly:=f)
    Fp:=GF(p);
    Rp<z>:=PolynomialRing(Fp);
    nf := -Evaluate(minpoly,-x); // minpoly of -a, we use this to match Rubin-Silverberg
    fp:=Rp!nf;
    S<t>:=quo<Rp|fp>;
    if p mod 3 eq 1 then
        _,u,v := NormEquation(3,p);
        if u mod 3 eq 2 then u := -u; end if; // make u = 1 mod 3
        if v mod 3 eq 1 then v := -v; end if; // make v = 2 mod 3
        uu := 3*v-u;  vv := 3*v+u;
        h := t^((p-1) div 6);
        A := [[2*u, 1, 1], [-2*u, -1, 1], [-vv, 2*u, uu], [-uu, 2*u, vv], [vv, -2*u, uu], [uu, -2*u, vv]];
        L := &*[R| (p*T^2 - a[1]*T + 1)^Degree(GCD(fp,a[3]*Evaluate(h,z)-a[2])) : a in A];
        if Degree(L) eq 6 then return L; end if;
        if Degree(L) eq 0 then
            h := t^((p^3-1) div 6);  assert Degree(h) eq 0;
            s := u^3 - 9*u*v^2;  v := 3*u^2*v - 3*v^3; u := s;
            uu := 3*v-u;  vv := 3*v+u;
            A := [[2*u, 1, 1], [-2*u, -1, 1], [-vv, 2*u, uu], [-uu, 2*u, vv], [vv, -2*u, uu], [uu, -2*u, vv]];
            a := [a:a in A|a[3]*h eq a[2]][1];
            return p^3*T^6 - a[1]*T^3 + 1;
        else
            assert Degree(L) eq 2;
            fp := ExactQuotient(fp,GCD(fp,Evaluate(h^6,z)-1));  assert Degree(fp) eq 2;
            S<t>:=quo<Rp|fp>;
            s := u^2-3*v^2; v := 2*u*v; u := s;
            uu := 3*v-u;  vv := 3*v+u;
            A := [[2*u, 1, 1], [-2*u, -1, 1], [-vv, 2*u, uu], [-uu, 2*u, vv], [vv, -2*u, uu], [uu, -2*u, vv]];
            h := t^((p^2-1) div 6);  assert Degree(h) eq 0;
            a := [a:a in A|a[3]*h eq a[2]][1];
            return L*(p^2*T^4 - a[1]*T^2 + 1);
        end if;
    else
        gp := GCD(fp,Evaluate(t^(p-1)-1,z));
        if Degree(gp) eq 0 then return p^3*T^6+1; end if;
        if Degree(gp) eq 3 then return (p*T^2+1)^3; end if;
        fp := ExactQuotient(fp,gp);  assert Degree(fp) eq 2;
        S<t>:=quo<Rp|fp>;
        h := t^((p^2-1) div 6);
        if h eq 1 then return (p*T^2+1)^3; end if;
        if -h eq 1 then return (p*T^2+1)*(p*T^2-1)^2; end if;
        h := h^3;
        if h eq 1 then return (p*T^2+1)*(p^2*T^4-p*T^2+1); end if;
        if -h eq 1 then return (p*T^2+1)*(p^2*T^4+p*T^2+1); end if;
        assert false;
    end if;
end function;

twists := [MinimalPolynomial(c):c in [a,a^5,2^2*a,2^2*a^5,2^4*a,2^4*a^5,3^3*83^3*a,3^3*83^3*a^5,3^3*83^3*2^2*a,3^3*83^3*2^2*a^5,3^3*83^3*2^4*a,3^3*83^3*2^4*a^5]];
L13r := [ Coefficients(EulerRes(13:minpoly:=twists[i]))[2..4] : i in [1..12]];
assert Set(L13t) eq Set(L13r);
indexes := [Index(L13t,L13r[i]):i in [1..12]];

P cat:= [19,83];

for p in PrimesInInterval(1,2^10) do if not p in P then assert &and[EulerRes(p:minpoly:=twists[i]) eq EulerTwist(p,indexes[i]): i in [1..12]]; end if; end for;

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
