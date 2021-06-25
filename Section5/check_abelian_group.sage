#
# This SageMath code performs a consistency check described in Remark 5.42,
# confirming that Table 4 and the lists in Proposition 5.23 and Proposition 5.25 are complete.
#

def f(u,v,w):
   return u/v + u/w + v/u + v/w + w/u + w/v

def reorder(i,j,k):
   l = [i,j,k]
   l.sort()
   return tuple([lcm([i.denom(), j.denom(), k.denom()])] + l)

def check_subgroup(ord,i,j,k,z,n):
   return all(f(z^(m*n*i),z^(m*n*j),z^(m*n*k)).is_integer() for m in range(ord))

## Find all groups of the form A(1, n\phi).

def find_gens(n):
   K.<z> = CyclotomicField(n)
   l = [reorder(i/n, j/n, ((-i-j)%n)/n) for i in range(n) for j in range(n)]
   l = list(set(l))
   l2 = []
   for (ord,i,j,k) in l:
      ans = (ord,i,j,k)
      for m in range(ord):
         if gcd(m, ord) == 1:
            ans2 = reorder((i*n*m%n)/n, (j*n*m%n)/n, (k*n*m%n)/n)
            if ans2 < ans:
               ans = ans2
      if not ans in l2:
         l2.append(ans)
   l2.sort()
   l3 = []
   for (ord,i,j,k) in l2:
      if ord%3 != 0: 
         continue
      if (i*ord/3).is_integer() or (j*ord/3).is_integer() or (k*ord/3).is_integer(): 
         continue
      if check_subgroup(ord,i,j,k,z,n):
         l3.append((ord,i,j,k))
   return l3

assert all(not (ord%7==0 and ord>7*3) for (ord,i,j,k) in find_gens(3*49))
assert all(not (ord%7==0 and ord>7*3) for (ord,i,j,k) in find_gens(3*42))
l = find_gens(3*48); print(l)

## Identify groups of the form A(m, n phi) for m in {2,3,4,6}, as well as those with action by \alt{3}.
n = 72
K.<z> = CyclotomicField(n)
for m in [2,3,4,6]:
   l2 = [(ord,i,j,k) for (ord,i,j,k) in l if ord%m==0 and 
         all(check_subgroup(ord,i+i1/m,j+j1/m,k+(-i1-j1)/m,z,72) for i1 in range(m) for j1 in range(m))]
   print(m, l2)
   for ord in divisors(36):
      l3 = [(i,j,k) for (ord1,i,j,k) in l2 if ord==ord1]
      if len(l3) == 0: continue
      s = []
      for (i,j,k) in l3: 
         s1 = [((t*i+i1/m)*n%n/n, (t*j+j1/m)*n%n/n, (t*k-i1/m-j1/m)*n%n/n) for i1 in range(m) for j1 in range(m) for t in range(ord)]
         s1.sort()
         if all((j1,k1,i1) in s1 for (i1,j1,k1) in s1):
            print("alt3: ", m, (ord,i,j,k))
         s.append(tuple(s1))
      print(m,ord,len(s), len(set(s)))



