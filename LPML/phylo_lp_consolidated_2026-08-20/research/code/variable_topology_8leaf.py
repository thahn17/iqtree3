import itertools, numpy as np
from scipy.optimize import linprog
pi=np.full(4,.25)
def jc(t):
 e=np.exp(-4*t/3);P=np.full((4,4),.25-.25*e);np.fill_diagonal(P,.25+.75*e);return P
Pt=jc(.10); P2=[jc(t) for t in (.12,.18,.24,.30)]; P4=[jc(.20),jc(.40)]
alpha,beta=Pt.min(),Pt.max()
def cb(n):return alpha**n,beta**n
def qb(n):return alpha**n,pi.max()*beta**(n-1)
class LP:
 def __init__(s):s.names=[];s.bounds=[];s.Au=[];s.bu=[];s.Ae=[];s.be=[]
 def var(s,n,lo=None,hi=None):s.names.append(n);s.bounds.append((lo,hi));return len(s.names)-1
 def row(s,c):
  r=np.zeros(len(s.names));
  for i,v in c.items():r[i]=v
  return r
 def le(s,c,b):s.Au.append(s.row(c));s.bu.append(b)
 def ge(s,c,b):s.le({i:-v for i,v in c.items()},-b)
 def eq(s,c,b):s.Ae.append(s.row(c));s.be.append(b)
 def mat(s,R):
  if not R:return None
  n=len(s.names);return np.array([np.pad(r,(0,n-len(r))) for r in R])
 def solve(s,obj,maxi=False):
  c=np.zeros(len(s.names));
  for i,v in obj.items():c[i]=-v if maxi else v
  a=linprog(c,A_ub=s.mat(s.Au),b_ub=np.array(s.bu) if s.bu else None,A_eq=s.mat(s.Ae),b_eq=np.array(s.be) if s.be else None,bounds=s.bounds,method='highs')
  if not a.success:raise RuntimeError(a.message)
  return (-a.fun if maxi else a.fun),a

def mcc(lp,z,x,y,lx,ux,ly,uy):
 lp.ge({z:1,x:-ly,y:-lx},-lx*ly);lp.ge({z:1,x:-uy,y:-ux},-ux*uy)
 lp.le({z:1,x:-ly,y:-ux},-ux*ly);lp.le({z:1,x:-uy,y:-lx},-lx*uy)
def routeprod(lp,w,y,x,l,u):
 # w = y*x, 0<=y<=1, l<=x<=u
 lp.ge({w:1,y:-l},0);lp.le({w:1,y:-u},0)
 lp.ge({w:1,x:-1,y:-u},-u);lp.le({w:1,x:-1,y:-l},-l)
def tighten_bounds(lp,vs):
 out=[]
 for v in vs:
  lo,_=lp.solve({v:1});hi,_=lp.solve({v:1},True);out.append((lo,hi))
 return out
# leaf e candidates among 4 count2 parents
cand=[[0,1],[0,1],[1,2],[1,2],[2,3],[2,3],[3,0],[3,0]]
def exact(pat):
 lm=[Pt[:,b] for b in pat]; best=-1;nt=0
 for bits in itertools.product((0,1),repeat=8):
  pars=[cand[e][bits[e]] for e in range(8)]
  if any(pars.count(p)!=2 for p in range(4)):continue
  nt+=1; m2=[]
  for p in range(4):
   ch=[e for e in range(8) if pars[e]==p]
   z=lm[ch[0]]*lm[ch[1]];m2.append(P2[p]@z)
  for left in itertools.combinations(range(4),2):
   L=set(left);R=[p for p in range(4) if p not in L]
   z0=m2[list(L)[0]]*m2[list(L)[1]];z1=m2[R[0]]*m2[R[1]]
   val=pi@((P4[0]@z0)*(P4[1]@z1))
   if val>best:best=val
 return best,nt

def relax(pat,count=False,obbt=False):
 lm=np.array([Pt[:,b] for b in pat]);lp=LP()
 # tier 1: 8 count1 edges -> 4 count2 parents, 2 slots each
 slots=[(p,s) for p in range(4) for s in range(2)]
 y={}
 for e in range(8):
  for si,(p,s) in enumerate(slots):
   if p in cand[e]:y[e,si]=lp.var(f'y1_{e}_{si}',0,1)
 for e in range(8):lp.eq({v:1 for (ee,_),v in y.items() if ee==e},1)
 for si in range(8):lp.eq({v:1 for (e,ss),v in y.items() if ss==si},1)
 qslot=[[None]*4 for _ in range(8)];qbnd=[[None]*4 for _ in range(8)]
 for si,(p,s) in enumerate(slots):
  elig=[e for e in range(8) if (e,si) in y]
  for i in range(4):
   vals=lm[elig,i];lo,hi=float(vals.min()),float(vals.max());qbnd[si][i]=(lo,hi)
   qslot[si][i]=lp.var(f'q1_{si}_{i}',lo,hi);c={qslot[si][i]:1}
   for e in elig:c[y[e,si]]=-lm[e,i]
   lp.eq(c,0)
 z2=[];m2=[]
 for p in range(4):
  a,b=2*p,2*p+1;z=[]
  for i in range(4):
   lx,ux=qbnd[a][i];ly,uy=qbnd[b][i];lo,hi=lx*ly,ux*uy
   if count:
    cl,ch=cb(2);lo=max(lo,cl);hi=min(hi,ch)
   v=lp.var(f'z2_{p}_{i}',lo,hi);mcc(lp,v,qslot[a][i],qslot[b][i],lx,ux,ly,uy);z.append(v)
  if count:
   ql,qh=qb(2);lp.ge({z[i]:pi[i] for i in range(4)},ql);lp.le({z[i]:pi[i] for i in range(4)},qh)
  z2.append(z)
  mm=[]
  for i in range(4):
   lo,hi=cb(2) if count else (0,1);v=lp.var(f'm2_{p}_{i}',lo,hi);c={v:1}
   for j in range(4):c[z[j]]=-P2[p][i,j]
   lp.eq(c,0);mm.append(v)
  m2.append(mm)
 # optionally tighten count2 message boxes before variable routing
 m2b=[]
 for p in range(4):m2b.append(tighten_bounds(lp,m2[p]) if obbt else [lp.bounds[v] for v in m2[p]])
 for p in range(4):
  for i,v in enumerate(m2[p]):lp.bounds[v]=m2b[p][i]
 # tier 2: four count2 edges -> two count4 parents, every child has 2 possible parents
 slots2=[(q,s) for q in range(2) for s in range(2)]
 y2={(p,si):lp.var(f'y2_{p}_{si}',0,1) for p in range(4) for si in range(4)}
 for p in range(4):lp.eq({y2[p,si]:1 for si in range(4)},1)
 for si in range(4):lp.eq({y2[p,si]:1 for p in range(4)},1)
 # routed flows w=y*m2 and slot messages q2
 q2=[[lp.var(f'q2_{si}_{i}',*(cb(2) if count else (0,1))) for i in range(4)] for si in range(4)]
 for p in range(4):
  for i in range(4):
   ws=[];l,u=m2b[p][i]
   for si in range(4):
    w=lp.var(f'w_{p}_{si}_{i}',min(0,l),max(0,u));routeprod(lp,w,y2[p,si],m2[p][i],l,u);ws.append(w)
   # RLT conservation: exactly one upward attachment
   c={w:1 for w in ws};c[m2[p][i]]=-1;lp.eq(c,0)
 for si in range(4):
  for i in range(4):
   c={q2[si][i]:1}
   for p in range(4):
    # find w by name index (small model; dictionary would be cleaner)
    idx=lp.names.index(f'w_{p}_{si}_{i}');c[idx]=-1
   lp.eq(c,0)
 z4=[];m4=[]
 for q in range(2):
  a,b=2*q,2*q+1;z=[]
  for i in range(4):
   lx,ux=lp.bounds[q2[a][i]];ly,uy=lp.bounds[q2[b][i]];lo,hi=lx*ly,ux*uy
   if count:
    cl,ch=cb(4);lo=max(lo,cl);hi=min(hi,ch)
   v=lp.var(f'z4_{q}_{i}',lo,hi);mcc(lp,v,q2[a][i],q2[b][i],lx,ux,ly,uy);z.append(v)
  if count:
   ql,qh=qb(4);lp.ge({z[i]:pi[i] for i in range(4)},ql);lp.le({z[i]:pi[i] for i in range(4)},qh)
  z4.append(z);mm=[]
  for i in range(4):
   lo,hi=cb(4) if count else (0,1);v=lp.var(f'm4_{q}_{i}',lo,hi);c={v:1}
   for j in range(4):c[z[j]]=-P4[q][i,j]
   lp.eq(c,0);mm.append(v)
  m4.append(mm)
 m4b=[tighten_bounds(lp,m4[q]) if obbt else [lp.bounds[v] for v in m4[q]] for q in range(2)]
 for q in range(2):
  for i,v in enumerate(m4[q]):lp.bounds[v]=m4b[q][i]
 # root count8
 r=[]
 for i in range(4):
  lx,ux=m4b[0][i];ly,uy=m4b[1][i];lo,hi=lx*ly,ux*uy
  if count:
   cl,ch=cb(8);lo=max(lo,cl);hi=min(hi,ch)
  v=lp.var(f'r_{i}',lo,hi);mcc(lp,v,m4[0][i],m4[1][i],lx,ux,ly,uy);r.append(v)
 if count:
  ql,qh=qb(8);lp.ge({r[i]:pi[i] for i in range(4)},ql);lp.le({r[i]:pi[i] for i in range(4)},qh)
 val,a=lp.solve({r[i]:pi[i] for i in range(4)},True)
 # routing fractionality diagnostic for y2 (tier 2)
 frac=max(min(a.x[v],1-a.x[v]) for v in y2.values())
 return val,frac

patterns={'AAAAAAAA':(0,)*8,'AAAACCCC':(0,0,0,0,1,1,1,1),'AACCGGTT':(0,0,1,1,2,2,3,3),'AAAACCGT':(0,0,0,0,1,1,2,3)}
print('alpha,beta',alpha,beta)
for n in (2,4,8):print('n',n,'component',cb(n),'mass',qb(n))
print('pattern,exact,valid_leaf_assignments,plain,count,count_obbt,plain_ratio,count_ratio,count_obbt_ratio,frac')
for name,p in patterns.items():
 ex,nt=exact(p);pl,_=relax(p,False,False);co,_=relax(p,True,False);cobbt,fr=relax(p,True,True)
 print(name,ex,nt,pl,co,cobbt,pl/ex,co/ex,cobbt/ex,fr,sep=',')
