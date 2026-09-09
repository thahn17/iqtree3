import itertools,math,importlib.util,sys
import numpy as np
from scipy.linalg import expm
spec=importlib.util.spec_from_file_location('la','/mnt/data/multisite_lowerarm_taxonflow.py');la=importlib.util.module_from_spec(spec);sys.modules['la']=la;spec.loader.exec_module(la)
PI,Q,BASE_INDEX=la.PI,la.Q,la.BASE_INDEX

# Canonical rooted binary trees as nested tuples of leaf ints.
def key(t): return str(t)
def trees(leaves):
    leaves=tuple(sorted(leaves))
    if len(leaves)==1:return [leaves[0]]
    first=leaves[0];rest=leaves[1:];out=[]
    # choose left subset containing first, nonempty proper; canonical by string key
    for mask in range(1<<(len(rest))):
        A=(first,)+tuple(rest[i] for i in range(len(rest)) if mask>>i&1)
        if len(A)==len(leaves):continue
        B=tuple(x for x in leaves if x not in A)
        for x in trees(A):
            for y in trees(B):
                if key(x)<=key(y):out.append((x,y))
                else:out.append((y,x))
    # dedup
    d={key(t):t for t in out};return list(d.values())

def count(t): return 1 if isinstance(t,int) else count(t[0])+count(t[1])
def edges(t):
    if isinstance(t,int):return []
    out=[(id(t),0,count(t[0])),(id(t),1,count(t[1]))]
    return out+edges(t[0])+edges(t[1])

def eval_site(t,pat,choice,emap):
    def rec(node):
        if isinstance(node,int):
            v=np.zeros(4);v[BASE_INDEX[pat[node]]]=1;return v
        a,b=node;za=rec(a);zb=rec(b)
        ma=count(a);mb=count(b)
        ta=la.branch_times(ma)*(2.0 if choice[emap[(id(node),0,ma)]] else 1.0)
        tb=la.branch_times(mb)*(2.0 if choice[emap[(id(node),1,mb)]] else 1.0)
        return (expm(Q*ta)@za)*(expm(Q*tb)@zb)
    return float(PI@rec(t))

cols=[tuple('AACC'),tuple('ACGT'),tuple('AGCT')]
best=-1e300;bestinfo=None;ts=trees(range(4));print('trees',len(ts))
for t in ts:
    ee=edges(t);emap={e:k for k,e in enumerate(ee)}
    for bits in itertools.product([0,1],repeat=len(ee)):
        val=sum(math.log(eval_site(t,p,bits,emap)) for p in cols)
        if val>best:best=val;bestinfo=(t,bits)
print('exact',best,bestinfo)
M=la.BoundaryTaxonModel('/mnt/data/lowerarm_units4/library.json',cols,1.05,2.0);v,r,dt=M.solve();print('LP',v,'gap',v-best,'vars',len(M.lp.names),'eq',len(M.lp.eq),'le',len(M.lp.le),'solve',dt,'root',M.root_taxon_sums(r))
