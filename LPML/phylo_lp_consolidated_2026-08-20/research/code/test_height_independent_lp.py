import itertools, math, time
import numpy as np
import pandas as pd
from scipy.optimize import linprog
from scipy.sparse import csr_matrix

PI=np.full(4,0.25)

def jc69(t):
    e=np.exp(-4*t/3)
    P=np.full((4,4),0.25-0.25*e)
    np.fill_diagonal(P,0.25+0.75*e)
    return P
PT=jc69(0.10)
ALPHA=float(PT.min()); BETA=float(PT.max())

def cb(n): return ALPHA**n, BETA**n
def mb(n): return ALPHA**n, .25*BETA**(n-1)
def dob(P): return max(.5*np.abs(P[i]-P[k]).sum() for i in range(4) for k in range(4))

def edge_matrix(level,e):
    # level 1 means count-2 edge, etc.; deliberately heterogeneous local lengths
    bases=[0.14,0.22,0.31,0.40,0.48]
    t=bases[min(level-1,len(bases)-1)] + 0.025*((3*e+level)%4)
    return jc69(t)

def candidate_parents(nlower):
    m=nlower//2
    return {e:(e//2,(e//2+1)%m) for e in range(nlower)}

def candidate_edges_for_parent(nlower,p):
    cp=candidate_parents(nlower)
    return [e for e in range(nlower) if p in cp[e]]

def valid_assignments(nlower):
    m=nlower//2; cp=candidate_parents(nlower); out=[]; caps=[0]*m; arr=[0]*nlower
    def rec(e):
        if e==nlower:
            if all(c==2 for c in caps): out.append(tuple(arr))
            return
        for p in cp[e]:
            if caps[p]<2:
                caps[p]+=1; arr[e]=p; rec(e+1); caps[p]-=1
    rec(0); return out

ASSIGN_CACHE={}
def assigns(n):
    if n not in ASSIGN_CACHE: ASSIGN_CACHE[n]=valid_assignments(n)
    return ASSIGN_CACHE[n]

def exact_best(pattern):
    N=len(pattern); assert N&(N-1)==0 and N>=4
    msgs=[PT[:,b] for b in pattern]
    levels=int(math.log2(N))
    best=-1.0
    def rec(cur, level):
        nonlocal best
        nlower=len(cur)
        if nlower==2:
            val=float(PI @ (cur[0]*cur[1]));
            if val>best: best=val
            return
        m=nlower//2
        Ps=[edge_matrix(level,p) for p in range(m)]
        for ass in assigns(nlower):
            kids=[[] for _ in range(m)]
            for e,p in enumerate(ass): kids[p].append(e)
            nxt=[]
            for p in range(m):
                a,b=kids[p]
                z=cur[a]*cur[b]
                nxt.append(Ps[p]@z)
            rec(nxt,level+1)
    rec(msgs,1)
    return best

class SparseLP:
    def __init__(self):
        self.names=[];self.bounds=[];self.ub=[];self.ubrhs=[];self.eq=[];self.eqrhs=[]
    def var(self,name,lb=0,ub=None):
        j=len(self.names);self.names.append(name);self.bounds.append((lb,ub));return j
    def le(self,d,rhs): self.ub.append(dict(d));self.ubrhs.append(rhs)
    def ge(self,d,rhs): self.ub.append({j:-v for j,v in d.items()});self.ubrhs.append(-rhs)
    def eql(self,d,rhs): self.eq.append(dict(d));self.eqrhs.append(rhs)
    def mat(self,rows):
        data=[];ri=[];ci=[]
        for r,d in enumerate(rows):
            for j,v in d.items():
                if abs(v)>0: ri.append(r);ci.append(j);data.append(v)
        return csr_matrix((data,(ri,ci)),shape=(len(rows),len(self.names))) if rows else None
    def solve(self,obj,maximize=False,scale=1.0):
        c=np.zeros(len(self.names))
        for j,v in obj.items(): c[j]=(-v if maximize else v)*scale
        res=linprog(c,A_ub=self.mat(self.ub),b_ub=np.array(self.ubrhs) if self.ub else None,
                    A_eq=self.mat(self.eq),b_eq=np.array(self.eqrhs) if self.eq else None,
                    bounds=self.bounds,method='highs',options={'dual_feasibility_tolerance':1e-9,'primal_feasibility_tolerance':1e-9})
        if not res.success: raise RuntimeError(res.message)
        val=sum(v*res.x[j] for j,v in obj.items())
        return val,res

# Polytope descriptor for a message after an edge with descendant count n.
def desc(n,P):
    lo,hi=cb(n); qlo,qhi=mb(n); D=dob(P)*(hi-lo)
    ratio={}
    for i in range(4):
        for k in range(i+1,4):
            rr=P[i]/P[k]; ratio[i,k]=(float(rr.min()),float(rr.max()))
    return {'lo':lo,'hi':hi,'qlo':qlo,'qhi':qhi,'D':D,'ratio':ratio}

def add_scaled_poly(lp,vec,scale,D):
    lo,hi=D['lo'],D['hi']
    for v in vec:
        lp.ge({v:1,scale:-lo},0); lp.le({v:1,scale:-hi},0)
    d={v:.25 for v in vec}; d[scale]=-D['qlo']; lp.ge(d,0)
    d={v:.25 for v in vec}; d[scale]=-D['qhi']; lp.le(d,0)
    for (i,k),(rl,rh) in D['ratio'].items():
        lp.ge({vec[i]:1,vec[k]:-rl},0); lp.le({vec[i]:1,vec[k]:-rh},0)
        lp.le({vec[i]:1,vec[k]:-1,scale:-D['D']},0)
        lp.le({vec[k]:1,vec[i]:-1,scale:-D['D']},0)

def add_poly(lp,vec,D):
    one=lp.var(f'_one_{len(lp.names)}',1,1)
    add_scaled_poly(lp,vec,one,D)

def add_route_perspective(lp,x,W,a,D):
    add_scaled_poly(lp,W,a,D)
    # remainder belongs to (1-a)P
    b=lp.var(f'remscale_{len(lp.names)}',0,1); lp.eql({b:1,a:1},1)
    R=[]
    for i in range(4):
        r=lp.var(f'rem_{len(lp.names)}_{i}',0,D['hi']);lp.eql({r:1,x[i]:-1,W[i]:1},0);R.append(r)
    add_scaled_poly(lp,R,b,D)

def add_persp_mcc(lp,t,A,B,h,Dx,Dy):
    lx,ux=Dx['lo'],Dx['hi'];ly,uy=Dy['lo'],Dy['hi']
    # t = h*x*y, A=h*x, B=h*y
    lp.ge({t:1,B:-lx,A:-ly,h:lx*ly},0)
    lp.ge({t:1,B:-ux,A:-uy,h:ux*uy},0)
    lp.le({t:1,B:-ux,A:-ly,h:ux*ly},0)
    lp.le({t:1,B:-lx,A:-uy,h:lx*uy},0)

def add_scaled_row_poly(lp,row,scale,D):
    # same as scaled poly but scale may be a signal component rather than [0,1].
    # inequalities are homogeneous after multiplying P constraints by scale.
    lo,hi=D['lo'],D['hi']
    for v in row:
        lp.ge({v:1,scale:-lo},0);lp.le({v:1,scale:-hi},0)
    d={v:.25 for v in row};d[scale]=-D['qlo'];lp.ge(d,0)
    d={v:.25 for v in row};d[scale]=-D['qhi'];lp.le(d,0)
    for (i,k),(rl,rh) in D['ratio'].items():
        lp.ge({row[i]:1,row[k]:-rl},0);lp.le({row[i]:1,row[k]:-rh},0)
        lp.le({row[i]:1,row[k]:-1,scale:-D['D']},0)
        lp.le({row[k]:1,row[i]:-1,scale:-D['D']},0)

def build_lp(pattern,strategy):
    # strategies: diag, vector, tensor, tensor_mass
    N=len(pattern); lp=SparseLP()
    msgs=[[None]*4 for _ in range(N)]
    descs=[None]*N
    # leaf messages are constants represented by fixed variables
    for e,b in enumerate(pattern):
        vals=PT[:,b]
        msgs[e]=[lp.var(f'leaf_{e}_{i}',float(vals[i]),float(vals[i])) for i in range(4)]
        # not used for leaf product relaxation
    ndesc=1; level=1
    while len(msgs)>2:
        lower=msgs; nlower=len(lower); m=nlower//2; cp=candidate_parents(nlower)
        # descriptors for lower messages: leaves fixed on first tier, otherwise stored
        lowerD=descs
        # assignment vars
        a={}
        for e in range(nlower):
            for p in cp[e]: a[e,p]=lp.var(f'a_L{level}_{e}_{p}',0,1)
            lp.eql({a[e,p]:1 for p in cp[e]},1)
        # pair marginals at each parent
        H={}; pairs_by_p={}
        for p in range(m):
            cand=candidate_edges_for_parent(nlower,p)
            pairs=[(e,f) for e,f in itertools.combinations(cand,2)]
            pairs_by_p[p]=pairs
            for e,f in pairs:H[p,e,f]=lp.var(f'H_L{level}_{p}_{e}_{f}',0,1)
            lp.eql({H[p,e,f]:1 for e,f in pairs},1)
            for e in cand:
                d={a[e,p]:-1}
                for u,v in pairs:
                    if e==u or e==v:d[H[p,u,v]]=d.get(H[p,u,v],0)+1
                lp.eql(d,0)
        # if not first tier, construct arc-routed vectors W_{e,p}; perspective only for vector+ strategies
        W={}
        if level>1:
            for e in range(nlower):
                D=lowerD[e]
                for p in cp[e]:
                    vec=[lp.var(f'W_L{level}_{e}_{p}_{i}',0,D['hi']) for i in range(4)]
                    W[e,p]=vec
                    if strategy in ('vector','tensor','tensor_mass'):
                        add_route_perspective(lp,lower[e],vec,a[e,p],D)
                    else:
                        # component McCormick a*x plus conservation
                        for i in range(4):
                            x=lower[e][i];w=vec[i];lo,hi=D['lo'],D['hi'];aa=a[e,p]
                            lp.ge({w:1,aa:-lo},0);lp.le({w:1,aa:-hi},0)
                            lp.ge({w:1,x:-1,aa:-hi},-hi);lp.le({w:1,x:-1,aa:-lo},-lo)
            for e in range(nlower):
                for i in range(4):lp.eql({W[e,p][i]:1 for p in cp[e]}|{lower[e][i]:-1},0)
        zparents=[]
        for p in range(m):
            pair_contrib=[]
            # disaggregate each child arc among local pairs
            pairA={};pairB={}
            for e,f in pairs_by_p[p]:
                h=H[p,e,f]
                if level==1:
                    # exact local product: children fixed leaf messages
                    vals=[]
                    for i in range(4):
                        coeff=float(PT[i,pattern[e]]*PT[i,pattern[f]])
                        t=lp.var(f'T_L{level}_{p}_{e}_{f}_{i}',0,coeff)
                        lp.eql({t:1,h:-coeff},0);vals.append(t)
                    pair_contrib.append(vals)
                    continue
                De,Df=lowerD[e],lowerD[f]
                Ae=[lp.var(f'A_L{level}_{p}_{e}_{f}_{i}',0,De['hi']) for i in range(4)]
                Bf=[lp.var(f'B_L{level}_{p}_{e}_{f}_{i}',0,Df['hi']) for i in range(4)]
                pairA[e,f]=Ae;pairB[e,f]=Bf
                if strategy in ('vector','tensor','tensor_mass'):
                    add_scaled_poly(lp,Ae,h,De);add_scaled_poly(lp,Bf,h,Df)
                else:
                    for i in range(4):
                        # only component scaled boxes; linkage through marginals below
                        lp.ge({Ae[i]:1,h:-De['lo']},0);lp.le({Ae[i]:1,h:-De['hi']},0)
                        lp.ge({Bf[i]:1,h:-Df['lo']},0);lp.le({Bf[i]:1,h:-Df['hi']},0)
                if strategy in ('tensor','tensor_mass'):
                    C=[[lp.var(f'C_L{level}_{p}_{e}_{f}_{i}_{j}',0,De['hi']*Df['hi']) for j in range(4)] for i in range(4)]
                    for i in range(4):
                        for j in range(4):add_persp_mcc(lp,C[i][j],Ae[i],Bf[j],h,De,Df)
                    # row = Ae_i * feasible f-vector; col = Bf_j * feasible e-vector
                    for i in range(4):add_scaled_row_poly(lp,C[i],Ae[i],Df)
                    for j in range(4):add_scaled_row_poly(lp,[C[i][j] for i in range(4)],Bf[j],De)
                    vals=[C[i][i] for i in range(4)]
                    if strategy=='tensor_mass':
                        qpair=lp.var(f'qp_L{level}_{p}_{e}_{f}',0,1)
                        lp.eql({qpair:1,**{vals[i]:-.25 for i in range(4)}},0)
                        massA=lp.var(f'mA_L{level}_{p}_{e}_{f}',0,1);lp.eql({massA:1,**{Ae[i]:-.25 for i in range(4)}},0)
                        massB=lp.var(f'mB_L{level}_{p}_{e}_{f}',0,1);lp.eql({massB:1,**{Bf[i]:-.25 for i in range(4)}},0)
                        lp.le({qpair:1,massA:-Df['hi']},0);lp.le({qpair:1,massB:-De['hi']},0)
                else:
                    vals=[]
                    for i in range(4):
                        t=lp.var(f'T_L{level}_{p}_{e}_{f}_{i}',0,De['hi']*Df['hi'])
                        add_persp_mcc(lp,t,Ae[i],Bf[i],h,De,Df);vals.append(t)
                pair_contrib.append(vals)
            # pair A/B marginals equal arc-routed W
            if level>1:
                for e in candidate_edges_for_parent(nlower,p):
                    for i in range(4):
                        d={W[e,p][i]:-1}
                        for u,v in pairs_by_p[p]:
                            if e==u:d[pairA[u,v][i]]=d.get(pairA[u,v][i],0)+1
                            elif e==v:d[pairB[u,v][i]]=d.get(pairB[u,v][i],0)+1
                        lp.eql(d,0)
            # aggregate parent partial z
            nnew=2*ndesc; lz,uz=cb(nnew)
            z=[]
            for i in range(4):
                v=lp.var(f'z_L{level}_{p}_{i}',lz,uz);d={v:-1}
                for vals in pair_contrib:d[vals[i]]=d.get(vals[i],0)+1
                lp.eql(d,0);z.append(v)
            q=lp.var(f'qz_L{level}_{p}',*mb(nnew));lp.eql({q:1,**{z[i]:-.25 for i in range(4)}},0)
            # transition over the new edge unless this will still be below root (always here since while lower>2)
            P=edge_matrix(level,p);x=[]
            for i in range(4):
                v=lp.var(f'x_L{level}_{p}_{i}',lz,uz);d={v:1}
                for j in range(4):d[z[j]]=d.get(z[j],0)-P[i,j]
                lp.eql(d,0);x.append(v)
            Dnew=desc(nnew,P);add_poly(lp,x,Dnew)
            zparents.append((x,Dnew))
        msgs=[x for x,D in zparents];descs=[D for x,D in zparents]
        ndesc*=2;level+=1
    # root product of final two edge messages. Use tensor if selected, otherwise diagonal McCormick.
    x,y=msgs;Dx,Dy=descs
    z=[]
    if strategy in ('tensor','tensor_mass'):
        # Here H=1, so A=x B=y and cross matrix directly.
        one=lp.var('root_one',1,1)
        C=[[lp.var(f'CR_{i}_{j}',0,Dx['hi']*Dy['hi']) for j in range(4)] for i in range(4)]
        for i in range(4):
            for j in range(4):add_persp_mcc(lp,C[i][j],x[i],y[j],one,Dx,Dy)
        for i in range(4):add_scaled_row_poly(lp,C[i],x[i],Dy)
        for j in range(4):add_scaled_row_poly(lp,[C[i][j] for i in range(4)],y[j],Dx)
        z=[C[i][i] for i in range(4)]
    else:
        one=lp.var('root_one',1,1)
        for i in range(4):
            t=lp.var(f'rootz_{i}',0,Dx['hi']*Dy['hi']);add_persp_mcc(lp,t,x[i],y[i],one,Dx,Dy);z.append(t)
    q=lp.var('root_like',0,1);lp.eql({q:1,**{z[i]:-.25 for i in range(4)}},0)
    # scale objective based on leaf count to keep tiny heterogeneous cases numerically visible
    val,res=lp.solve({q:1},maximize=True,scale=1e6)
    return val,len(lp.names),len(lp.ub)+len(lp.eq)

def patterns_for(N):
    if N==8:
        return {
          'grouped_AC':(0,0,0,0,1,1,1,1),
          'fourway_blocks':(0,0,1,1,2,2,3,3),
          'interleaved':(0,1,2,3,0,1,2,3),
          'skew':(0,0,0,0,1,1,2,3)}
    if N==16:
        return {
          'grouped_AC':tuple([0]*8+[1]*8),
          'fourway_blocks':tuple([0]*4+[1]*4+[2]*4+[3]*4),
          'interleaved':tuple([0,1,2,3]*4),
          'skew':tuple([0]*8+[1]*4+[2]*2+[3]*2)}
    raise ValueError

if __name__=='__main__':
    rows=[]
    for N in (8,16):
        for pname,pat in patterns_for(N).items():
            t=time.time();ex=exact_best(pat);tex=time.time()-t
            for strat in ('diag','vector','tensor','tensor_mass'):
                t=time.time();v,nv,nc=build_lp(pat,strat);dt=time.time()-t
                rows.append({'N':N,'pattern':pname,'strategy':strat,'exact':ex,'upper':v,'ratio':v/ex,'vars':nv,'constraints':nc,'seconds':dt,'exact_seconds':tex})
                print(rows[-1],flush=True)
    df=pd.DataFrame(rows);df.to_csv('/mnt/data/height_independent_results.csv',index=False)
