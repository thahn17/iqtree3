import itertools
import numpy as np
from scipy.optimize import linprog

pi=np.full(4,0.25)

def jc69(t):
    e=np.exp(-4*t/3)
    P=np.full((4,4),0.25-0.25*e)
    np.fill_diagonal(P,0.25+0.75*e)
    return P

Pt=jc69(0.10)
P2=[jc69(t) for t in (0.12,0.20,0.32,0.45)]
P4=[jc69(t) for t in (0.18,0.38)]
alpha,beta=float(Pt.min()),float(Pt.max())

def cb(n): return alpha**n,beta**n
def mb(n): return alpha**n,0.25*beta**(n-1)
def dob(P): return max(0.5*np.abs(P[i]-P[k]).sum() for i in range(4) for k in range(4))

class LP:
    def __init__(self): self.names=[];self.bounds=[];self.A=[];self.b=[];self.E=[];self.f=[]
    def var(self,name,lb=0,ub=None):
        j=len(self.names);self.names.append(name);self.bounds.append((lb,ub));return j
    def _row(self,d):
        r=np.zeros(len(self.names));
        for j,v in d.items():r[j]=v
        return r
    def le(self,d,rhs):self.A.append(self._row(d));self.b.append(rhs)
    def ge(self,d,rhs):self.le({j:-v for j,v in d.items()},-rhs)
    def eq(self,d,rhs):self.E.append(self._row(d));self.f.append(rhs)
    def _pad(self,R):
        if not R:return None
        n=len(self.names);return np.array([np.pad(r,(0,n-len(r))) for r in R])
    def solve(self,obj,maximize=False):
        c=np.zeros(len(self.names))
        for j,v in obj.items():c[j]= -v if maximize else v
        res=linprog(c,A_ub=self._pad(self.A),b_ub=np.array(self.b) if self.b else None,
                    A_eq=self._pad(self.E),b_eq=np.array(self.f) if self.f else None,
                    bounds=self.bounds,method='highs')
        if not res.success: raise RuntimeError(res.message)
        return (-res.fun if maximize else res.fun),res

def mcc(lp,z,x,y,lx,ux,ly,uy):
    lp.ge({z:1,x:-ly,y:-lx},-lx*ly)
    lp.ge({z:1,x:-uy,y:-ux},-ux*uy)
    lp.le({z:1,x:-ly,y:-ux},-ux*ly)
    lp.le({z:1,x:-uy,y:-lx},-lx*uy)

def routeprod(lp,w,a,x,lx,ux):
    lp.ge({w:1},0)
    lp.ge({w:1,x:-1,a:-ux},-ux)
    lp.le({w:1,a:-ux},0)
    lp.le({w:1,x:-1,a:-lx},-lx)

def mass(lp,vec,name):
    q=lp.var(name,0,1);d={q:1}
    for i in range(4):d[vec[i]]=d.get(vec[i],0)-0.25
    lp.eq(d,0);return q

def obbt(lp,vec):
    out=[]
    for v in vec:
        lo,_=lp.solve({v:1});hi,_=lp.solve({v:1},maximize=True);out.append((lo,hi))
    return out

def osc_bound(lp,vec):
    D=0.0
    for i in range(4):
        for k in range(i+1,4):
            lo,_=lp.solve({vec[i]:1,vec[k]:-1})
            hi,_=lp.solve({vec[i]:1,vec[k]:-1},maximize=True)
            D=max(D,abs(lo),abs(hi))
    return D

leaf_parents={0:(0,1),1:(0,1),2:(1,2),3:(1,2),4:(2,3),5:(2,3),6:(3,0),7:(3,0)}

def exact_best(pat):
    lm=[Pt[:,b] for b in pat];best=-1
    for bits in itertools.product((0,1),repeat=8):
        par=[leaf_parents[e][bits[e]] for e in range(8)]
        if [par.count(p) for p in range(4)] != [2]*4:continue
        z2=[]
        for p in range(4):
            kids=[e for e in range(8) if par[e]==p]
            z2.append(lm[kids[0]]*lm[kids[1]])
        x2=[P2[p]@z2[p] for p in range(4)]
        for up in itertools.product((0,1),repeat=4):
            if sum(up)!=2:continue
            z4=[]
            for q in range(2):
                kids=[p for p in range(4) if up[p]==q]
                z4.append(x2[kids[0]]*x2[kids[1]])
            x4=[P4[q]@z4[q] for q in range(2)]
            best=max(best,0.25*np.dot(x4[0],x4[1]))
    return best

def build(pat,mode):
    use_count=mode!='plain';use_edge=mode in ('edge','joint');use_joint=mode=='joint'
    lm=np.array([Pt[:,b] for b in pat]);lp=LP()
    slots2=[(p,s) for p in range(4) for s in range(2)]
    y={}
    for e in range(8):
        ks=[k for k,(p,s) in enumerate(slots2) if p in leaf_parents[e]]
        for k in ks:y[e,k]=lp.var(f'y1_{e}_{k}',0,1)
        lp.eq({y[e,k]:1 for k in ks},1)
    for k in range(8):lp.eq({y[e,k]:1 for e in range(8) if (e,k) in y},1)
    s2=[]
    for k,(p,s) in enumerate(slots2):
        vec=[]
        for i in range(4):
            es=[e for e in range(8) if (e,k) in y];vals=[lm[e,i] for e in es]
            v=lp.var(f's2_{k}_{i}',min(vals),max(vals));d={v:1}
            for e in es:d[y[e,k]]=d.get(y[e,k],0)-lm[e,i]
            lp.eq(d,0);vec.append(v)
        s2.append(vec)
    z2=[];q2=[]
    for p in range(4):
        A,B=s2[2*p],s2[2*p+1];z=[]
        for i in range(4):
            lx,ux=lp.bounds[A[i]];ly,uy=lp.bounds[B[i]];cl,cu=cb(2)
            v=lp.var(f'z2_{p}_{i}',max(lx*ly,cl) if use_count else lx*ly,min(ux*uy,cu) if use_count else ux*uy)
            mcc(lp,v,A[i],B[i],lx,ux,ly,uy);z.append(v)
        z2.append(z);q=mass(lp,z,f'q2_{p}');q2.append(q)
        if use_count:
            l,u=mb(2);lp.bounds[q]=(l,u);lp.ge({q:1},l);lp.le({q:1},u)
        if use_edge:
            _,cu=cb(2);_,qu=mb(2)
            for e in range(8):
                if p not in leaf_parents[e]:continue
                ks=[k for k,(pp,s) in enumerate(slots2) if pp==p and (e,k) in y]
                # Stronger: parent <= child message times the best possible sibling envelope.
                sibling_candidates=[f for f in range(8) if f!=e and p in leaf_parents[f]]
                Usib=np.array([max(lm[f,i] for f in sibling_candidates) for i in range(4)])
                for i in range(4):
                    target=Usib[i]*lm[e,i]
                    M=max(0,cu-target);d={z[i]:1}
                    for k in ks:d[y[e,k]]=d.get(y[e,k],0)+M
                    lp.le(d,target+M)
                targetq=0.25*float(np.dot(Usib,lm[e]))
                M=max(0,qu-targetq);d={q:1}
                for k in ks:d[y[e,k]]=d.get(y[e,k],0)+M
                lp.le(d,targetq+M)
    x2=[];D2=[]
    for p in range(4):
        vec=[]
        for i in range(4):
            v=lp.var(f'x2_{p}_{i}',0,1);d={v:1}
            for j in range(4):d[z2[p][j]]=d.get(z2[p][j],0)-P2[p][i,j]
            lp.eq(d,0);vec.append(v)
        x2.append(vec);qx=mass(lp,vec,f'qx2_{p}');lp.eq({qx:1,q2[p]:-1},0)
        if use_count:
            for i in range(4):
                for k in range(i+1,4):
                    rr=P2[p][i]/P2[p][k];lp.ge({vec[i]:1,vec[k]:-rr.min()},0);lp.le({vec[i]:1,vec[k]:-rr.max()},0)
        l,u=cb(2);D2.append(dob(P2[p])*(u-l))
    bx2=[]
    for p in range(4):
        b=obbt(lp,x2[p]) if use_count else [(0,1)]*4;bx2.append(b)
        if use_count:
            for i in range(4):lp.bounds[x2[p][i]]=b[i]
    slots4=[(q,s) for q in range(2) for s in range(2)];y2={}
    for p in range(4):
        for k in range(4):y2[p,k]=lp.var(f'y2_{p}_{k}',0,1)
        lp.eq({y2[p,k]:1 for k in range(4)},1)
    for k in range(4):lp.eq({y2[p,k]:1 for p in range(4)},1)
    s4=[];qs4=[];Ds=[];W={}
    for k,(q,s) in enumerate(slots4):
        vec=[]
        for i in range(4):
            lb=min(bx2[p][i][0] for p in range(4));ub=max(bx2[p][i][1] for p in range(4));v=lp.var(f's4_{k}_{i}',lb,ub);d={v:1}
            for p in range(4):
                w=lp.var(f'w_{p}_{k}_{i}',0,ub);W[p,k,i]=w;routeprod(lp,w,y2[p,k],x2[p][i],bx2[p][i][0],bx2[p][i][1]);d[w]=-1
            lp.eq(d,0);vec.append(v)
        s4.append(vec);qs4.append(mass(lp,vec,f'qs4_{k}'));Ds.append(max(D2))
    for p in range(4):
        for i in range(4):lp.eq({**{W[p,k,i]:1 for k in range(4)},x2[p][i]:-1},0)
    if use_joint:
        for k in range(4):
            b=obbt(lp,s4[k])
            for i in range(4): lp.bounds[s4[k][i]]=b[i]
            qlo,_=lp.solve({qs4[k]:1}); qhi,_=lp.solve({qs4[k]:1},maximize=True)
            lp.bounds[qs4[k]]=(qlo,qhi)
            Ds[k]=osc_bound(lp,s4[k])
    z4=[];q4=[]
    for q in range(2):
        A,B=s4[2*q],s4[2*q+1];qa,qb=qs4[2*q],qs4[2*q+1];z=[]
        for i in range(4):
            lx,ux=lp.bounds[A[i]];ly,uy=lp.bounds[B[i]];cl,cu=cb(4);v=lp.var(f'z4_{q}_{i}',max(lx*ly,cl) if use_count else lx*ly,min(ux*uy,cu) if use_count else ux*uy);mcc(lp,v,A[i],B[i],lx,ux,ly,uy);z.append(v)
        z4.append(z);qp=mass(lp,z,f'q4_{q}');q4.append(qp)
        if use_count:
            l,u=mb(4);lp.bounds[qp]=(l,u);lp.ge({qp:1},l);lp.le({qp:1},u)
        if use_edge:
            lp.le({qp:1,qa:-1},0);lp.le({qp:1,qb:-1},0)
            UA=max(lp.bounds[v][1] for v in A);UB=max(lp.bounds[v][1] for v in B)
            lp.le({qp:1,qa:-UB},0);lp.le({qp:1,qb:-UA},0)
            _,qu=mb(4);cl,cu=cb(4)
            for p in range(4):
                ks=[2*q,2*q+1]
                # If p is attached here, bound using the best possible *other* count-2 edge.
                others=[r for r in range(4) if r!=p]
                Usib=np.array([max(bx2[r][i][1] for r in others) for i in range(4)])
                rhs={x2[p][i]:0.25*Usib[i] for i in range(4)}
                rhs_lb=sum(0.25*Usib[i]*bx2[p][i][0] for i in range(4))
                M=max(0,qu-rhs_lb);d={qp:1}
                for i,c in rhs.items():d[i]=d.get(i,0)-c
                for k in ks:d[y2[p,k]]=d.get(y2[p,k],0)+M
                lp.le(d,M)
                for i in range(4):
                    M=max(0,cu-Usib[i]*bx2[p][i][0]);d={z[i]:1,x2[p][i]:-Usib[i]}
                    for k in ks:d[y2[p,k]]=d.get(y2[p,k],0)+M
                    lp.le(d,M)
        if use_joint:
            lqa,uqa=lp.bounds[qa];lqb,uqb=lp.bounds[qb];r=lp.var(f'r4_{q}',lqa*lqb,uqa*uqb);mcc(lp,r,qa,qb,lqa,uqa,lqb,uqb);lp.le({qp:1,r:-1},Ds[2*q]*Ds[2*q+1]/4)
    x4=[];D4=[]
    for q in range(2):
        vec=[]
        for i in range(4):
            v=lp.var(f'x4_{q}_{i}',0,1);d={v:1}
            for j in range(4):d[z4[q][j]]=d.get(z4[q][j],0)-P4[q][i,j]
            lp.eq(d,0);vec.append(v)
        x4.append(vec);qx=mass(lp,vec,f'qx4_{q}');lp.eq({qx:1,q4[q]:-1},0)
        if use_count:
            for i in range(4):
                for k in range(i+1,4):
                    rr=P4[q][i]/P4[q][k];lp.ge({vec[i]:1,vec[k]:-rr.min()},0);lp.le({vec[i]:1,vec[k]:-rr.max()},0)
        l,u=cb(4);D4.append(dob(P4[q])*(u-l))
    bx4=[]
    for q in range(2):
        b=obbt(lp,x4[q]) if use_count else [(0,1)]*4;bx4.append(b)
        if use_count:
            for i in range(4):lp.bounds[x4[q][i]]=b[i]
        if use_joint:
            D4[q]=osc_bound(lp,x4[q])
            qlo,_=lp.solve({q4[q]:1});qhi,_=lp.solve({q4[q]:1},maximize=True)
            lp.bounds[q4[q]]=(qlo,qhi)
    z8=[]
    for i in range(4):
        lx,ux=lp.bounds[x4[0][i]];ly,uy=lp.bounds[x4[1][i]];cl,cu=cb(8);v=lp.var(f'z8_{i}',max(lx*ly,cl) if use_count else lx*ly,min(ux*uy,cu) if use_count else ux*uy);mcc(lp,v,x4[0][i],x4[1][i],lx,ux,ly,uy);z8.append(v)
    q8=mass(lp,z8,'q8')
    if use_count:
        l,u=mb(8);lp.bounds[q8]=(l,u);lp.ge({q8:1},l);lp.le({q8:1},u)
    if use_edge:
        lp.le({q8:1,q4[0]:-1},0);lp.le({q8:1,q4[1]:-1},0)
        U0=max(lp.bounds[v][1] for v in x4[0]);U1=max(lp.bounds[v][1] for v in x4[1]);lp.le({q8:1,q4[0]:-U1},0);lp.le({q8:1,q4[1]:-U0},0)
    if use_joint:
        l0,u0=lp.bounds[q4[0]];l1,u1=lp.bounds[q4[1]];r=lp.var('r8',l0*l1,u0*u1);mcc(lp,r,q4[0],q4[1],l0,u0,l1,u1);lp.le({q8:1,r:-1},D4[0]*D4[1]/4)
    return lp.solve({q8:1},maximize=True)[0]

patterns={'AAAAAAAA':(0,)*8,'AAAACCCC':(0,0,0,0,1,1,1,1),'AACCGGTT':(0,0,1,1,2,2,3,3),'AAAACCGT':(0,0,0,0,1,1,2,3)}
print('pattern,exact,plain,count,edge,joint,joint/exact,gap_removed')
for name,pat in patterns.items():
    ex=exact_best(pat);vals={m:build(pat,m) for m in ('plain','count','edge','joint')};gr=(vals['plain']-vals['joint'])/(vals['plain']-ex) if vals['plain']>ex else float('nan')
    print(name,ex,vals['plain'],vals['count'],vals['edge'],vals['joint'],vals['joint']/ex,gr,sep=',')
print('alpha,beta',alpha,beta)
for n in (2,4,8):print('n',n,'cb',cb(n),'mb',mb(n))
print('D2',[dob(P) for P in P2]);print('D4',[dob(P) for P in P4])
