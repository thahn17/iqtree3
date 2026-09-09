import importlib.util, itertools, time, math
import numpy as np, pandas as pd
spec=importlib.util.spec_from_file_location('g','/mnt/data/test_height_independent_lp.py')
g=importlib.util.module_from_spec(spec);spec.loader.exec_module(g)

# Fixed 12-direction family: 4 coordinates, total mass, and 7 nontrivial full-sign contrasts.
DIRS=[]
for i in range(4):
    d=[0]*4;d[i]=1;DIRS.append(tuple(d))
DIRS.append((1,1,1,1))
for tail in itertools.product((-1,1),repeat=3):
    d=(1,)+tail
    if d!=(1,1,1,1):DIRS.append(d)
DIRS=list(dict.fromkeys(DIRS))
print('direction_count',len(DIRS))

def generic_desc(n,P,support=None):
    D=g.desc(n,P)
    D['support']=support or {}
    if support:
        # tighten coordinate and mass bounds from support
        los=[];his=[]
        for i in range(4):
            key=tuple(1 if k==i else 0 for k in range(4));lo,hi=support[key];los.append(lo);his.append(hi)
        # descriptor currently assumes common lo/hi. Use safe common envelope.
        D['lo']=min(los);D['hi']=max(his)
        loT,hiT=support[(1,1,1,1)];D['qlo']=loT/4;D['qhi']=hiT/4
        # tighten common oscillation
        DD=0; have=False
        for i,j in itertools.combinations(range(4),2):
            key=tuple(1 if k==i else -1 if k==j else 0 for k in range(4))
            if key in support:
                lo,hi=support[key]; DD=max(DD,abs(lo),abs(hi)); have=True
        if have: D['D']=min(D['D'],DD)
    return D

def add_scaled_poly_support(lp,vec,scale,D):
    g.add_scaled_poly(lp,vec,scale,D)
    for d,(lo,hi) in D.get('support',{}).items():
        row={vec[i]:d[i] for i in range(4) if d[i]}
        row[scale]=row.get(scale,0)-lo;lp.ge(row,0)
        row={vec[i]:d[i] for i in range(4) if d[i]}
        row[scale]=row.get(scale,0)-hi;lp.le(row,0)

def local_pair_support(De,Df,P):
    # Outer hull for one possible child pair, using fixed tensor constraints.
    lp=g.SparseLP();one=lp.var('one',1,1)
    x=[lp.var(f'x{i}',De['lo'],De['hi']) for i in range(4)]
    y=[lp.var(f'y{i}',Df['lo'],Df['hi']) for i in range(4)]
    add_scaled_poly_support(lp,x,one,De);add_scaled_poly_support(lp,y,one,Df)
    C=[[lp.var(f'C{i}{j}',0,De['hi']*Df['hi']) for j in range(4)] for i in range(4)]
    for i in range(4):
        for j in range(4):g.add_persp_mcc(lp,C[i][j],x[i],y[j],one,De,Df)
    for i in range(4):add_scaled_poly_support(lp,C[i],x[i],Df)
    for j in range(4):add_scaled_poly_support(lp,[C[i][j] for i in range(4)],y[j],De)
    z=[C[i][i] for i in range(4)]
    out=[]
    for i in range(4):
        v=lp.var(f'o{i}',0,1);dd={v:1}
        for j in range(4):dd[z[j]]=dd.get(z[j],0)-P[i,j]
        lp.eql(dd,0);out.append(v)
    supp={}
    for d in DIRS:
        obj={out[i]:d[i] for i in range(4) if d[i]}
        lo,_=lp.solve(obj,maximize=False,scale=1e4);hi,_=lp.solve(obj,maximize=True,scale=1e4)
        supp[d]=(lo,hi)
    return supp

def build_descriptors(pattern):
    N=len(pattern);lowerD=None;ndesc=1;level=1; allD=[]
    while N//(2**(level-1))>2:
        nlower=N//(2**(level-1));m=nlower//2;newD=[]
        for p in range(m):
            P=g.edge_matrix(level,p);cand=g.candidate_edges_for_parent(nlower,p)
            pairs=list(itertools.combinations(cand,2))
            if level==1:
                vals=[]
                for e,f in pairs:
                    z=g.PT[:,pattern[e]]*g.PT[:,pattern[f]];vals.append(P@z)
                vals=np.array(vals);supp={d:(float((vals@np.array(d)).min()),float((vals@np.array(d)).max())) for d in DIRS}
            else:
                # union over possible local edge pairs: support extrema are extrema over pair outer hulls
                per=[local_pair_support(lowerD[e],lowerD[f],P) for e,f in pairs]
                supp={d:(min(s[d][0] for s in per),max(s[d][1] for s in per)) for d in DIRS}
            newD.append(generic_desc(2*ndesc,P,supp))
        allD.append(newD);lowerD=newD;ndesc*=2;level+=1
    return allD

def build_global(pattern,Ds):
    # Same fixed local tensor model as before, but every perspective uses recursive support facets.
    N=len(pattern);lp=g.SparseLP();msgs=[]
    for e,b in enumerate(pattern):
        vals=g.PT[:,b];msgs.append([lp.var(f'leaf{e}_{i}',float(vals[i]),float(vals[i])) for i in range(4)])
    prevD=None;level=1;ndesc=1
    while len(msgs)>2:
        lower=msgs;nlower=len(lower);m=nlower//2;cp=g.candidate_parents(nlower)
        a={}
        for e in range(nlower):
            for p in cp[e]:a[e,p]=lp.var(f'a{level}_{e}_{p}',0,1)
            lp.eql({a[e,p]:1 for p in cp[e]},1)
        H={};pairsP={}
        for p in range(m):
            pairs=list(itertools.combinations(g.candidate_edges_for_parent(nlower,p),2));pairsP[p]=pairs
            for e,f in pairs:H[p,e,f]=lp.var(f'H{level}_{p}_{e}_{f}',0,1)
            lp.eql({H[p,e,f]:1 for e,f in pairs},1)
            for e in g.candidate_edges_for_parent(nlower,p):
                dd={a[e,p]:-1}
                for u,v in pairs:
                    if e in (u,v):dd[H[p,u,v]]=dd.get(H[p,u,v],0)+1
                lp.eql(dd,0)
        W={}
        if level>1:
            for e in range(nlower):
                D=prevD[e]
                for p in cp[e]:
                    v=[lp.var(f'W{level}_{e}_{p}_{i}',0,D['hi']) for i in range(4)];W[e,p]=v
                    add_scaled_poly_support(lp,v,a[e,p],D)
                    b=lp.var(f'b{level}_{e}_{p}',0,1);lp.eql({b:1,a[e,p]:1},1)
                    R=[]
                    for i in range(4):
                        r=lp.var(f'R{level}_{e}_{p}_{i}',0,D['hi']);lp.eql({r:1,lower[e][i]:-1,v[i]:1},0);R.append(r)
                    add_scaled_poly_support(lp,R,b,D)
            for e in range(nlower):
                for i in range(4):lp.eql({**{W[e,p][i]:1 for p in cp[e]},lower[e][i]:-1},0)
        new=[]
        for p in range(m):
            pairAB={};diag=[]
            for e,f in pairsP[p]:
                h=H[p,e,f]
                if level==1:
                    vals=[]
                    for i in range(4):
                        co=float(g.PT[i,pattern[e]]*g.PT[i,pattern[f]]);v=lp.var(f't{level}_{p}_{e}_{f}_{i}',0,co);lp.eql({v:1,h:-co},0);vals.append(v)
                    diag.append(vals);continue
                De,Df=prevD[e],prevD[f]
                A=[lp.var(f'A{level}_{p}_{e}_{f}_{i}',0,De['hi']) for i in range(4)]
                B=[lp.var(f'B{level}_{p}_{e}_{f}_{i}',0,Df['hi']) for i in range(4)]
                add_scaled_poly_support(lp,A,h,De);add_scaled_poly_support(lp,B,h,Df);pairAB[e,f]=(A,B)
                C=[[lp.var(f'C{level}_{p}_{e}_{f}_{i}_{j}',0,De['hi']*Df['hi']) for j in range(4)] for i in range(4)]
                for i in range(4):
                    for j in range(4):g.add_persp_mcc(lp,C[i][j],A[i],B[j],h,De,Df)
                for i in range(4):add_scaled_poly_support(lp,C[i],A[i],Df)
                for j in range(4):add_scaled_poly_support(lp,[C[i][j] for i in range(4)],B[j],De)
                diag.append([C[i][i] for i in range(4)])
            if level>1:
                for e in g.candidate_edges_for_parent(nlower,p):
                    for i in range(4):
                        dd={W[e,p][i]:-1}
                        for u,v in pairsP[p]:
                            A,B=pairAB[u,v]
                            if e==u:dd[A[i]]=dd.get(A[i],0)+1
                            elif e==v:dd[B[i]]=dd.get(B[i],0)+1
                        lp.eql(dd,0)
            lz,uz=g.cb(2*ndesc);z=[]
            for i in range(4):
                v=lp.var(f'z{level}_{p}_{i}',lz,uz);dd={v:-1}
                for vv in diag:dd[vv[i]]=dd.get(vv[i],0)+1
                lp.eql(dd,0);z.append(v)
            P=g.edge_matrix(level,p);x=[]
            for i in range(4):
                v=lp.var(f'x{level}_{p}_{i}',lz,uz);dd={v:1}
                for j in range(4):dd[z[j]]=dd.get(z[j],0)-P[i,j]
                lp.eql(dd,0);x.append(v)
            D=Ds[level-1][p];one=lp.var(f'one{level}_{p}',1,1);add_scaled_poly_support(lp,x,one,D);new.append(x)
        msgs=new;prevD=Ds[level-1];ndesc*=2;level+=1
    x,y=msgs;Dx,Dy=prevD
    one=lp.var('rootone',1,1);C=[[lp.var(f'cr{i}{j}',0,Dx['hi']*Dy['hi']) for j in range(4)] for i in range(4)]
    for i in range(4):
        for j in range(4):g.add_persp_mcc(lp,C[i][j],x[i],y[j],one,Dx,Dy)
    for i in range(4):add_scaled_poly_support(lp,C[i],x[i],Dy)
    for j in range(4):add_scaled_poly_support(lp,[C[i][j] for i in range(4)],y[j],Dx)
    q=lp.var('q',0,1);lp.eql({q:1,**{C[i][i]:-.25 for i in range(4)}},0)
    val,res=lp.solve({q:1},maximize=True,scale=1e6)
    return val,len(lp.names),len(lp.ub)+len(lp.eq)

if __name__=='__main__':
    cases=[(16,'grouped_AC'),(16,'fourway_blocks'),(32,'grouped_AC'),(32,'fourway_blocks')]
    rows=[]
    for N,name in cases:
        if N==16:pat=g.patterns_for(16)[name];ex=g.exact_best(pat)
        else:
            pat=tuple([0]*16+[1]*16) if name=='grouped_AC' else tuple([0]*8+[1]*8+[2]*8+[3]*8);ex=float('nan')
        t=time.time();Ds=build_descriptors(pat);prep=time.time()-t
        t=time.time();v,nv,nc=build_global(pat,Ds);dt=time.time()-t
        row={'N':N,'pattern':name,'exact':ex,'upper':v,'ratio':v/ex if ex==ex else float('nan'),'vars':nv,'constraints':nc,'prep':prep,'solve':dt};rows.append(row);print(row,flush=True)
    pd.DataFrame(rows).to_csv('/mnt/data/support12_results.csv',index=False)
