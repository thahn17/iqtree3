import itertools, importlib.util, contextlib, io, time
import numpy as np, pandas as pd
buf=io.StringIO()
with contextlib.redirect_stdout(buf):
    spec=importlib.util.spec_from_file_location('b','/mnt/data/test_varphylo2.py')
    b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b)
Pt=b.Pt;P2=b.P2;P4=b.P4;LP=b.LP;leaf_parents=b.leaf_parents

# Eligible unordered lower-edge pairs for each count-2 parent.
elig={}
for p in range(4):
    es=[e for e in range(8) if p in leaf_parents[e]]
    elig[p]=list(itertools.combinations(es,2))

def build(pat):
    lp=LP()
    # h[p,pair] = continuous marginal that this parent receives exactly this lower-edge pair.
    h={}
    for p in range(4):
        for ef in elig[p]:h[p,ef]=lp.var(f'h_{p}_{ef[0]}_{ef[1]}',0,1)
        lp.eq({h[p,ef]:1 for ef in elig[p]},1)
    # Every lower edge is used by exactly one of its two possible parents.
    for e in range(8):
        d={}
        for p in leaf_parents[e]:
            for ef in elig[p]:
                if e in ef:d[h[p,ef]]=d.get(h[p,ef],0)+1
        lp.eq(d,1)

    # Exact signal template for each parent/pair.
    temp={}
    for p in range(4):
        for ef in elig[p]:
            e,f=ef
            temp[p,ef]=P2[p]@(Pt[:,pat[e]]*Pt[:,pat[f]])

    # Route each pair-marginal to one of four upper child slots.
    # J splits h exactly; no a*x multiplication exists.
    J={};a={}
    for p in range(4):
        for k in range(4):
            a[p,k]=lp.var(f'a_{p}_{k}',0,1)
            for ef in elig[p]:J[p,k,ef]=lp.var(f'J_{p}_{k}_{ef[0]}_{ef[1]}',0,1)
            lp.eq({**{J[p,k,ef]:1 for ef in elig[p]},a[p,k]:-1},0)
        lp.eq({a[p,k]:1 for k in range(4)},1)
        for ef in elig[p]:
            d={h[p,ef]:-1}
            for k in range(4):d[J[p,k,ef]]=1
            lp.eq(d,0)
    for k in range(4):lp.eq({a[p,k]:1 for p in range(4)},1)

    # K at each count-4 parent: joint selection of two distinct count-2 children and their lower-edge pairs.
    K=[];cfgsig=[]
    for q in range(2):
        k0,k1=2*q,2*q+1;kd={};sig={}
        for p in range(4):
            for r in range(4):
                if p==r:continue
                for ef in elig[p]:
                    for gh in elig[r]:
                        # In a discrete scenario, the two count2 blocks cannot share a lower edge.
                        if set(ef)&set(gh):continue
                        cfg=(p,r,ef,gh);v=lp.var(f'K_{q}_{p}_{r}_{ef[0]}_{ef[1]}_{gh[0]}_{gh[1]}',0,1);kd[cfg]=v
                        z=temp[p,ef]*temp[r,gh]
                        sig[cfg]=P4[q]@z
        for p in range(4):
            for ef in elig[p]:
                d={J[p,k0,ef]:-1}
                for cfg,v in kd.items():
                    if cfg[0]==p and cfg[2]==ef:d[v]=d.get(v,0)+1
                lp.eq(d,0)
        for r in range(4):
            for gh in elig[r]:
                d={J[r,k1,gh]:-1}
                for cfg,v in kd.items():
                    if cfg[1]==r and cfg[3]==gh:d[v]=d.get(v,0)+1
                lp.eq(d,0)
        K.append(kd);cfgsig.append(sig)

    # Global root joint marginal.  It preserves both count2-child identity and the underlying lower-edge pairing.
    M={};C0={c:[] for c in K[0]};C1={c:[] for c in K[1]}
    for c0 in K[0]:
        p,r,ef,gh=c0
        leaves0=set(ef)|set(gh)
        for c1 in K[1]:
            s,v,ij,kl=c1
            if len({p,r,s,v})!=4:continue
            leaves1=set(ij)|set(kl)
            if leaves0 & leaves1:continue
            # Four disjoint pairs imply all 8 lower edges are represented once.
            if len(leaves0|leaves1)!=8:continue
            m=lp.var('M_'+str(len(M)),0,1);M[c0,c1]=m;C0[c0].append((c1,m));C1[c1].append((c0,m))
    for c0,v0 in K[0].items():
        d={v0:-1}
        for c1,m in C0[c0]:d[m]=1
        lp.eq(d,0)
    for c1,v1 in K[1].items():
        d={v1:-1}
        for c0,m in C1[c1]:d[m]=1
        lp.eq(d,0)

    scale=1e10
    obj={m:scale*0.25*float(np.dot(cfgsig[0][c0],cfgsig[1][c1])) for (c0,c1),m in M.items()}
    val,sol=lp.solve(obj,maximize=True)
    val/=scale
    return val,len(lp.names),len(lp.A)+len(lp.E),len(K[0]),len(M)

if __name__=='__main__':
    cases=[
      ('AAAACCCC',(0,0,0,0,1,1,1,1),0.006220127561334139),
      ('AACCGGTT',(0,0,1,1,2,2,3,3),0.00032032190163812385),
      ('AAAACCGT',(0,0,0,0,1,1,2,3),0.00016375505988646638),
      ('ACGTACGT',(0,1,2,3,0,1,2,3),5.5485524598073915e-08),
    ]
    rows=[]
    for name,pat,ex in cases:
        t=time.time();v,n,c,k,m=build(pat);row={'pattern':name,'exact benchmark':ex,'edge-pair marginal LP':v,'LP/exact':v/ex,'variables':n,'constraints':c,'K0 configs':k,'root M configs':m,'seconds':time.time()-t};rows.append(row);print(row,flush=True)
    pd.DataFrame(rows).to_csv('/mnt/data/edgepair_marginal_results.csv',index=False)
