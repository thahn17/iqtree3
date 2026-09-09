import itertools, importlib.util, contextlib, io, time
import numpy as np, pandas as pd
buf=io.StringIO()
with contextlib.redirect_stdout(buf):
    spec=importlib.util.spec_from_file_location('b','/mnt/data/test_varphylo2.py')
    b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b)
Pt,P2,P4=b.Pt,b.P2,b.P4; LP=b.LP; mass=b.mass; mcc=b.mcc; obbt=b.obbt
cb,mb=b.cb,b.mb; leaf_parents=b.leaf_parents; exact_best=b.exact_best
CLASSES=[(a,c) for a in range(4) for c in range(4)]
PAIRS=[(p,r) for p in range(4) for r in range(4) if p!=r]

def desc_from_points(points):
    pts=np.asarray(points)
    comp=[(float(pts[:,i].min()),float(pts[:,i].max())) for i in range(4)]
    vals=.25*pts.sum(axis=1); massb=(float(vals.min()),float(vals.max()))
    diff={}
    for i in range(4):
        for j in range(i+1,4):
            v=pts[:,i]-pts[:,j];diff[i,j]=(float(v.min()),float(v.max()))
    dirs=[]
    for s in itertools.product((-1,1),repeat=4):
        if s[0]!=1 or all(x==1 for x in s):continue
        v=pts@np.array(s);dirs.append((s,float(v.min()),float(v.max())))
    return {'comp':comp,'mass':massb,'diff':diff,'dirs':dirs}

def scaled_only(lp,w,a,D):
    for i,(lo,hi) in enumerate(D['comp']):
        lp.ge({w[i]:1,a:-lo},0);lp.le({w[i]:1,a:-hi},0)
    lo,hi=D['mass'];d={w[i]:.25 for i in range(4)};d[a]=-lo;lp.ge(d,0);d={w[i]:.25 for i in range(4)};d[a]=-hi;lp.le(d,0)
    for (i,j),(lo,hi) in D['diff'].items():
        lp.ge({w[i]:1,w[j]:-1,a:-lo},0);lp.le({w[i]:1,w[j]:-1,a:-hi},0)
    for s,lo,hi in D['dirs']:
        d={w[i]:s[i] for i in range(4)};d[a]=-lo;lp.ge(d,0)
        d={w[i]:s[i] for i in range(4)};d[a]=-hi;lp.le(d,0)

def persp_prod(lp,t,A,B,h,DA,DB):
    for i in range(4):
        lx,ux=DA['comp'][i];ly,uy=DB['comp'][i]
        lp.ge({t[i]:1,B[i]:-lx,A[i]:-ly,h:lx*ly},0)
        lp.ge({t[i]:1,B[i]:-ux,A[i]:-uy,h:ux*uy},0)
        lp.le({t[i]:1,B[i]:-ux,A[i]:-ly,h:ux*ly},0)
        lp.le({t[i]:1,B[i]:-lx,A[i]:-uy,h:lx*uy},0)

def build(pat):
    lp=LP();pat=tuple(pat)
    slots2=[(p,s) for p in range(4) for s in range(2)];y={}
    for e in range(8):
        ks=[k for k,(p,s) in enumerate(slots2) if p in leaf_parents[e]]
        for k in ks:y[e,k]=lp.var(f'y1_{e}_{k}',0,1)
        lp.eq({y[e,k]:1 for k in ks},1)
    for k in range(8):lp.eq({y[e,k]:1 for e in range(8) if (e,k) in y},1)
    Hbase={};classmsg={}
    for p in range(4):
        k0,k1=2*p,2*p+1;E0=[e for e in range(8) if (e,k0) in y];E1=[e for e in range(8) if (e,k1) in y]
        rb=[];cbv=[]
        for a in range(4):
            r=lp.var(f'r_{p}_{a}',0,1);d={r:1}
            for e in E0:
                if pat[e]==a:d[y[e,k0]]=d.get(y[e,k0],0)-1
            lp.eq(d,0);rb.append(r)
            c=lp.var(f'c_{p}_{a}',0,1);d={c:1}
            for e in E1:
                if pat[e]==a:d[y[e,k1]]=d.get(y[e,k1],0)-1
            lp.eq(d,0);cbv.append(c)
        elig=set(E0)|set(E1)
        for a,c in CLASSES:
            ub=0 if (a==c and sum(1 for e in elig if pat[e]==a)<2) else 1
            Hbase[p,a,c]=lp.var(f'Hb_{p}_{a}_{c}',0,ub)
            classmsg[p,a,c]=P2[p]@(Pt[:,a]*Pt[:,c])
        for a in range(4):
            d={rb[a]:-1}
            for c in range(4):d[Hbase[p,a,c]]=1
            lp.eq(d,0)
        for c in range(4):
            d={cbv[c]:-1}
            for a in range(4):d[Hbase[p,a,c]]=1
            lp.eq(d,0)
    # count2 routing latent classes
    y2={};J={}
    for p in range(4):
        for k in range(4):y2[p,k]=lp.var(f'y2_{p}_{k}',0,1)
        lp.eq({y2[p,k]:1 for k in range(4)},1)
    for k in range(4):lp.eq({y2[p,k]:1 for p in range(4)},1)
    for p in range(4):
        for k in range(4):
            for a,c in CLASSES:J[p,k,a,c]=lp.var(f'J_{p}_{k}_{a}_{c}',0,1)
            lp.eq({**{J[p,k,a,c]:1 for a,c in CLASSES},y2[p,k]:-1},0)
        for a,c in CLASSES:
            d={Hbase[p,a,c]:-1}
            for k in range(4):d[J[p,k,a,c]]=1
            lp.eq(d,0)
    # count4 local coupling; retain pair identity Hpair and weighted root-facing vector U
    Hpair={};Upair={};Dpair={};x4=[];q4=[]
    for q in range(2):
        k0,k1=2*q,2*q+1;K={}
        for p,r in PAIRS:
            # conditional point set for this local child-edge pair
            pts=[]
            for a,c in CLASSES:
                for d,e in CLASSES:
                    pts.append(P4[q]@(classmsg[p,a,c]*classmsg[r,d,e]))
            Dpair[q,p,r]=desc_from_points(pts)
            Hpair[q,p,r]=lp.var(f'Hpair_{q}_{p}_{r}',0,1)
            for a,c in CLASSES:
                for d,e in CLASSES:K[p,r,a,c,d,e]=lp.var(f'K_{q}_{p}_{r}_{a}_{c}_{d}_{e}',0,1)
            # Hpair = all local class mass for this edge pair
            dd={Hpair[q,p,r]:-1}
            for a,c in CLASSES:
                for d,e in CLASSES:dd[K[p,r,a,c,d,e]]=1
            lp.eq(dd,0)
        # J marginals
        for p in range(4):
            for a,c in CLASSES:
                dd={J[p,k0,a,c]:-1}
                for r in range(4):
                    if r==p:continue
                    for d,e in CLASSES:dd[K[p,r,a,c,d,e]]=1
                lp.eq(dd,0)
        for r in range(4):
            for d,e in CLASSES:
                dd={J[r,k1,d,e]:-1}
                for p in range(4):
                    if p==r:continue
                    for a,c in CLASSES:dd[K[p,r,a,c,d,e]]=1
                lp.eq(dd,0)
        # root-facing weighted vector for each pair identity
        for p,r in PAIRS:
            vec=[]
            for i in range(4):
                hi=Dpair[q,p,r]['comp'][i][1];u=lp.var(f'U_{q}_{p}_{r}_{i}',0,hi);dd={u:1}
                for a,c in CLASSES:
                    vp=classmsg[p,a,c]
                    for d,e in CLASSES:
                        coeff=(P4[q]@(vp*classmsg[r,d,e]))[i]
                        dd[K[p,r,a,c,d,e]]=dd.get(K[p,r,a,c,d,e],0)-coeff
                lp.eq(dd,0);vec.append(u)
            Upair[q,p,r]=vec
        # aggregate count4/root-facing vector
        vec=[]
        for i in range(4):
            v=lp.var(f'x4_{q}_{i}',0,1);dd={v:1}
            for p,r in PAIRS:dd[Upair[q,p,r][i]]=-1
            lp.eq(dd,0);vec.append(v)
        x4.append(vec);qq=mass(lp,vec,f'q4_{q}');q4.append(qq);l,u=mb(4);lp.bounds[qq]=(l,u);lp.ge({qq:1},l);lp.le({qq:1},u)
    # Tight aggregate count-4 message boxes before the root coupling.
    bx4=[]
    for q in range(2):
        bd=obbt(lp,x4[q]);bx4.append(bd)
        for i in range(4):lp.bounds[x4[q][i]]=bd[i]

    # root joint pair identity coupling, only complementary/disjoint count2 pairs allowed
    M={};A={};B={};T={}
    compatible=[]
    for p,r in PAIRS:
        for s,t in PAIRS:
            if len({p,r,s,t})<4:continue
            compatible.append((p,r,s,t));h=lp.var(f'M_{p}_{r}_{s}_{t}',0,1);M[p,r,s,t]=h
            avec=[];bvec=[];tvec=[]
            for i in range(4):
                ai=lp.var(f'A_{p}_{r}_{s}_{t}_{i}',0,Dpair[0,p,r]['comp'][i][1]);bi=lp.var(f'B_{p}_{r}_{s}_{t}_{i}',0,Dpair[1,s,t]['comp'][i][1]);ti=lp.var(f'T_{p}_{r}_{s}_{t}_{i}',0,Dpair[0,p,r]['comp'][i][1]*Dpair[1,s,t]['comp'][i][1]);avec.append(ai);bvec.append(bi);tvec.append(ti)
            scaled_only(lp,avec,h,Dpair[0,p,r]);scaled_only(lp,bvec,h,Dpair[1,s,t]);persp_prod(lp,tvec,avec,bvec,h,Dpair[0,p,r],Dpair[1,s,t])
            A[p,r,s,t]=avec;B[p,r,s,t]=bvec;T[p,r,s,t]=tvec
    # marginals M = Hpair
    for p,r in PAIRS:
        dd={Hpair[0,p,r]:-1}
        for s,t in PAIRS:
            if (p,r,s,t) in M:dd[M[p,r,s,t]]=1
        lp.eq(dd,0)
        for i in range(4):
            dd={Upair[0,p,r][i]:-1}
            for s,t in PAIRS:
                if (p,r,s,t) in M:dd[A[p,r,s,t][i]]=1
            lp.eq(dd,0)
    for s,t in PAIRS:
        dd={Hpair[1,s,t]:-1}
        for p,r in PAIRS:
            if (p,r,s,t) in M:dd[M[p,r,s,t]]=1
        lp.eq(dd,0)
        for i in range(4):
            dd={Upair[1,s,t][i]:-1}
            for p,r in PAIRS:
                if (p,r,s,t) in M:dd[B[p,r,s,t][i]]=1
            lp.eq(dd,0)
    z=[]
    for i in range(4):
        lo,hi=cb(8);v=lp.var(f'z8_{i}',lo,hi);dd={v:-1}
        for key in compatible:dd[T[key][i]]=1
        lp.eq(dd,0);z.append(v)
    # Intersect the disaggregated pair-identity root representation with the
    # ordinary aggregate McCormick hull. This prevents the disaggregation from
    # weakening easy cases such as strongly interleaved patterns.
    for i in range(4):
        lx,ux=bx4[0][i];ly,uy=bx4[1][i]
        mcc(lp,z[i],x4[0][i],x4[1][i],lx,ux,ly,uy)
    q8=mass(lp,z,'q8');l,u=mb(8);lp.bounds[q8]=(l,u);lp.ge({q8:1},l);lp.le({q8:1},u)
    lp.le({q8:1,q4[0]:-1},0);lp.le({q8:1,q4[1]:-1},0)
    val,res=lp.solve({q8:1e8},maximize=True)
    return val/1e8,len(lp.names),len(lp.A)+len(lp.E)

if __name__=='__main__':
    pat=(0,1,2,3,0,1,2,3)
    ex=exact_best(pat);t=time.time();v,nv,nc=build(pat);dt=time.time()-t
    print('ACGTACGT',ex,v,v/ex,nv,nc,dt)
