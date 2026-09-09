#!/usr/bin/env python3
"""Research LP: lower-arm branch lengths + shared named-taxon flow tightening.

Key ideas:
- branch length selector belongs to each split-unit lower port (left/right), shared across sites;
- routers carry raw child partials; transition is applied at the requesting parent port;
- optional shared named-taxon membership is propagated through the SAME physical Beneš switches;
- every taxon's memberships across root candidates sum exactly to 1;
- count-1 likelihood requests are the linear image of that same taxon membership;
- optional root-partition likelihood facets condition each site's root contribution on the
  actual shared fractional taxa assigned to the root's left side.

Pure continuous LP.
"""
from __future__ import annotations
import importlib.util,sys,json,math,time
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np
from scipy.linalg import expm
from scipy.optimize import linprog

spec=importlib.util.spec_from_file_location('vl','/mnt/data/multisite_variable_length_compressed.py')
vl=importlib.util.module_from_spec(spec);sys.modules['vl']=vl;spec.loader.exec_module(vl)
SparseLP=vl.SparseLP; add_switch_hull=vl.add_switch_hull
PI=vl.PI;Q=vl.Q;BASES=vl.BASES;BASE_INDEX=vl.BASE_INDEX


def branch_times(m): return vl.short_time(m)
def trans_pair(m,long_ratio):
    t=branch_times(m);return expm(Q*t),expm(Q*t*long_ratio)

def raw_boxes(global_counts,long_ratio=2.0,variable_lengths=True):
    """Raw (pre-parent-edge) partial boxes by (count,composition)."""
    g=tuple(map(int,global_counts));L=sum(g)
    raw={};by=defaultdict(list)
    for b in range(4):
        if g[b]<=0:continue
        c=tuple(1 if k==b else 0 for k in range(4));v=np.zeros(4);v[b]=1.0
        raw[(1,c)]=np.stack([v,v],axis=1);by[1].append(c)
    for m in range(2,L+1):
        acc={}
        for a in range(1,m//2+1):
            b=m-a
            if not by[a] or not by[b]:continue
            PaS,PaL=trans_pair(a,long_ratio);PbS,PbL=trans_pair(b,long_ratio)
            for c1 in by[a]:
                A=raw[(a,c1)]; alo,ahi=A[:,0],A[:,1]
                amlos=[PaS@alo];amhis=[PaS@ahi]
                if variable_lengths:amlos.append(PaL@alo);amhis.append(PaL@ahi)
                amlo=np.min(np.stack(amlos),axis=0);amhi=np.max(np.stack(amhis),axis=0)
                for c2 in by[b]:
                    c=tuple(c1[k]+c2[k] for k in range(4))
                    if any(c[k]>g[k] for k in range(4)):continue
                    B=raw[(b,c2)];blo,bhi=B[:,0],B[:,1]
                    bmlos=[PbS@blo];bmhis=[PbS@bhi]
                    if variable_lengths:bmlos.append(PbL@blo);bmhis.append(PbL@bhi)
                    bmlo=np.min(np.stack(bmlos),axis=0);bmhi=np.max(np.stack(bmhis),axis=0)
                    zlo=amlo*bmlo;zhi=amhi*bmhi
                    if c not in acc:acc[c]=[zlo.copy(),zhi.copy()]
                    else:
                        acc[c][0]=np.minimum(acc[c][0],zlo);acc[c][1]=np.maximum(acc[c][1],zhi)
        by[m]=list(acc)
        for c,(lo,hi) in acc.items():raw[(m,c)]=np.stack([lo,hi],axis=1)
    return raw

def root_partition_values(counts,raw,a,long_ratio,variable_lengths):
    L=sum(counts);b=L-a; vals={}
    PaS,PaL=trans_pair(a,long_ratio);PbS,PbL=trans_pair(b,long_ratio)
    for (mm,c1),A in raw.items():
        if mm!=a:continue
        c2=tuple(counts[k]-c1[k] for k in range(4))
        if min(c2)<0 or (b,c2) not in raw:continue
        ahi=A[:,1];bhi=raw[(b,c2)][:,1]
        au=[PaS@ahi];bu=[PbS@bhi]
        if variable_lengths:au.append(PaL@ahi);bu.append(PbL@bhi)
        muA=np.max(np.stack(au),axis=0);muB=np.max(np.stack(bu),axis=0)
        vals[c1]=float(PI@(muA*muB))
    return vals

def affine_majorants(comps,vals,max_facets=20):
    # y <= a0+dA*A+dC*C+dG*G; exact/tight supports at representative compositions.
    X=np.array([[1,c[0],c[1],c[2]] for c in comps],float);y=np.array(vals,float)
    if len(X)==0:return []
    idx=list(range(len(X))) if len(X)<=max_facets else sorted(set(np.linspace(0,len(X)-1,max_facets,dtype=int).tolist()))
    fs=[]
    for ii in idx:
        res=linprog(X[ii],A_ub=-X,b_ub=-y,bounds=[(None,None)]*4,method='highs')
        if res.success:
            f=tuple(np.round(res.x,11))
            if f not in fs and np.all(X@res.x+1e-8>=y):fs.append(f)
    return fs

class Model:
    def __init__(self,libpath,columns,log_ratio=1.1,long_ratio=2.0,variable_lengths=True,taxon_flow=True,root_partition=True):
        cnt=Counter(tuple(c) for c in columns);self.patterns=list(cnt);self.weights=np.array([cnt[p] for p in self.patterns],float)
        self.lib=json.loads(Path(libpath).read_text());self.L=int(self.lib['summary']['leaf_count']);self.units=self.lib['units'];self.roots=set(self.lib['root_candidates'])
        self.splits=[u for u,x in self.units.items() if x['kind']=='split'];self.tax=[u for u,x in self.units.items() if x['kind']=='taxon']
        self.by=defaultdict(list)
        for u in self.splits:self.by[int(self.units[u]['descendant_count'])].append(u)
        self.reqdef=defaultdict(list)
        for u in self.splits:
            for port in ('left','right'):
                f=int(self.units[u][port+'_flow']);self.reqdef[f].append((u,port))
        self.log_ratio=log_ratio;self.long_ratio=long_ratio;self.variable_lengths=variable_lengths;self.taxon_flow=taxon_flow;self.root_partition=root_partition
        self.lp=SparseLP();self.act={u:self.lp.var('act:'+u,0,1) for u in self.splits};self.lp.add_eq({self.act[u]:1 for u in self.roots},1)
        self.wireact={};self.sw={};self.routermeta={};self._build_topology()
        # lower-arm long activities, one per split port
        self.ell={}
        if variable_lengths:
            for u in self.splits:
                for port in ('left','right'):
                    e=self.lp.var(f'long:{u}:{port}',0,1);self.ell[u,port]=e;self.lp.add_le({e:1,self.act[u]:-1},0)
        self.tmem={}
        if taxon_flow:self._build_taxon_flow()
        self.info=[];self._prep()
        self.req={};self.hlong={};self.z={};self.q={};self.lam={};self.obj={};self.log_error_bound=0
        self._build_sites()
        if taxon_flow and root_partition:self._add_root_partition_facets()
    def _build_topology(self):
        for mstr,r in self.lib['routers'].items():
            m=int(mstr);W=int(r['padded_width']);dims=list(r['dimensions']);D=len(dims);self.routermeta[m]=(r,W,D)
            if D==0:
                rr=r['input_requests'][0];oo=r['output_targets'][0];a=self.act[rr['parent_unit']]
                if m==1:self.lp.add_eq({a:1},1)
                else:self.lp.add_eq({a:1,self.act[oo['unit_id']]:-1},0)
                continue
            for l in range(D+1):
                for w in range(W):self.wireact[m,l,w]=self.lp.var(f'wa:{m}:{l}:{w}',0,1)
            rin={int(x['slot']):x for x in r['input_requests']};rout={int(x['slot']):x for x in r['output_targets']}
            for w in range(W):
                v=self.wireact[m,0,w]
                if w in rin:self.lp.add_eq({v:1,self.act[rin[w]['parent_unit']]:-1},0)
                else:self.lp.add_eq({v:1},0)
                v=self.wireact[m,D,w]
                if w in rout:
                    if m==1:self.lp.add_eq({v:1},1)
                    else:self.lp.add_eq({v:1,self.act[rout[w]['unit_id']]:-1},0)
                else:self.lp.add_eq({v:1},0)
            for st,dim in enumerate(dims):
                seen=set()
                for w in range(W):
                    b=min(w,w^(1<<dim))
                    if b in seen:continue
                    seen.add(b);o=b^(1<<dim);s=self.lp.var(f'sw:{m}:{st}:{b}',0,1);self.sw[m,st,b]=s
                    add_switch_hull(self.lp,self.wireact[m,st,b],self.wireact[m,st,o],self.wireact[m,st+1,b],self.wireact[m,st+1,o],s)
    def _req_slot(self,u,port):
        m=int(self.units[u][port+'_flow']);r,W,D=self.routermeta[m]
        slot=next(int(x['slot']) for x in r['input_requests'] if x['parent_unit']==u and x['parent_port']==port)
        return m,slot,D
    def _build_taxon_flow(self):
        L=self.L; lp=self.lp
        # create memberships on all router wires, including one wire for D=0
        for m,(r,W,D) in self.routermeta.items():
            layers=range(D+1) if D>0 else [0]
            for l in layers:
                for w in range(W):
                    if D==0:a=self.act[r['input_requests'][0]['parent_unit']]
                    else:a=self.wireact[m,l,w]
                    row={a:-float(m)}
                    for t in range(L):
                        v=lp.var(f'tm:{m}:{l}:{w}:{t}',0,1);self.tmem[m,l,w,t]=v;lp.add_le({v:1,a:-1},0);row[v]=1
                    lp.add_eq(row,0)
            if D>0:
                for st,dim in enumerate(r['dimensions']):
                    seen=set()
                    for w in range(W):
                        b=min(w,w^(1<<dim))
                        if b in seen:continue
                        seen.add(b);o=b^(1<<dim);s=self.sw[m,st,b]
                        for t in range(L):add_switch_hull(lp,self.tmem[m,st,b,t],self.tmem[m,st,o,t],self.tmem[m,st+1,b,t],self.tmem[m,st+1,o,t],s)
        taxidx={uid:int(self.units[uid]['label'].split('_')[-1])-1 for uid in self.tax}
        # output endpoints
        for m,(r,W,D) in self.routermeta.items():
            outlayer=D if D>0 else 0
            for x in r['output_targets']:
                w=int(x['slot']);uid=x['unit_id']
                if m==1:
                    tt=taxidx[uid]
                    for t in range(L):lp.add_eq({self.tmem[m,outlayer,w,t]:1},1 if t==tt else 0)
                else:
                    for t in range(L):
                        row={self.tmem[m,outlayer,w,t]:1}
                        for port in ('left','right'):
                            fm,slot,Df=self._req_slot(uid,port);v=self.tmem[fm,0,slot,t];row[v]=row.get(v,0)-1
                        lp.add_eq(row,0)
        # Split disjointness: a named taxon can occur at most once inside any active subtree.
        # For roots the relation is equality because every full root contains every taxon.
        for u in self.splits:
            for t in range(L):
                row={}
                for port in ('left','right'):
                    m,slot,D=self._req_slot(u,port);v=self.tmem[m,0,slot,t];row[v]=row.get(v,0)+1
                if u in self.roots:
                    row[self.act[u]]=row.get(self.act[u],0)-1;lp.add_eq(row,0)
                else:
                    row[self.act[u]]=row.get(self.act[u],0)-1;lp.add_le(row,0)
        # Redundant global certificate: every named taxon reaches the root mixture with total mass exactly 1.
        for t in range(L):
            row={}
            for u in self.roots:
                for port in ('left','right'):
                    m,slot,D=self._req_slot(u,port);v=self.tmem[m,0,slot,t];row[v]=row.get(v,0)+1
            lp.add_eq(row,1)
    def _prep(self):
        for pat in self.patterns:
            counts=tuple(pat.count(b) for b in BASES);raw=raw_boxes(counts,self.long_ratio,self.variable_lengths)
            scale={};ubraw={};ubmsg={}
            for m in range(1,self.L+1):
                Bs=[B for (mm,c),B in raw.items() if mm==m]
                hi=np.max(np.stack([B[:,1] for B in Bs]),axis=0);sc=max(float(hi.max()),1e-300);scale[m]=sc;ubraw[m]=hi/sc
                if m<self.L:
                    Ps,Pl=trans_pair(m,self.long_ratio); his=[Ps@(hi/sc)]
                    if self.variable_lengths:his.append(Pl@(hi/sc))
                    ubmsg[m]=np.max(np.stack(his),axis=0)
            rootvals={}
            for u in self.roots:
                a=int(self.units[u]['left_flow']); vals=root_partition_values(counts,raw,a,self.long_ratio,self.variable_lengths);rootvals[u]=vals
            self.info.append({'counts':counts,'raw':raw,'scale':scale,'ubraw':ubraw,'ubmsg':ubmsg,'rootvals':rootvals})
    def _request_message_expr(self,p,u,port,i):
        m=int(self.units[u][port+'_flow']);R=self.req[p,u,port];Ps,Pl=trans_pair(m,self.long_ratio)
        d={R[j]:float(Ps[i,j]) for j in range(4)}
        if self.variable_lengths:
            H=self.hlong[p,u,port]
            for j in range(4):d[H[j]]=d.get(H[j],0)+float(Pl[i,j]-Ps[i,j])
        return d
    def _build_sites(self):
        for p,pat in enumerate(self.patterns):
            I=self.info[p]
            # raw request vectors and lower-arm length split
            for m in range(1,self.L):
                ub=I['ubraw'][m]
                for u,port in self.reqdef[m]:
                    a=self.act[u];rv=[];hh=[]
                    for j in range(4):
                        v=self.lp.var(f'req:{p}:{u}:{port}:{j}',0,float(ub[j]));rv.append(v);self.lp.add_le({v:1,a:-float(ub[j])},0)
                    self.req[p,u,port]=rv
                    if self.variable_lengths:
                        ell=self.ell[u,port]
                        for j,v in enumerate(rv):
                            U=float(ub[j]);h=self.lp.var(f'h:{p}:{u}:{port}:{j}',0,U);hh.append(h)
                            self.lp.add_le({h:1,v:-1},0);self.lp.add_le({h:1,ell:-U},0);self.lp.add_le({v:1,h:-1,a:-U,ell:U},0)
                        self.hlong[p,u,port]=hh
            # if taxon flow, count-1 raw vector is exact composition of the named taxa on that request
            if self.taxon_flow:
                for u,port in self.reqdef[1]:
                    m,slot,D=self._req_slot(u,port);R=self.req[p,u,port]
                    for j in range(4):
                        row={R[j]:1}
                        for t,b in enumerate(pat):
                            if BASE_INDEX[b]==j:
                                v=self.tmem[1,0,slot,t];row[v]=row.get(v,0)-1
                        self.lp.add_eq(row,0)
            # parent raw product
            for u in self.splits:
                aN=int(self.units[u]['left_flow']);bN=int(self.units[u]['right_flow']);m=int(self.units[u]['descendant_count']);act=self.act[u]
                k=I['scale'][aN]*I['scale'][bN]/I['scale'][m];uA=I['ubmsg'][aN];uB=I['ubmsg'][bN];zz=[]
                for i in range(4):
                    z=self.lp.var(f'z:{p}:{u}:{i}',0,float(I['ubraw'][m][i]));zz.append(z)
                    A=self._request_message_expr(p,u,'left',i);B=self._request_message_expr(p,u,'right',i)
                    # z/k perspective McCormick, [0,uA]x[0,uB]
                    row={z:1}
                    for v,c in A.items():row[v]=row.get(v,0)-k*float(uB[i])*c
                    for v,c in B.items():row[v]=row.get(v,0)-k*float(uA[i])*c
                    row[act]=row.get(act,0)+k*float(uA[i]*uB[i]);self.lp.add_ge(row,0)
                    row={z:1}
                    for v,c in A.items():row[v]=row.get(v,0)-k*float(uB[i])*c
                    self.lp.add_le(row,0)
                    row={z:1}
                    for v,c in B.items():row[v]=row.get(v,0)-k*float(uA[i])*c
                    self.lp.add_le(row,0)
                self.z[p,u]=zz
            # raw count conservation through each router
            for m in range(1,self.L):
                for j in range(4):
                    row={self.req[p,u,port][j]:1 for u,port in self.reqdef[m]}
                    if m==1:
                        if not self.taxon_flow:
                            rhs=sum(1.0 for b in pat if BASE_INDEX[b]==j);self.lp.add_eq(row,rhs)
                        else:self.lp.add_eq(row,sum(1.0 for b in pat if BASE_INDEX[b]==j)) # redundant, useful numerically
                    else:
                        for u in self.by[m]:row[self.z[p,u][j]]=row.get(self.z[p,u][j],0)-1
                        self.lp.add_eq(row,0)
            # root q
            q=self.lp.var(f'q:{p}',0,1);self.q[p]=q;row={q:1}
            for u in self.roots:
                for i in range(4):row[self.z[p,u][i]]=row.get(self.z[p,u][i],0)-float(PI[i])
                # old constant root cut retained as baseline
                vals=I['rootvals'][u];hi=max(vals.values())/I['scale'][self.L]
                expr={self.z[p,u][i]:float(PI[i]) for i in range(4)};expr[self.act[u]]=expr.get(self.act[u],0)-hi;self.lp.add_le(expr,0)
            self.lp.add_eq(row,0)
            # one-pass log chords
            # conservative lower floor for interpolation; keep numerical floor moderate in test harness
            pts=vl.geometric_breaks(1e-12,self.log_ratio);lam=[]
            for kk,x in enumerate(pts):
                v=self.lp.var(f'lam:{p}:{kk}',0,1);lam.append(v);self.obj[v]=self.weights[p]*math.log(x)
            self.lam[p]=lam;self.lp.add_eq({v:1 for v in lam},1);self.lp.add_eq({q:1,**{lam[k]:-pts[k] for k in range(len(pts))}},0)
            self.log_error_bound+=self.weights[p]*vl.chord_gap(self.log_ratio)
    def _add_root_partition_facets(self):
        self.root_facets=0
        for p,pat in enumerate(self.patterns):
            I=self.info[p]
            for u in self.roots:
                aN=int(self.units[u]['left_flow']);vals=I['rootvals'][u]
                comps=list(vals);y=[vals[c]/I['scale'][self.L] for c in comps];fac=affine_majorants(comps,y,30)
                m,slot,D=self._req_slot(u,'left');act=self.act[u];lexpr={self.z[p,u][i]:float(PI[i]) for i in range(4)}
                for a0,dA,dC,dG in fac:
                    row=dict(lexpr);row[act]=row.get(act,0)-a0
                    for t,b in enumerate(pat):
                        coeff={'A':dA,'C':dC,'G':dG,'T':0.0}[b]
                        if coeff:
                            v=self.tmem[m,0,slot,t];row[v]=row.get(v,0)-coeff
                    self.lp.add_le(row,0);self.root_facets+=1
    def solve(self):
        val,res,dt=self.lp.solve(self.obj,True)
        const=sum(self.weights[p]*math.log(self.info[p]['scale'][self.L]) for p in range(len(self.patterns)))
        return val+const,res,dt
    def root_taxon_sums(self,res):
        if not self.taxon_flow:return None
        out=[]
        for t in range(self.L):
            z=0
            for u in self.roots:
                for port in ('left','right'):
                    m,slot,D=self._req_slot(u,port);z+=res.x[self.tmem[m,0,slot,t]]
            out.append(z)
        return out

if __name__=='__main__':
    cols=[tuple('AAAACCCC'),tuple('ACGTACGT'),tuple('AACCGGTT'),tuple('AGCTAGCT'),tuple('AAAACCGT'),tuple('AACCGGTT')]
    lib='/mnt/data/varlen_units8/library.json'
    for tf,rp in [(False,False),(True,False),(True,True)]:
        M=Model(lib,cols,1.1,2.0,True,tf,rp);v,r,t=M.solve();print('tf',tf,'rp',rp,'obj',v,'vars',len(M.lp.names),'eq',len(M.lp.eq),'le',len(M.lp.le),'solve',t,'rootfac',getattr(M,'root_facets',0))
        if tf:print('rootsum',min(M.root_taxon_sums(r)),max(M.root_taxon_sums(r)))


def affine_majorants_nd(features,vals,max_facets=40):
    X=np.column_stack([np.ones(len(features)),np.asarray(features,float)]);y=np.asarray(vals,float)
    if not len(X):return []
    # target every point when small, otherwise evenly spaced plus extrema by value
    if len(X)<=max_facets:idx=list(range(len(X)))
    else:
        idx=set(np.linspace(0,len(X)-1,max_facets-2,dtype=int).tolist());idx.add(int(np.argmin(y)));idx.add(int(np.argmax(y)));idx=sorted(idx)
    fs=[]
    for ii in idx:
        res=linprog(X[ii],A_ub=-X,b_ub=-y,bounds=[(None,None)]*X.shape[1],method='highs')
        if res.success:
            f=tuple(np.round(res.x,11))
            if f not in fs and np.all(X@res.x+1e-8>=y):fs.append(f)
    return fs

class RootLengthPartitionModel(Model):
    def __init__(self,libpath,columns,log_ratio=1.1,long_ratio=2.0):
        # Build taxon flow but skip composition-only root facets; add length-aware facets instead.
        Model.__init__(self,libpath,columns,log_ratio,long_ratio,True,True,False)
        self._add_root_length_partition_facets()
    def _add_root_length_partition_facets(self):
        self.root_length_facets=0
        for p,pat in enumerate(self.patterns):
            I=self.info[p];counts=I['counts'];raw=I['raw']
            for u in self.roots:
                aN=int(self.units[u]['left_flow']);bN=int(self.units[u]['right_flow']);PsA,PlA=trans_pair(aN,self.long_ratio);PsB,PlB=trans_pair(bN,self.long_ratio)
                feats=[];vals=[]
                for (mm,c1),A in raw.items():
                    if mm!=aN:continue
                    c2=tuple(counts[k]-c1[k] for k in range(4))
                    if min(c2)<0 or (bN,c2) not in raw:continue
                    B=raw[(bN,c2)];
                    for la,PA in enumerate((PsA,PlA)):
                        for lb,PB in enumerate((PsB,PlB)):
                            val=float(PI@((PA@A[:,1])*(PB@B[:,1])))/I['scale'][self.L]
                            feats.append((c1[0],c1[1],c1[2],la,lb));vals.append(val)
                fac=affine_majorants_nd(feats,vals,60)
                m,slot,D=self._req_slot(u,'left');act=self.act[u];eL=self.ell[u,'left'];eR=self.ell[u,'right'];lexpr={self.z[p,u][i]:float(PI[i]) for i in range(4)}
                for f in fac:
                    a0,dA,dC,dG,dL,dR=f;row=dict(lexpr);row[act]=row.get(act,0)-a0;row[eL]=row.get(eL,0)-dL;row[eR]=row.get(eR,0)-dR
                    for t,b in enumerate(pat):
                        coeff={'A':dA,'C':dC,'G':dG,'T':0.0}[b]
                        if coeff:
                            v=self.tmem[m,0,slot,t];row[v]=row.get(v,0)-coeff
                    self.lp.add_le(row,0);self.root_length_facets+=1

if __name__=='__main__':
    print('--- length aware root facets ---')
    cols=[tuple('AAAACCCC'),tuple('ACGTACGT'),tuple('AACCGGTT'),tuple('AGCTAGCT'),tuple('AAAACCGT'),tuple('AACCGGTT')]
    M=RootLengthPartitionModel('/mnt/data/varlen_units8/library.json',cols,1.1,2.0);v,r,t=M.solve();print('obj',v,'vars',len(M.lp.names),'eq',len(M.lp.eq),'le',len(M.lp.le),'fac',M.root_length_facets,'solve',t,'rootsum',min(M.root_taxon_sums(r)),max(M.root_taxon_sums(r)))

class BoundaryTaxonModel(RootLengthPartitionModel):
    def _build_taxon_flow(self):
        # Projection: named-taxon membership only at real router request boundaries.
        # Physical Beneš activity/switch topology remains shared, but taxon identities are conserved by count globally.
        L=self.L;lp=self.lp
        self.tmem={}
        # request membership vectors
        for m,(r,W,D) in self.routermeta.items():
            for x in r['input_requests']:
                slot=int(x['slot']);u=x['parent_unit'];port=x['parent_port'];a=self.act[u]
                row={a:-float(m)}
                for t in range(L):
                    v=lp.var(f'bm:{m}:{slot}:{t}',0,1);self.tmem[m,0,slot,t]=v;lp.add_le({v:1,a:-1},0);row[v]=1
                lp.add_eq(row,0)
        # local split disjointness; root equality
        for u in self.splits:
            for t in range(L):
                row={}
                for port in ('left','right'):
                    m,slot,D=self._req_slot(u,port);v=self.tmem[m,0,slot,t];row[v]=row.get(v,0)+1
                row[self.act[u]]=row.get(self.act[u],0)-1
                if u in self.roots:lp.add_eq(row,0)
                else:lp.add_le(row,0)
        # count-wise taxon conservation through the universal router projection
        for m in range(1,self.L):
            for t in range(L):
                row={}
                for u,port in self.reqdef[m]:
                    mm,slot,D=self._req_slot(u,port);v=self.tmem[m,0,slot,t];row[v]=row.get(v,0)+1
                if m==1:
                    lp.add_eq(row,1)
                else:
                    for child in self.by[m]:
                        for port in ('left','right'):
                            mm,slot,D=self._req_slot(child,port);v=self.tmem[mm,0,slot,t];row[v]=row.get(v,0)-1
                    lp.add_eq(row,0)
        # explicit redundant root certificate
        for t in range(L):
            row={}
            for u in self.roots:
                for port in ('left','right'):
                    m,slot,D=self._req_slot(u,port);v=self.tmem[m,0,slot,t];row[v]=row.get(v,0)+1
            lp.add_eq(row,1)

if __name__=='__main__':
    print('--- boundary taxon projection ---')
    cols=[tuple('AAAACCCC'),tuple('ACGTACGT'),tuple('AACCGGTT'),tuple('AGCTAGCT'),tuple('AAAACCGT'),tuple('AACCGGTT')]
    M=BoundaryTaxonModel('/mnt/data/varlen_units8/library.json',cols,1.1,2.0);v,r,t=M.solve();print('obj',v,'vars',len(M.lp.names),'eq',len(M.lp.eq),'le',len(M.lp.le),'fac',M.root_length_facets,'solve',t,'rootsum',min(M.root_taxon_sums(r)),max(M.root_taxon_sums(r)))

class BoundaryShapeModel(BoundaryTaxonModel):
    def __init__(self,libpath,columns,log_ratio=1.1,long_ratio=2.0,max_shape_facets=8):
        self.max_shape_facets=max_shape_facets
        BoundaryTaxonModel.__init__(self,libpath,columns,log_ratio,long_ratio)
        self._add_membership_shape_facets()
    def _add_membership_shape_facets(self):
        self.shape_facets=0
        for p,pat in enumerate(self.patterns):
            I=self.info[p];raw=I['raw'];counts=I['counts']
            cache={}
            for m in range(1,self.L+1):
                comps=[c for (mm,c) in raw if mm==m]
                if not comps:continue
                for j in range(4):
                    vals=[float(raw[(m,c)][j,1]/I['scale'][m]) for c in comps]
                    cache[m,j]=affine_majorants(comps,vals,self.max_shape_facets)
            # request raw vectors: composition = membership on that parent port
            for m in range(2,self.L):
                for u,port in self.reqdef[m]:
                    _,slot,_=self._req_slot(u,port);act=self.act[u]
                    for j in range(4):
                        for a0,dA,dC,dG in cache.get((m,j),[]):
                            row={self.req[p,u,port][j]:1,act:-a0}
                            for t,b in enumerate(pat):
                                coeff={'A':dA,'C':dC,'G':dG,'T':0.0}[b]
                                if coeff:
                                    v=self.tmem[m,0,slot,t];row[v]=row.get(v,0)-coeff
                            self.lp.add_le(row,0);self.shape_facets+=1
            # split raw output: composition = union of its two lower-port memberships
            for u in self.splits:
                m=int(self.units[u]['descendant_count']);act=self.act[u]
                slots=[]
                for port in ('left','right'):
                    fm,slot,D=self._req_slot(u,port);slots.append((fm,slot))
                for j in range(4):
                    for a0,dA,dC,dG in cache.get((m,j),[]):
                        row={self.z[p,u][j]:1,act:-a0}
                        for fm,slot in slots:
                            for t,b in enumerate(pat):
                                coeff={'A':dA,'C':dC,'G':dG,'T':0.0}[b]
                                if coeff:
                                    v=self.tmem[fm,0,slot,t];row[v]=row.get(v,0)-coeff
                        self.lp.add_le(row,0);self.shape_facets+=1

if __name__=='__main__':
    print('--- boundary shape facets ---')
    cols=[tuple('AAAACCCC'),tuple('ACGTACGT'),tuple('AACCGGTT'),tuple('AGCTAGCT'),tuple('AAAACCGT'),tuple('AACCGGTT')]
    M=BoundaryShapeModel('/mnt/data/lowerarm_units8/library.json',cols,1.1,2.0,6);v,r,t=M.solve();print('obj',v,'vars',len(M.lp.names),'eq',len(M.lp.eq),'le',len(M.lp.le),'rootfac',M.root_length_facets,'shapefac',M.shape_facets,'solve',t)

class BoundaryLongMembershipModel(BoundaryTaxonModel):
    def _build_taxon_flow(self):
        BoundaryTaxonModel._build_taxon_flow(self)
        # Joint membership of the subtree with the long-path scenario on each parent lower arm.
        self.lmem={}
        for m in range(1,self.L):
            for u,port in self.reqdef[m]:
                act=self.act[u];ell=self.ell[u,port];row={ell:-float(m)}
                mm,slot,D=self._req_slot(u,port)
                for t in range(self.L):
                    x=self.tmem[m,0,slot,t];y=self.lp.var(f'lm:{u}:{port}:{t}',0,1);self.lmem[u,port,t]=y
                    self.lp.add_le({y:1,x:-1},0)
                    self.lp.add_le({y:1,ell:-1},0)
                    self.lp.add_le({x:1,y:-1,act:-1,ell:1},0)
                    row[y]=1
                self.lp.add_eq(row,0)
    def _build_sites(self):
        Model._build_sites(self)
        # For a leaf request, raw partial is one-hot, so its long part is exactly the
        # long-scenario taxon/base membership.  This removes state-specific length cheating at count 1.
        for p,pat in enumerate(self.patterns):
            for u,port in self.reqdef[1]:
                H=self.hlong[p,u,port]
                for j in range(4):
                    row={H[j]:1}
                    for t,b in enumerate(pat):
                        if BASE_INDEX[b]==j:row[self.lmem[u,port,t]]=row.get(self.lmem[u,port,t],0)-1
                    self.lp.add_eq(row,0)

class RootLongMembershipFacetModel(BoundaryLongMembershipModel):
    def __init__(self,libpath,columns,log_ratio=1.1,long_ratio=2.0):
        # Need BoundaryLongMembership dynamic methods, but skip scalar-length root facets.
        Model.__init__(self,libpath,columns,log_ratio,long_ratio,True,True,False)
        self._add_root_long_membership_facets()
    def _add_root_long_membership_facets(self):
        self.root_longmem_facets=0
        for p,pat in enumerate(self.patterns):
            I=self.info[p];counts=I['counts'];raw=I['raw']
            for u in self.roots:
                aN=int(self.units[u]['left_flow']);bN=int(self.units[u]['right_flow']);PsA,PlA=trans_pair(aN,self.long_ratio);PsB,PlB=trans_pair(bN,self.long_ratio)
                feats=[];vals=[]
                for (mm,c1),A in raw.items():
                    if mm!=aN:continue
                    c2=tuple(counts[k]-c1[k] for k in range(4))
                    if min(c2)<0 or (bN,c2) not in raw:continue
                    B=raw[(bN,c2)]
                    for la,PA in enumerate((PsA,PlA)):
                        for lb,PB in enumerate((PsB,PlB)):
                            val=float(PI@((PA@A[:,1])*(PB@B[:,1])))/I['scale'][self.L]
                            feats.append((c1[0],c1[1],c1[2],la*c1[0],la*c1[1],la*c1[2],lb*c2[0],lb*c2[1],lb*c2[2],la,lb));vals.append(val)
                fac=affine_majorants_nd(feats,vals,60)
                mL,slotL,_=self._req_slot(u,'left');mR,slotR,_=self._req_slot(u,'right');act=self.act[u];eL=self.ell[u,'left'];eR=self.ell[u,'right'];lexpr={self.z[p,u][i]:float(PI[i]) for i in range(4)}
                for f in fac:
                    a0,*d=f;row=dict(lexpr);row[act]=row.get(act,0)-a0;row[eL]=row.get(eL,0)-d[9];row[eR]=row.get(eR,0)-d[10]
                    for t,bch in enumerate(pat):
                        bi={'A':0,'C':1,'G':2,'T':3}[bch]
                        if bi<3:
                            x=self.tmem[mL,0,slotL,t];row[x]=row.get(x,0)-d[bi]
                            yl=self.lmem[u,'left',t];yr=self.lmem[u,'right',t]
                            row[yl]=row.get(yl,0)-d[3+bi];row[yr]=row.get(yr,0)-d[6+bi]
                    self.lp.add_le(row,0);self.root_longmem_facets+=1

if __name__=='__main__':
    print('--- long membership ---')
    cols=[tuple('AAAACCCC'),tuple('ACGTACGT'),tuple('AACCGGTT'),tuple('AGCTAGCT'),tuple('AAAACCGT'),tuple('AACCGGTT')]
    for Cls in (BoundaryLongMembershipModel,RootLongMembershipFacetModel):
        M=Cls('/mnt/data/lowerarm_units8/library.json',cols,1.1,2.0);v,r,t=M.solve();print(Cls.__name__,v,'vars',len(M.lp.names),'eq',len(M.lp.eq),'le',len(M.lp.le),'fac',getattr(M,'root_longmem_facets',getattr(M,'root_length_facets',0)),'solve',t)

class CombinedLongMembershipFacetModel(BoundaryLongMembershipModel):
    def __init__(self,libpath,columns,log_ratio=1.1,long_ratio=2.0):
        BoundaryLongMembershipModel.__init__(self,libpath,columns,log_ratio,long_ratio)
        RootLongMembershipFacetModel._add_root_long_membership_facets(self)

class LongShapeModel(BoundaryLongMembershipModel):
    def __init__(self,libpath,columns,log_ratio=1.1,long_ratio=2.0,max_facets=2):
        self.max_long_shape_facets=max_facets
        BoundaryLongMembershipModel.__init__(self,libpath,columns,log_ratio,long_ratio)
        self._add_long_shape_facets()
    def _add_long_shape_facets(self):
        self.long_shape_facets=0
        for p,pat in enumerate(self.patterns):
            I=self.info[p];raw=I['raw'];cache={}
            for m in range(2,self.L):
                comps=[c for (mm,c) in raw if mm==m]
                for j in range(4):
                    vals=[float(raw[(m,c)][j,1]/I['scale'][m]) for c in comps]
                    cache[m,j]=affine_majorants(comps,vals,self.max_long_shape_facets)
            for m in range(2,self.L):
                for u,port in self.reqdef[m]:
                    act=self.act[u];ell=self.ell[u,port];R=self.req[p,u,port];H=self.hlong[p,u,port];_,slot,_=self._req_slot(u,port)
                    for j in range(4):
                        for a0,dA,dC,dG in cache.get((m,j),[]):
                            # long H <= a0*ell + d*C_long
                            row={H[j]:1,ell:-a0}
                            # short R-H <= a0*(act-ell)+d*(C-C_long)
                            row2={R[j]:1,H[j]:-1,act:-a0,ell:a0}
                            for t,b in enumerate(pat):
                                coeff={'A':dA,'C':dC,'G':dG,'T':0.0}[b]
                                if coeff:
                                    y=self.lmem[u,port,t];x=self.tmem[m,0,slot,t]
                                    row[y]=row.get(y,0)-coeff
                                    row2[x]=row2.get(x,0)-coeff;row2[y]=row2.get(y,0)+coeff
                            self.lp.add_le(row,0);self.lp.add_le(row2,0);self.long_shape_facets+=2
