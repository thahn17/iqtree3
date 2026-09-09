import importlib.util,sys,math,time,itertools
from collections import defaultdict
import numpy as np
from scipy.optimize import linprog
spec=importlib.util.spec_from_file_location('base','/mnt/data/multisite_lowerarm_taxonflow.py')
base=importlib.util.module_from_spec(spec);sys.modules['base']=base;spec.loader.exec_module(base)
BASES=base.BASES; BASE_INDEX=base.BASE_INDEX; PI=base.PI

class FlowPackingModel(base.BoundaryTaxonModel):
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)
        self.pack_cuts=0
        self._add_flow_packing()
    def _add_flow_packing(self):
        # Same taxon can occupy at most one active split subtree of a fixed descendant count.
        for m in range(2,self.L):
            for t in range(self.L):
                row={}
                for u in self.by[m]:
                    for port in ('left','right'):
                        fm,slot,D=self._req_slot(u,port)
                        v=self.tmem[fm,0,slot,t]
                        row[v]=row.get(v,0)+1
                if row:
                    self.lp.add_le(row,1.0); self.pack_cuts+=1

# fixed support directions, normalized so coefficient magnitude <=1
DIRS=[]
for i in range(4):
    d=np.zeros(4);d[i]=1;DIRS.append(d)
for i in range(4):
    for j in range(i+1,4):
        d=np.zeros(4);d[i]=d[j]=0.5;DIRS.append(d)
DIRS.append(np.ones(4)/4)
for i in range(1,4):
    d=np.zeros(4);d[0]=0.5;d[i]=-0.5;DIRS.append(d)
# unique
_u=[]
for d in DIRS:
    if not any(np.allclose(d,e) for e in _u):_u.append(d)
DIRS=_u

def affine_majorants_features(X,y,max_facets=10):
    X=np.asarray(X,float); y=np.asarray(y,float)
    A=np.column_stack([np.ones(len(X)),X])
    if len(A)==0:return []
    if len(A)<=max_facets: idx=list(range(len(A)))
    else:
        # include extrema and evenly spaced after sorting by y
        order=np.argsort(y)
        idx=set([int(order[0]),int(order[-1])])
        for z in np.linspace(0,len(order)-1,max_facets,dtype=int):idx.add(int(order[z]))
        idx=sorted(idx)
    fs=[]
    for ii in idx:
        res=linprog(A[ii],A_ub=-A,b_ub=-y,bounds=[(None,None)]*A.shape[1],method='highs')
        if res.success and np.all(A@res.x+1e-8>=y):
            key=tuple(np.round(res.x,10))
            if key not in fs:fs.append(key)
    return fs

def box_support(box,d,scale=1.0):
    lo=box[:,0]/scale; hi=box[:,1]/scale
    return float(np.sum(np.where(d>=0,d*hi,d*lo)))

class MultiStateShapeModel(FlowPackingModel):
    def __init__(self,*args,max_facets=6,dirs=None,**kwargs):
        self.ms_max_facets=max_facets; self.ms_dirs=DIRS if dirs is None else dirs
        super().__init__(*args,**kwargs)
        self.ms_cuts=0
        self._add_multistate_shape()
    def _add_multistate_shape(self):
        for p,pat in enumerate(self.patterns):
            I=self.info[p]; raw=I['raw']
            cache={}
            for m in range(2,self.L+1):
                comps=[c for (mm,c) in raw if mm==m]
                X=np.array([[c[0],c[1],c[2]] for c in comps],float)
                for di,d in enumerate(self.ms_dirs):
                    y=[box_support(raw[(m,c)],d,I['scale'][m]) for c in comps]
                    cache[m,di]=affine_majorants_features(X,y,self.ms_max_facets)
            # raw request vectors tied to exact shared membership composition
            for m in range(2,self.L):
                for u,port in self.reqdef[m]:
                    _,slot,_=self._req_slot(u,port); act=self.act[u]; R=self.req[p,u,port]
                    for di,d in enumerate(self.ms_dirs):
                        for f in cache.get((m,di),[]):
                            a0,dA,dC,dG=f; row={act:-a0}
                            for j in range(4):
                                if d[j]: row[R[j]]=row.get(R[j],0)+float(d[j])
                            for t,b in enumerate(pat):
                                coeff={'A':dA,'C':dC,'G':dG,'T':0.0}[b]
                                if coeff:
                                    x=self.tmem[m,0,slot,t];row[x]=row.get(x,0)-coeff
                            self.lp.add_le(row,0);self.ms_cuts+=1
            # split raw output tied to union membership
            for u in self.splits:
                m=int(self.units[u]['descendant_count']); act=self.act[u]; Z=self.z[p,u]
                slots=[self._req_slot(u,port)[:2] for port in ('left','right')]
                for di,d in enumerate(self.ms_dirs):
                    for f in cache.get((m,di),[]):
                        a0,dA,dC,dG=f;row={act:-a0}
                        for j in range(4):
                            if d[j]:row[Z[j]]=row.get(Z[j],0)+float(d[j])
                        for fm,slot in slots:
                            for t,b in enumerate(pat):
                                coeff={'A':dA,'C':dC,'G':dG,'T':0.0}[b]
                                if coeff:
                                    x=self.tmem[fm,0,slot,t];row[x]=row.get(x,0)-coeff
                        self.lp.add_le(row,0);self.ms_cuts+=1

# Joint split facets: state support of parent z conditioned jointly on left/right composition and length endpoints.
class JointSplitSupportModel(FlowPackingModel):
    def __init__(self,*args,max_facets=6,dirs=None,**kwargs):
        self.js_max_facets=max_facets; self.js_dirs=DIRS if dirs is None else dirs
        super().__init__(*args,**kwargs)
        self.js_cuts=0
        self._add_joint_split_support()
    def _add_joint_split_support(self):
        for p,pat in enumerate(self.patterns):
            I=self.info[p];raw=I['raw']
            # cache by (a,b,direction) features [cL A,C,G, cR A,C,G, lenL,lenR]
            cache={}
            splittypes=sorted(set((int(self.units[u]['left_flow']),int(self.units[u]['right_flow'])) for u in self.splits))
            for a,b in splittypes:
                compsA=[c for (mm,c) in raw if mm==a]; compsB=[c for (mm,c) in raw if mm==b]
                if not compsA or not compsB:continue
                PsA,PlA=base.trans_pair(a,self.long_ratio);PsB,PlB=base.trans_pair(b,self.long_ratio)
                for di,d in enumerate(self.js_dirs):
                    X=[];y=[]
                    for cA in compsA:
                        BA=raw[(a,cA)]
                        for cB in compsB:
                            # global site counts bound union
                            if any(cA[k]+cB[k]>I['counts'][k] for k in range(4)):continue
                            BB=raw[(b,cB)]
                            for la,PA in enumerate((PsA,PlA)):
                                Alo=PA@BA[:,0];Ahi=PA@BA[:,1]
                                for lb,PB in enumerate((PsB,PlB)):
                                    Blo=PB@BB[:,0];Bhi=PB@BB[:,1]
                                    lo=Alo*Blo/I['scale'][a]/I['scale'][b]
                                    hi=Ahi*Bhi/I['scale'][a]/I['scale'][b]
                                    # Model z scaling k=scale[a]*scale[b]/scale[m], so physical product/scale[m]
                                    # above /scale[a]/scale[b], multiply k
                                    m=a+b;k=I['scale'][a]*I['scale'][b]/I['scale'][m]
                                    lo*=k;hi*=k
                                    sup=float(np.sum(np.where(d>=0,d*hi,d*lo)))
                                    X.append([cA[0],cA[1],cA[2],cB[0],cB[1],cB[2],la,lb]);y.append(sup)
                    if X:cache[a,b,di]=affine_majorants_features(X,y,self.js_max_facets)
            for u in self.splits:
                a=int(self.units[u]['left_flow']);b=int(self.units[u]['right_flow']);act=self.act[u];eL=self.ell[u,'left'];eR=self.ell[u,'right'];Z=self.z[p,u]
                mL,sL,_=self._req_slot(u,'left');mR,sR,_=self._req_slot(u,'right')
                for di,d in enumerate(self.js_dirs):
                    for f in cache.get((a,b,di),[]):
                        a0,*coef=f;row={act:-a0,eL:-coef[6],eR:-coef[7]}
                        for j in range(4):
                            if d[j]:row[Z[j]]=row.get(Z[j],0)+float(d[j])
                        for t,ch in enumerate(pat):
                            bi={'A':0,'C':1,'G':2,'T':3}[ch]
                            if bi<3:
                                xL=self.tmem[mL,0,sL,t];xR=self.tmem[mR,0,sR,t]
                                row[xL]=row.get(xL,0)-coef[bi]; row[xR]=row.get(xR,0)-coef[3+bi]
                        self.lp.add_le(row,0);self.js_cuts+=1

if __name__=='__main__':
    cols8=[tuple('AAAACCCC'),tuple('ACGTACGT'),tuple('AACCGGTT'),tuple('AGCTAGCT'),tuple('AAAACCGT'),tuple('AACCGGTT')]
    lib8='/mnt/data/test_match_units8/library.json'
    tests=[('baseline',base.BoundaryTaxonModel,{}),('packing',FlowPackingModel,{}),('multistate',MultiStateShapeModel,{'max_facets':4}),('jointsplit',JointSplitSupportModel,{'max_facets':4})]
    for name,Cls,kw in tests:
        st=time.time();M=Cls(lib8,cols8,1.1,2.0,**kw);build=time.time()-st;v,r,sol=M.solve()
        print(name,'obj',v,'vars',len(M.lp.names),'eq',len(M.lp.eq),'le',len(M.lp.le),'build',build,'solve',sol,'pack',getattr(M,'pack_cuts',0),'ms',getattr(M,'ms_cuts',0),'js',getattr(M,'js_cuts',0))

class JointSiteRootModel(FlowPackingModel):
    def __init__(self,*args,max_facets=4,pairwise=True,**kwargs):
        self.jr_max_facets=max_facets; self.jr_pairwise=pairwise
        super().__init__(*args,**kwargs)
        self.jr_cuts=0
        self._add_joint_site_root()
    def _site_root_endpoint_val(self,p,u,subset,la,lb):
        I=self.info[p];pat=self.patterns[p];a=int(self.units[u]['left_flow']);b=int(self.units[u]['right_flow'])
        c1=[0,0,0,0]
        for t in subset:c1[BASE_INDEX[pat[t]]]+=1
        c1=tuple(c1);c2=tuple(I['counts'][k]-c1[k] for k in range(4))
        if (a,c1) not in I['raw'] or (b,c2) not in I['raw']:return None
        A=I['raw'][(a,c1)][:,1];B=I['raw'][(b,c2)][:,1]
        PsA,PlA=base.trans_pair(a,self.long_ratio);PsB,PlB=base.trans_pair(b,self.long_ratio)
        PA=(PsA,PlA)[la];PB=(PsB,PlB)[lb]
        return float(PI@((PA@A)*(PB@B)))/I['scale'][self.L]
    def _add_joint_site_root(self):
        P=len(self.patterns)
        pairs=list(itertools.combinations(range(P),2)) if self.jr_pairwise else []
        for u in self.roots:
            a=int(self.units[u]['left_flow']); act=self.act[u];eL=self.ell[u,'left'];eR=self.ell[u,'right']
            mL,sL,_=self._req_slot(u,'left')
            subsets=list(itertools.combinations(range(self.L),a))
            # precompute q endpoint vector per point
            X=[]; Qv=[]
            for S in subsets:
                ind=np.zeros(self.L);ind[list(S)]=1
                for la in (0,1):
                    for lb in (0,1):
                        vals=[self._site_root_endpoint_val(p,u,S,la,lb) for p in range(P)]
                        if all(v is not None for v in vals):
                            X.append(list(ind)+[la,lb]);Qv.append(vals)
            X=np.asarray(X,float);Qv=np.asarray(Qv,float)
            for p,q in pairs:
                # normalized equal-weight support so neither pattern dominates numerically
                up=max(Qv[:,p].max(),1e-12); uq=max(Qv[:,q].max(),1e-12)
                weights=(1.0/up,1.0/uq)
                y=Qv[:,p]*weights[0]+Qv[:,q]*weights[1]
                fac=affine_majorants_features(X,y,self.jr_max_facets)
                for f in fac:
                    a0,*coef=f;row={act:-a0,eL:-coef[self.L],eR:-coef[self.L+1]}
                    for i in range(4):
                        row[self.z[p,u][i]]=row.get(self.z[p,u][i],0)+weights[0]*float(PI[i])
                        row[self.z[q,u][i]]=row.get(self.z[q,u][i],0)+weights[1]*float(PI[i])
                    for t in range(self.L):
                        c=coef[t]
                        if c:
                            x=self.tmem[mL,0,sL,t];row[x]=row.get(x,0)-c
                    self.lp.add_le(row,0);self.jr_cuts+=1

class JointAllSiteRootModel(FlowPackingModel):
    def __init__(self,*args,max_facets=8,**kwargs):
        self.ja_max_facets=max_facets
        super().__init__(*args,**kwargs);self.ja_cuts=0;self._add_joint_all()
    def _site_val(self,p,u,S,la,lb):
        return JointSiteRootModel._site_root_endpoint_val(self,p,u,S,la,lb)
    def _add_joint_all(self):
        P=len(self.patterns)
        for u in self.roots:
            a=int(self.units[u]['left_flow']);act=self.act[u];eL=self.ell[u,'left'];eR=self.ell[u,'right'];mL,sL,_=self._req_slot(u,'left')
            X=[];Qv=[]
            for S in itertools.combinations(range(self.L),a):
                ind=np.zeros(self.L);ind[list(S)]=1
                for la in (0,1):
                    for lb in (0,1):
                        vals=[self._site_val(p,u,S,la,lb) for p in range(P)]
                        if all(v is not None for v in vals):X.append(list(ind)+[la,lb]);Qv.append(vals)
            X=np.asarray(X,float);Qv=np.asarray(Qv,float)
            # fixed normalized directions: all sites together, and two alternating groups
            umax=np.maximum(Qv.max(axis=0),1e-12)
            dirs=[1/umax]
            if P>=3:
                d=np.zeros(P);d[::2]=1/umax[::2];dirs.append(d)
                d=np.zeros(P);d[1::2]=1/umax[1::2];dirs.append(d)
            for w in dirs:
                y=Qv@w;fac=affine_majorants_features(X,y,self.ja_max_facets)
                for f in fac:
                    a0,*coef=f;row={act:-a0,eL:-coef[self.L],eR:-coef[self.L+1]}
                    for p in range(P):
                        if not w[p]:continue
                        for i in range(4):row[self.z[p,u][i]]=row.get(self.z[p,u][i],0)+w[p]*float(PI[i])
                    for t in range(self.L):
                        if coef[t]:
                            x=self.tmem[mL,0,sL,t];row[x]=row.get(x,0)-coef[t]
                    self.lp.add_le(row,0);self.ja_cuts+=1

def packing_dp_values(global_counts,k,comp_values,qmax):
    g=tuple(global_counts)
    # dp[h][usedtuple]=max value
    dp=[{(0,0,0,0):0.0}]+[{} for _ in range(qmax)]
    comps=list(comp_values)
    for h in range(1,qmax+1):
        cur={}
        for used,val0 in dp[h-1].items():
            for c in comps:
                nu=tuple(used[i]+c[i] for i in range(4))
                if any(nu[i]>g[i] for i in range(4)):continue
                v=val0+comp_values[c]
                if v>cur.get(nu,-1e300):cur[nu]=v
        dp[h]=cur
    B=[0.0]
    for h in range(1,qmax+1):B.append(max(dp[h].values()) if dp[h] else -1e300)
    return B

def upper_lines_1d(xs,ys):
    X=np.column_stack([np.ones(len(xs)),np.asarray(xs,float)]);y=np.asarray(ys,float)
    fs=[]
    for ii in range(len(xs)):
        res=linprog(X[ii],A_ub=-X,b_ub=-y,bounds=[(None,None),(None,None)],method='highs')
        if res.success and np.all(X@res.x+1e-9>=y):
            f=tuple(np.round(res.x,12))
            if f not in fs:fs.append(f)
    return fs

class AntichainLikelihoodPackingModel(base.BoundaryTaxonModel):
    def __init__(self,*args,dirs=None,**kwargs):
        self.ap_dirs=dirs if dirs is not None else [np.eye(4)[i] for i in range(4)]+[np.ones(4)/4]
        super().__init__(*args,**kwargs);self.ap_cuts=0;self._add_antichain_likelihood_packing()
    def _add_antichain_likelihood_packing(self):
        for p,pat in enumerate(self.patterns):
            I=self.info[p];g=I['counts'];raw=I['raw']
            for k in range(2,self.L):
                units=self.by[k]
                if not units:continue
                q=min(len(units),self.L//k)
                if q<=0:continue
                comps=[c for (mm,c) in raw if mm==k]
                for d in self.ap_dirs:
                    cv={c:box_support(raw[(k,c)],d,I['scale'][k]) for c in comps}
                    B=packing_dp_values(g,k,cv,q)
                    finite=[h for h,b in enumerate(B) if b>-1e200]
                    xs=finite;ys=[B[h] for h in finite]
                    for a0,beta in upper_lines_1d(xs,ys):
                        row={}
                        for u in units:
                            for j in range(4):
                                if d[j]:row[self.z[p,u][j]]=row.get(self.z[p,u][j],0)+float(d[j])
                            row[self.act[u]]=row.get(self.act[u],0)-beta
                        self.lp.add_le(row,a0);self.ap_cuts+=1

class RichJointRootModel(FlowPackingModel):
    def __init__(self,*args,max_facets=4,mode='allsubsets',**kwargs):
        self.rj_max=max_facets;self.rj_mode=mode
        super().__init__(*args,**kwargs);self.rj_cuts=0;self._add_rich()
    def _site_val(self,p,u,S,la,lb):return JointSiteRootModel._site_root_endpoint_val(self,p,u,S,la,lb)
    def _add_rich(self):
        P=len(self.patterns)
        for u in self.roots:
            a=int(self.units[u]['left_flow']);act=self.act[u];eL=self.ell[u,'left'];eR=self.ell[u,'right'];mL,sL,_=self._req_slot(u,'left')
            X=[];Qv=[]
            for S in itertools.combinations(range(self.L),a):
                ind=np.zeros(self.L);ind[list(S)]=1
                for la in (0,1):
                    for lb in (0,1):
                        vals=[self._site_val(p,u,S,la,lb) for p in range(P)]
                        if all(v is not None for v in vals):X.append(list(ind)+[la,lb]);Qv.append(vals)
            X=np.asarray(X,float);Qv=np.asarray(Qv,float);um=np.maximum(Qv.max(axis=0),1e-12)
            dirs=[]
            # always include objective-like all-site direction
            dirs.append(self.weights/um)
            dirs.append(np.ones(P)/um)
            # pair directions
            for p,q in itertools.combinations(range(P),2):
                d=np.zeros(P);d[p]=1/um[p];d[q]=1/um[q];dirs.append(d)
            if self.rj_mode=='allsubsets' and P<=8:
                for mask in range(1,1<<P):
                    d=np.zeros(P)
                    for p in range(P):
                        if mask>>p&1:d[p]=1/um[p]
                    dirs.append(d)
            # dedup normalized patterns
            uniq=[];seen=set()
            for d in dirs:
                key=tuple(np.round(d/np.max(d),8)) if np.max(d)>0 else tuple(d)
                if key not in seen:seen.add(key);uniq.append(d)
            for w in uniq:
                y=Qv@w
                for f in affine_majorants_features(X,y,self.rj_max):
                    a0,*coef=f;row={act:-a0,eL:-coef[self.L],eR:-coef[self.L+1]}
                    for p in range(P):
                        if w[p]:
                            for i in range(4):row[self.z[p,u][i]]=row.get(self.z[p,u][i],0)+w[p]*float(PI[i])
                    for t in range(self.L):
                        if coef[t]:
                            x=self.tmem[mL,0,sL,t];row[x]=row.get(x,0)-coef[t]
                    self.lp.add_le(row,0);self.rj_cuts+=1

class JointSiteInternalModel(JointAllSiteRootModel):
    def __init__(self,*args,max_internal_facets=2,**kwargs):
        self.ji_max=max_internal_facets
        super().__init__(*args,**kwargs);self.ji_cuts=0;self._add_joint_internal()
    def _parent_state_upper(self,p,a,b,SL,SR,la,lb,i):
        I=self.info[p];pat=self.patterns[p]
        cA=[0,0,0,0];cB=[0,0,0,0]
        for t in SL:cA[BASE_INDEX[pat[t]]]+=1
        for t in SR:cB[BASE_INDEX[pat[t]]]+=1
        cA=tuple(cA);cB=tuple(cB)
        if (a,cA) not in I['raw'] or (b,cB) not in I['raw']:return None
        A=I['raw'][(a,cA)][:,1];B=I['raw'][(b,cB)][:,1]
        PsA,PlA=base.trans_pair(a,self.long_ratio);PsB,PlB=base.trans_pair(b,self.long_ratio)
        PA=(PsA,PlA)[la];PB=(PsB,PlB)[lb]
        return float((PA@A)[i]*(PB@B)[i]/I['scale'][a+b])
    def _add_joint_internal(self):
        P=len(self.patterns);cache={}
        splittypes=sorted(set((int(self.units[u]['left_flow']),int(self.units[u]['right_flow'])) for u in self.splits))
        for a,b in splittypes:
            # enumerate disjoint left/right taxa; unused taxa allowed for m<L
            X=[];Vals=[]
            taxa=range(self.L)
            for SL in itertools.combinations(taxa,a):
                rem=[t for t in taxa if t not in SL]
                for SR in itertools.combinations(rem,b):
                    indL=np.zeros(self.L);indR=np.zeros(self.L);indL[list(SL)]=1;indR[list(SR)]=1
                    for la in (0,1):
                        for lb in (0,1):
                            vv=np.empty((P,4));ok=True
                            for p in range(P):
                                for i in range(4):
                                    z=self._parent_state_upper(p,a,b,SL,SR,la,lb,i)
                                    if z is None:ok=False;break
                                    vv[p,i]=z
                                if not ok:break
                            if ok:X.append(list(indL)+list(indR)+[la,lb]);Vals.append(vv)
            if not X:continue
            X=np.asarray(X,float);Vals=np.asarray(Vals,float)
            for i in range(4):
                umax=np.maximum(Vals[:,:,i].max(axis=0),1e-12)
                w=self.weights/umax
                y=Vals[:,:,i]@w
                cache[a,b,i]=(w,affine_majorants_features(X,y,self.ji_max))
        for u in self.splits:
            a=int(self.units[u]['left_flow']);b=int(self.units[u]['right_flow']);act=self.act[u];eL=self.ell[u,'left'];eR=self.ell[u,'right']
            mL,sL,_=self._req_slot(u,'left');mR,sR,_=self._req_slot(u,'right')
            for i in range(4):
                if (a,b,i) not in cache:continue
                w,fac=cache[a,b,i]
                for f in fac:
                    a0,*coef=f;row={act:-a0,eL:-coef[2*self.L],eR:-coef[2*self.L+1]}
                    for p in range(P):row[self.z[p,u][i]]=row.get(self.z[p,u][i],0)+w[p]
                    for t in range(self.L):
                        if coef[t]:
                            x=self.tmem[mL,0,sL,t];row[x]=row.get(x,0)-coef[t]
                        if coef[self.L+t]:
                            x=self.tmem[mR,0,sR,t];row[x]=row.get(x,0)-coef[self.L+t]
                    self.lp.add_le(row,0);self.ji_cuts+=1

def comp_moment_features_from_subset(pat,S):
    c=np.zeros(4,dtype=int)
    for t in S:c[BASE_INDEX[pat[t]]]+=1
    feats=[c[0],c[1],c[2]]
    # 10 second raw moments upper triangle b<=c: n_b*n_c (including squares)
    for b in range(4):
        for c2 in range(b,4):feats.append(c[b]*c[c2])
    return feats

class RootPairMomentModel(FlowPackingModel):
    def __init__(self,*args,max_facets=5,joint=True,**kwargs):
        self.rpm_max=max_facets;self.rpm_joint=joint
        super().__init__(*args,**kwargs)
        self.rpair={};self._add_root_pair_moments();self.rpm_cuts=0;self._add_pairmoment_root_facets()
    def _add_root_pair_moments(self):
        for u in self.roots:
            a=int(self.units[u]['left_flow']);A=self.act[u];m,s,_=self._req_slot(u,'left')
            if a<=1:continue
            for t in range(self.L):
                for v in range(t+1,self.L):
                    y=self.lp.var(f'rpair:{u}:{t}:{v}',0,1);self.rpair[u,t,v]=y
                    xt=self.tmem[m,0,s,t];xv=self.tmem[m,0,s,v]
                    self.lp.add_le({y:1,xt:-1},0);self.lp.add_le({y:1,xv:-1},0)
                    self.lp.add_ge({y:1,xt:-1,xv:-1,A:1},0)
            for t in range(self.L):
                row={}
                for v in range(self.L):
                    if v==t:continue
                    key=(u,min(t,v),max(t,v));row[self.rpair[key]]=row.get(self.rpair[key],0)+1
                xt=self.tmem[m,0,s,t];row[xt]=row.get(xt,0)-(a-1);self.lp.add_eq(row,0)
    def _moment_expr_coeffs(self,u,pat):
        # Return list of dicts for 3 first moments + 10 second raw moments.
        m,s,_=self._req_slot(u,'left');expr=[]
        # first A,C,G
        for b in range(3):
            d={}
            for t,ch in enumerate(pat):
                if BASE_INDEX[ch]==b:d[self.tmem[m,0,s,t]]=1
            expr.append(d)
        # second moments
        for b in range(4):
            for c in range(b,4):
                d={}
                if b==c:
                    # n_b^2 = n_b + 2 sum same-base pairs
                    for t,ch in enumerate(pat):
                        if BASE_INDEX[ch]==b:d[self.tmem[m,0,s,t]]=d.get(self.tmem[m,0,s,t],0)+1
                    for t in range(self.L):
                        if BASE_INDEX[pat[t]]!=b:continue
                        for v in range(t+1,self.L):
                            if BASE_INDEX[pat[v]]==b and (u,t,v) in self.rpair:d[self.rpair[u,t,v]]=d.get(self.rpair[u,t,v],0)+2
                else:
                    for t in range(self.L):
                        bt=BASE_INDEX[pat[t]]
                        if bt not in (b,c):continue
                        for v in range(t+1,self.L):
                            bv=BASE_INDEX[pat[v]]
                            if {bt,bv}=={b,c} and (u,t,v) in self.rpair:d[self.rpair[u,t,v]]=d.get(self.rpair[u,t,v],0)+1
                expr.append(d)
        return expr
    def _site_root_val(self,p,u,S,la,lb):return JointSiteRootModel._site_root_endpoint_val(self,p,u,S,la,lb)
    def _add_pairmoment_root_facets(self):
        P=len(self.patterns)
        for u in self.roots:
            a=int(self.units[u]['left_flow']);A=self.act[u];eL=self.ell[u,'left'];eR=self.ell[u,'right']
            X=[];Qv=[]
            for S in itertools.combinations(range(self.L),a):
                # shared feature vector will be site-specific, so for joint facets concatenate? We instead
                # fit directly in shared taxon-pair y via site-dependent moment expressions only for single site;
                # for joint, use all sites' moment features concatenated, still expressions of same y.
                feats_all=[]
                for p,pat in enumerate(self.patterns):feats_all.extend(comp_moment_features_from_subset(pat,S))
                for la in (0,1):
                    for lb in (0,1):
                        vals=[self._site_root_val(p,u,S,la,lb) for p in range(P)]
                        if all(v is not None for v in vals):X.append(feats_all+[la,lb]);Qv.append(vals)
            X=np.asarray(X,float);Qv=np.asarray(Qv,float)
            # expressions corresponding to concatenated features
            exprs=[]
            for p,pat in enumerate(self.patterns):exprs.extend(self._moment_expr_coeffs(u,pat))
            dirs=[]
            umax=np.maximum(Qv.max(axis=0),1e-12)
            if self.rpm_joint:
                dirs=[self.weights/umax,np.ones(P)/umax]
            else:
                for p in range(P):
                    d=np.zeros(P);d[p]=1;dirs.append(d)
            for w in dirs:
                y=Qv@w
                for f in affine_majorants_features(X,y,self.rpm_max):
                    a0,*coef=f;row={A:-a0,eL:-coef[len(exprs)],eR:-coef[len(exprs)+1]}
                    for p in range(P):
                        if w[p]:
                            for i in range(4):row[self.z[p,u][i]]=row.get(self.z[p,u][i],0)+w[p]*float(PI[i])
                    for h,ex in enumerate(exprs):
                        c=coef[h]
                        if c:
                            for v,cc in ex.items():row[v]=row.get(v,0)-c*cc
                    self.lp.add_le(row,0);self.rpm_cuts+=1

class RootPairMomentCompositionModel(FlowPackingModel):
    def __init__(self,*args,max_facets=8,**kwargs):
        self.rpc_max=max_facets
        super().__init__(*args,**kwargs);self.rpair={};RootPairMomentModel._add_root_pair_moments(self);self.rpc_cuts=0;self._add_comp_quad_facets()
    def _moment_expr_coeffs(self,u,pat):return RootPairMomentModel._moment_expr_coeffs(self,u,pat)
    def _add_comp_quad_facets(self):
        for p,pat in enumerate(self.patterns):
            I=self.info[p];raw=I['raw'];counts=I['counts']
            for u in self.roots:
                a=int(self.units[u]['left_flow']);b=self.L-a;A=self.act[u];eL=self.ell[u,'left'];eR=self.ell[u,'right']
                PsA,PlA=base.trans_pair(a,self.long_ratio);PsB,PlB=base.trans_pair(b,self.long_ratio)
                X=[];y=[]
                for (mm,c1),BA in raw.items():
                    if mm!=a:continue
                    c2=tuple(counts[k]-c1[k] for k in range(4))
                    if min(c2)<0 or (b,c2) not in raw:continue
                    feats=[c1[0],c1[1],c1[2]]
                    for bb in range(4):
                        for cc in range(bb,4):feats.append(c1[bb]*c1[cc])
                    for la,PA in enumerate((PsA,PlA)):
                        for lb,PB in enumerate((PsB,PlB)):
                            val=float(PI@((PA@BA[:,1])*(PB@raw[(b,c2)][:,1])))/I['scale'][self.L]
                            X.append(feats+[la,lb]);y.append(val)
                exprs=self._moment_expr_coeffs(u,pat);lexpr={self.z[p,u][i]:float(PI[i]) for i in range(4)}
                for f in affine_majorants_features(X,y,self.rpc_max):
                    a0,*coef=f;row=dict(lexpr);row[A]=row.get(A,0)-a0;row[eL]=row.get(eL,0)-coef[len(exprs)];row[eR]=row.get(eR,0)-coef[len(exprs)+1]
                    for h,ex in enumerate(exprs):
                        c=coef[h]
                        if c:
                            for v,cc in ex.items():row[v]=row.get(v,0)-c*cc
                    self.lp.add_le(row,0);self.rpc_cuts+=1

class RootPairMomentCompositionPruned(RootPairMomentCompositionModel):
    def __init__(self,*args,pair_bounds='upper',**kwargs):
        self.pair_bounds=pair_bounds
        super().__init__(*args,**kwargs)
    def _add_root_pair_moments(self):
        for u in self.roots:
            a=int(self.units[u]['left_flow']);A=self.act[u];m,s,_=self._req_slot(u,'left')
            if a<=1:continue
            for t in range(self.L):
                for v in range(t+1,self.L):
                    y=self.lp.var(f'rpair:{u}:{t}:{v}',0,1);self.rpair[u,t,v]=y
                    xt=self.tmem[m,0,s,t];xv=self.tmem[m,0,s,v]
                    if self.pair_bounds in ('upper','full'):
                        self.lp.add_le({y:1,xt:-1},0);self.lp.add_le({y:1,xv:-1},0)
                    if self.pair_bounds=='full':self.lp.add_ge({y:1,xt:-1,xv:-1,A:1},0)
            for t in range(self.L):
                row={}
                for v in range(self.L):
                    if v==t:continue
                    row[self.rpair[u,min(t,v),max(t,v)]]=row.get(self.rpair[u,min(t,v),max(t,v)],0)+1
                xt=self.tmem[m,0,s,t];row[xt]=row.get(xt,0)-(a-1);self.lp.add_eq(row,0)
class RootPairMomentCompositionPruned2(FlowPackingModel):
    def __init__(self,*args,max_facets=3,pair_bounds='upper',**kwargs):
        self.rpc_max=max_facets;self.pair_bounds=pair_bounds
        super().__init__(*args,**kwargs);self.rpair={};self._add_pairs();self.rpc_cuts=0;RootPairMomentCompositionModel._add_comp_quad_facets(self)
    def _add_pairs(self):
        for u in self.roots:
            a=int(self.units[u]['left_flow']);A=self.act[u];m,s,_=self._req_slot(u,'left')
            if a<=1:continue
            for t in range(self.L):
                for v in range(t+1,self.L):
                    y=self.lp.var(f'rpair:{u}:{t}:{v}',0,1);self.rpair[u,t,v]=y
                    xt=self.tmem[m,0,s,t];xv=self.tmem[m,0,s,v]
                    if self.pair_bounds in ('upper','full'):
                        self.lp.add_le({y:1,xt:-1},0);self.lp.add_le({y:1,xv:-1},0)
                    if self.pair_bounds=='full':self.lp.add_ge({y:1,xt:-1,xv:-1,A:1},0)
            for t in range(self.L):
                row={}
                for v in range(self.L):
                    if v==t:continue
                    row[self.rpair[u,min(t,v),max(t,v)]]=row.get(self.rpair[u,min(t,v),max(t,v)],0)+1
                xt=self.tmem[m,0,s,t];row[xt]=row.get(xt,0)-(a-1);self.lp.add_eq(row,0)
    def _moment_expr_coeffs(self,u,pat):return RootPairMomentModel._moment_expr_coeffs(self,u,pat)
class SelectiveRootPairMoment(FlowPackingModel):
    def __init__(self,*args,root_left_sizes=None,max_facets=3,**kwargs):
        self.sel=set(root_left_sizes or []);self.rpc_max=max_facets
        super().__init__(*args,**kwargs);self.rpair={};self._add_sel_pairs();self.rpc_cuts=0;self._add_sel_facets()
    def _add_sel_pairs(self):
        for r in self.roots:
            a=int(self.units[r]['left_flow'])
            if a<=1 or a not in self.sel:continue
            A=self.act[r];m,s,_=self._req_slot(r,'left')
            for t in range(self.L):
                for v in range(t+1,self.L):
                    y=self.lp.var(f'rpair:{r}:{t}:{v}',0,1);self.rpair[r,t,v]=y
                    xt=self.tmem[m,0,s,t];xv=self.tmem[m,0,s,v]
                    self.lp.add_le({y:1,xt:-1},0);self.lp.add_le({y:1,xv:-1},0);self.lp.add_ge({y:1,xt:-1,xv:-1,A:1},0)
            for t in range(self.L):
                row={}
                for v in range(self.L):
                    if v==t:continue
                    row[self.rpair[r,min(t,v),max(t,v)]]=row.get(self.rpair[r,min(t,v),max(t,v)],0)+1
                row[self.tmem[m,0,s,t]]=-(a-1);self.lp.add_eq(row,0)
    def _moment_expr(self,r,pat):return RootPairMomentModel._moment_expr_coeffs(self,r,pat)
    def _add_sel_facets(self):
        for p,pat in enumerate(self.patterns):
            I=self.info[p];raw=I['raw'];counts=I['counts']
            for r in self.roots:
                a=int(self.units[r]['left_flow'])
                if a not in self.sel:continue
                b=self.L-a;A=self.act[r];eL=self.ell[r,'left'];eR=self.ell[r,'right'];PsA,PlA=base.trans_pair(a,self.long_ratio);PsB,PlB=base.trans_pair(b,self.long_ratio)
                X=[];y=[]
                for (mm,c1),BA in raw.items():
                    if mm!=a:continue
                    c2=tuple(counts[k]-c1[k] for k in range(4))
                    if min(c2)<0 or (b,c2) not in raw:continue
                    feat=[c1[0],c1[1],c1[2]]+[c1[bb]*c1[cc] for bb in range(4) for cc in range(bb,4)]
                    for la,PA in enumerate((PsA,PlA)):
                        for lb,PB in enumerate((PsB,PlB)):
                            X.append(feat+[la,lb]);y.append(float(PI@((PA@BA[:,1])*(PB@raw[(b,c2)][:,1])))/I['scale'][self.L])
                ex=self._moment_expr(r,pat) if a>1 else []
                # a=1 needs no pair moment; skip selective moment facet here
                if a<=1:continue
                qrow={self.z[p,r][i]:float(PI[i]) for i in range(4)}
                for f in affine_majorants_features(X,y,self.rpc_max):
                    a0,*coef=f;row=dict(qrow);row[A]=row.get(A,0)-a0;row[eL]=row.get(eL,0)-coef[len(ex)];row[eR]=row.get(eR,0)-coef[len(ex)+1]
                    for h,d in enumerate(ex):
                        if coef[h]:
                            for v,c in d.items():row[v]=row.get(v,0)-coef[h]*c
                    self.lp.add_le(row,0);self.rpc_cuts+=1
