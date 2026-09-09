#!/usr/bin/env python3
"""One-pass pruned multisite relaxation.

Research prototype combining only pruning rules that are structural / build-time certified:
- no per-column internal pruning network by default; exact root composition+root-arm-length upper hull replaces it;
- shared named-taxon boundary memberships and fixed-flow packing remain;
- shared root taxon-pair moments for second-order cross-site matching;
- pair bounds are omitted for root side size 2 (provably redundant with nonnegativity+degree equations);
- site/root second-order facets are omitted only when the feasible composition set is fully enumerated and
  forms an affine simplex, so every quadratic feature is exactly affine in first-order counts;
- an entire root-size pair-moment layer is omitted only if that certificate holds for every site pattern;
- explicit q_s variable is eliminated: root-contribution sum is equated directly to log-chord interpolation;
- optional equal-child-flow orientation row is added only when both arms have identical allowed transition families.

All decisions are made before the LP solve and consider the full candidate-tree library, not a selected tree.
"""
from __future__ import annotations
import importlib.util,sys,math
from collections import defaultdict
import numpy as np

spec=importlib.util.spec_from_file_location('t','/mnt/data/test_rootonly_pruning.py')
t=importlib.util.module_from_spec(spec);sys.modules['t']=t;spec.loader.exec_module(t)
spec2=importlib.util.spec_from_file_location('h','/mnt/data/test_exact_upper_hull_facets.py')
h=importlib.util.module_from_spec(spec2);sys.modules['h']=h;spec2.loader.exec_module(h)
base=t.base;PI=base.PI;BASE_INDEX=base.BASE_INDEX


def feasible_compositions_limited(g,a,limit=5):
    """Return up to limit feasible 4-base compositions, plus whether enumeration completed."""
    g=tuple(map(int,g));out=[]
    loA=max(0,a-g[1]-g[2]-g[3]);hiA=min(g[0],a)
    for A in range(loA,hiA+1):
        loC=max(0,a-A-g[2]-g[3]);hiC=min(g[1],a-A)
        for C in range(loC,hiC+1):
            loG=max(0,a-A-C-g[3]);hiG=min(g[2],a-A-C)
            for G in range(loG,hiG+1):
                T=a-A-C-G
                if 0<=T<=g[3]:
                    out.append((A,C,G,T))
                    if len(out)>=limit:return out,False
    return out,True


def certified_quadratic_redundant(g,a):
    """Cheap sufficient *exact* certificate.

    If all feasible compositions are fully enumerated (<=4 here) and are affinely independent
    up to their count, then any scalar function on those composition points -- hence every quadratic
    composition feature -- has an affine representation in (1,nA,nC,nG).
    False means 'do not prune'; it is not a claim that quadratic moments are necessary.
    """
    C,complete=feasible_compositions_limited(g,a,5)
    if not complete:return False
    if len(C)<=1:return True
    F=np.asarray([[1,c[0],c[1],c[2]] for c in C],float)
    r=int(np.linalg.matrix_rank(F,tol=1e-10))
    return len(C)<=r


class PrunedGeneralModel(t.RootOnly):
    def __init__(self,libpath,columns,log_ratio=1.1,long_ratio=2.0,
                 max_quad_facets=3, add_flow_packing=True, add_equal_split_symmetry=True):
        self._pg_quad=max_quad_facets
        self._pg_add_pack=add_flow_packing
        self._pg_add_sym=add_equal_split_symmetry
        self._first_hull_cache={};self._quad_cache={};self._quad_need_cache={}
        # RootOnly invokes our build sites / first hull, then calls our pair and quadratic methods.
        super().__init__(libpath,columns,log_ratio,long_ratio,
                         pair_moments=True,max_first_facets=0,max_quad_facets=max_quad_facets,
                         pair_bounds='full')
        if add_flow_packing:self._add_flow_pack()
        if add_equal_split_symmetry:self._add_equal_split_symmetry()

    # ---------- per-column network pruning + q elimination ----------
    def _build_sites(self):
        self.req={};self.hlong={};self.z={};self.q={};self.lam={};self.obj={};self.log_error_bound=0;self.qr={}
        pts=base.vl.geometric_breaks(1e-12,self.log_ratio)
        for p,pat in enumerate(self.patterns):
            I=self.info[p]
            for r in self.roots:
                # upper bound only for variable box; activity scaling is supplied by exact hull facets
                vals=I['rootvals'][r];hi=max(vals.values())/I['scale'][self.L]
                self.qr[p,r]=self.lp.var(f'qr:{p}:{r}',0,float(hi))
            lam=[]
            for kk,x in enumerate(pts):
                v=self.lp.var(f'lam:{p}:{kk}',0,1);lam.append(v);self.obj[v]=self.weights[p]*math.log(x)
            self.lam[p]=lam
            self.lp.add_eq({v:1 for v in lam},1)
            # eliminate q_s: sum_r qr == sum_k lambda_k q_k
            row={self.qr[p,r]:1 for r in self.roots}
            for k,v in enumerate(lam):row[v]=row.get(v,0)-pts[k]
            self.lp.add_eq(row,0)
            self.log_error_bound+=self.weights[p]*base.vl.chord_gap(self.log_ratio)

    # ---------- exact first-order root hull, cached by (global counts, split size) ----------
    def _first_hull_for(self,p,r):
        I=self.info[p];counts=tuple(I['counts']);a=int(self.units[r]['left_flow']);key=(counts,a)
        if key in self._first_hull_cache:return self._first_hull_cache[key]
        b=self.L-a;raw=I['raw'];PsA,PlA=base.trans_pair(a,self.long_ratio);PsB,PlB=base.trans_pair(b,self.long_ratio)
        X=[];vals=[]
        for (mm,c1),A in raw.items():
            if mm!=a:continue
            c2=tuple(counts[k]-c1[k] for k in range(4))
            if min(c2)<0 or (b,c2) not in raw:continue
            B=raw[(b,c2)]
            for la,PA in enumerate((PsA,PlA)):
                for lb,PB in enumerate((PsB,PlB)):
                    X.append((c1[0],c1[1],c1[2],la,lb))
                    vals.append(float(PI@((PA@A[:,1])*(PB@B[:,1])))/I['scale'][self.L])
        fs=h.upper_hull_facets(X,vals)
        self._first_hull_cache[key]=fs
        return fs

    def _add_root_length_partition_facets(self):
        self.root_length_facets=0
        for p,pat in enumerate(self.patterns):
            for r in self.roots:
                m,s,_=self._req_slot(r,'left');Aact=self.act[r];eL=self.ell[r,'left'];eR=self.ell[r,'right']
                for f in self._first_hull_for(p,r):
                    a0,dA,dC,dG,dL,dR=f
                    row={self.qr[p,r]:1,Aact:-a0,eL:-dL,eR:-dR}
                    for tax,ch in enumerate(pat):
                        c={'A':dA,'C':dC,'G':dG,'T':0.0}[ch]
                        if c:
                            x=self.tmem[m,0,s,tax];row[x]=row.get(x,0)-c
                    self.lp.add_le(row,0);self.root_length_facets+=1

    # ---------- certified second-order pruning ----------
    def _needs_quad(self,p,a):
        key=(p,a)
        if key not in self._quad_need_cache:
            self._quad_need_cache[key]=not certified_quadratic_redundant(tuple(self.info[p]['counts']),a)
        return self._quad_need_cache[key]

    def _add_pairs(self):
        self.quad_needed_sizes={a:any(self._needs_quad(p,a) for p in range(len(self.patterns))) for a in range(2,self.L//2+1)}
        self.pair_bound_rows=0;self.pair_degree_rows=0;self.pair_pruned_sizes=[]
        for r in self.roots:
            a=int(self.units[r]['left_flow']);A=self.act[r];m,s,_=self._req_slot(r,'left')
            if a<=1 or not self.quad_needed_sizes.get(a,False):
                if a>1:self.pair_pruned_sizes.append(a)
                continue
            for u in range(self.L):
                for v in range(u+1,self.L):
                    y=self.lp.var(f'rpair:{r}:{u}:{v}',0,1);self.rpair[r,u,v]=y
                    # a=2: nonnegativity + degree equations already imply y<=x_u,x_v and Frechet lower bound.
                    if a>=3:
                        xt=self.tmem[m,0,s,u];xv=self.tmem[m,0,s,v]
                        self.lp.add_le({y:1,xt:-1},0);self.lp.add_le({y:1,xv:-1},0)
                        self.lp.add_ge({y:1,xt:-1,xv:-1,A:1},0);self.pair_bound_rows+=3
            for u in range(self.L):
                row={}
                for v in range(self.L):
                    if v==u:continue
                    row[self.rpair[r,min(u,v),max(u,v)]]=row.get(self.rpair[r,min(u,v),max(u,v)],0)+1
                row[self.tmem[m,0,s,u]]=-(a-1);self.lp.add_eq(row,0);self.pair_degree_rows+=1

    def _moment_exprs(self,r,pat):
        # same as RootOnly, but called only where pair variables exist
        return t.RootOnly._moment_exprs(self,r,pat)

    def _quad_facets_for(self,p,r):
        I=self.info[p];counts=tuple(I['counts']);a=int(self.units[r]['left_flow']);key=(counts,a,self._pg_quad)
        if key in self._quad_cache:return self._quad_cache[key]
        b=self.L-a;raw=I['raw'];PsA,PlA=base.trans_pair(a,self.long_ratio);PsB,PlB=base.trans_pair(b,self.long_ratio)
        X=[];vals=[]
        for (mm,c1),BA in raw.items():
            if mm!=a:continue
            c2=tuple(counts[k]-c1[k] for k in range(4))
            if min(c2)<0 or (b,c2) not in raw:continue
            feat=[c1[0],c1[1],c1[2]]+[c1[bb]*c1[cc] for bb in range(4) for cc in range(bb,4)]
            for la,PA in enumerate((PsA,PlA)):
                for lb,PB in enumerate((PsB,PlB)):
                    X.append(feat+[la,lb]);vals.append(float(PI@((PA@BA[:,1])*(PB@raw[(b,c2)][:,1])))/I['scale'][self.L])
        fs=t.facetfun(X,vals,self._pg_quad)
        self._quad_cache[key]=fs
        return fs

    def _add_quad_facets(self):
        self.root_moment_facets=0;self.quad_site_root_pruned=0
        for p,pat in enumerate(self.patterns):
            for r in self.roots:
                a=int(self.units[r]['left_flow'])
                if a<=1 or not self._needs_quad(p,a):
                    if a>1:self.quad_site_root_pruned+=1
                    continue
                ex=self._moment_exprs(r,pat);Aact=self.act[r];eL=self.ell[r,'left'];eR=self.ell[r,'right']
                for f in self._quad_facets_for(p,r):
                    a0,*coef=f;row={self.qr[p,r]:1,Aact:-a0,eL:-coef[len(ex)],eR:-coef[len(ex)+1]}
                    for hh,d in enumerate(ex):
                        c=coef[hh]
                        if c:
                            for v,cc in d.items():row[v]=row.get(v,0)-c*cc
                    self.lp.add_le(row,0);self.root_moment_facets+=1

    # ---------- flow-k packing certificate ----------
    def _add_flow_pack(self):
        self.flow_packing_cuts=0
        for k in range(2,self.L):
            for tax in range(self.L):
                row={}
                for u in self.by[k]:
                    for port in ('left','right'):
                        m,s,_=self._req_slot(u,port);x=self.tmem[m,0,s,tax];row[x]=row.get(x,0)+1
                if row:self.lp.add_le(row,1);self.flow_packing_cuts+=1

    # ---------- exact symmetry break when equal arms are exchangeable ----------
    def _add_equal_split_symmetry(self):
        self.equal_split_symmetry_rows=0
        # Current model uses trans_pair(flow,long_ratio) for both ports, so equal flows have identical families.
        weights=[(i+1)/self.L for i in range(self.L)]
        for u in self.splits:
            a=int(self.units[u]['left_flow']);b=int(self.units[u]['right_flow'])
            if a!=b:continue
            ml,sl,_=self._req_slot(u,'left');mr,sr,_=self._req_slot(u,'right');row={}
            for tax,w in enumerate(weights):
                xl=self.tmem[ml,0,sl,tax];xr=self.tmem[mr,0,sr,tax];row[xl]=row.get(xl,0)+w;row[xr]=row.get(xr,0)-w
            self.lp.add_ge(row,0);self.equal_split_symmetry_rows+=1

    def solve(self):
        val,res,dt=self.lp.solve(self.obj,True)
        const=sum(self.weights[p]*math.log(self.info[p]['scale'][self.L]) for p in range(len(self.patterns)))
        return val+const,res,dt

if __name__=='__main__':
    cols=[tuple('AAAACCCC'),tuple('ACGTACGT'),tuple('AACCGGTT'),tuple('AGCTAGCT'),tuple('AAAACCGT'),tuple('AACCGGTT')]
    for lib in ['/mnt/data/test_match_units8/library.json','/mnt/data/test_match_units12/library.json']:
        L=8 if '8/' in lib else 12
        cc=cols if L==8 else t.random_columns(12,4,7000)
        M=PrunedGeneralModel(lib,cc,1.1,2.0);v,_,dt=M.solve()
        print(L,v,len(M.lp.names),len(M.lp.eq),len(M.lp.le),'first',M.root_length_facets,'quad',M.root_moment_facets,'qpruned',M.quad_site_root_pruned,'pairpruned',M.pair_pruned_sizes,'solve',dt)
