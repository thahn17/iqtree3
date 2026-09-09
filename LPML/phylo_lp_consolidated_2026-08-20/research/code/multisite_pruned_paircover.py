#!/usr/bin/env python3
"""Pruned general root model plus one-pass taxon-pair-cover anchor networks.
Anchors are chosen from the alignment only, before optimization, to distinguish every taxon pair
that is distinguishable anywhere in the full alignment. Internal anchor propagation uses a 3-contrast
branch-length lift (one state common mode removed because (P_long-P_short)1=0).
"""
from __future__ import annotations
import importlib.util,sys
spec=importlib.util.spec_from_file_location('p','/mnt/data/multisite_pruned_general.py');p=importlib.util.module_from_spec(spec);sys.modules['p']=p;spec.loader.exec_module(p)
base=p.base;BASE_INDEX=base.BASE_INDEX;PI=base.PI

def pair_cover_anchors(patterns):
    if not patterns:return []
    L=len(patterns[0]);U={(i,j) for i in range(L) for j in range(i+1,L) if any(col[i]!=col[j] for col in patterns)}
    cover=[{q for q in U if col[q[0]]!=col[q[1]]} for col in patterns];sel=[]
    while U:
        best=max(range(len(cover)),key=lambda k:(len(cover[k]&U),-k));gain=cover[best]&U
        if not gain:break
        sel.append(best);U-=gain
    return sel

class PairCoverHybrid(p.PrunedGeneralModel):
    def __init__(self,*args,anchors='pair-cover',**kw):
        super().__init__(*args,**kw)
        if anchors=='pair-cover':self.anchor_ids=pair_cover_anchors(self.patterns)
        else:self.anchor_ids=list(anchors or [])
        self.areq={};self.ag={};self.az={}
        for q in self.anchor_ids:self._add_anchor(q)

    def _msg(self,q,u,port,i):
        m=int(self.units[u][port+'_flow']);R=self.areq[q,u,port];Ps,Pl=base.trans_pair(m,self.long_ratio);D=Pl-Ps
        row={R[j]:float(Ps[i,j]) for j in range(4)};G=self.ag[q,u,port]
        for j in range(3):row[G[j]]=row.get(G[j],0)+float(D[i,j])
        return row

    def _add_anchor(self,q):
        pat=self.patterns[q];I=self.info[q]
        for m in range(1,self.L):
            ub=I['ubraw'][m]
            for u,port in self.reqdef[m]:
                a=self.act[u];ell=self.ell[u,port];R=[];G=[]
                for j in range(4):
                    v=self.lp.var(f'areq:{q}:{u}:{port}:{j}',0,float(ub[j]));R.append(v);self.lp.add_le({v:1,a:-float(ub[j])},0)
                for j in range(3):
                    lo=-float(ub[3]);hi=float(ub[j]);g=self.lp.var(f'ag:{q}:{u}:{port}:{j}',min(0.,lo),max(0.,hi));G.append(g)
                    self.lp.add_le({g:1,ell:-hi},0);self.lp.add_le({g:-1,ell:lo},0)
                    self.lp.add_le({R[j]:1,R[3]:-1,g:-1,a:-hi,ell:hi},0)
                    self.lp.add_le({R[j]:-1,R[3]:1,g:1,a:lo,ell:-lo},0)
                self.areq[q,u,port]=R;self.ag[q,u,port]=G
        # exact flow-1 state composition from shared named taxa
        for u,port in self.reqdef[1]:
            m,s,_=self._req_slot(u,port);R=self.areq[q,u,port]
            for j in range(4):
                row={R[j]:1}
                for tax,ch in enumerate(pat):
                    if BASE_INDEX[ch]==j:
                        x=self.tmem[m,0,s,tax];row[x]=row.get(x,0)-1
                self.lp.add_eq(row,0)
        # products
        for u in self.splits:
            aN=int(self.units[u]['left_flow']);bN=int(self.units[u]['right_flow']);m=int(self.units[u]['descendant_count']);act=self.act[u]
            k=I['scale'][aN]*I['scale'][bN]/I['scale'][m];uA=I['ubmsg'][aN];uB=I['ubmsg'][bN];Z=[]
            for i in range(4):
                z=self.lp.var(f'az:{q}:{u}:{i}',0,float(I['ubraw'][m][i]));Z.append(z);A=self._msg(q,u,'left',i);B=self._msg(q,u,'right',i)
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
            self.az[q,u]=Z
        # count conservation
        for m in range(1,self.L):
            for j in range(4):
                row={self.areq[q,u,port][j]:1 for u,port in self.reqdef[m]}
                if m==1:self.lp.add_eq(row,sum(1. for ch in pat if BASE_INDEX[ch]==j))
                else:
                    for u in self.by[m]:row[self.az[q,u][j]]=row.get(self.az[q,u][j],0)-1
                    self.lp.add_eq(row,0)
        # anchor root contribution equals propagated likelihood; general root hull stays in force as extra valid cuts
        for r in self.roots:
            row={self.qr[q,r]:1}
            for i in range(4):row[self.az[q,r][i]]=row.get(self.az[q,r][i],0)-float(PI[i])
            self.lp.add_eq(row,0)

if __name__=='__main__':
 cols=[tuple('AAAACCCC'),tuple('ACGTACGT'),tuple('AACCGGTT'),tuple('AGCTAGCT'),tuple('AAAACCGT'),tuple('AACCGGTT')]
 M=PairCoverHybrid('/mnt/data/test_match_units8/library.json',cols,1.1,2.0);v,_,dt=M.solve();print('anchors',M.anchor_ids,'obj',v,'vars',len(M.lp.names),'eq',len(M.lp.eq),'le',len(M.lp.le),'solve',dt)
