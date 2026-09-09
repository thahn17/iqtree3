#!/usr/bin/env python3
from __future__ import annotations
import importlib.util,sys,time,csv,math
import numpy as np
spec=importlib.util.spec_from_file_location('t','/mnt/data/test_rootonly_pruning.py');t=importlib.util.module_from_spec(spec);sys.modules['t']=t;spec.loader.exec_module(t)
spec2=importlib.util.spec_from_file_location('h','/mnt/data/test_exact_upper_hull_facets.py');h=importlib.util.module_from_spec(spec2);sys.modules['h']=h;spec2.loader.exec_module(h)
base=t.base;PI=base.PI

class RootOnlyLambda(t.RootOnly):
    def _add_root_length_partition_facets(self):
        self.root_length_facets=0;self.root_lambda_vars=0;self.root_point_total=0
        for p,pat in enumerate(self.patterns):
            I=self.info[p];counts=I['counts'];raw=I['raw']
            for r in self.roots:
                a=int(self.units[r]['left_flow']);b=self.L-a;PsA,PlA=base.trans_pair(a,self.long_ratio);PsB,PlB=base.trans_pair(b,self.long_ratio)
                points=[]
                for (mm,c1),A in raw.items():
                    if mm!=a:continue
                    c2=tuple(counts[k]-c1[k] for k in range(4))
                    if min(c2)<0 or (b,c2) not in raw:continue
                    B=raw[(b,c2)]
                    for la,PA in enumerate((PsA,PlA)):
                        for lb,PB in enumerate((PsB,PlB)):
                            v=float(PI@((PA@A[:,1])*(PB@B[:,1])))/I['scale'][self.L]
                            points.append((c1[0],c1[1],c1[2],la,lb,v))
                self.root_point_total+=len(points)
                lams=[]
                for j,pt in enumerate(points):
                    z=self.lp.var(f'rlam:{p}:{r}:{j}',0,1);lams.append(z);self.root_lambda_vars+=1
                Aact=self.act[r];eL=self.ell[r,'left'];eR=self.ell[r,'right'];m,s,_=self._req_slot(r,'left')
                # perspective convex-combination mass
                row={z:1 for z in lams};row[Aact]=-1;self.lp.add_eq(row,0)
                # composition features
                for bi,ch in enumerate('ACG'):
                    row={}
                    for j,pt in enumerate(points):
                        if pt[bi]:row[lams[j]]=float(pt[bi])
                    for tax,bch in enumerate(pat):
                        if bch==ch:
                            x=self.tmem[m,0,s,tax];row[x]=row.get(x,0)-1
                    self.lp.add_eq(row,0)
                row={lams[j]:float(points[j][3]) for j in range(len(points))};row[eL]=-1;self.lp.add_eq(row,0)
                row={lams[j]:float(points[j][4]) for j in range(len(points))};row[eR]=-1;self.lp.add_eq(row,0)
                # q <= weighted upper value
                row={self.qr[p,r]:1}
                for j,z in enumerate(lams):row[z]=row.get(z,0)-float(points[j][5])
                self.lp.add_le(row,0);self.root_length_facets+=1

class RootOnlyExactHull(h.RootOnlyHull):
    pass

def bench(label,lib,cols):
    rows=[]
    variants=[
      ('K5',t.RootOnly,dict(max_first_facets=5,max_quad_facets=1)),
      ('K20',t.RootOnly,dict(max_first_facets=20,max_quad_facets=1)),
      ('K60',t.RootOnly,dict(max_first_facets=60,max_quad_facets=1)),
      ('lambda',RootOnlyLambda,dict(max_first_facets=0,max_quad_facets=1)),
      ('exact_hull',RootOnlyExactHull,dict(max_quad_facets=1)),
    ]
    for name,Cls,kw in variants:
        st=time.perf_counter()
        try:
            M=Cls(lib,cols,1.1,2.0,**kw);build=time.perf_counter()-st
            st2=time.perf_counter();v,_,solve=M.solve();wall=time.perf_counter()-st2
            row=dict(case=label,method=name,objective=v,build_s=build,solve_s=solve,solve_wall_s=wall,total_s=build+wall,
                     vars=len(M.lp.names),eq=len(M.lp.eq),le=len(M.lp.le),root_rows=getattr(M,'root_length_facets',None),root_lambda_vars=getattr(M,'root_lambda_vars',0),root_points=getattr(M,'root_point_total',0))
            print(row,flush=True);rows.append(row)
        except Exception as e:
            row=dict(case=label,method=name,error=repr(e));print(row,flush=True);rows.append(row)
    return rows

if __name__=='__main__':
    allr=[]
    allr+=bench('L8_1pat','/mnt/data/test_match_units8/library.json',[tuple('ACGTACGT')])
    allr+=bench('L12_1pat','/mnt/data/test_match_units12/library.json',[t.random_columns(12,1,702)[0]])
    allr+=bench('L16_1pat','/mnt/data/test_match_units16/library.json',[t.random_columns(16,1,703)[0]])
    keys=sorted({k for r in allr for k in r})
    with open('/mnt/data/root_envelope_tradeoff.csv','w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(allr)
