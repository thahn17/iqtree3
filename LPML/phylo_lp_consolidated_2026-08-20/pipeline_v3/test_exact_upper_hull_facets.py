#!/usr/bin/env python3
"""Small exact-upper-hull helper retained for oracle/small-L preprocessing."""
from __future__ import annotations
import numpy as np
from scipy.spatial import ConvexHull, QhullError
import test_rootonly_pruning as t


def upper_hull_facets(features,vals):
    X=np.asarray(features,float);y=np.asarray(vals,float)
    if len(X)==0:return []
    P=np.column_stack([X,y])
    try:
        hull=ConvexHull(P,qhull_options='QJ')
        fs=[]
        for eq in hull.equations:
            n=eq[:-1];off=eq[-1];nq=n[-1]
            if nq<=1e-10:continue
            f=np.concatenate(([-off/nq],-n[:-1]/nq))
            if np.all(np.column_stack([np.ones(len(X)),X])@f+1e-7>=y):
                k=tuple(np.round(f,11))
                if k not in fs:fs.append(k)
        if fs:return fs
    except (QhullError,ValueError):pass
    return t.facetfun(X,y,min(max(len(X),1),60))

class RootOnlyHull(t.RootOnly):
    def _first_facets(self,p,r):
        I=self.info[p];counts=I['counts'];a=int(self.units[r]['left_flow']);b=self.L-a;raw=I['raw'];X=[];y=[]
        PsA,PlA=t.base.trans_pair(a,self.long_ratio);PsB,PlB=t.base.trans_pair(b,self.long_ratio)
        for (mm,c1),A in raw.items():
            if mm!=a:continue
            c2=tuple(counts[k]-c1[k] for k in range(4))
            if min(c2)<0 or (b,c2) not in raw:continue
            for la,PA in enumerate((PsA,PlA)):
                for lb,PB in enumerate((PsB,PlB)):
                    X.append([c1[0],c1[1],c1[2],la,lb]);y.append(float(t.base.PI@((PA@A[:,1])*(PB@raw[(b,c2)][:,1])))/I['scale'][self.L])
        return upper_hull_facets(X,y)
