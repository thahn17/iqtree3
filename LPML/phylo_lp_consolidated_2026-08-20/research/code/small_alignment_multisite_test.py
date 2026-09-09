#!/usr/bin/env python3
"""Small shared-tree alignment experiment.

- one shared continuous topology for all columns
- exact duplicate-column compression
- count-adaptive likelihood scaling on top of beta^m scaling
- split-aware composition root cuts
- adaptive tangent outer approximation of log
- optional exact 8-taxon benchmark

This is a research harness, not a production solver.
"""
import importlib.util, math, numpy as np
from collections import Counter

spec=importlib.util.spec_from_file_location('ms','/mnt/data/multisite_block_scaled.py')
ms=importlib.util.module_from_spec(spec);spec.loader.exec_module(ms)

LIB='/mnt/data/multisite_units8/library.json'
# Columns are listed across taxa leaf_1 ... leaf_8.
COLUMNS=[
    tuple('AAAACCCC'),
    tuple('ACGTACGT'),
    tuple('AACCGGTT'),
    tuple('AGCTAGCT'),
    tuple('AAAACCGT'),
    tuple('AACCGGTT'),  # exact duplicate; compressed to weight 2
]


def main():
    counts=Counter(COLUMNS)
    patterns=list(counts)
    weights=np.array([counts[p] for p in patterns],dtype=float)
    exact,info,_=ms.base.exact_shared_best(patterns,weights)

    M=ms.Model(LIB,[[p] for p in patterns])
    obj={M.logv[i]:float(weights[i]) for i in range(len(patterns))}
    hist=[]; total=0.0
    for it in range(15):
        val,res,dt=M.lp.solve(obj,True); total+=dt
        if not res.success:
            raise RuntimeError(res.message)
        gaps=[]
        for i in range(len(patterns)):
            q=max(res.x[M.q[i]],1e-12)
            gaps.append(res.x[M.logv[i]]-math.log(q))
        hist.append((it,max(gaps),dt))
        if max(gaps)<1e-4:
            break
        for i,g in enumerate(gaps):
            if g>1e-4:
                M.addtan(i,max(res.x[M.q[i]],1e-12))

    const=sum(weights[i]*(math.log(M.info[i]['sroot'])+8*math.log(ms.BETA)) for i in range(len(patterns)))
    upper=val+const
    print('columns:',len(COLUMNS),'unique:',len(patterns),'weights:',weights.tolist())
    print('exact shared-tree logL:',exact)
    print('LP logL upper:',upper)
    print('log gap:',upper-exact)
    print('variables:',len(M.lp.names),'equalities:',len(M.lp.eq),'inequalities:',len(M.lp.le))
    print('solve seconds:',total,'log OA history:',hist)
    print('fractional root activities:',[(u,res.x[M.act[u]]) for u in M.roots if res.x[M.act[u]]>1e-7])

if __name__=='__main__':
    main()
