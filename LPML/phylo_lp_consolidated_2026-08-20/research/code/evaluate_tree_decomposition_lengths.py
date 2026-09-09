#!/usr/bin/env python3
from __future__ import annotations
import importlib.util,sys,random,csv,time

def load(name,path):
 s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);sys.modules[name]=m;s.loader.exec_module(m);return m
pg=load('pg_eval','/mnt/data/multisite_pruned_general.py')
dec=load('dec_eval','/mnt/data/fractional_tree_decomposition.py')
rows=[];comp_rows=[]
cases=[]
cols8=[tuple('AAAACCCC'),tuple('ACGTACGT'),tuple('AACCGGTT'),tuple('AGCTAGCT'),tuple('AAAACCGT'),tuple('AACCGGTT')]
cases.append((8,'/mnt/data/test_match_units8/library.json',cols8,800,.8,10))
rng=random.Random(7);cols12=[tuple(rng.choice('ACGT') for _ in range(12)) for _ in range(4)]
cases.append((12,'/mnt/data/test_match_units12/library.json',cols12,400,.65,8))
for L,lib,cols,ns,temp,pool in cases:
 M=pg.PrunedGeneralModel(lib,cols,1.1,2.0);obj,res,_=M.solve();D=dec.FractionalTreeDecomposer(M,res)
 for strat in ['after','binary-before']:
  t=time.time();out=D.decompose(n_samples=ns,seed=17,temperature=temp,candidate_pool=pool,include_pairs=True,length_strategy=strat,min_proportion=1e-3,max_fit_candidates=160);dt=time.time()-t
  rows.append({'L':L,'strategy':strat,'candidates_fit':out.candidate_count,'components_ge_0.001':len(out.components),'topology_explained':out.topology_explained_fraction,'topology_rmse':out.topology_rmse,'length_explained_postfit':out.length_explained_fraction,'length_rmse_postfit':out.length_postfit_rmse,'seconds':dt,'top_share':out.components[0].proportion if out.components else 0})
  if strat=='after':
   for rank,c in enumerate(out.components[:15],1):comp_rows.append({'L':L,'rank':rank,'proportion':c.proportion,'root_unit':c.root_unit,'newick':c.newick})
with open('/mnt/data/tree_decomposition_length_strategy.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
with open('/mnt/data/fractional_tree_components_example.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=comp_rows[0].keys());w.writeheader();w.writerows(comp_rows)
