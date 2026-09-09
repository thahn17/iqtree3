#!/usr/bin/env python3
"""End-to-end FASTA -> fractional phylogenetic LP -> top-K integral trees -> ML branch lengths.

Main phases (all timed):
  1. FASTA read / column sample / substitution-model setup
  2. universal fixed-flow library + LP build/preprocessing
  3. continuous LP solve
  4. non-destructive top-K integral-tree identification
  5. continuous branch-length optimization and final likelihood ranking

The LP keeps variable lower-arm branch lengths.  Integral-tree identification uses topology
frequency/taxon structure, then each top-K topology is rescored with full continuous branch
length optimization using Felsenstein pruning + analytic gradients + L-BFGS-B.
"""
from __future__ import annotations
import argparse,csv,json,math,random,sys,time
from dataclasses import asdict
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parent
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))

from phylo_substitution import BASES,BASE_INDEX,jc69,gtr_f,set_active_model,build_spectral_relaxation,is_jc_like
from phylo_branch_opt import IUPAC,TreeLikelihood


def read_fasta(path):
    records=[];name=None;buf=[]
    with open(path,encoding='utf-8') as f:
        for raw in f:
            line=raw.strip()
            if not line:continue
            if line.startswith('>'):
                if name is not None:records.append((name,''.join(buf).upper()))
                name=line[1:].strip()
                if not name:raise ValueError('Empty FASTA header')
                # Use full header as taxon ID; it is quoted in Newick when needed.
                buf=[]
            else:
                if name is None:raise ValueError('FASTA sequence before first header')
                buf.append(''.join(line.split()))
    if name is not None:records.append((name,''.join(buf).upper()))
    if not records:raise ValueError('No FASTA records found')
    if len({n for n,s in records})!=len(records):raise ValueError('Duplicate FASTA taxon IDs')
    n=len(records[0][1])
    if n==0 or any(len(s)!=n for _,s in records):raise ValueError('FASTA is not a rectangular nonempty alignment')
    return records


def select_columns(records,max_columns=None,seed=1):
    n=len(records[0][1])
    if max_columns is None or max_columns<=0 or max_columns>=n:idx=list(range(n))
    else:idx=sorted(random.Random(seed).sample(range(n),int(max_columns)))
    out=[(name,''.join(seq[j] for j in idx)) for name,seq in records]
    return out,idx


def empirical_frequencies(records,pseudocount=1e-8):
    c=np.full(4,float(pseudocount))
    for _,seq in records:
        for ch in seq.upper():
            v=np.asarray(IUPAC.get(ch,IUPAC['N']),float);z=v.sum()
            if z<=0 or z>=4:continue  # gaps/N are uninformative for +F
            c+=v/z
    if c.sum()<=0:raise ValueError('Cannot estimate +F frequencies: no informative DNA characters')
    return c/c.sum()


def lp_columns(records):
    """Columns for the LP: exact A/C/G/T; any ambiguity is safely relaxed to null N."""
    seqs=[]
    for _,s in records:
        seqs.append(''.join(('T' if c=='U' else c) if ('T' if c=='U' else c) in BASE_INDEX else 'N' for c in s.upper()))
    n=len(seqs[0]);return [tuple(seqs[t][j] for t in range(len(seqs))) for j in range(n)]


def write_fasta(records,path,width=80):
    with open(path,'w',encoding='utf-8') as f:
        for name,seq in records:
            f.write('>'+name+'\n')
            for i in range(0,len(seq),width):f.write(seq[i:i+width]+'\n')


def configure_lp_model(model):
    # Configure retained research modules as well as newly imported ones.
    import multisite_variable_length_compressed as vl;vl.configure_substitution_model(model)
    import multisite_lowerarm_taxonflow as low;low.configure_substitution_model(model)
    import test_rootonly_pruning as tr;tr.base.configure_substitution_model(model);tr.PI=model.pi
    import multisite_pruned_general as pg;pg.base.configure_substitution_model(model);pg.PI=model.pi
    return pg


def build_library(labels,out_dir):
    from fixed_flow_unit_generator_variable_lengths import LibraryBuilder
    out_dir=Path(out_dir);out_dir.mkdir(parents=True,exist_ok=True)
    b=LibraryBuilder(labels,variable_length_paths=True,long_path_edges=2,variable_length_scope='branches')
    b.build();paths=b.write_library(out_dir)
    return str(out_dir/'library.json'),paths


def parse_rates(s):
    vals=[float(x) for x in s.replace('/',',').replace(':',',').split(',') if x.strip()]
    if len(vals)!=6:raise argparse.ArgumentTypeError('Need six GTR rates: AC,AG,AT,CG,CT,GT')
    return vals


def parse_branch_max(s):
    if isinstance(s,(int,float)):return float(s)
    q=str(s).strip().lower()
    if q in ('auto','none','inf','infinity'):return None
    try:v=float(q)
    except ValueError as e:raise argparse.ArgumentTypeError('branch max must be a positive number or auto') from e
    if not math.isfinite(v) or v<=0:raise argparse.ArgumentTypeError('branch max must be positive or auto')
    return v


def timing_print(name,seconds):print(f'[timing] {name}: {seconds:.6f} s',flush=True)


def run(args):
    out=Path(args.out_prefix);out.parent.mkdir(parents=True,exist_ok=True)
    times={};total0=time.perf_counter()

    # 1. input/model
    st=time.perf_counter();records=read_fasta(args.fasta);original_sites=len(records[0][1])
    records,indices=select_columns(records,args.max_columns,args.seed)
    if args.library:
        libmeta=json.loads(Path(args.library).read_text());taxa=list(libmeta.get('taxa',[]));rd=dict(records)
        if set(taxa)!=set(rd):raise ValueError('FASTA taxon IDs must exactly match prebuilt library taxa')
        records=[(t,rd[t]) for t in taxa]
    labels=[n for n,s in records];seqs=[s for n,s in records]
    if args.model=='jc':model=jc69()
    else:
        pi=empirical_frequencies(records,args.freq_pseudocount);model=gtr_f(pi,args.gtr_rates)
    set_active_model(model);pg=configure_lp_model(model)
    cols=lp_columns(records)
    # LP duplicate compression occurs inside the model after ambiguity->N mapping.
    lp_unique_patterns=len(set(cols))
    exact_unique_patterns=len(set(tuple(seq[j].upper() for seq in seqs) for j in range(len(seqs[0]))))
    sampled_path=None
    if len(indices)<original_sites:
        sampled_path=str(out)+'.sampled.fa';write_fasta(records,sampled_path)
    times['01_input_model_s']=time.perf_counter()-st;timing_print('01 input/sample/model',times['01_input_model_s'])

    # 2. library + LP build/preprocessing
    st=time.perf_counter()
    if args.library:libpath=args.library
    else:libpath,_=build_library(labels,str(out)+'.units')
    M=pg.PrunedGeneralModel(libpath,cols,log_ratio=args.log_ratio,long_ratio=args.long_ratio,
                            max_quad_facets=args.max_quad_facets,preprocess_mode=args.preprocess_mode,
                            substitution_model=model,branch_relaxation=args.branch_relaxation,
                            spectral_support_points=args.spectral_support_points)
    times['02_library_lp_build_s']=time.perf_counter()-st;timing_print('02 library + LP build/preprocess',times['02_library_lp_build_s'])

    # 3. LP solve
    st=time.perf_counter();lp_score,res,solver_reported=M.solve();times['03_lp_solve_s']=time.perf_counter()-st
    timing_print('03 LP solve',times['03_lp_solve_s'])

    # 4. top-K structural candidates, non-destructive
    st=time.perf_counter()
    import fractional_tree_hierarchical_units as hd
    cand,frac=hd.rank_integral_tree_candidates(M,res,topk=args.top_k,root_beam=args.root_beam,
        leaf_beam=args.leaf_beam,child_branching=args.child_branching,leaf_branching=args.leaf_branching,
        merge_branching=args.merge_branching,max_expansions=args.max_expansions)
    if not cand:raise RuntimeError('No complete integral tree candidate was recovered from the fractional LP')
    times['04_tree_identification_s']=time.perf_counter()-st;timing_print('04 top-K tree identification',times['04_tree_identification_s'])

    # 5. optimize continuous branch lengths and rank by actual alignment likelihood
    st=time.perf_counter();scored=[]
    short_time=pg.base.branch_times
    for srank,z in enumerate(cand,1):
        TL=TreeLikelihood(z.topology,labels,seqs,model,decoded=z,frac=frac,short_time=short_time,long_ratio=args.long_ratio)
        opt=TL.optimize(args.branch_min,args.branch_max,args.branch_maxiter,args.branch_ftol,multistart=args.branch_multistart)
        scored.append(dict(structural_rank=srank,bottleneck_support=float(z.proportion),secondary_score=float(z.secondary_score),
                           direction=z.direction,log_likelihood=float(opt.log_likelihood),optimization_success=bool(opt.success),
                           optimization_iterations=opt.iterations,likelihood_evaluations=opt.evaluations,newick=opt.newick,
                           optimizer_message=opt.message,branch_min_bound=opt.min_length_bound,branch_max_bound=opt.max_length_bound))
    scored.sort(key=lambda r:r['log_likelihood'],reverse=True)
    for i,r in enumerate(scored,1):r['likelihood_rank']=i
    best=scored[0]
    times['05_branch_opt_rescore_s']=time.perf_counter()-st;timing_print('05 branch optimization + rescoring',times['05_branch_opt_rescore_s'])
    times['total_s']=time.perf_counter()-total0;timing_print('TOTAL',times['total_s'])

    # Outputs
    topk_path=str(out)+'.topk.tsv'
    fields=['likelihood_rank','structural_rank','log_likelihood','bottleneck_support','secondary_score','direction',
            'optimization_success','optimization_iterations','likelihood_evaluations','branch_min_bound','branch_max_bound','newick','optimizer_message']
    with open(topk_path,'w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(scored)
    best_path=str(out)+'.best.tree';Path(best_path).write_text(best['newick']+'\n',encoding='utf-8')
    timing_path=str(out)+'.timings.json';Path(timing_path).write_text(json.dumps(times,indent=2)+'\n')
    summary={
      'fasta':str(args.fasta),'taxa':len(labels),'original_columns':original_sites,'used_columns':len(indices),
      'sampled_column_indices_1based':[i+1 for i in indices] if len(indices)<original_sites else 'all',
      'sampled_fasta':sampled_path,'library':libpath,'substitution_model':model.name,
      'unique_exact_site_patterns':exact_unique_patterns,'unique_lp_site_patterns':lp_unique_patterns,
      'duplicate_columns_merged_lp':len(indices)-lp_unique_patterns,
      'stationary_frequencies':model.pi.tolist(),'gtr_rates_AC_AG_AT_CG_CT_GT':list(model.rates) if model.rates else None,
      'preprocess_mode':args.preprocess_mode,'branch_relaxation_requested':args.branch_relaxation,
      'branch_relaxation_effective':getattr(M,'branch_relaxation','two-point'),
      'spectral_support_points':args.spectral_support_points if getattr(M,'branch_relaxation','two-point')=='spectral' else None,
      'spectral_decay_rates':getattr(M,'spectral_mu',np.array([])).tolist() if getattr(M,'branch_relaxation','two-point')=='spectral' else None,
      'spectral_mode_count':len(getattr(M,'spectral_mu',[])) if getattr(M,'branch_relaxation','two-point')=='spectral' else None,
      'spectral_polytope_constraints':len(pg.base.spectral_relaxation().b) if getattr(M,'branch_relaxation','two-point')=='spectral' else None,
      'spectral_polytope_vertices':len(pg.base.spectral_relaxation().vertices) if getattr(M,'branch_relaxation','two-point')=='spectral' else None,
      'lp_variables':len(M.lp.names),'lp_equalities':len(M.lp.eq),'lp_inequalities':len(M.lp.le),
      'lp_relaxation_log_score':float(lp_score),'lp_solver_reported_s':float(solver_reported),
      'top_k_requested':args.top_k,'top_k_recovered':len(cand),'best_topk_log_likelihood':best['log_likelihood'],
      'best_tree_file':best_path,'topk_table':topk_path,'timings_file':timing_path,
      'final_branch_min_bound':best['branch_min_bound'],'final_branch_max_bound':best['branch_max_bound'],
      'final_branch_max_mode':'auto_saturation' if args.branch_max is None else 'explicit','branch_multistart':args.branch_multistart,
      'timings':times,
    }
    summary_path=str(out)+'.summary.json';Path(summary_path).write_text(json.dumps(summary,indent=2)+'\n')
    print(f'[result] LP relaxation log score: {lp_score:.10f}')
    print(f'[result] best top-K optimized logL: {best["log_likelihood"]:.10f}')
    print(f'[result] best Newick: {best_path}')
    print(f'[result] top-K table: {topk_path}')
    print(f'[result] timing log: {timing_path}')
    return summary


def make_parser():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('fasta',help='Input FASTA alignment (.fa/.fasta); all sequences must have equal length')
    p.add_argument('--out-prefix',default='phylo_lp_run',help='Output prefix')
    p.add_argument('--library',default=None,help='Optional prebuilt library.json. If omitted, build one using FASTA taxon IDs.')
    p.add_argument('--max-columns',type=int,default=None,help='Randomly sample at most this many alignment columns; default uses all columns')
    p.add_argument('--seed',type=int,default=1,help='Seed for column sampling')
    p.add_argument('--model',choices=['jc','gtr+f'],default='jc',help='DNA substitution model for both LP and final rescoring')
    p.add_argument('--gtr-rates',type=parse_rates,default=None,metavar='AC,AG,AT,CG,CT,GT',help='Six positive GTR exchangeabilities; required with --model gtr+f')
    p.add_argument('--branch-relaxation',choices=['auto','two-point','spectral'],default='auto',
                   help='LP branch-effect relaxation. auto = two-point for JC, 3-mode spectral for non-JC (default).')
    p.add_argument('--spectral-support-points',type=int,default=7,metavar='N',
                   help='Sample orientations used to build the full-range non-JC spectral branch hull (default: 7; minimum 5).')
    p.add_argument('--freq-pseudocount',type=float,default=1e-8,help='Tiny +F pseudocount used only if bases are absent/rare')
    p.add_argument('--preprocess-mode',choices=['fast','balanced','accuracy','oracle'],default='balanced',help='LP preprocessing/root-envelope policy; balanced is default')
    p.add_argument('--log-ratio',type=float,default=1.1,help='Geometric log-chord spacing')
    p.add_argument('--long-ratio',type=float,default=2.0,help='LP long-path / short-path branch-length ratio')
    p.add_argument('--max-quad-facets',type=int,default=None,help='Override mode default for quadratic root facets')
    p.add_argument('--top-k',type=int,default=10,help='Number of non-destructive structural candidates to branch-optimize and rescore')
    p.add_argument('--root-beam',type=int,default=256);p.add_argument('--leaf-beam',type=int,default=8)
    p.add_argument('--child-branching',type=int,default=16);p.add_argument('--leaf-branching',type=int,default=8);p.add_argument('--merge-branching',type=int,default=128)
    p.add_argument('--max-expansions',type=int,default=500000)
    p.add_argument('--branch-min',type=float,default=1e-12,help='Minimum final integral-tree branch length; effectively zero')
    p.add_argument('--branch-max',type=parse_branch_max,default=None,help='Maximum final integral-tree branch length, or auto (default). Auto uses the substitution-model spectral gap so the upper bound is effectively infinite at 1e-12 transition tolerance.')
    p.add_argument('--branch-maxiter',type=int,default=200,help='L-BFGS-B iterations per top-K tree')
    p.add_argument('--branch-ftol',type=float,default=1e-10)
    p.add_argument('--branch-multistart',type=int,default=4,help='Deterministic L-BFGS-B starts per recovered topology (default: 4)')
    return p


def main(argv=None):
    p=make_parser();args=p.parse_args(argv)
    if args.model=='gtr+f' and args.gtr_rates is None:p.error('--gtr-rates is required with --model gtr+f')
    if args.model=='jc' and args.gtr_rates is not None:p.error('--gtr-rates only applies to --model gtr+f')
    run(args)

if __name__=='__main__':main()
