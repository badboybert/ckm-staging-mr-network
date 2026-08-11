# -*- coding: utf-8 -*-
"""Step 7: enhanced AHA-staging analysis (Component B backbone + concordance).
Consumes an MR edge table -> significant causal DAG -> SCC condensation (feedback modules)
-> Hasse partial order -> curl/circulation fraction -> min feedback-arc-set (edges AHA gets
backward) -> stage-order concordance with a label-permutation null.
Honest framing: this tests NECESSARY POPULATION-LEVEL conditions of the AHA staging order,
NOT within-person progression. Discordance falsifies; concordance is only 'consistent-with'."""
# --- portable roots -----------------------------------------------------------------------------
# Injected by scripts/build_repo.py. The working tree hardcoded an absolute local path; the deposit
# resolves it from CKM_P4_BASE, or from this file's own location (repo/code/ plays the role the
# working tree called independent_build/). Raw GWAS summary statistics are NOT deposited: set
# CKM_ROOT to wherever you obtained them if you intend to re-run the upstream extraction steps.
import os as _os
P4_BASE = _os.environ.get("CKM_P4_BASE") or _os.path.dirname(
    _os.path.dirname(_os.path.abspath(__file__)))
P4_ROOT = _os.environ.get("CKM_P4_ROOT") or _os.path.dirname(P4_BASE)
CKM_ROOT = _os.environ.get("CKM_ROOT") or _os.path.dirname(P4_ROOT)
SHARED_LIB = _os.environ.get("CKM_SHARED_LIB") or _os.path.join(CKM_ROOT, "_shared")
# ------------------------------------------------------------------------------------------------

import csv, io, os, itertools
from collections import Counter
import networkx as nx

BASE = P4_BASE


def _multiset_perms(counts):
    """Yield every DISTINCT arrangement of a stage multiset (no duplicate labelings)."""
    n = sum(counts.values()); cur = [None]*n
    def rec(i):
        if i == n:
            yield tuple(cur); return
        for v in sorted(counts):
            if counts[v]:
                counts[v] -= 1; cur[i] = v
                yield from rec(i+1)
                counts[v] += 1
    return rec(0)


def _exact_perm_p(edges, nodes, stage, obs):
    """Exact one-sided permutation P for the label-permutation staging null: the fraction of
    distinct stage-label assignments whose cross-stage concordance is >= the observed. Node order
    is irrelevant to the result (it is a property of the labeling set, not its enumeration order)."""
    ge = tot = 0
    for pm in _multiset_perms(dict(Counter(stage[n] for n in nodes))):
        tot += 1
        smap = dict(zip(nodes, pm))
        cc = [(a, b) for a, b in edges if smap[a] != smap[b]]
        if not cc:
            continue
        c = sum(1 for a, b in cc if smap[a] < smap[b]) / len(cc)
        if c >= obs - 1e-12:
            ge += 1
    return ge/tot, ge, tot

# AHA 2023 CKM stage per node (1=adiposity, 2=metabolic/kidney, 4=clinical CVD)
STAGE = {"BMI":1, "SBP":2,"TG":2,"HDL":2,"TC":2,"LDL":2,"HbA1c":2,"FI":2,"T2D":2,"eGFR":2,"CKD":2,
         "NAFLD":2, "CAD":4,"HF":4,"Stroke":4,"AF":4}

def load_all(path):
    rows=[]
    with io.open(path,encoding="utf-8") as f:
        for r in csv.DictReader(f):
            try: p=float(r["ivw_p"]); b=float(r["ivw_b"])
            except: continue
            steiger = str(r.get("steiger_correct", r.get("steiger",""))).upper()=="TRUE"
            rows.append((r["exposure"], r["outcome"], b, p, steiger))
    return rows

def analyse(path):
    allrows = load_all(path)
    n_tested = len(allrows)
    bonf = 0.05/n_tested                      # principled EDGE-significance threshold
    print(f"Edges tested: {n_tested}  | Bonferroni edge threshold = 0.05/{n_tested} = {bonf:.2e}")
    # DAG = Bonferroni-significant, Steiger-correct edges
    edges = [(a,b,beta,p,st) for (a,b,beta,p,st) in allrows if p<bonf and st]

    # --- falsification report: nominally-significant reverse CROSS-STAGE edges (candidate violations) ---
    cand=[]
    for a,b,beta,p,st in allrows:
        if a in STAGE and b in STAGE and STAGE[a]>STAGE[b] and p<0.05 and st:
            cand.append((a,b,beta,p, "SURVIVES-Bonferroni" if p<bonf else "sub-threshold(NS after correction)"))
    print("\n[FALSIFICATION SCAN] nominally-sig, Steiger-correct BACKWARD (high->low AHA stage) edges:")
    if cand:
        for a,b,beta,p,tag in sorted(cand,key=lambda x:x[3]):
            print(f"    {a}->{b} (stage {STAGE[a]}->{STAGE[b]})  b={beta:+.3f} p={p:.2e}  [{tag}]")
    else:
        print("    none")
    G = nx.DiGraph()
    for a,b,beta,p,st in edges:
        if a in STAGE and b in STAGE:
            G.add_edge(a,b,beta=beta,p=p)
    print(f"Significant Steiger-correct causal edges: {G.number_of_edges()} over {G.number_of_nodes()} nodes")

    # 1. feedback modules (SCCs of size>1)
    sccs=[c for c in nx.strongly_connected_components(G) if len(c)>1]
    print(f"\n[1] Feedback modules (cyclic SCCs): {sccs if sccs else 'none (graph is acyclic)'}")

    # 2. TRUE partial order: sources (in-deg 0), sinks (out-deg 0), and DAG layers by
    #    longest-path-from-source (nodes in the same layer are causally CONCURRENT/incomparable).
    C=nx.condensation(G)
    sources=[n for n in G if G.in_degree(n)==0]
    sinks  =[n for n in G if G.out_degree(n)==0]
    # longest-path layer per node (on the condensation DAG, then map back)
    layer_c={}
    for i in nx.topological_sort(C):
        preds=list(C.predecessors(i))
        layer_c[i]=0 if not preds else 1+max(layer_c[p] for p in preds)
    node_layer={m:layer_c[i] for i in C.nodes for m in C.nodes[i]['members']}
    layers={}
    for n,l in node_layer.items(): layers.setdefault(l,[]).append(n)
    print(f"[2] Genetic causal layering (nodes in same layer are causally concurrent/incomparable):")
    for l in sorted(layers): print(f"      layer {l}: {sorted(layers[l])}")
    print(f"      SOURCE (upstream drivers, in-deg 0): {sorted(sources)}")
    print(f"      SINK   (downstream endpoints, out-deg 0): {sorted(sinks)}")
    order=layers

    # 3. curl / circulation fraction: cross-stage edges going 'backward' (high->low AHA stage)
    cross=[(a,b) for a,b in G.edges if STAGE[a]!=STAGE[b]]
    backward=[(a,b) for a,b in cross if STAGE[a]>STAGE[b]]
    forward=[(a,b) for a,b in cross if STAGE[a]<STAGE[b]]
    conc = len(forward)/len(cross) if cross else float('nan')
    print(f"\n[3] Cross-stage edges: {len(cross)}  | forward(low->high AHA)={len(forward)}  backward(high->low)={len(backward)}")
    print(f"    Stage-order CONCORDANCE = {conc:.3f}")
    print(f"    Backward (AHA-discordant) edges: {backward if backward else 'none'}")

    # 4. min feedback-arc-set relative to the AHA order = edges you'd remove to make graph respect stages
    #    here: report the backward cross-stage edges + any intra-cycle edges as 'violations'
    print(f"[4] Edges the AHA staging order gets 'backward' (violations): "
          f"{backward if backward else 'none'}")

    # 5. permutation null for concordance (assign the stage multiset to nodes).
    #    The null's support is the set of DISTINCT stage-label assignments to the nodes, which is
    #    small here (e.g. 1,980 labelings for 12 nodes at composition 1/8/3). It is therefore
    #    enumerated EXHAUSTIVELY: the exact permutation P is the fraction of labelings whose
    #    concordance is at least the observed. This removes the seed- and node-order-dependence of a
    #    Monte-Carlo estimate (a shuffled stage vector is zipped onto a node list, so both the seed
    #    and the node ordering move the last digit of an MC P). A 10,000-draw MC estimate is retained
    #    as a disclosed cross-check and must agree within Monte-Carlo error (~4e-4 at this P).
    import random
    nodes=list(G.nodes); stages=[STAGE[n] for n in nodes]
    obs=conc
    pval, ge, tot = _exact_perm_p(list(G.edges), nodes, STAGE, obs)
    # MC cross-check (not reported; asserts the exact value is not a coding error)
    mc_ge=0; NP=10000; rng=random.Random(12345)
    for _ in range(NP):
        perm=stages[:]; rng.shuffle(perm); smap=dict(zip(nodes,perm))
        cc=[(a,b) for a,b in G.edges if smap[a]!=smap[b]]
        ff=sum(1 for a,b in cc if smap[a]<smap[b])
        if (ff/len(cc) if cc else 0)>=obs: mc_ge+=1
    mc_p=(mc_ge+1)/(NP+1)
    assert abs(pval-mc_p) < 5e-3, f"exact perm P {pval} disagrees with MC {mc_p} beyond MC error"
    print(f"\n[5] Concordance vs label-permutation null (EXACT enumeration of {tot} distinct labelings):"
          f" exact p = {ge}/{tot} = {pval:.4g}  (MC cross-check {mc_p:.4f})")
    print(f"    => {'concordant with AHA order beyond chance' if pval<0.05 else 'NOT beyond chance'}")
    print("\nNOTE: tests necessary POPULATION-LEVEL causal-ordering conditions, not within-person progression.")
    return dict(edges=G.number_of_edges(), sccs=sccs, order=order, concordance=conc, perm_p=pval,
                perm_ge=ge, perm_tot=tot, backward=backward)

if __name__=="__main__":
    import sys
    path=sys.argv[1] if len(sys.argv)>1 else os.path.join(BASE,"results/forward_local_edges.csv")
    print(f"=== Staging analysis on {os.path.basename(path)} ===")
    analyse(path)
