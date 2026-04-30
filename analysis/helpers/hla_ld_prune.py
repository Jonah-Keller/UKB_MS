#!/usr/bin/env python3
"""
Compute pairwise r² between all FDR<0.05 HLA alleles and LD-prune to
independent signals.  Outputs:
  data/ukb/genetics/hla_ld_clusters.csv   — every allele + cluster ID + keep flag
  data/ukb/genetics/hla_allele_enrichment_pruned.csv — pruned enrichment table
"""
import csv, math, sys
from pathlib import Path

PROJ   = Path(__file__).resolve().parents[2]
RAW    = PROJ / "app178779_20241219160805.dataset.csv"
ENR    = PROJ / "data/ukb/genetics/hla_allele_enrichment.csv"
OUT_LD = PROJ / "data/ukb/genetics/hla_ld_clusters.csv"
OUT_PR = PROJ / "data/ukb/genetics/hla_allele_enrichment_pruned.csv"

R2_THRESHOLD = 0.64   # r > 0.8 → same haplotype

# ── load enrichment table ───────────────────────────────────────────────────
enr_rows = []
with open(ENR) as f:
    reader = csv.DictReader(f)
    for row in reader:
        enr_rows.append(row)

sig_idx = [int(r["idx"]) for r in enr_rows if float(r["fdr"]) < 0.05]
print(f"Significant alleles (FDR<0.05): {len(sig_idx)}")
print(f"Indices: {sorted(sig_idx)}")

# ── stream full HIBAG array, extract only needed indices ────────────────────
print("Streaming dataset...", flush=True)
dosages = {i: [] for i in sig_idx}
n_rows  = 0

with open(RAW) as f:
    reader = csv.DictReader(f)
    for row in reader:
        arr_str = row["participant.p22182"].strip('"')
        vals    = arr_str.split(",")
        for i in sig_idx:
            try:
                dosages[i].append(float(vals[i]))
            except (IndexError, ValueError):
                dosages[i].append(float("nan"))
        n_rows += 1
        if n_rows % 100_000 == 0:
            print(f"  {n_rows:,} rows read", flush=True)

print(f"Done — {n_rows:,} participants")

# ── compute pairwise r² (skip NaN pairs) ───────────────────────────────────
def pearson_r(xs, ys):
    pairs = [(x, y) for x, y in zip(xs, ys)
             if not (math.isnan(x) or math.isnan(y))]
    n = len(pairs)
    if n < 100:
        return float("nan")
    mx = sum(p[0] for p in pairs) / n
    my = sum(p[1] for p in pairs) / n
    num = sum((p[0]-mx)*(p[1]-my) for p in pairs)
    vx  = sum((p[0]-mx)**2 for p in pairs)
    vy  = sum((p[1]-my)**2 for p in pairs)
    denom = math.sqrt(vx * vy)
    return num / denom if denom > 0 else float("nan")

print("Computing pairwise r²...", flush=True)
r2_mat = {}
idxs   = sorted(sig_idx)
for i, a in enumerate(idxs):
    for b in idxs[i+1:]:
        r  = pearson_r(dosages[a], dosages[b])
        r2 = r*r if not math.isnan(r) else float("nan")
        r2_mat[(a, b)] = r2

# ── greedy LD clustering (sort by abs enrichment, most significant first) ──
enr_map = {int(r["idx"]): r for r in enr_rows}

def abs_enr(idx):
    return abs(float(enr_map[idx]["enrichment"])) if idx in enr_map else 0

ordered = sorted(sig_idx, key=abs_enr, reverse=True)

cluster_id = {}
cluster    = 0
keep       = {}

for a in ordered:
    if a in cluster_id:
        continue
    cluster_id[a] = cluster
    keep[a]        = True
    # find all alleles in LD with a
    for b in ordered:
        if b in cluster_id:
            continue
        pair = (min(a,b), max(a,b))
        r2   = r2_mat.get(pair, float("nan"))
        if not math.isnan(r2) and r2 >= R2_THRESHOLD:
            cluster_id[b] = cluster
            keep[b]        = False
    cluster += 1

# ── print summary ─────────────────────────────────────────────────────────
print("\n=== LD clusters (r² threshold =", R2_THRESHOLD, ") ===")
from collections import defaultdict
by_cluster = defaultdict(list)
for idx, cid in cluster_id.items():
    by_cluster[cid].append(idx)

for cid, members in sorted(by_cluster.items()):
    canonical = [i for i in members if keep[i]][0]
    others    = [i for i in members if not keep[i]]
    canon_row = enr_map.get(canonical, {})
    print(f"\nCluster {cid}  [canonical idx {canonical} = {canon_row.get('label','?')} "
          f"enr={float(canon_row.get('enrichment',0)):+.1f}%  "
          f"gene={canon_row.get('gene_inferred','?')}]")
    for m in members:
        row = enr_map.get(m, {})
        flag = "KEEP" if keep[m] else "DROP"
        print(f"  {flag}  idx {m:3d}  freq={float(row.get('allele_freq',0)):.3f}  "
              f"ms={float(row.get('ms_carrier',0)):.1f}%  "
              f"hc={float(row.get('hc_carrier',0)):.1f}%  "
              f"enr={float(row.get('enrichment',0)):+.1f}%  "
              f"label={row.get('label','')}")

# ── write outputs ─────────────────────────────────────────────────────────
print(f"\nKept: {sum(keep.values())} / {len(sig_idx)} alleles")

with open(OUT_LD, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["idx","cluster_id","keep","label","gene","allele_freq",
                "enrichment","fdr","ms_carrier","hc_carrier"])
    for idx in sorted(cluster_id):
        row = enr_map.get(idx, {})
        # infer gene
        if idx <= 17:    gene = "HLA-A"
        elif idx <= 51:  gene = "HLA-B"
        elif idx <= 65:  gene = "HLA-C"
        elif idx <= 100: gene = "HLA-DP/DQA"
        elif idx <= 230: gene = "HLA-DQB1"
        else:            gene = "HLA-DRB1"
        w.writerow([idx, cluster_id[idx], int(keep[idx]),
                    row.get("label",""), gene,
                    row.get("allele_freq",""), row.get("enrichment",""),
                    row.get("fdr",""), row.get("ms_carrier",""),
                    row.get("hc_carrier","")])

# pruned enrichment = all rows, with keep column merged in
with open(ENR) as fin, open(OUT_PR, "w", newline="") as fout:
    reader = csv.DictReader(fin)
    fields = reader.fieldnames + ["keep_ld"]
    writer = csv.DictWriter(fout, fieldnames=fields)
    writer.writeheader()
    for row in reader:
        idx = int(row["idx"])
        row["keep_ld"] = keep.get(idx, None)  # None for non-sig alleles
        writer.writerow(row)

print(f"Wrote {OUT_LD}")
print(f"Wrote {OUT_PR}")
