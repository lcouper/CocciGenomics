#!/usr/bin/env python3
"""Collapse twisst weights for one bootstrap replicate. Works for any reference-group spec.

Replaces the K=2-specific inline block in twisst_boot.sbatch and k3b_summarise.py.
Bin naming. Each topology is labelled by which reference groups are sister to Focal:
  one group           -> that group's short name
  every group         -> own_branch
  a proper subset     -> "pair" + the short names joined
Short name is the group name up to the first "_", so C1_KernRiver -> C1 while
Bakersfield stays Bakersfield. That reproduces the bin labels both the K=2 and K=3
downstream readers already expect.

  python3 twisst_summarise.py -r 7 -g groups/refs_K3b.tsv \
      -s summaries/rep7.csv -w summaries/rep7.windows.csv.gz w_*.csv.gz
"""
import argparse, gzip, csv, os
from ete3 import Tree

p = argparse.ArgumentParser()
p.add_argument("-r", "--rep", required=True)
p.add_argument("-g", "--groups", required=True, help="refs tsv used for this spec")
p.add_argument("-s", "--summary", required=True)
p.add_argument("-w", "--windows", required=True)
p.add_argument("weights", nargs="+")
a = p.parse_args()

refs = sorted({f[1] for f in (l.split() for l in open(a.groups)) if len(f) == 2 and f[1] != "Outgroup"})
if not refs:
    raise SystemExit("no reference groups parsed from %s" % a.groups)
short = {g: g.split("_")[0] for g in refs}


def classify(nwk):
    t = Tree(nwk.strip())
    f = t & "Focal"
    names = []
    for c in f.up.children:
        if c is f:
            continue
        names += [c.name] if c.is_leaf() else c.get_leaf_names()
    cl = sorted(short[n] for n in names if n in short)
    if len(cl) == 0:
        raise SystemExit("topology has no reference group sister to Focal -- "
                         "group names in %s do not match the tree: %s" % (a.groups, nwk.strip()))
    if len(cl) == 1:
        return cl[0]
    if len(cl) == len(refs):
        return "own_branch"
    return "pair" + "".join(cl)


def read_topologies(fh):
    """Leading '[#]topoN <newick>' lines, then the tab-separated column header.

    Detected by newick content rather than a leading '#', since twisst writes the
    marker inconsistently across versions.
    """
    topos = []
    for line in fh:
        if "(" in line and ";" in line:
            topos.append(line.split(None, 1)[1])
        else:
            break                                  # this line is the column header
    if not topos:
        raise SystemExit("no topology lines found -- unexpected weights file header")
    return topos


sf = open(a.summary, "w", newline="")
wf = gzip.open(a.windows, "wt", newline="")
sw, ww = csv.writer(sf), csv.writer(wf)

nt = None
for path in sorted(a.weights):
    iso = os.path.basename(path)[2:-7]                 # w_<isolate>.csv.gz
    with gzip.open(path, "rt") as fh:
        tb = [classify(t) for t in read_topologies(fh)]
        if nt is None:                                 # header depends on topology count
            nt = len(tb)
            sw.writerow(["rep", "isolate", "topo", "bin", "weight", "n_windows",
                         "modal_topo", "modal_topo_frac", "modal_bin", "modal_bin_frac"])
            ww.writerow(["rep", "isolate", "window"] + ["topo%d" % (k + 1) for k in range(nt)])
        elif len(tb) != nt:
            raise SystemExit("%s has %d topologies, expected %d" % (path, len(tb), nt))

        bins = sorted(set(tb))
        acc, domT = [0.0] * nt, [0] * nt
        domB = dict.fromkeys(bins, 0)
        wins = 0
        for i, line in enumerate(fh, start=1):
            try:
                v = [float(x) for x in line.split()]
            except ValueError:
                ww.writerow([a.rep, iso, i] + [""] * nt)
                continue
            s = sum(v)
            if s <= 0 or s != s:
                ww.writerow([a.rep, iso, i] + [""] * nt)
                continue
            v = [x / s for x in v]
            for k in range(nt):
                acc[k] += v[k]
            domT[v.index(max(v))] += 1
            b = dict.fromkeys(bins, 0.0)
            for k in range(nt):
                b[tb[k]] += v[k]
            domB[max(b, key=b.get)] += 1
            wins += 1
            ww.writerow([a.rep, iso, i] + ["%.6f" % x for x in v])

    mT = domT.index(max(domT)) + 1 if wins else ""
    mTf = max(domT) / wins if wins else float("nan")
    mB = max(domB, key=domB.get) if wins else ""
    mBf = domB[mB] / wins if wins else float("nan")
    for k in range(nt):
        sw.writerow([a.rep, iso, k + 1, tb[k],
                     "%.6f" % (acc[k] / wins) if wins else "", wins,
                     mT, "%.6f" % mTf, mB, "%.6f" % mBf])

sf.close()
wf.close()
