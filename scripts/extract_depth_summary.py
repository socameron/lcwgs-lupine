#!/usr/bin/env python3
import argparse
import gzip
import re
import sys
from pathlib import Path
import numpy as np
import pandas as pd

POP_IN_PATH = re.compile(r"/all_by_popln/([^/]+)/", flags=re.IGNORECASE)

def infer_pop(filepath: str) -> str:
    """Infer population from path or filename."""
    m = POP_IN_PATH.search(filepath)
    if m:
        return m.group(1)
    # Fallback: basename starts with POP_...
    base = Path(filepath).name
    return base.split("_", 1)[0]

def iter_depth_chunks(path, chunksize=1_000_000):
    """Yield numpy arrays of depth values from a gzipped TSV."""
    with gzip.open(path, "rt") as f:
        head = f.readline()
        f.seek(0)
        has_header = ("depth" in head.split())
        if has_header:
            usecols, header = ["depth"], 0
        else:
            first = f.readline().rstrip("\n")
            ncols = len(first.split("\t"))
            f.seek(0)
            usecols, header = [ncols - 1], None

        for chunk in pd.read_csv(
            f, sep="\t", chunksize=chunksize, usecols=usecols,
            dtype="int32", header=header
        ):
            arr = chunk.iloc[:, 0].to_numpy(copy=False)
            if arr.size:
                yield arr[arr >= 0]

def add_to_hist(hist, arr):
    """Accumulate bincounts of arr into hist (auto-resizing)."""
    if arr.size == 0:
        return hist
    maxv = int(arr.max())
    if hist is None:
        hist = np.zeros(maxv + 1, dtype=np.int64)
    elif maxv >= hist.size:
        new = np.zeros(maxv + 1, dtype=np.int64)
        new[: hist.size] = hist
        hist = new
    hist += np.bincount(arr, minlength=hist.size)
    return hist

def quantiles_from_hist(hist, qs=(5, 10, 98)):
    """Return quantiles dict, total sites, mean, median."""
    if hist is None or hist.sum() == 0:
        return {q: np.nan for q in qs}, 0, np.nan, np.nan
    total = int(hist.sum())
    cdf = np.cumsum(hist)
    mean = (np.arange(hist.size, dtype=np.float64) * hist).sum() / total
    median_idx = int(np.searchsorted(cdf, (total - 1) / 2))
    median = float(median_idx)
    out = {}
    for q in qs:
        k = max(1, int(np.ceil(q / 100.0 * total)))
        idx = int(np.searchsorted(cdf, k))
        out[q] = float(idx)
    return out, total, mean, median

def main():
    ap = argparse.ArgumentParser(
        description="Summarize pooled depth cutoffs (p05, p10, p98) per population into one CSV."
    )
    ap.add_argument("--out", required=True, help="Output CSV path.")
    ap.add_argument("files", nargs="+", help="All per-scaffold depth files (*.txt.gz) for all populations.")
    args = ap.parse_args()

    # histograms per population
    pop_hist = {}
    pop_scaf_counts = {}

    for fp in args.files:
        pop = infer_pop(fp)
        pop_scaf_counts[pop] = pop_scaf_counts.get(pop, 0) + 1
        for arr in iter_depth_chunks(fp):
            pop_hist[pop] = add_to_hist(pop_hist.get(pop), arr)

    rows = []
    for pop, hist in sorted(pop_hist.items()):
        qvals, total, mean, median = quantiles_from_hist(hist, qs=(5, 10, 98))
        rows.append({
            "population": pop,
            "n_scaffolds": pop_scaf_counts.get(pop, 0),
            "n_sites": total,
            "mean_depth": mean,
            "median_depth": median,
            "p05": qvals[5],
            "p10": qvals[10],
            "p98": qvals[98],
        })

    df = pd.DataFrame(rows, columns=[
        "population","n_scaffolds","n_sites","mean_depth","median_depth","p05","p10","p98"
    ]).sort_values("population")

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, index=False)

if __name__ == "__main__":
    sys.exit(main())
