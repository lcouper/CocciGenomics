#!/usr/bin/env python3
"""Bootstrap SNP columns within each fixed-size window of a .geno file.

For each contiguous block of W sites within a scaffold, draw W sites WITH
replacement from that same block. Output rows carry the block's original
sorted positions (in order), so positions stay monotonic and windows line
up 1:1 with the unbootstrapped run.
"""
import argparse, gzip, random

p = argparse.ArgumentParser()
p.add_argument("-i", "--input",  required=True)
p.add_argument("-o", "--output", required=True)
p.add_argument("-w", "--window", type=int, default=100)
p.add_argument("-s", "--seed",   type=int, required=True)
a = p.parse_args()

rng = random.Random(a.seed)

def emit(block, out):
    n = len(block)
    pos = [r[1] for r in block]                       # original positions, sorted
    picks = [block[rng.randrange(n)] for _ in range(n)]
    for newpos, row in zip(pos, picks):
        out.write("\t".join([row[0], newpos] + row[2:]) + "\n")

with gzip.open(a.input, "rt") as fh, gzip.open(a.output, "wt") as out:
    first = fh.readline()
    carry = None
    if first.startswith("#"):
        out.write(first)                              # pass the header through
    elif first.strip():
        carry = first.rstrip("\n").split("\t")

    block, cur = ([carry], carry[0]) if carry else ([], None)

    for line in fh:
        line = line.rstrip("\n")
        if not line:
            continue
        f = line.split("\t")
        if f[0] != cur or len(block) == a.window:      # new scaffold, or window full
            if block:
                emit(block, out)
            block, cur = [], f[0]
        block.append(f)
    if block:
        emit(block, out)
