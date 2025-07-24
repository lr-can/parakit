#!/usr/bin/env python3
"""
Build two tables for the CYP2D6/CYP2D7 locus from GFF3 files:

1) hprc.cyp2d6.coords.tsv       # 3 columns  (copy‑level, c1/c2)
2) hprc.cyp2d6.haps.coords.tsv  # 2 columns  (hap‑level union)

Input: cat_cyp/*.cyp2d6.gff3   (#1 and #2 hap files per sample)
 └─ each GFF must contain exactly one CYP2D6 and/or CYP2D7 gene feature

Rules:
  • column 1  = sample‑hap  (e.g. HG00621.1)
  • column 2  = sample#hap#contig:start-end
  • column 3  = c1_…  for CYP2D6   |   c2_…  for CYP2D7
  • hap‑level line spans the min(start)‑max(end) of the two copies
"""

import re
from pathlib import Path

gff_dir      = Path("cat_cyp")
coords_path  = "hprc.cyp2d6.coords.tsv"
haps_path    = "hprc.cyp2d6.haps.coords.tsv"

# RegEx to catch the gene name in the 9th GFF column
rx_name = re.compile(r"Name=(CYP2D[67])\b")

with open(coords_path, "w") as fout_c, open(haps_path, "w") as fout_h:
    for gff in sorted(gff_dir.glob("*.cyp2d6.gff3")):
        sample_hap = gff.stem.replace(".cyp2d6", "")           # HG00438.1
        sample, hap = sample_hap.split(".")                    # HG00438 , 1

        genes = {}  # { 'CYP2D6': (chrom, start, end), 'CYP2D7': (...) }

        with gff.open() as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                cols = line.rstrip().split("\t")
                if len(cols) < 9 or cols[2] != "gene":
                    continue
                m = rx_name.search(cols[8])
                if not m:
                    continue
                gene    = m.group(1)                           # CYP2D6 / CYP2D7
                chrom   = cols[0]
                start   = int(cols[3])
                end     = int(cols[4])
                genes[gene] = (chrom, start, end)

        # --- copy‑level output ------------------------------------------------
        if "CYP2D6" in genes:
            chrom, s, e = genes["CYP2D6"]
            coord = f"{sample}#{hap}#{chrom}:{s}-{e}"
            fout_c.write(f"{sample_hap}\t{coord}\tc1_{sample_hap}\n")

        if "CYP2D7" in genes:
            chrom, s, e = genes["CYP2D7"]
            coord = f"{sample}#{hap}#{chrom}:{s}-{e}"
            fout_c.write(f"{sample_hap}\t{coord}\tc2_{sample_hap}\n")

        # --- hap‑level output (union span) ------------------------------------
        if genes:
            chroms = {v[0] for v in genes.values()}
            if len(chroms) != 1:
                raise ValueError(f"{gff.name}: CYP2D6 and CYP2D7 on different contigs")
            chrom = chroms.pop()
            start = min(v[1] for v in genes.values())
            end   = max(v[2] for v in genes.values())
            coord = f"{sample}#{hap}#{chrom}:{start}-{end}"
            fout_h.write(f"{sample_hap}\t{coord}\n")
