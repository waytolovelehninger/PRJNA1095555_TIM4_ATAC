#!/usr/bin/env python3
# Input : GENCODE vM25 GTF (.gtf.gz), mm10 chromosome sizes (2 columns: name, length; e.g. UCSC mm10.chrom.sizes)
# Output: TSS BED (stdout), one row per unique protein-coding TSS
# Rule (re-implementation of prepare_tss() in the pilot QC script): GTF 'gene' features with
# gene_type "protein_coding" on chr1-19/X/Y; TSS = gene start (+ strand) or gene end - 1 (- strand), 0-based;
# drop a TSS whose +-2 kb window leaves the chromosome; genes sharing chrom/position/strand are collapsed
# into one row (gene IDs joined by ';'). Rows sorted by (chrom as string, position, strand).
import collections, gzip, re, sys

gtf, sizes = sys.argv[1:3]
CANON = {f'chr{i}' for i in range(1, 20)} | {'chrX', 'chrY'}
length = dict((l.split()[0], int(l.split()[1])) for l in open(sizes))
sites = collections.defaultdict(list)
with gzip.open(gtf, 'rt') as f:
    for line in f:
        if line.startswith('#'): continue
        x = line.rstrip('\n').split('\t')
        if x[2] != 'gene' or x[0] not in CANON: continue
        a = dict(re.findall(r'(\w+) "([^"]*)"', x[8]))
        if a.get('gene_type') != 'protein_coding': continue
        s, e, strand = int(x[3]) - 1, int(x[4]), x[6]
        tss = s if strand == '+' else e - 1
        if tss - 2000 < 0 or tss + 2000 >= length[x[0]]: continue
        sites[(x[0], tss, strand)].append(a['gene_id'])
for (c, p, strand), genes in sorted(sites.items()):
    sys.stdout.write(f'{c}\t{p}\t{p + 1}\t{";".join(sorted(genes))}\t0\t{strand}\n')
