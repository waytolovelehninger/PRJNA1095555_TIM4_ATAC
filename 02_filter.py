#!/usr/bin/env python3
# Input : $WORK/<run>/work/markdup.bam (01), $BLACKLIST (ENCODE mm10 blacklist v2); run after `source config.sh`
# Output: $WORK/<run>/filtered/<run>.primary.filtered.bam (+ .bai); fragments left after each rule printed to stderr
import bisect, collections, os, re, subprocess, sys

SAM, WORK = os.environ['SAMTOOLS'], os.environ['WORK']
CANON = {f'chr{i}' for i in range(1, 20)} | {'chrX', 'chrY'}
STEPS = ['primary_pair', 'proper_pair', 'not_qcfail', 'not_duplicate', 'mapq', 'not_chrM', 'canonical', 'not_blacklist']
CIGAR = re.compile(r'(\d+)([MIDNSHP=X])')

blacklist = collections.defaultdict(list)          # merged half-open intervals per chromosome
for line in open(os.environ['BLACKLIST']):
    c, s, e = line.split()[:3]
    blacklist[c].append([int(s), int(e)])
for c, iv in blacklist.items():
    merged = []
    for s, e in sorted(iv):
        if merged and s <= merged[-1][1]: merged[-1][1] = max(e, merged[-1][1])
        else: merged.append([s, e])
    blacklist[c] = ([s for s, _ in merged], [e for _, e in merged])

def in_blacklist(c, s, e):
    starts, ends = blacklist.get(c, ([], []))
    i = bisect.bisect_left(starts, e) - 1
    return i >= 0 and ends[i] > s

def ref_end(pos, cigar):
    parts = CIGAR.findall(cigar)
    if not parts or ''.join(n + op for n, op in parts) != cigar or any(int(n) == 0 for n, _ in parts): return -1
    return pos + sum(int(n) for n, op in parts if op in 'MDN=X')

def rules_passed(pair, chrom_len):
    """Number of consecutive rules (STEPS) passed by one read name's primary records."""
    if len(pair) != 2: return 0
    (af, ac, ap, aq, ag, ar, an, at), (bf, bc, bp, bq, bg, br, bn, bt) = [
        (int(f[1]), f[2], int(f[3]) - 1, int(f[4]), f[5], f[2] if f[6] == '=' else f[6], int(f[7]) - 1, int(f[8]))
        for f in (r.split('\t') for r in pair)]
    if (af & 192, bf & 192) not in ((64, 128), (128, 64)): return 0
    ae, be = ref_end(ap, ag), ref_end(bp, bg)
    s, e = min(ap, bp), max(ae, be)
    proper = (af & bf & 1 and not (af | bf) & 4 and af & bf & 2 and ac == bc == ar == br and an == bp and bn == ap
              and bool(af & 8) == bool(bf & 4) and bool(bf & 8) == bool(af & 4)
              and bool(af & 32) == bool(bf & 16) and bool(bf & 32) == bool(af & 16) and bool(af & 16) != bool(bf & 16)
              and ae > ap and be > bp and 0 <= s < e <= chrom_len.get(ac, -1) and at == -bt and abs(at) == e - s)
    checks = [proper, not (af | bf) & 512, not (af | bf) & 1024,
              all(q == 255 or 30 <= q <= 254 for q in (aq, bq)), ac != 'chrM', ac in CANON, not in_blacklist(ac, s, e)]
    return 1 + next((i for i, ok in enumerate(checks) if not ok), len(checks))

def filter_run(run):
    inp, out = f'{WORK}/{run}/work/markdup.bam', f'{WORK}/{run}/filtered/{run}.primary.filtered.bam'
    tmp = f'{WORK}/{run}/work/tmp'
    os.makedirs(os.path.dirname(out), exist_ok=True); os.makedirs(tmp, exist_ok=True)
    collate = subprocess.Popen([SAM, 'collate', '-@', '1', '-O', '--output-fmt', 'SAM', inp, f'{tmp}/filter_collate'],
                               stdout=subprocess.PIPE, text=True)
    sort = subprocess.Popen([SAM, 'sort', '-@', '2', '-m', '768M', '-T', f'{tmp}/final_sort', '-O', 'BAM', '-o', out, '-'],
                            stdin=subprocess.PIPE, text=True)
    chrom_len, kept, name, pair = {}, collections.Counter(), None, []
    def finish():
        n = rules_passed(pair, chrom_len)
        for step in STEPS[:n]: kept[step] += 1
        if n == len(STEPS): sort.stdin.writelines(pair)
    for line in collate.stdout:
        if line.startswith('@'):
            if line.startswith('@SQ'):
                tags = dict(t.split(':', 1) for t in line.rstrip('\n').split('\t')[1:])
                chrom_len[tags['SN']] = int(tags['LN'])
            sort.stdin.write(line); continue
        qname, flag = line.split('\t', 2)[:2]
        if qname != name:
            if name is not None: finish()
            name, pair = qname, []
        if not int(flag) & 2304: pair.append(line)    # primary records only
    if name is not None: finish()
    sort.stdin.close()
    if collate.wait() or sort.wait(): sys.exit(f'samtools failed for {run}')
    subprocess.run([SAM, 'index', '-@', '2', '-b', out], check=True)
    print(run, *(f'{s}={kept[s]}' for s in STEPS), sep='\t', file=sys.stderr)

for run in (os.environ.get('RUNS') or '').split() or [l.split('\t')[0] for l in open(os.environ['SAMPLES']).readlines()[1:]]:
    filter_run(run)
