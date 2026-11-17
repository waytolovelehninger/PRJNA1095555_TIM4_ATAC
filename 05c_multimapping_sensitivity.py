#!/usr/bin/env python3
# Input : filtered BAMs of male runs (02), $RESULTS/04_consensus/common_peaks.bed (04); run after `source config.sh`
# Output: $RESULTS/05c_multimapping/counts_{noXS,ASminusXS_gt5}.tsv (runs x peaks, same layout as 04) and
#         fragment_ambiguity.tsv (fragments per run by MAPQ type and smallest AS - XS of the two mates)
# Fragments are counted as in 04 but dropped when either mate has an XS tag (noXS)
# or when either mate has AS - XS <= 5 (ASminusXS_gt5). Peaks are unchanged.
import bisect, collections, os, subprocess

SAM, WORK, OUT = os.environ['SAMTOOLS'], os.environ['WORK'], os.environ['RESULTS'] + '/05c_multimapping'
rows = [l.rstrip('\n').split('\t') for l in open(os.environ['SAMPLES'])]
col = {k: i for i, k in enumerate(rows[0])}
runs = [r[0] for r in rows[1:] if r[col['sex_verified']] == 'Male']
os.makedirs(OUT, exist_ok=True)

starts, ends, index, ids = collections.defaultdict(list), collections.defaultdict(list), collections.defaultdict(list), []
for i, line in enumerate(open(os.environ['RESULTS'] + '/04_consensus/common_peaks.bed')):
    c, s, e, n = line.split()[:4]
    starts[c].append(int(s)); ends[c].append(int(e)); index[c].append(i); ids.append(n)

def ab_bin(d):
    return 'no_XS' if d is None else 'le5' if d <= 5 else '6_15' if d <= 15 else 'gt15'

counts = {k: {} for k in ('noXS', 'ASminusXS_gt5')}
amb = collections.Counter()
for r in runs:
    cnt = {k: [0] * len(ids) for k in counts}
    pending = {}
    view = subprocess.Popen([SAM, 'view', '-F', '2304', f'{WORK}/{r}/filtered/{r}.primary.filtered.bam'],
                            stdout=subprocess.PIPE, text=True)
    for line in view.stdout:
        x = line.rstrip('\n').split('\t')
        t = dict((f[:2], int(f[5:])) for f in x[11:] if f[:2] in ('AS', 'XS'))
        rec = (int(x[1]), x[2], int(x[3]), int(x[7]), int(x[8]), int(x[4]), t['AS'] - t['XS'] if 'XS' in t else None)
        m = pending.pop(x[0], None)
        if m is None: pending[x[0]] = rec; continue
        a = rec if rec[0] & 64 else m                               # read 1 defines the fragment (as in 04)
        ds = [v[6] for v in (rec, m) if v[6] is not None]
        d = min(ds) if ds else None
        amb[(r, 'MAPQ255' if a[5] == 255 else 'numeric', ab_bin(d))] += 1
        c, s = a[1], min(a[2], a[3]) - 1
        e = s + abs(a[4])
        for j in range(bisect.bisect_right(ends[c], s), bisect.bisect_left(starts[c], e)) if c in ends else ():
            if d is None: cnt['noXS'][index[c][j]] += 1
            if d is None or d > 5: cnt['ASminusXS_gt5'][index[c][j]] += 1
    if view.wait() or pending: raise SystemExit(f'samtools failed or unpaired records for {r}')
    for k in counts: counts[k][r] = cnt[k]

for k, m in counts.items():
    with open(f'{OUT}/counts_{k}.tsv', 'w', newline='\n') as f:
        f.write('sample\t' + '\t'.join(ids) + '\n')
        f.writelines(r + '\t' + '\t'.join(map(str, m[r])) + '\n' for r in runs)
with open(f'{OUT}/fragment_ambiguity.tsv', 'w', newline='\n') as f:
    f.write('run\tmapq\tmin_AS_minus_XS\tfragments\n')
    f.writelines(f'{r}\t{q}\t{b}\t{v}\n' for (r, q, b), v in sorted(amb.items()))
