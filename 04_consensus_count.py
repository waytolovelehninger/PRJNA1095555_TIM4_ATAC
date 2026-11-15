#!/usr/bin/env python3
# Input : $WORK/<run>/peaks/<run>_peaks.broadPeak and filtered BAMs for all runs in $SAMPLES; run after `source config.sh`
# Output: $RESULTS/04_consensus/common_peaks.bed and sample_by_common_peak.raw_fragment_counts.tsv (runs x peaks)
import bisect, collections, os, subprocess

SAM, WORK, OUT = os.environ['SAMTOOLS'], os.environ['WORK'], os.environ['RESULTS'] + '/04_consensus'
runs = [l.split('\t')[0] for l in open(os.environ['SAMPLES']).readlines()[1:]]
bam = lambda r: f'{WORK}/{r}/filtered/{r}.primary.filtered.bam'
os.makedirs(OUT, exist_ok=True)

# union of overlapping or adjacent peaks (0-based half-open) from all runs, chromosomes in BAM header order
peaks = collections.defaultdict(list)
for r in runs:
    for line in open(f'{WORK}/{r}/peaks/{r}_peaks.broadPeak'):
        c, s, e = line.split('\t')[:3]
        peaks[c].append((int(s), int(e)))
header = subprocess.run([SAM, 'view', '-H', bam(runs[0])], capture_output=True, text=True, check=True).stdout
chroms = [l.split('\t')[1][3:] for l in header.splitlines() if l.startswith('@SQ')]
common = []
for c in chroms:
    merged = []
    for s, e in sorted(peaks[c]):
        if merged and s <= merged[-1][1]: merged[-1][1] = max(merged[-1][1], e)
        else: merged.append([s, e])
    common += [(c, s, e) for s, e in merged]
ids = [f'common_{i + 1:07d}' for i in range(len(common))]
with open(f'{OUT}/common_peaks.bed', 'w', newline='\n') as f:
    f.writelines(f'{c}\t{s}\t{e}\t{i}\n' for (c, s, e), i in zip(common, ids))

# each fragment (read 1 of a retained pair) counted once in every peak it overlaps by >= 1 bp
starts, ends, index = collections.defaultdict(list), collections.defaultdict(list), collections.defaultdict(list)
for i, (c, s, e) in enumerate(common):
    starts[c].append(s); ends[c].append(e); index[c].append(i)
with open(f'{OUT}/sample_by_common_peak.raw_fragment_counts.tsv', 'w', newline='\n') as f:
    f.write('sample\t' + '\t'.join(ids) + '\n')
    for r in runs:
        counts = [0] * len(common)
        view = subprocess.Popen([SAM, 'view', '-f', '64', '-F', '2304', bam(r)], stdout=subprocess.PIPE, text=True)
        for line in view.stdout:
            x = line.split('\t')
            c, a = x[2], min(int(x[3]), int(x[7])) - 1
            b = a + abs(int(x[8]))
            for j in range(bisect.bisect_right(ends[c], a), bisect.bisect_left(starts[c], b)):
                counts[index[c][j]] += 1
        if view.wait(): raise SystemExit(f'samtools failed for {r}')
        f.write(r + '\t' + '\t'.join(map(str, counts)) + '\n')
