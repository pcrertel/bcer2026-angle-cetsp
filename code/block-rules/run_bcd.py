#!/usr/bin/env python3
"""run_bcd.py -- comparison of the block selection rules of the BCD (Section 5.1 of the paper).

Runs ./bcd_compare (built by run.sh) on the three sets of 100 random instances in instances/
(radii-x1, radii-x4 and radii-x10: r_max multiplied by 1, 4 and 10), writes radii-x1.csv,
radii-x4.csv and radii-x10.csv, and prints the summary reported in the paper.

Usage: python3 run_bcd.py
"""
import csv, glob, os, statistics, subprocess

SETS = ['radii-x1', 'radii-x4', 'radii-x10']
FIELDS = ['instance', 'm', 'psi_cyclic', 'psi_queue', 'updates_cyclic', 'updates_queue',
          'cycles_cyclic', 'cycles_queue', 'time_cyclic_s', 'time_queue_s']
TOL = 1e-9          # tie: relative difference at most TOL


def run(folder):
    rows = []
    for fn in sorted(glob.glob(os.path.join('instances', folder, '*.txt'))):
        a = subprocess.run(['./bcd_compare', fn], capture_output=True, text=True,
                           timeout=3600).stdout.split()
        v = [os.path.basename(fn), int(a[0]), float(a[1]), float(a[2])] + \
            [int(t) for t in a[3:7]] + [float(t) for t in a[7:9]]
        rows.append(dict(zip(FIELDS, v)))
    return rows


for s in SETS:
    rows = run(s)
    with open(s + '.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(rows)
    rel = [100.0 * (r['psi_queue'] - r['psi_cyclic']) / r['psi_cyclic'] for r in rows]
    queue_lower = sum(1 for d in rel if d < -100 * TOL)
    cyclic_lower = sum(1 for d in rel if d > 100 * TOL)
    upd = sum(r['updates_queue'] for r in rows) / sum(r['updates_cyclic'] for r in rows)
    tim = sum(r['time_queue_s'] for r in rows) / sum(r['time_cyclic_s'] for r in rows)
    print('== %s (%d instances) ==' % (s, len(rows)))
    print('  cycles of the cyclic rule: mean %.2f, max %d'
          % (statistics.mean(r['cycles_cyclic'] for r in rows), max(r['cycles_cyclic'] for r in rows)))
    print('  queue lower in %d, cyclic lower in %d, equal in %d; mean relative difference %+.4f%%, '
          'largest %.4f%%' % (queue_lower, cyclic_lower, len(rows) - queue_lower - cyclic_lower,
                              statistics.mean(rel), max(abs(d) for d in rel)))
    print('  block updates of the queue: %.1f%% of the cyclic rule; running time: %.1f%%'
          % (100 * upd, 100 * tim))
