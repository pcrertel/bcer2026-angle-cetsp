#!/usr/bin/env python3
"""feasibility.py -- checks the feasibility of solution files and recomputes their total turning angle.

Usage:
  python3 feasibility.py INSTANCE_FILE SOLUTION_FILE        one solution
  python3 feasibility.py --dirs INSTANCE_DIR SOLUTION_DIR   every solution of SOLUTION_DIR, each checked against
                                                            the instance file with the same name in INSTANCE_DIR
Example, from the root of the repository:
  python3 scripts/feasibility.py --dirs instances/angle-tsp-disk solutions/angle-tsp-disk > report.csv

Instance file: first line m; then m lines "x y r" (a depot, if any, has r = 0).
Solution file: first line "name seed label cost ..."; second line m; then m lines "i x y", in the order of the
route, with i the index of the disk (1..m, in the order of the instance file) and (x, y) its visit point.

A solution is feasible when
  (1) the route has m points and the indices form a permutation of 1..m;
  (2) each visit point lies in its disk, ||p_i - c_i|| <= r_i, with a relative tolerance of 1e-9 of the scale
      of the instance (for r_i = 0, the point is the center).
The recomputed cost is the total turning angle of the closed route, the sum of atan2(|u x v|, <u,v>) at the
vertices, u and v being the incoming and outgoing edges. Consecutive coincident points (edge shorter than 1e-12
of the scale) count as a single vertex, and a route with one or two distinct vertices costs 2*pi.
Output: one line per solution, "status;solution;m;max_distance_outside;reported_cost;recomputed_cost;
relative_difference", where status is OK or INFEASIBLE:<reason>, and a summary on stderr.
"""
import math, os, sys


def read_numbers(path):
    with open(path) as f:
        return [line.split() for line in f if line.strip()]


def turning(points):
    sc = max((abs(x) + abs(y) for x, y in points), default=0.0) or 1.0
    q = []
    for x, y in points:                                   # merges consecutive coincident points
        if q and abs(x - q[-1][0]) + abs(y - q[-1][1]) <= 1e-12 * sc:
            continue
        q.append((x, y))
    while len(q) > 1 and abs(q[-1][0] - q[0][0]) + abs(q[-1][1] - q[0][1]) <= 1e-12 * sc:
        q.pop()
    if len(q) < 3:
        return 2.0 * math.pi
    s = 0.0
    n = len(q)
    for k in range(n):
        u = (q[k][0] - q[k - 1][0], q[k][1] - q[k - 1][1])
        v = (q[(k + 1) % n][0] - q[k][0], q[(k + 1) % n][1] - q[k][1])
        s += math.atan2(abs(u[0] * v[1] - u[1] * v[0]), u[0] * v[0] + u[1] * v[1])
    return s


def check(instance_file, solution_file):
    rows = read_numbers(instance_file)
    m = int(float(rows[0][0]))
    disks = [(float(r[0]), float(r[1]), float(r[2])) for r in rows[1:]]
    scale = max((abs(x) + abs(y) + r for x, y, r in disks), default=0.0) or 1.0
    status = 'OK'
    rows = read_numbers(solution_file)
    reported = float(rows[0][3])
    mr = int(float(rows[1][0]))
    points, seen, out = [], set(), 0.0
    for r in rows[2:]:
        i, x, y = int(float(r[0])), float(r[1]), float(r[2])
        points.append((x, y))
        if not 1 <= i <= m:
            status = 'INFEASIBLE:index %s' % r[0]
            continue
        if i in seen:
            status = 'INFEASIBLE:disk %d repeated' % i
            continue
        seen.add(i)
        cx, cy, cr = disks[i - 1]
        d = math.hypot(x - cx, y - cy) - cr
        out = max(out, d)
        if d > 1e-9 * scale and status == 'OK':
            status = 'INFEASIBLE:point outside disk %d' % i
    if len(disks) != m:
        status = 'INFEASIBLE:instance has %d of %d disks' % (len(disks), m)
    elif status == 'OK' and (mr != m or len(points) != m):
        status = 'INFEASIBLE:route with %d points for %d disks' % (len(points), m)
    rec = turning(points) if points else 0.0
    rel = abs(rec - reported) / abs(reported) if reported != 0 else abs(rec)
    return status, m, out / scale, reported, rec, rel


def main():
    a = sys.argv[1:]
    if len(a) == 3 and a[0] == '--dirs':
        pairs = [(os.path.join(a[1], f), os.path.join(a[2], f)) for f in sorted(os.listdir(a[2])) if f.endswith('.txt')]
    elif len(a) == 2:
        pairs = [(a[0], a[1])]
    else:
        sys.exit(__doc__)
    print('status;solution;m;max_distance_outside;reported_cost;recomputed_cost;relative_difference')
    ok, worst, which = 0, 0.0, ''
    for fi, fs in pairs:
        st, m, out, rep, rec, rel = check(fi, fs)
        print('%s;%s;%d;%.3e;%.17g;%.15g;%.3e' % (st, fs, m, out, rep, rec, rel))
        ok += st == 'OK'
        if rel > worst:
            worst, which = rel, fs
    print('# %d solutions: %d feasible, %d infeasible; largest relative difference between reported and recomputed '
          'cost %.3e (%s)' % (len(pairs), ok, len(pairs) - ok, worst, which), file=sys.stderr)


if __name__ == '__main__':
    main()
