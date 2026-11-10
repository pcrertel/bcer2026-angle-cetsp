#!/usr/bin/env python3
"""run_alns.py -- runs the ALNS-BCD on one instance with the configuration used in the paper.

Usage:  python3 run_alns.py INSTANCE_FILE SEED [--iterations N] [--time-limit SECONDS] [--outdir DIR]
Example: python3 run_alns.py ../../instances/angle-tsp-radii/angletsp-n050-01.txt 123456

Instance file: the number m of disks in the first line, followed by one line "x y r" per disk (a depot,
if any, is a disk of radius zero).

The run is executed in OUTDIR (default: run-<instance>-s<seed>), which receives
  log.txt               the output of the program: the cost of the initial solution, one line
                        "iter = ... Best sol = ... Time = ..." per new best solution, the statistics
                        of the operators, and a last line "Cost = <value> Time = <seconds>";
  solutions/best-*.txt  the best solution: a header line (name, seed, label, total turning angle,
                        iterations, time), the number m of disks and m lines "i x y", in the order of
                        the route, with the index i of the disk (1..m, in the order of the instance
                        file) and its visit point (x, y) in the coordinates of the instance.

Configuration (Section 5.3 of the paper): N = 200000 iterations, three temperature cycles of ceil(N/3)
iterations, beta = 3, chi = 0.05, T_end = 0.01, lambda = 0.1, weight segment of 100 iterations, scores
(10, 5, 3), largest number of removed disks 100, p_split = 0.75, p_blink = 0.01, the eleven removal
operators and the insertion operators Greedy, Regret-2, Regret-3 and SISR greedy with blinks, the
closed-form block update and the queue block selection rule. The ALNS-BCD is deterministic for a given
seed. The seeds used in the paper are 123456, 42, 34567, 99 and 2026, and, for the Angle-TSP instances,
also 7, 2718, 31415, 1618 and 8128.
"""
import argparse, os, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(HERE, 'alns_cetsp')


def input_lines(instance, seed, n_iter):
    # one value per line, in the order read by inputdata (instance.f90):
    # instance, seed, N, beta, chi, T_end, lambda, weight segment, sigma1, sigma2, sigma3,
    # largest number of removed disks, p_split, p_blink
    lines = ['"%s"' % instance, seed, n_iter, '3.0', '0.05', '0.01', '0.10', '100',
             '10.0', '5.0', '3.0', '100', '0.75', '0.01']
    return '\n'.join(str(x) for x in lines) + '\n'


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('instance')
    ap.add_argument('seed')
    ap.add_argument('--iterations', type=int, default=200000)
    ap.add_argument('--time-limit', type=int, default=3600)
    ap.add_argument('--outdir')
    a = ap.parse_args()
    if not os.path.exists(BIN):
        sys.exit('alns_cetsp not found: run build.sh first')
    inst = os.path.abspath(a.instance)
    name = os.path.splitext(os.path.basename(inst))[0]
    out = a.outdir or 'run-%s-s%s' % (name, a.seed)
    os.makedirs(out, exist_ok=True)
    env = dict(os.environ, LC_NUMERIC='C')
    with open(os.path.join(out, 'log.txt'), 'w') as log:
        p = subprocess.run(['timeout', str(a.time_limit), BIN], cwd=out, env=env, text=True,
                           input=input_lines(inst, a.seed, a.iterations), stdout=log, stderr=subprocess.STDOUT)
    print('exit status %d; output in %s' % (p.returncode, out))


if __name__ == '__main__':
    main()
