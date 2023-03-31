#!/usr/bin/env python3

import argparse
import pysam
from intervaltree import IntervalTree
from collections import Counter


argparser = argparse.ArgumentParser(description = 'Computes overall DP across provided chromosomal regions.')
argparser.add_argument('-b', '--bed', metavar = 'file', dest = 'in_bed', type = str, required = True, help = 'Input BED file with regions. Overlapping regions will be automatically merged.')
argparser.add_argument('-d', '--depth', metavar = 'file', dest = 'in_depth', type = str, required = True, help = 'Tabix indexed tab-delimited file with three columns: chromosome, position, depth. No header required. If position is missing, the DP=0 is assumed.')
argparser.add_argument('-c', '--chromosome', metavar = 'name', dest = 'in_chrom', type = str, required = False, help = 'Chromosome name. To include only autosomals use `auto`.')
argparser.add_argument('-o', '--output', metavar = 'file', dest = 'out_filename', type = str, required = True, help = 'Output file name.')


CHROMOSOMES = [ str(i) for i in range(1, 23) ] + ['X', 'Y']


# This merges overlappin BED intervals, so that bases are not counted multiple times
def load_intervals(filename, chromosomes):
    chromosome_intervals = {}
    with open(filename, 'rt') as ifile:
        for line in ifile:
            fields = line.rstrip().split()
            chrom = fields[0]
            if chrom.startswith('chr'):
                chrom = chrom[3:]
            if chrom not in chromosomes:
                continue
            start = int(fields[1])
            stop = int(fields[2])
            intervals = chromosome_intervals.setdefault(chrom, IntervalTree())
            intervals.addi(start, stop) # 0-based, last position is not included
    for chrom, intervals in chromosome_intervals.items():
        intervals.merge_overlaps()
    return chromosome_intervals


if __name__ == '__main__':
    args = argparser.parse_args()
    if args.in_chrom is not None:
        if args.in_chrom.lower() == 'auto':
            chromosomes = set([ str(i) for i in range(1, 23) ])
        else:
            chromosomes = set([ args.in_chrom ])
    else:
        chromosomes = set(CHROMOSOMES)
    chromosome_intervals = load_intervals(args.in_bed, chromosomes)
    N = 0
    DP_SUM = 0
    DP_BINS = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110]
    DP_BINS_COUNTS = Counter()
    with pysam.TabixFile(args.in_depth) as itabix:
        tabix_chroms = set(list(itabix.contigs))
        for chrom, intervals in chromosome_intervals.items():
            if chrom not in tabix_chroms:
                chrom = chrom[3:] if chrom.startswith('chr') else f'chr{chrom}'
                if chrom not in tabix_chroms:
                    continue
            for interval in sorted(intervals):
                N += interval.end - interval.begin
                for i, row in enumerate(itabix.fetch(chrom, interval.begin, interval.end, parser = pysam.asTuple()), 0): # fetch is 0-based and bed is 0-based, so everything should be good
                    dp = int(row[2])
                    DP_SUM += dp
                    for min_dp in DP_BINS:
                        if dp > min_dp:
                            DP_BINS_COUNTS[min_dp] += 1
                        else:
                            break
    with open(args.out_filename, 'wt', buffering = 1) as ofile:
        ofile.write('# Chromosomes: {}\n'.format(','.join(sorted(list(chromosomes)))))
        ofile.write('N\tAVG_DP\t{}\n'.format('\t'.join([f'DP{min_dp}' for min_dp in DP_BINS])))
        ofile.write('{}\t{}\t{}\n'.format(N, DP_SUM / N, '\t'.join([ f'{DP_BINS_COUNTS[min_dp] / N}' for min_dp in DP_BINS])))

