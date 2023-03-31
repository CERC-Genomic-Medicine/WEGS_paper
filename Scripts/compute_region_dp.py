#!/usr/bin/env python3

import argparse
import pysam
from intervaltree import IntervalTree
from collections import Counter
import numpy as np

argparser = argparse.ArgumentParser(description = 'Computes DP and %GC for regions in BED file.')
argparser.add_argument('-b', '--bed', metavar = 'file', dest = 'in_bed', type = str, required = True, help = 'Input BED file with regions.')
argparser.add_argument('-r', '--reference', metavar = 'file', dest = 'in_ref_fasta', type = str, required = True, help = 'Input reference FASTA file (indexed).')
argparser.add_argument('-d', '--depth', metavar = 'file', dest = 'in_depth', type = str, required = True, help = 'Tabix indexed tab-delimited file with three columns: chromosome, position, depth. No header required. If position is missing, the DP=0 is assumed.')
argparser.add_argument('-c', '--chromosome', metavar = 'name', dest = 'in_chrom', type = str, required = False, help = 'Chromosome name. To include only autosomals use `auto`.')
argparser.add_argument('-o', '--output', metavar = 'file', dest = 'out_filename', type = str, required = True, help = 'Output file name.')


CHROMOSOMES = [ str(i) for i in range(1, 23) ] + ['X', 'Y']

DP_BINS = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110]

# We load regions into the interval tree just to get them sorted when processing later
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
        n_before = len(intervals)
        intervals.merge_overlaps()
        assert len(intervals) == n_before # we don't expect any overlap in target regions
    return chromosome_intervals


def get_gc_content(fasta_filename, chrom, begin, end):
    bases = Counter()
    with pysam.FastaFile(fasta_filename) as ifasta:
        fasta_chroms = set(list(ifasta.references))
        if chrom not in fasta_chroms:
            chrom = chrom[3:] if chrom.startswith('chr') else f'chr{chrom}'
            if chrom not in fasta_chroms:
                return { 'N_ACGT' : None, 'N_GC': None } 
        for base in ifasta.fetch(chrom, begin, end):
            bases[base] += 1
        n_acgt = sum(bases[x] for x in ['A', 'C', 'G', 'T'])
        n_gc = bases['C'] + bases['G']
    return { 'N_ACGT': n_acgt, 'N_GC': n_gc }


def get_dp(depth_file, chrom, begin, end):
    n_total = end - begin
    bases_dp = np.zeros(n_total, dtype = np.uint16)
    bases_dp_bins = Counter()
    result =  { 'SUM_DP': None, 'AVG_DP': None, 'MEDIAN_DP': None }
    with pysam.TabixFile(depth_file) as itabix:
        tabix_chroms = set(list(itabix.contigs))
        if chrom not in tabix_chroms:
            chrom = chrom[3:] if chrom.startswith('chr') else f'chr{chrom}'
            if chrom not in tabix_chroms:
                return result
        for i, row in enumerate(itabix.fetch(chrom, begin, end, parser = pysam.asTuple()), 0): # fetch is 0-based and bed is 0-based, so everything should be good
            bases_dp[i] = np.uint16(row[2])
            for min_dp in DP_BINS:
                if bases_dp[i] > min_dp:
                    bases_dp_bins[min_dp] += 1
                else:
                    break
        result =  { 'SUM_DP': np.sum(bases_dp), 'AVG_DP': np.mean(bases_dp), 'MEDIAN_DP': np.median(bases_dp) }
        result.update(bases_dp_bins)
        return result


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
    columns = ['N_ACGT', 'N_GC', 'SUM_DP', 'AVG_DP', 'MEDIAN_DP'] + DP_BINS
    with open(args.out_filename, 'wt', buffering = 1) as ofile:
        ofile.write('CHROM\tSTART\tSTOP\t{}\n'.format('\t'.join([f'MIN_DP_{c}' for c in columns])))
        for chrom, intervals in chromosome_intervals.items():
            for interval in sorted(intervals):
                stats = get_gc_content(args.in_ref_fasta, chrom, interval.begin, interval.end)
                stats.update(get_dp(args.in_depth, chrom, interval.begin, interval.end))
                ofile.write('{}\t{}\t{}\t{}\n'.format(chrom, interval.begin, interval.end, '\t'.join([str(stats.get(c, 0)) for  c in columns ])))

