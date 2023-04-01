#!/usr/bin/env python3

import argparse
import os
import pysam
#from collections import Counter
import numpy as np

argparser = argparse.ArgumentParser(description = 'Compare DP and BQ in false negative calls in two samples.')
argparser.add_argument('-v1', '--vcf-happy1', metavar = 'file', dest = 'in_happy1', type = str, nargs = '*', required = True, help = 'The VCF output from hap.py for the first sample. Files will be matched by prefix with `-d1`.')
argparser.add_argument('-d1', '--dp-bq1', metavar = 'file', dest = 'in_dp_bq1', type = str, nargs = '*', required = True, help = 'The tabix-indexed file with DP and BQ for each position. Files will be matched by prefix with `-h1`.')


argparser.add_argument('-v2', '--vcf-happy2', metavar = 'file', dest = 'in_happy2', type = str, nargs = '*', required = True, help = 'The VCF output from hap.py for the first sample. Files will be matched by prefix with `-d2`.')
argparser.add_argument('-d2', '--dp-bq2', metavar = 'file', dest = 'in_dp_bq2', type = str, nargs = '*', required = True, help = 'The tabix-indexed file with DP and BQ for each position. Files will be matched by prefix with `-h2`.')


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


def ascii2bq(bq_string):
    return [ord(ascii_char) - 33 for ascii_char in bq_string]


def count_ref_bases(bases_string):
    n_ref = 0
    for i, code in enumerate(bases_string):
        if code == '.' or code == ',':
            if i > 0 and bases_string[i - 1] != '^':
                n_ref += 1
    return n_ref


def bq_test(bases_string, bq_string):
    #if '-' not in bases_string:
    #    return ([], [])
    #print(bases_string, bq_string)  
    #assert '*' not in bases_string
    #if '*' in bases_string:
    #    print(bases_string, bq_string)
    #assert '#' not in bases_string
    i = 0
    bq_i = 0
    ref_bq = []
    nonref_bq = []
    while i < len(bases_string):
        if bases_string[i] == '$':
            i += 1
        elif bases_string[i] == '^':
            i += 2
        elif bases_string[i] == '-' or bases_string[i] == '+':
            i += 1
            indel_length = ''
            while bases_string[i].isnumeric():
                indel_length += bases_string[i]
                i += 1
            indel_length = int(indel_length)
            #print('\t', indel_length, i, bases_string[i], bases_string[i + indel_length - 1])
            #base = bases_string[i:i + del_length]
            i += indel_length
        else:
            base = bases_string[i]
            #print(bq_i, bq_string[bq_i], base)
            if base == '.' or base == ',':
                ref_bq.append(ord(bq_string[bq_i]) - 33)
            else:
                nonref_bq.append(ord(bq_string[bq_i]) - 33)
            bq_i += 1
            i += 1
    return (ref_bq, nonref_bq)




def load_false_negatives(happy_filename, dp_bq_filename, chromosomes, fn_variants):
    with pysam.VariantFile(happy_filename) as ivcf, pysam.TabixFile(dp_bq_filename) as itabix:
        for record in ivcf:
            chrom = record.chrom[:3] if record.chrom.startswith('chr') else record.chrom
            if not chrom in chromosomes:
                continue
            if record.samples['TRUTH']['BD'] != 'FN':
                continue
            if len(record.alts) > 1: # ignore multi-allelic single row entries for simplicity in comparison
                continue
            if len(record.ref) != 1 or len(record.alts[0]) != 1: # ignore InDels for simplicity in comparison
                continue
            dp = None
            n_ref = None
            mean_ref_bq = None
            median_ref_bq = None
            mean_nonref_bq = None
            median_nonref_bq = None
            for row in itabix.fetch(chrom, record.pos - 1, record.pos, parser = pysam.asTuple()): # fetch is 0-based, so we substract 1
                assert dp is None, 'The second iteration was not supposed to happen!'
                dp = np.uint16(row[2])
                if dp == 0:
                    break
                n_ref = count_ref_bases(row[3])
                ref_bq, nonref_bq = bq_test(row[3], row[4])
                assert len(ref_bq) + len(nonref_bq) == dp
                if len(ref_bq) > 0:
                    mean_ref_bq = np.mean(ref_bq)
                    median_ref_bq = np.median(ref_bq)
                if len(nonref_bq) > 0:
                    mean_nonref_bq = np.mean(nonref_bq)
                    median_nonref_bq = np.median(nonref_bq)
            if dp is None:
                dp = 0
            variant_name = f'{chrom}_{record.pos}'
            variant_info = fn_variants.get(variant_name, None)
            if variant_info is None:
                fn_variants[variant_name] = {
                    'N_1': 1,
                    'DP_1': dp,
                    'N_REF_1': n_ref,
                    'MEAN_REF_BQ_1': mean_ref_bq,
                    'MEDIAN_REF_BQ_1': median_ref_bq,
                    'MEAN_NONREF_BQ_1': mean_nonref_bq,
                    'MEDIAN_NONREF_BQ_1': median_nonref_bq
                }
            else:
                fn_variants['N_1'] += 1
                if dp > variant_info['DP_1']: # keep only largest DP to be very conservative
                    variant_info['DP_1'] = dp
                    variant_info['N_REF_1'] = n_ref
                    variant_info['MEAN_REF_BQ_1'] = mean_ref_bq,
                    variant_info['MEDIAN_REF_BQ_1'] =  median_ref_bq,
                    variant_info['MEAN_NONREF_BQ_1'] =  mean_nonref_bq,
                    variant_info['MEDIAN_NONREF_BQ_1'] =  median_nonref_bq


def compare_false_negatives(happy_filename, dp_bq_filename, fn_variants):
    with pysam.VariantFile(happy_filename) as ivcf, pysam.TabixFile(dp_bq_filename) as itabix:
        for record in ivcf:
            variant_name = f'{record.chrom}_{record.pos}'
            if variant_name not in fn_variants:
                continue
            fn_variant = fn_variants[variant_name]
            if record.samples['TRUTH']['BD'] != 'TP' or record.samples['QUERY']['BD'] != 'TP':
                if 'N_2' not in fn_variant:
                    fn_variant['N_2'] = 0
                    fn_variant['DP_2'] = None
                    fn_variant['N_REF_2'] = None
                    fn_variant['MEAN_REF_BQ_2'] = None
                    fn_variant['MEDIAN_REF_BQ_2'] = None
                    fn_variant['MEAN_NONREF_BQ_2'] = None
                    fn_variant['MEDIAN_NONREF_BQ_2'] = None
                    #fn_variant['MEAN_BQ_2'] = None
                    #fn_variant['MEDIAN_BQ_2'] = None
                    #fn_variant['BQ10_2'] = None
                    #fn_variant['BQ20_2'] = None
                    #fn_variant['BQ30_2'] = None
                    #fn_variant['BQ36_2'] = None
                    #fn_variant['BQ38_2'] = None
                    #fn_variant['BQ40_2'] = None
                continue 
            dp = None
            n_ref = None
            mean_ref_bq = None
            median_ref_bq = None
            mean_nonref_bq = None
            median_nonref_bq = None
            #mean_bq = None
            #median_bq = None
            #bq10 = bq20 = bq30 = bq36 = bq38  = bq40  = None
            for row in itabix.fetch(record.chrom, record.pos - 1, record.pos, parser = pysam.asTuple()): # fetch is 0-based, so we substract 1
                assert dp is None, 'The second iteration was not supposed to happen!'
                dp = np.uint16(row[2])
                n_ref = count_ref_bases(row[3])
                ref_bq, nonref_bq = bq_test(row[3], row[4])
                assert len(ref_bq) + len(nonref_bq) == dp
                if len(ref_bq) > 0:
                    mean_ref_bq = np.mean(ref_bq)
                    median_ref_bq = np.median(ref_bq)
                if len(nonref_bq) > 0:
                    mean_nonref_bq = np.mean(nonref_bq)
                    median_nonref_bq = np.median(nonref_bq)
                #bq = ascii2bq(row[4])
                #assert len(bq) == dp
                #mean_bq = np.mean(bq)
                #median_bq = np.median(bq)
                #bq10 = bq20 = bq30 = bq36 = bq38  = bq40  = 0
                #for v in bq:
                #    if v < 10:
                #        bq10 += 1
                #    elif v < 20:
                #        bq20 += 1
                #    elif v < 30:
                #        bq30 += 1
                #    elif v < 36:
                #        bq36 += 1
                #    elif v < 38:
                #        bq38 += 1
                #    elif v < 40:
                #        bq40 += 1
            #assert dp is not None, record  # at this point we expect that there will be always DP>0 because the variant is TP. 
            if dp is None: # workaround for GATK
                dp = 2 
            if 'N_2' not in fn_variant:
                #print('here', fn_variant)
                fn_variant['N_2'] = 1
                fn_variant['DP_2'] = dp
                fn_variant['N_REF_2'] = n_ref
                fn_variant['MEAN_REF_BQ_2'] = mean_ref_bq
                fn_variant['MEDIAN_REF_BQ_2'] = median_ref_bq
                fn_variant['MEAN_NONREF_BQ_2'] =  mean_nonref_bq
                fn_variant['MEDIAN_NONREF_BQ_2'] =  median_nonref_bq
                #fn_variant['MEAN_BQ_2'] = mean_bq
                #fn_variant['MEDIAN_BQ_2'] = median_bq
                #fn_variant['BQ10_2'] = bq10
                #fn_variant['BQ20_2'] = bq20
                #fn_variant['BQ30_2'] = bq30
                #fn_variant['BQ36_2'] = bq36
                #fn_variant['BQ38_2'] = bq38
                #fn_variant['BQ40_2'] = bq40
                #print(fn_variant)
            else:
                fn_variant['N_2'] += 1
                if fn_variant['DP_2'] is None or  dp < fn_variant['DP_2']: # take minimal DP to by very conservative
                    fn_variant['DP_2'] = dp
                    fn_variant['N_REF_2'] = n_ref
                    fn_variant['MEAN_REF_BQ_2'] = mean_ref_bq
                    fn_variant['MEDIAN_REF_BQ_2'] = median_ref_bq
                    fn_variant['MEAN_NONREF_BQ_2'] =  mean_nonref_bq
                    fn_variant['MEDIAN_NONREF_BQ_2'] =  median_nonref_bq
                    #fn_variant['MEAN_BQ_2'] = mean_bq
                    #fn_variant['MEDIAN_BQ_2'] = median_bq
                    #fn_variant['BQ10_2'] = bq10
                    #fn_variant['BQ20_2'] = bq20
                    #fn_variant['BQ30_2'] = bq30
                    #fn_variant['BQ36_2'] = bq36
                    #fn_variant['BQ38_2'] = bq38
                    #fn_variant['BQ40_2'] = bq40


if __name__ == '__main__':
    args = argparser.parse_args()

    if args.in_chrom is not None:
        if args.in_chrom.lower() == 'auto':
            chromosomes = set([ str(i) for i in range(1, 23) ])
        else:
            chromosomes = set([ args.in_chrom ])
    else:
        chromosomes = set(CHROMOSOMES)

    in_files1 = []
    fn_variants1 = dict()
    in_files2 = []

    assert len(args.in_happy1) == len(args.in_dp_bq1), f'{len(args.in_happy1)} vs {len(args.in_dp_bq1)}. Number of files in `-h1` and `-d1` must match!'
    for happy_file, dp_bq_file in zip(sorted(args.in_happy1), sorted(args.in_dp_bq1)):
        happy_file_prefix = os.path.basename(happy_file).split('.')[0]
        dp_bq_file_prefix = os.path.basename(dp_bq_file).split('.')[0]
        assert happy_file_prefix == dp_bq_file_prefix, f'Prefixes for files {happy_file} and {dp_bq_file} don\'t match!'
        in_files1.append((happy_file, dp_bq_file))


    assert len(args.in_happy2) == len(args.in_dp_bq2), f'{len(args.in_happy2)} vs {len(args.in_dp_bq2)}. Number of files in `-h2` and `-d2` must match!'
    for happy_file, dp_bq_file in zip(sorted(args.in_happy2), sorted(args.in_dp_bq2)):
        happy_file_prefix = os.path.basename(happy_file).split('.')[0]
        dp_bq_file_prefix = os.path.basename(dp_bq_file).split('.')[0]
        assert happy_file_prefix == dp_bq_file_prefix, f'Prefixes for files {happy_file} and {dp_bq_file} don\'t match!'
        in_files2.append((happy_file, dp_bq_file))

    print("Loading FN calls")
    for happy_file, dp_bq_file in in_files1:
        print(f'\t{happy_file}')
        load_false_negatives(happy_file, dp_bq_file, chromosomes, fn_variants1)

    print("Comparing FN calls")
    for happy_file, dp_bq_file in in_files2:
        #print(happy_file, dp_bq_file)
        compare_false_negatives(happy_file, dp_bq_file, fn_variants1)
        #break

    #columns = ['N_1', 'DP_1', 'N_REF_1', 'MEAN_BQ_1', 'MEDIAN_BQ_1', 'BQ10_1', 'BQ20_1', 'BQ30_1', 'BQ36_1', 'BQ38_1', 'BQ40_1', 'N_2', 'DP_2', 'N_REF_2', 'MEAN_BQ_2', 'MEDIAN_BQ_2', 'BQ10_2', 'BQ20_2', 'BQ30_2', 'BQ36_2', 'BQ38_2', 'BQ40_2']
    columns = ['N_1', 'DP_1', 'N_REF_1', 'MEAN_REF_BQ_1', 'MEDIAN_REF_BQ_1', 'MEAN_NONREF_BQ_1', 'MEDIAN_NONREF_BQ_1', 'N_2', 'DP_2', 'N_REF_2', 'MEAN_REF_BQ_2', 'MEDIAN_REF_BQ_2', 'MEAN_NONREF_BQ_2', 'MEDIAN_NONREF_BQ_2']

    format_template = '\t'.join('{}' for c in columns) + '\n'

    with open(args.out_filename, 'wt') as ofile:
        ofile.write('\t'.join(columns) + '\n')
        for name, stats in fn_variants1.items():
            ofile.write('\t'.join(str(stats.get(c)) for c in columns) + '\n')
            
    '''    
    fn1 = []
    for filename in args.in_happy1:
        fn1.append(load_false_negatives(filename, chromosomes))

    for fn in fn1:
        print(len(fn))

    fn1 = set.intersection(*fn1)
    print(len(fn1))


    fn2 = []
    for filename in args.in_happy2:
        fn2.append(load_false_negatives(filename, chromosomes))

    for fn in fn2:
        print(len(fn))

    fn2 = set.intersection(*fn2)
    print(len(fn2))

    #fn1 = load_false_negatives(args.in_happy1, chromosomes)
    #fn2 = load_false_negatives(args.in_happy2, chromosomes)

    print(len(fn1.intersection(fn2)))
    print(len(fn1 - fn2))
    print(len(fn2 - fn1))
    
    #print(len(fn1), len(fn2), len(fn1.intersection(fn2)))
    '''
