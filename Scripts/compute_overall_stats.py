#!/usr/bin/env python3

import argparse
import pysam


argparser = argparse.ArgumentParser(description = 'Computes basic stats for BAM file.')
argparser.add_argument('-s', '--sequences', metavar = 'file', dest = 'in_bam', type = str, required = True, help = 'Input BAM file.')
argparser.add_argument('-F', '--exclude-flag', metavar = 'integer', dest = 'exclude_flag', type = int, required = True, help = 'Samtools exclusion flag (see -F option in `samtools view`). Only integer numbers are accepted. Specify 0 for no filters.')
argparser.add_argument('-q', '--min-read-map-quality', metavar = 'integer', dest = 'min_read_map_quality', type = int, required = True, help = 'Minimal read mapping quality.')
argparser.add_argument('-c', '--chromosome', metavar = 'name', dest = 'in_chrom', type = str, required = False, help = 'Chromosome name. To include only autosomals use `auto`.')
argparser.add_argument('-o', '--output', metavar = 'file', dest = 'out_filename', type = str, required = True, help = 'Output file name.')


CHROMOSOMES = [ str(i) for i in range(1, 23) ] + ['X', 'Y']


def get_reads_stats(bam_filename, chromosomes, exclude_flag, min_map_quality):
    n_nofilter_reads = 0
    n_nofilter_bases = 0
    n_nofilter_bqsum = 0
    n_qcfail_reads = 0
    n_qcfail_bases = 0
    n_qcfail_bqsum = 0
    n_reads = 0
    n_bases = 0
    n_bqsum = 0
    n_unmapped_reads = 0
    n_duplicated_reads = 0
    n_minq_reads = 0
    n_minq_bases = 0
    n_minq_bqsum = 0
    with pysam.AlignmentFile(bam_filename, 'rb') as ibam:
        bam_chroms = set(list(ibam.references))
        for chrom in chromosomes:
            if chrom not in bam_chroms:
                chrom = chrom[3:] if chrom.startswith('chr') else f'chr{chrom}'
            for read in ibam.fetch(chrom):
                if not read.is_paired:
                    continue
                if read.is_secondary or read.is_supplementary:
                    continue
                bases = len(read.query_alignment_qualities)
                bqsum = sum(read.query_alignment_qualities)
                n_nofilter_reads += 1
                n_nofilter_bases += bases
                n_nofilter_bqsum += bqsum
                if read.flag & 512 != 0:
                    n_qcfail_reads += 1
                    n_qcfail_bases += bases
                    n_qcfail_bqsum += bqsum
                if read.flag & exclude_flag != 0:
                    continue
                n_reads += 1
                n_bases += bases
                n_bqsum += bqsum
                if read.flag & 4 != 0:
                    n_unmapped_reads += 1
                if read.flag & 1024 != 0:
                    n_duplicated_reads += 1
                if read.mapping_quality < min_map_quality:
                    continue
                n_minq_reads += 1
                n_minq_bases += bases
                n_minq_bqsum += bqsum
    return { 
            'N_NOFILTER_READS': n_nofilter_reads, 
            'N_NOFILTER_BASES': n_nofilter_bases,
            'N_NOFILTER_BQSUM': n_nofilter_bqsum,
            'N_QCFAIL_READS': n_qcfail_reads,
            'N_QCFAIL_BASES': n_qcfail_bases,
            'N_QCFAIL_BQSUM': n_qcfail_bqsum,
            'N_READS': n_reads,
            'N_BASES': n_bases,
            'N_BQSUM': n_bqsum,
            'N_UNMAPPED_READS': n_unmapped_reads,
            'N_DUPLICATED_READS': n_duplicated_reads,
            'N_MINQ_READS': n_minq_reads,
            'N_MINQ_BASES': n_minq_bases, 
            'N_MINQ_BQSUM': n_minq_bqsum 
        }


if __name__ == '__main__':
    args = argparser.parse_args()
    if args.in_chrom is not None:
        if args.in_chrom.lower() == 'auto':
            chromosomes = set([ str(i) for i in range(1, 23) ])
        else:
            chromosomes = set([ args.in_chrom ])
    else:
        chromosomes = set(CHROMOSOMES)
    columns = ['N_NOFILTER_READS', 'N_NOFILTER_BASES', 'N_NOFILTER_BQSUM', 'N_QCFAIL_READS', 'N_QCFAIL_BASES', 'N_QCFAIL_BQSUM', 'N_READS', 'N_BASES', 'N_BQSUM', 'N_UNMAPPED_READS', 'N_DUPLICATED_READS', 'N_MINQ_READS', 'N_MINQ_BASES', 'N_MINQ_BQSUM']
    with open(args.out_filename, 'wt', buffering = 1) as ofile:
        ofile.write('# Chromosomes: {}\n'.format(','.join(sorted(list(chromosomes)))))
        ofile.write('{}\n'.format('\t'.join(columns)))
        stats = get_reads_stats(args.in_bam, chromosomes, args.exclude_flag, args.min_read_map_quality)
        ofile.write('{}\n'.format('\t'.join([str(stats[c]) for  c in columns ])))
