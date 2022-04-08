#!/usr/bin/env python3

################
# vcf2align.py #
################
#
# returns a multiple sequence alignment (by default in fasta format) of samples in a VCF/BCF file for a given genomic region
# heterozygous and missing genotypes are coded as ambiguous nucleotides (R, Y, S, W, K, M, N)
# the alignment can include only the positions in the VCF or fill in the invariant positions from a reference genome
# for options and command line syntax, run:
# $ python vcf2align.py -h


import sys
from argparse import ArgumentParser

from pysam import VariantFile
from pysam import FastaFile

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment


def get_gt_nucleotide(GT, alleles):
    '''
    returns a nucleotide given a genotype (GT) and possible alleles
    GT is a tuple as returned by record.gt, e.g. (0,0), (0,1) or (1,1)
    alleles is a tuple as returned by record.alleles, e.g. ('A','T')

    nucleotide is denoted by IUPAC ambiguity code if heterozygous or missing (N),
    X if unrecognized combination of alleles (in case of indels)
    '''

    ambiguity_codes = {'A/T': 'W',
                       'T/A': 'W',
                       'C/G': 'S',
                       'G/C': 'S',
                       'A/C': 'M',
                       'C/A': 'M',
                       'G/T': 'K',
                       'T/G': 'K',
                       'C/T': 'Y',
                       'T/C': 'Y',
                       'A/G': 'R',
                       'G/A': 'R'
                       }

    if None in GT:
        genotype = 'N'*len(alleles[0])
    else:
        bases = (alleles[GT[0]], alleles[GT[1]])

        if bases[0] == bases[1]:
            genotype = bases[0]
        else:
            genotype = ambiguity_codes.get(f"{bases[0]}/{bases[1]}", 'X'*len(bases[0]))

    return genotype


def pad_alleles(alleles, length):
    '''pads each element of alleles to length by adding '-' '''

    return tuple(map(lambda s: s.ljust(length, '-'), alleles))


def get_region_align(vcf, region, reference=None):
    '''
    returns a sequence alignment (MultipleSeqAlignment object) for a region in the VCF
    if a reference genome is given, (invariant) positions not in the VCF are filled
    vcf and reference are expected to be the file paths
    region is expected to be a samtools formated region string, i.e.
    'contig:start-end' with 1-based inclusive coordinates
    '''

    variant_file = VariantFile(vcf)

    # initializing empty dict for sequence records for each sample
    sample_list = [s for s in variant_file.header.samples]
    record_dict = {sample: SeqRecord(Seq(''), id=sample, description='')
                   for sample in sample_list}

    if reference:
        ref_fasta = FastaFile(reference)
        contig, start, end = region.replace(':','-').split('-')
        start, end = int(start), int(end)
        expected_position = start
        record_dict['reference'] = SeqRecord(Seq(''), id='reference', description='')

    # fetching region from VCF
    vcf_region = variant_file.fetch(region=region)

    # looping over each position in the VCF
    for record in vcf_region:

        # fill positions from reference BEFORE the position in the VCF
        if reference:
            if record.pos > expected_position:
                fill_ref_region = f"{contig}:{expected_position}-{record.pos-1}"
                fill_ref_seq = ref_fasta.fetch(region=fill_ref_region)
                for sample in record_dict:
                    record_dict[sample].seq += fill_ref_seq

        # process position
        alleles = record.alleles

        # check for alleles with size >1
        # if yes, warn the user and pad the alleles to the same size
        max_length_alleles = max(map(len, alleles))
        if max_length_alleles > 1:
            warning = f"WARNING: alleles with length >1 detected at {record.chrom}:{record.pos}, possible indels\n"
            sys.stderr.write(warning)
            alleles = pad_alleles(alleles, max_length_alleles)

        for sample in record.samples:
            GT = record.samples[sample].get('GT')
            genotype = get_gt_nucleotide(GT, alleles)
            record_dict[sample].seq += genotype

        if reference:
            record_dict['reference'].seq += alleles[0]
            expected_position = record.pos + 1

    # fill positions from reference AFTER last position in the VCF
    if reference:
        if expected_position <= end:
            fill_ref_region = f"{contig}:{expected_position}-{end}"
            fill_ref_seq = ref_fasta.fetch(region=fill_ref_region)
            for sample in record_dict:
                record_dict[sample].seq += fill_ref_seq

    # turning seq records dict into multiple sequence alignment object
    align = MultipleSeqAlignment(list(record_dict.values()))

    return align


if __name__ == "__main__":

    # creating and parsing command line arguments
    arg_parser = ArgumentParser(
        description='Creates alignments from VCF files')
    arg_parser.add_argument('vcf', help='an indexed VCF/BCF file')
    arg_parser.add_argument('region', help='a genomic region')
    arg_parser.add_argument('-f', '--format',
                            help='file format of the output alignment (default: fasta)',
                            type=str,
                            default='fasta',
                            choices=[
                                'clustal',
                                'emboss',
                                'fasta',
                                'fasta-m10',
                                'ig',
                                'maf',
                                'mauve',
                                'msf',
                                'nexus',
                                'phylip',
                                'phylip-sequential',
                                'phylip-relaxed',
                                'stockholm'
                                ]
                            )
    arg_parser.add_argument('-ref', '--reference',
                            help='indexed reference genome fasta',
                            type=str,
                            default=None
                            )
    arg_parser.add_argument('-o', '--output', help='output file name', type=str, default=None)
    args = arg_parser.parse_args()

    # run
    align = get_region_align(args.vcf, args.region, reference=args.reference)
    formatted_align = format(align, args.format)

    if args.output:
        with open(args.output, 'w') as outf:
            outf.write(formatted_align)
    else:
        sys.stdout.write(formatted_align)
