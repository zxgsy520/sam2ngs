#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import sys
import json
import logging
import argparse

LOG = logging.getLogger(__name__)

__version__ = "1.1.0"
__author__ = ("Xingguo Zhang",)
__email__ = "113178210@qq.com"
__all__ = []


def read_tsv(file, sep=None):
    """
    read tsv joined with sep
    :param file: file name
    :param sep: separator
    :return: list
    """
    LOG.info("reading message from %r" % file)

    for line in open(file):
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        yield line.split(sep)


def read_stdin(sep=None):

    for line in sys.stdin:
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        yield line.split(sep)


def read_fasta(file):
    '''Read fasta file'''

    if file.endswith(".gz"):
        fp = gzip.open(file)
    elif file.endswith(".fasta") or file.endswith(".fa"):
        fp = open(file)
    else:
        raise Exception("%r file format error" % file)

    seq = ''

    for line in fp:
        if type(line) == type(b''):
            line = line.decode('utf-8')
        line = line.strip()

        if not line:
            continue
        if not seq:
            seq += "%s\n" % line.strip(">").split()[0]
            continue
        if line.startswith(">"):
            line = line.strip(">").split()[0]
            seq = seq.split('\n')

            yield [seq[0], seq[1]]
            seq = ''
            seq += "%s\n" % line
        else:
            seq += line

    seq = seq.split('\n')
    if len(seq)==2:
        yield [seq[0], seq[1]]
    fp.close()


def read_fastq(file):
    '''Read fastq file'''

    if file.endswith(".gz"):
        fp = gzip.open(file)
    elif file.endswith(".fastq") or file.endswith(".fq"):
        fp = open(file)
    else:
        raise Exception("%r file format error" % file)

    seq = []

    for line in fp:
        line = line.strip()

        if not line:
            continue
        if line.startswith('@') and (len(seq)==0 or len(seq)>=5):
            seq.append(line.strip("@").split()[0])
            continue
        if line.startswith('@') and len(seq)==4:
            yield seq
            seq = []
            seq.append(line.strip("@").split()[0])
        else:
            seq.append(line)

    if len(seq)==4:
        yield seq
    fp.close()


def reverse_complement(string):

    rule = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}

    string = string.upper()
    string = "".join(map(lambda x:rule[x], string[::-1]))

    return string


def read_sam(file):

    if file == "":
        fs = read_stdin('\t')
    elif file.endswith(".sam"):
        fs = read_tsv(file, '\t')
    else:
        raise Exception("%r file format error" % file)
    
    reads = []
    for line in fs:
        if line[0].startswith("@") or line[2]=="*":
            continue

        seqid = line[0].split('/')[0]
        tigid = line[2]
        seq = line[9]
        quality = line[10]
        if line[0].startswith("-"):
            seq = reverse_complement(line[9])
            quality = line[10][::-1]

        if len(reads)==0 or len(reads)>2:
            reads = []
            reads.append([seqid, seq, quality, tigid])
            continue

        if len(reads)==2:
            if reads[0][0]==reads[1][0]:
                yield reads
            reads = []

        reads.append([seqid, seq, quality, tigid])

    if len(reads)==2:
        if reads[0][0]==reads[1][0]:
            yield reads


def choose_reads(file, genome, prefix, depth='all'):

    if depth=='all':
        depth = float("inf")
    else:
        depth = int(depth)
    if genome.endswith(".fastq") or genome.endswith(".fq") or genome.endswith(".fastq.gz") or genome.endswith(".fq.gz"):
        fh = read_fastq(genome)
    elif genome.endswith(".fasta") or genome.endswith(".fa") or genome.endswith(".fasta.gz") or genome.endswith(".fa.gz"):
        fh = read_fasta(genome)
    else:
        raise Exception("%r file format error" % genome)

    genome = {}

    for line in fh:
        genome[line[0]] = len(line[1])

    data = {}
    connection = {}
    fr1 = open('%s.sam_R1.fq' % prefix, 'w')
    fr2 = open('%s.sam_R2.fq' % prefix, 'w')

    for r1, r2 in read_sam(file):
        if r1[3]!=r2[3]:
            conn = "%s %s" % (r1[3], r2[3])
            if conn not in connection:
                connection[conn] = 0
            connection[conn] += 1
            #LOG.info("Contig %s assembly may be wrong or broken" % r1[3])
            continue
        tigid = r1[3]

        if tigid not in data:
            data[tigid] = 0
        if data[tigid]>=genome[tigid]*depth:
            continue
        fr1.write("@%s 1:N:0:TGATAACG\n%s\n+\n%s\n" % (r1[0], r1[1], r1[2]))
        fr2.write("@%s 2:N:0:TGATAACG\n%s\n+\n%s\n" % (r2[0], r2[1], r2[2]))

        data[tigid] += len(r1[1])
        data[tigid] += len(r1[2])

    fr1.close()
    fr2.close()

    with open("%s.link.json" % prefix, "w") as fh:
        json.dump(connection, fh, indent=2)


def add_hlep_args(parser):

    parser.add_argument("-i", "--input", metavar='FILE', type=str, default="",
        help="Input files in sam.")
    parser.add_argument("-g", "--genome", metavar='FILE', type=str, required=True,
        help="Inport genome files.")
    parser.add_argument("-d", "--depth", metavar='STR', type=str, default=200,
        help="Set the reads depth selected by each contig(Select all=all), default=200")
    parser.add_argument("-p", "--prefix", metavar='STR', type=str, default='out',
        help="Inport the file prefix.")

    return parser


def main():

    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
        description='''
name:
    sam2ngs Select no reads on the map

attention:
    sam2ngs -i input.sam -g genome.fa -p name
    minimap2 -t 8 -ax sr genome.fa r1.fq r2.fq|samblaster -a |sam2ngs -g genome.fa -p name

version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_hlep_args(parser).parse_args()

    choose_reads(args.input, args.genome, args.prefix, args.depth)


if __name__ == "__main__":

    main()
