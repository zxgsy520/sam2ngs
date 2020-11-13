#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import sys
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


def sam2unmap(file, minmapq=20):

    data = {}
    mapid = set()

    if file == "":
        fs = read_stdin('\t')
    elif file.endswith(".sam"):
        fs = read_tsv(file, '\t')
    else:
        raise Exception("%r file format error" % file)

    for line in fs:
        if "@PG" in line[0]:
            continue
        if line[0] in mapid:
            continue
        if line[2]!="*" or int(line[4]) > minmapq:
            mapid.add(line[0])
            continue
        if line[0] not in data:
            data[line[0]] = []

        seq = [line[9], line[10]]
        if seq in data[line[0]]:
            continue
        else:
            data[line[0]].append(seq)

    for seqid in data:
        if seqid in mapid:
            continue
        n = 0
        for seq, quality in data[line[0]]:
            n += 1
            print("@%s_%s\n%s\n+\n%s" % (seqid, n, seq, quality))

    return 0


def add_hlep_args(parser):

    parser.add_argument("-i", "--input", metavar='FILE', type=str, default="",
        help="Input files in sam.")
    parser.add_argument("-mq", "--minmapq", metavar='INT', type=int, default=20,
        help="Set the minimum quality value for reads alignment, default=20")

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
    sam2unmap Recover unmap sequence

attention:
    sam2unmap -i input.sam >unmap.fastq
    minimap2 -t 8 -ax sr genome.fa r1.fq r2.fq|sam2unmap >unmap.fastq

version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_hlep_args(parser).parse_args()

    sam2unmap(args.input, args.minmapq)


if __name__ == "__main__":

    main()
