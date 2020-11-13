#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import sys
import gzip
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


def filter_paf(file, identities=50, mapq=1):

    if file == "":
        fh = read_stdin('\t')
    elif file.endswith(".paf"):
        fh = read_tsv(file, '\t')
    else:
        raise Exception("%r file format error" % file)

    for line in fh:
        match = int(line[9])
        blen = int(line[10])
        match = min(blen, match)

        if (match*100.0/blen) < identities or int(line[11]) < mapq:
            continue
        print('\t'.join(line))


def add_hlep_args(parser):

    parser.add_argument("-i", "--input", metavar='FILE', type=str, default='',
        help="Input files in paf")
    parser.add_argument("-id", "--identities", metavar='FLOAT', type=float, default=50.0,
        help="Set the minimum identities of reads comparison, default=50.0")
    parser.add_argument("-mq", "--mapq", metavar='INT', type=int, default=2,
        help="Set the minimum quality value for reads alignment, default=2")

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
    filter_paf.py  Filter the reads comparison result file

attention:
    filter_paf.py -i file.paf >file_clean.paf
    cat file.paf | filter_paf.py >file_clean.paf

version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_hlep_args(parser).parse_args()

    filter_paf(args.input, args.identities, args.mapq)


if __name__ == "__main__":

    main()
