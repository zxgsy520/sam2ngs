#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import sys
import logging
import argparse

LOG = logging.getLogger(__name__)

__version__ = "1.0.0"
__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__all__ = []


def stdin_read_tsv(file, sep=None):

    if file == "":
        fh = sys.stdin
    else:
        fh = open(file)

    for line in fh:
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        yield line.split(sep)


def ngs_sam2unmapid(file, dtype="mgi"):

    unmapids = set()

    for line in stdin_read_tsv(file, '\t'):
        if line[2] != "*" or line[5] != "*":
            continue
        unmapids.add(line[0])

    for i in unmapids:
        if dtype=="mgi":
            print("%s/1" % i)
            print("%s/2" % i)
        else:
            print(i)

    return 0


def add_hlep_args(parser):

    parser.add_argument("-s", "--sam", metavar='FILE', type=str, default='',
        help="Input the compared sam file.")
    parser.add_argument("--dtype", metavar='STR', type=str, choices=["mgi", "illumina", "other"], default="illumina",
        help='Set up the sequencing platform of the data, default=illumina.')

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
    ngs_sam2unmapid.py  Extract the unmap sequence id from the sam file.

attention:
    ngs_sam2unmapid.py -s file.sam >unmap.id
    samtools view file.sam |ngs_sam2unmapid.py --dtype mgi >unmap.id

version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_hlep_args(parser).parse_args()

    ngs_sam2unmapid(args.sam, args.dtype)


if __name__ == "__main__":

    main()
