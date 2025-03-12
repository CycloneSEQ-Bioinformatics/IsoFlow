#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
from os import path
import pandas as pd
import numpy as np
from functools import reduce

def argparser():
    """Argument parser for entrypoint."""
    parser = argparse.ArgumentParser(
    description="""Merge tab separated files on a given field using pandas.""")
    parser.add_argument(
        '-j', metavar='join', help="Join type (outer).", default="outer")
    parser.add_argument(
        '-f', metavar='field',
        help="Join on this field (Reference).", default="Reference")
    parser.add_argument(
        '-o', metavar='out_tsv',
        help="Output tsv (merge_tsvs.tsv).", default="merge_tsvs.tsv")
    parser.add_argument(
        '-z', action="store_true",
        help="Fill NA values with zero.", default=False)
    parser.add_argument(
        '-tpm', type=bool, default=False,
        help="TPM instead of counts")
    parser.add_argument(
        '-tsvs', metavar='input_tsvs', nargs='*',
        help="Input tab separated files.")

    return parser


def main(args):
    """Run entry point."""
    dfs = {x: pd.read_csv(x, sep="\t") for x in args.tsvs}

    ndfs = []
    for x, df in dfs.items():
        # Transform counts to integers:
        if args.tpm:
            df = df.rename(columns={'TPM': 'Count', 'Name': args.f})
        else:
            df = df.rename(columns={'NumReads': 'Count', 'Name': args.f})
            df.Count = np.array(df.Count, dtype=int)
        # Take only non-zero counts:
        df = df[df.Count > 0]
        df = df[[args.f, "Count"]]
        df = df.sort_values(by=["Count"], ascending=False)
        name = x.split('/')[1]
        df = df.rename(columns={'Count': name})
        ndfs.append(df)
    dfs = ndfs

    df_merged = reduce(lambda left, right: pd.merge(
        left, right, on=args.f, how=args.j), dfs)
    if args.z:
        df_merged = df_merged.fillna(0)

    df_merged.to_csv(args.o, sep="\t", index=False)

if __name__ == '__main__':
    parser = argparser()
    args = parser.parse_args()
    main(args)
