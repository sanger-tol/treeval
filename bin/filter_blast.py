#!/usr/bin/env python3

"""
Filter BLAST
------------
Author: dp24 / DLBPointon
Contact: dp24@sanger.ac.uk / dlb_pointon@hotmail.co.uk

A small script to take an input blast file
and filter it's contents to given percentage.

qstart & qend columns are then used to calculate strand direction.
pident is rounded up the next whole number.
sseqid is split to give transcript name and ensemble name columns.

input_file by default should be -outfrmt 6

Usage:
python3 filter_blast.py organism_name, data_type, file_location, percentage_to_filter

"""

import os
import sys
import pandas as pd
import numpy as np

pd.options.mode.chained_assignment = None
BLAST_HEADER = ['qseqid', 'sseqid','pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend','sstart','send','evalue','bitscore']

NEW_FORMAT = ['sseqid', 'sstart', 'send', 'qseqid_1', 'pident_2', 'strand', 'qseqid_2', 'qseqid_3']
NOTE = 'SORT BY COLUMN 0 AND 1'

def arg_check(arg_list):
    org = str(arg_list[1])
    dtype = str(arg_list[2])
    filt_percent = float(arg_list[4])

    if os.path.exists(arg_list[3]):
        input_file = arg_list[3]
    else:
        print("Can't find file")
        sys.exit()

    return org, dtype, input_file, filt_percent

def load_blast(input_file: str, blst_head: list):
    return pd.read_csv(input_file, '\t', names=blst_head)

def filter_blast(df, filt_percent: float):
    return df[(df['pident'] > filt_percent)]

def split_sseqid(df):
    return df.qseqid.str.extract('(?P<qseqid_1>\S+)\((?P<qseqid_2>\S+)\)', expand=True)

def rounded_df(df):
    return df['pident'].round(decimals = 0).astype(int)

def plus_neg(df):
    return np.where(df['sstart'].astype(int) < df['send'].astype(int), '+', '-')

def save_file(org: str, dtype: str, df, go_on: bool):
    if go_on:
        df.to_csv(f"{org}-{dtype}-filtered90.tsv", sep='\t', header=False, index=False)
    else:
        df.to_csv(f"{org}-{dtype}-EMPTY.tsv", sep='\t', header=False, index=False)

def main():
    arg_list = sys.argv
    print(arg_list)
    go_on = True

    org, dtype, input_file, filt_percent = arg_check(arg_list)

    if os.path.getsize(arg_list[3]) == 0:
        print("FILE SIZE: 0 - KILLING JOB")
        go_on = False
        final_blast = pd.DataFrame(columns=['EMPTY'])
    else:
        print(f"FILE SIZE: {os.path.getsize(arg_list[3])} - RUNNING JOB")

        blast_df = load_blast(input_file, BLAST_HEADER)
        # DROP COLUMN FUNCTION HERE?
        final_blast = filter_blast(blast_df, filt_percent)

        # Producing columns now
        pi_rounded = rounded_df(final_blast)
        stranded_df = plus_neg(final_blast)
        reorg_blast = split_sseqid(final_blast)

        # Stitch together the computed columns
        final_blast['pident_2'] = pi_rounded
        final_blast['strand'] = stranded_df
        final_blast[['qseqid_1','qseqid_2']] = reorg_blast
        final_blast['qseqid_3'] = final_blast['qseqid_2'].copy()

        final_blast = final_blast.reindex(columns=NEW_FORMAT)

        # Identifies where strand = '-' swap values in start and end sequence column.
        negative_strand = final_blast['strand'] == '-'
        final_blast.loc[negative_strand, ['sstart', 'send']] = (
            final_blast.loc[negative_strand, ['send', 'sstart']].values
            )
        final_blast[['sstart', 'send']] = final_blast[['sstart', 'send']].astype(int)

        final_blast = final_blast.sort_values(by=['sseqid','sstart'], ascending=True)


    print('Saving File')
    save_file(org, dtype, final_blast, go_on)

if __name__ == '__main__':
    main()
