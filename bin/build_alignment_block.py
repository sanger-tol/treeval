#!/usr/bin/env python

import pandas as pd
import os
import sys
import string
import random
import csv
import optparse
import pybedtools
from pybedtools import BedTool


def sort_blocks(df):
    return df.sort_values(["rstart", "qstart"], ascending=[True, True])


def get_block(df, index):
    block = pd.DataFrame([])
    if index < (len(small_cluster_sort.index) - 2):
        if df.iloc[index].loc["qstart"] > df.iloc[index + 1].loc["qstart"]:
            block = df[0 : index + 1]
            leftover = df[index + 1 : len(df.index) - 1]
            qmin = leftover[["qstart"]].min()
            print(qmin)
            index_list = list(range(0, index + 1))
            df.drop(df.index[index_list], inplace=True)

    return block, df


def arrange_fields(df):
    df["fragid"] = df["qchr"].str.cat(df["qstart"].astype(str), sep=":").str.cat(df["qend"].astype(str), sep=":")

    return df[["refchr", "rstart", "rend", "fragid", "qstrand"]]


def build_block(mylist):
    qmin = 0
    qlist = []
    nlist = []

    idx = 0
    while idx < len(mylist):
        x = mylist[idx]
        if idx < len(mylist) - 1:
            qcurrent = int(x[6])
            rcurrent = int(x[1])
            qnext = mylist[idx + 1][6]
            leftd = [int(y[6]) - qmin for y in mylist[idx:]]
            positives = list(filter(lambda x: x > 0, leftd))

            if positives:
                min_value = min(positives)
                indmin = leftd.index(min_value)
                rm = mylist[idx + indmin][1]

                if qcurrent > qmin and qcurrent < qnext and rm == rcurrent:
                    qmin = qcurrent
                    qlist.append(idx)
                elif qcurrent > qmin and qcurrent < qnext and rm > rcurrent:
                    nlist.append(idx)
                elif qcurrent > qmin and qcurrent > qnext or qcurrent < qmin and qcurrent > qnext:
                    nlist.append(idx)
            else:
                idx += 1
                continue

        if idx == len(mylist) - 1:
            if mylist[idx][6] > qmin:
                qlist.append(idx)
            else:
                nlist.append(idx)
        idx += 1

    alignment_chain = [mylist[i] for i in qlist]
    new_list = [mylist[i] for i in nlist]
    return alignment_chain, new_list



#########main##########

parser = optparse.OptionParser(version="%prog 1.0")
parser.add_option(
    "-i",
    "--input",
    dest="input_bedfile",
    default="default.input",
)
parser.add_option(
    "-o",
    "--output",
    dest="out_bedfile",
    default="default.output",
)
options, remainder = parser.parse_args()

inbed = options.input_bedfile
outbed = options.out_bedfile

fo = open(outbed, "a")

sc = pd.read_csv(inbed, sep="\t", comment="#", header=None)

sc.columns = ["refchr", "rstart", "rend", "qchr", "maplen", "qstrand", "qstart", "qend", "qlen"]

ans = [y for x, y in sc.groupby("refchr")]

for mycluster in ans:
    for small_cluster in [y for x, y in mycluster.groupby("qchr")]:
        small_cluster_sort = sort_blocks(small_cluster)

        newdf = small_cluster_sort.reset_index(drop=True)

        newlist = newdf.values.tolist()

        while newlist:
            blocks, newlist = build_block(newlist)

            newblocks = [
                [x if i != 3 else y[3] + ":" + str(y[6]) + ":" + str(y[7]) for i, x in enumerate(y)] for y in blocks
            ]
            a = pybedtools.BedTool(newblocks)
            merged = a.merge(d=100000, c="4,7,8", o="collapse,min,max", delim="|")
            fo.write(str(merged))

fo.close()
