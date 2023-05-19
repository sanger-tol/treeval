#!/usr/bin/env python

# Script originally developed by Yumi Sims (yy5@sanger.ac.uk)

import itertools, datetime, os, re, sys, time
import optparse
import sys
import io
import string
import random

def deal_headline (headline):
    direction = "+"
    qchrom = ""
    qlen = ""
    hl_comp = headline.strip().split(" ")
    if "Reverse" in headline:
        direction = "-"
        qchrom = hl_comp[1]
        qlen = hl_comp[6]
    else:
        qchrom = hl_comp[1]
        qlen = hl_comp[5]
    return qchrom + "\t" + str(qlen) + "\t" + direction


def deal_line(line):
    (refid, refpos, qpos, lmotif, qid, qlen, strand) = line.split("\t")
    refend = str(int(refpos) + int(lmotif) - 1)
    if strand == "+":
        qend = str(int(qpos) + int(lmotif) - 1)
    else:
        qend = str(int(qpos) - int(lmotif) + 1)
    return (
        refid
        + "\t"
        + refpos
        + "\t"
        + refend
        + "\t"
        + qid
        + "\t"
        + lmotif
        + "\t"
        + strand
        + "\t"
        + qpos
        + "\t"
        + qend
        + "\t"
        + qlen
    )   
            


parser = optparse.OptionParser(version="%prog 1.0")
parser.add_option(
    "-i",
    "--input",
    dest="input_mummerfile",
    default="default.input",
)
parser.add_option(
    "-l",
    "--motiflen",
    dest="motiflen",
    default="default.motiflen",
)

options, remainder = parser.parse_args()

motiflen = options.motiflen

inFile = open(options.input_mummerfile, "r")

keepCurrentSet = False
headline = ""
bedlist = []
for line in inFile:
    if line.startswith(">"):
        keepCurrentSet = False
    if keepCurrentSet:
        mypattern = re.compile("(\S+)\s+(\d+)\s+(\d+)\s+(\d+)$")
        num = mypattern.findall(line.strip())
        refid = num[0][0]
        refpos = num[0][1]
        qpos = num[0][2]
        mlen = num[0][3]
        if int(mlen) > int(motiflen):
            bedline = refid + "\t" + refpos + "\t" + qpos + "\t" + mlen + "\t" + deal_headline(headline)
            bedlist.append(bedline)
    if line.startswith(">"):
        keepCurrentSet = True
        headline = line

inFile.close()

for line in bedlist:
    print(deal_line(line))