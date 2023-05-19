#!/usr/bin/env python

# Script originally developed by Yumi Sims (yy5@sanger.ac.uk)

from hashlib import new
import optparse

def read_agp (agpfile):
     return {l.split("\t")[5]: l.split("\t")[0:3] for l in open(agpfile, "r")}


def get_pos(agp,mylist):
    myid = mylist[0]
    newpos = [agp[myid][0], int(mylist[1]) + int(agp[myid][1])-1, int(mylist[2]) + int(agp[myid][1]) - 1]
    return newpos


def makestring(myitem):
    return(str(myitem))


def main():
    parser = optparse.OptionParser(version="%prog 1.0")
    parser.add_option(
        "-i",
        "--bed",
        dest="bed",
        default="default.bed",
    )
    parser.add_option(
        "-r",
        "--agp",
        dest="agp",
        default="default.agp",
    )

    options, remainder = parser.parse_args()

    agp = read_agp(options.agp)

    for line in open(options.bed, "r"):
        linelist = line.strip().split("\t")
        newreflist = get_pos(agp, [x for x in linelist[0:3]])
        newqlist = get_pos(agp, [linelist[3], linelist[6], linelist[7]])

        linelist[0:3] = newreflist
        linelist[3] = newqlist[0]
        linelist[6:8] = newqlist[-2:]

        if newreflist != newqlist:
            print("\t".join(map(makestring, linelist)))

if __name__ == "__main__":
    main()