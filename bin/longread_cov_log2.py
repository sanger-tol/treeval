#!/usr/bin/env python

import optparse
import math

# Script originally developed by Will Eagles (we3@sanger.ac.uk)

def process_line(line):
    line_values = line.rsplit(None,1)

    try:
        cov_val = float(line_values[1])
    except:
        cov_val = 0

    if cov_val > 0:
        log_cov_val = math.log2(cov_val)
    else:
        log_cov_val = 0
    
    return line_values[0] + '\t' + str(log_cov_val)

def main():
    parser = optparse.OptionParser(version="%prog 1.0")
    parser.add_option(
        "-i",
        "--inputfile",
        dest="inputfile",
        default="default.input",
    )

    options, remainder = parser.parse_args()

    cov_bed = open(options.inputfile, "r")

    for line in cov_bed:
        print(process_line(line))

if __name__ == "__main__":
    main()