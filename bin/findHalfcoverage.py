#! /usr/bin/env python3

import re
import sys
from optparse import OptionParser


def load_scafsize(file) :
	# example is my.genome file, "scaffold\tsize"

	scafkey = {}
	scaffile = open(file, "r")
	for line in scaffile:
		line = line.replace("\n", "")
		name, size = re.split("\t", line)
		scafkey[name] = size

	scaffile.close()
	return scafkey


def getTotallength_undercov(file, cov, wiggleroom):
	# example is bed file of coverage,
	# scaffold_100_arrow      0       2       18

	coverage_cutoff = cov + wiggleroom

	myfile = open(file, "r")

	lowcoverage_sum = 0
	prev_scaf = ""
	scaf_lc = {}

	for line in myfile:
		line = line.replace("\n", "")
		objContents = re.split("\t", line)

		if (prev_scaf != objContents[0]):
			scaf_lc[prev_scaf] = lowcoverage_sum
			lowcoverage_sum = 0

		if (float(objContents[3]) < coverage_cutoff):
			length = float(objContents[2]) - float(objContents[1])
			lowcoverage_sum += length

		prev_scaf = objContents[0]

	scaf_lc[prev_scaf] = lowcoverage_sum
	myfile.close()

	return scaf_lc


def get_cov_peaks(file):
	# example is depthgraph.txt, "coverage\tbasepair count"

	myPeakFile = open(file, "r")

	rows = []
	for line in myPeakFile:
		line = line.replace("\n", "")
		items = re.split("\t", line)
		rows.append(items)

	myPeakFile.close()
	# print(rows[0])
	peakCov = sorted(rows, key=lambda cov: int(cov[1]), reverse=1)[0][0]

	if int(peakCov) == 0:
		peakCov = sorted(rows, key=lambda cov: int(cov[1]), reverse=1)[1][0]

	halfPeak = int(peakCov) / 2
	qrtPeak = int(peakCov) / 4

	print("#Coverage Peak is %s, HalfPeak is %s, QuarterPeak is %s " % (peakCov, halfPeak, qrtPeak) )
	
	return(peakCov, halfPeak, qrtPeak)


def calc_coverage(scafsize, totallowcov) :
	# calculate the % for lowcov coverage over entire scaffold.
	return totallowcov / scafsize*100


def getArguments():
	# get indivudual arguments from user

	parser = OptionParser(version="%prog 1.0")
	parser.add_option(
    	"-c",
    	"--coveragefile",
    	action="store", 
     	type="string",
      	dest="covfile",
        help="Scaffold Coverage filename"
    )
	parser.add_option(
    	"-m",
    	"--mygenome",
     	action="store",
      	type="string",
       	dest="mygenome",
        help="mygenome file, scaffold - size file"
    )
	parser.add_option(
    	"-d",
     	"--depthgraph", 
        action="store",
        type="string",
        dest="depth", 
        help="depthgraph file, bp count at each depth"
    )
	parser.add_option(
     	"-w",
    	"--wiggle",
     	action="store",
      	type="float",
       	dest="wig",
        default=5,
        help="wiggle room to add to depth cutoff ie 30X + wiggleroom.  Default is 5X"
    )
	parser.add_option(
     	"--cut", 
        action="store",
        type="float",
        dest="covcut",
        default=60,
        help="%Number for coverage cutoff to include in results.  ie 50% of scaffold needs to be under diploid peak etc.  Default is 60%"
	)	
	parser.add_option(
    	"-t",
     	"--totalsize",
		action="store",
  		type="int",
    	dest="totsize",
     	default=250000,
        help="total size that determines max coverage boundary."
    )

	(options, args) = parser.parse_args()

	if (options.covfile == None or options.mygenome == None or options.depth == None):
		print("Missing Options")
		exit()

	return options

def main():
	# main program	

	options = getArguments()		

	scaffold_sizes = load_scafsize(options.mygenome)
	(hapCov, dipCov, tetCov) = get_cov_peaks(options.depth)
	scaffold_lowcovsum = getTotallength_undercov(options.covfile, dipCov, options.wig)

	for scaffoldName in scaffold_lowcovsum:
		if (scaffoldName == ""):
			continue

		# print("==" + scaffoldName)
		totalSize = float(scaffold_sizes[scaffoldName])
		lowcovSize = float(scaffold_lowcovsum[scaffoldName])

		coverage = calc_coverage(totalSize, lowcovSize)

		if (coverage > options.covcut):

			if (totalSize > options.totsize):
				print( "**\t" + "\t".join([str(i) for i in [scaffoldName, int(totalSize), int(lowcovSize), "{:.1f}".format(coverage)]]))		
			else :
				print( "==\t" + "\t".join([str(i) for i in [scaffoldName, int(totalSize), int(lowcovSize), "{:.1f}".format(coverage)]]))


# -- script execuation -- #
if __name__ == "__main__":
	main()