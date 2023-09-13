#!/usr/local/bin/python

PRINT_ERROR = """Does not exist\n
					Get module installed before import attempt\n
					If running server side then contact your admin"""

try:
    import sys

    if sys.version_info[0] < 3 and sys.version_info[1] < 6:
        raise Exception(
            """Must be using Python 3.6 for the full
						functionality of this script"""
        )
    if sys.version_info[0] >= 3 and sys.version_info[1] >= 6:
        print("Your using at least Version 3.6, You are good to go...")
except ImportError:
    print(f"sys not imported \n {PRINT_ERROR}")
    sys.exit(0)

try:
    import os

    print("os imported")
except ImportError:
    print(f"os not imported \n {PRINT_ERROR}")
    sys.exit(0)

try:
    import argparse

    print("argparse imported")
except ImportError:
    print(f"argparse not imported \n {PRINT_ERROR}")
    sys.exit(0)

try:
    import re

    print("regex imported")
except ImportError:
    print(f"re not imported \n {PRINT_ERROR}")
    sys.exit(0)

DOCSTRING = """
----------------------------------------------------------------
					Gene Alignment Data Prep
					By Damon-Lee Pointon (dp24)
----------------------------------------------------------------
This script takes an input file and chunks it into 1000 sequence
files for cds and rna files or a user defined number (default 100)
for pep and cdna sequences.

----------------------------------------------------------------
Usage:

GA_data_prep.py {name}-{accession}.{datatype}.fasta {ncbi|ens} {chunk}

File name example is: ThalassiosiraPseudonana-ASM14940v2.rna.fasta
----------------------------------------------------------------
"""


def get_command_args(args=None):
    parser = argparse.ArgumentParser(
        prog="GA_data_prep.py (Python 3)", description=DOCSTRING, formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument("FASTA", action="store", help="Input unzipped fasta file", type=str)

    parser.add_argument("DB", action="store", help="database of origin", choices=["ncbi", "ens"], type=str)

    parser.add_argument("CHUNK", action="store", help="Chunk size for pep and cdna", default=100, type=int)

    parser.add_argument("-v", "--version", action="version", version="v7.0.0")

    options = parser.parse_args(args)
    return options


def entryfunction(file, dtype, org, entryper, filesavedto, options):
    """
    The entryfunction function splits a FASTA file into a defined
    number of entries per file, pep == 2000 enteries and everything
    else is split into 5000 enteries.
    :param seq_file:
    :param org:
    :param directory:
    :param entryper:
    :param option:
    """
    print("Entryfunction called")
    count = 0
    filecounter = 0
    entry = []
    print(file)
    if os.path.exists(file):
        print("File found at %s", file)
        with open(file, "r") as filetoparse:
            print("Renaming headers")
            for name, seq in read_fasta(filetoparse):
                new_name = massage(name, options)
                # print(new_name)  # Here as a manual check of headers
                nameseq = new_name, seq
                entry.append(nameseq)
                count += 1

                if count == entryper:
                    filecounter += 1
                    with open(f"{filesavedto}/{org}{filecounter}{dtype}.MOD.fa", "w") as done:
                        for head, body in entry:
                            done.write(f"{head}\n{body}\n")
                        count = 0
                        entry = []
                    print(f"File saved: -- {filesavedto}/{org}{filecounter}{dtype}.MOD.fa")

                filecounter += 1

            with open(f"{filesavedto}/{org}{filecounter}{dtype}.MOD.fa", "w") as done:
                for head, body in entry:
                    done.write(f"{head}\n{body}\n")
                entry = []

            print(f"File saved: -- {filesavedto}/{org}{filecounter}{dtype}.MOD.fa")


def massage(name, options):
    """
    A function to 'massage' the sequence headers into a more human
    readable style
    :param option:
    :param name:
    :return name:
    """

    if name.startswith(">"):
        if options.DB == "ncbi":
            gene_symbol = re.search(r"gene=([A-Z]\w+)", name)
            ens_code = re.search(r"GeneID:([1-9])\w+", name)
        else:
            gene_symbol = re.search(r"symbol:(\S+)", name)
            ens_code = re.search(r"ENS(\w+)T(\w+.\d+)", name)

        if gene_symbol:
            gene_symbol = gene_symbol.group(1)
        elif gene_symbol is None:
            gene_symbol = re.search(r"gene:(\S+)", name)

            if gene_symbol:
                gene_symbol = gene_symbol.group(1)

            elif gene_symbol is None:
                gene_symbol = re.search(r"PREDICTED: (.+) \[", name)
                if gene_symbol:
                    gene_symbol = gene_symbol.group(1)
                    gene_symbol = gene_symbol.split()
                    gene_symbol = "_".join(gene_symbol)
                else:
                    gene_symbol = "MissingInfo"

        if ens_code:
            ens_code = ens_code.group(0)

        elif ens_code is None:
            ens_code = re.search(r">(\S+)", name)
            if ens_code:
                ens_code = ens_code.group(1)
            elif ens_code is None:
                ens_code = "NoEnsCode"

        # print('Gene Symbol found as: %s', gene_symbol)
        # print('Ens Code found as: %s', ens_code)
        if gene_symbol == "MissingInfo":
            # print('MissingInfo replaced with %s', ens_code)
            gene_symbol = ens_code
        name = f">{gene_symbol}({ens_code})"

    else:
        print("Somethings gone wrongs, headers are wrong")
        sys.exit(0)

    return name


def read_fasta(filetoparse):
    """
    A function which opens and splits a fasta into name and seq.
    :param filetoparse:
    """
    print("Read_fasta called")
    counter = 0
    name, seq = None, []

    for line in filetoparse:
        line = line.rstrip()

        if line.startswith(">"):
            if name:
                yield name, "".join(seq)
            name, seq = line, []
        else:
            seq.append(line)

    if name:
        yield name, "".join(seq)
        counter += 1


def main():
    options = get_command_args()
    file = options.FASTA
    dtype = file.split(".")[1]
    org = file.split(".")[0]

    print(f"WORKING ON:\t\t{dtype}--{org}")

    directory = f"./{org.split('-')[0]}/{org.split('-')[0]}.{org.split('-')[1]}/{dtype}"

    try:
        os.makedirs(directory, mode=0o777)
    except:
        print("probably already exists")

    entryper = [1000 if dtype in ["cds", "rna"] else options.CHUNK]

    print(f"Records per file:\t{int(entryper[0])}")
    entryfunction(file, dtype, org.split("-")[0], int(entryper[0]), directory, options)


if __name__ == "__main__":
    main()
