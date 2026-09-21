#!/usr/bin/env python
import sys
import argparse


# Script originally developed by Charlotte Wright (cw22@sanger.ac.uk)
# -------------------
# Update for BUSCO 5.5.0 - by we3 (Will Eagles)
# Reorder start and end so smallest always second column. Also, trim range from scaffold name in first column.
# -------------------


def parse_table(table_file):
    with open(table_file, "r") as table:
        table_dict, chr_list = {}, []
        for line in table:
            if not line.startswith("#"):
                cols = line.rstrip("\n").split()
                buscoID, status = cols[0], cols[1]
                if status == "Complete":  # busco_type can either be "Complete" or "Duplicated"
                    chr = cols[2].split(":", 1)[0]
                    if int(cols[3]) > int(cols[4]):
                        start, stop = int(cols[4]), int(cols[3])
                    else:
                        start, stop = int(cols[3]), int(cols[4])

                    table_dict[buscoID] = [chr, start, stop]
                    if not chr in chr_list:
                        chr_list.append(chr)
    return table_dict, sorted(chr_list)


def parse_query_table(table_file):  # use this function only if interested in assiging duplicated buscos
    with open(table_file, "r") as table:
        table_dict, table_dict_dup, chr_list = {}, {}, []
        for line in table:
            if not line.startswith("#"):
                cols = line.rstrip("\n").split()
                buscoID, status = cols[0], cols[1]
                if status == "Duplicated":  # busco_type can either be "Complete" or "Duplicated"
                    chr = cols[2].split(":", 1)[0]
                    if int(cols[3]) > int(cols[4]):
                        start, stop = int(cols[4]), int(cols[3])
                    else:
                        start, stop = int(cols[3]), int(cols[4])

                    if buscoID in table_dict.keys():
                        table_dict_dup[buscoID] = [chr, start, stop]
                    else:
                        table_dict[buscoID] = [chr, start, stop]
                    if not chr in chr_list:
                        chr_list.append(chr)
    return table_dict, table_dict_dup, sorted(chr_list)


def print_summary_table(reference_table_dict, reference_chr_list, query_table_dict, query_chr_list, summary_table_file):
    with open(summary_table_file, "w") as summary_table:
        summary_table.write(
            "%s\t%s\t%s\t%s\t%s\t%s"
            % (
                "query_chr",
                "assigned_chr",
                "assigned_chr_busco_count",
                "total_query_chr_busco_chr",
                "perc",
                "\t".join("ref_" + (x) + "_count" for x in reference_chr_list),
            )
            + "\n"
        )
        assignment_dict = {}
        top_chr_dict = {}
        for chr in query_chr_list:
            assignment_dict[chr] = []
            for i in range(0, len(reference_chr_list)):
                assignment_dict[chr].append(0)
        for buscoID, position_list in query_table_dict.items():
            query_chr, query_start, query_stop = position_list
            try:
                reference_chr = reference_table_dict[buscoID][0]
            except KeyError:
                pass
            assignment_dict[query_chr][reference_chr_list.index(reference_chr)] += 1
        for chr, count_list in assignment_dict.items():
            top_chr, top_chr_count, total = "", 0, 0
            for reference_chr in reference_chr_list:
                index = reference_chr_list.index(reference_chr)
                count = count_list[index]
                total += count
                if not top_chr == "":
                    if count > top_chr_count:
                        top_chr, top_chr_count = reference_chr, count
                else:
                    top_chr, top_chr_count = reference_chr, count
            perc = round((top_chr_count / total), 2)
            summary_table.write(
                "%s\t%s\t%s\t%s\t%s\t%s"
                % (chr, top_chr, top_chr_count, total, perc, "\t".join(str(x) for x in count_list))
                + "\n"
            )
            top_chr_dict[chr] = top_chr
    return top_chr_dict


def print_location_table(reference_table_dict, query_table_dict, location_table_file, top_chr_dict):
    with open(location_table_file, "w") as location_table:
        location_table.write(
            "%s\t%s\t%s\t%s\t%s" % ("buscoID", "query_chr", "position", "assigned_chr", "status") + "\n"
        )
        for buscoID, position_list in query_table_dict.items():
            query_chr, query_start, query_stop = position_list
            position = (query_start + query_stop) / 2
            try:
                reference_chr = reference_table_dict[buscoID][0]
                if reference_chr != top_chr_dict[query_chr]:
                    status = reference_chr
                else:
                    status = "self"
                location_table.write(
                    "%s\t%s\t%s\t%s\t%s" % (buscoID, query_chr, position, reference_chr, status) + "\n"
                )
            except KeyError:
                pass


def print_dups_location_table(
    reference_table_dict, query_table_dict, query_table_dict2, location_table_file, top_chr_dict
):
    with open(location_table_file, "w") as location_table:
        location_table.write(
            "%s\t%s\t%s\t%s\t%s\t%s" % ("buscoID", "query_chr", "query_start", "query_end", "assigned_chr", "status")
            + "\n"
        )
        for buscoID, position_list in query_table_dict.items():
            query_chr, query_start, query_stop = position_list
            position = (query_start + query_stop) / 2
            try:
                reference_chr = reference_table_dict[buscoID][0]
                if reference_chr != top_chr_dict[query_chr]:
                    status = reference_chr
                else:
                    status = "self"
                location_table.write(
                    "%s\t%s\t%s\t%s\t%s\t%s" % (buscoID, query_chr, query_start, query_stop, reference_chr, status)
                    + "\n"
                )
            except KeyError:
                pass
        for buscoID, position_list in query_table_dict2.items():
            query_chr, query_start, query_stop = position_list
            position = (query_start + query_stop) / 2
            try:
                reference_chr = reference_table_dict[buscoID][0]
                if reference_chr != top_chr_dict[query_chr]:
                    status = reference_chr
                else:
                    status = "self"
                location_table.write(
                    "%s\t%s\t%s\t%s\t%s\t%s" % (buscoID, query_chr, query_start, query_stop, reference_chr, status)
                    + "\n"
                )
            except KeyError:
                pass


if __name__ == "__main__":
    SCRIPT = "buscopainter.py"
    # argument set up
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--version", action="version", version="1.0")
    parser.add_argument(
        "-r", "--reference_table", type=str, help="full_table.tsv file for reference species", required=True
    )
    parser.add_argument("-q", "--query_table", type=str, help="full_table.tsv for query species", required=True)
    parser.add_argument("-p", "--prefix", type=str, help="prefix for output file names", default="buscopainter")
    args = parser.parse_args()
    reference_table_file = args.reference_table
    query_table_file = args.query_table
    prefix = args.prefix
    summary_table_file = prefix + "_summary.tsv"
    location_table_file = prefix + "_complete_location.tsv"
    # run the functions
    reference_table_dict, reference_chr_list = parse_table(reference_table_file)  # get all complete buscos
    query_table_dict, query_chr_list = parse_table(query_table_file)  # get all complete buscos as usual
    top_chr_dict = print_summary_table(
        reference_table_dict, reference_chr_list, query_table_dict, query_chr_list, summary_table_file
    )  # work out top_chr dict using complete buscos as usual
    print("[+] Written output files successfully:")
    print("[+]	" + summary_table_file)
    print("[+]	" + location_table_file)
    query_table_dup1_dict, query_table_dup2_dict, query_chr_dup_list = parse_query_table(
        query_table_file
    )  # get all duplicated buscos
    print_location_table(reference_table_dict, query_table_dict, location_table_file, top_chr_dict)
    location_table_file = prefix + "_duplicated_location.tsv"
    print_dups_location_table(
        reference_table_dict, query_table_dup1_dict, query_table_dup2_dict, location_table_file, top_chr_dict
    )
    print("[+]	" + location_table_file)
