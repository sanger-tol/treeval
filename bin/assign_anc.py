#!/usr/bin/env python3
import pandas as pd
import optparse


# Script originally developed by Yumi Sims (yy5@sanger.ac.uk)

parser = optparse.OptionParser(version="%prog 1.0")
parser.add_option(
    "-l",
    "--locationfile",
    dest="locationfile",
    default="default.locationfile",
)

parser.add_option(
    "-f",
    "--fulltable",
    dest="fulltable",
    default="default.fulltable",
)

parser.add_option(
    "-c",
    "--csvfile",
    dest="csvfile",
    default="default.csvfile",
)

options, remainder = parser.parse_args()
locationfile = options.locationfile
fulltable = options.fulltable
csvfile = options.csvfile

location = pd.read_csv(locationfile, sep="\t", comment="#")

full_table = pd.read_csv(fulltable, sep="\t", header=None, comment="#")

fulltable_colnames = [
    "buscoID",
    "Status",
    "Sequence",
    "Gene Start",
    "Gene End",
    "Strand",
    "Score",
    "Length",
    "OrthoDB url",
    "Description",
]

full_table.columns = fulltable_colnames

df = location.merge(full_table, on="buscoID")

df_a = df.loc[:, "Sequence":"Gene End"]

df_new = df_a.join(df[["assigned_chr"]]).join(df[["Score"]]).join(df[["Strand"]]).join(df[["OrthoDB url"]])

df_new.fillna("NA", inplace=True)

dfnoNa = df_new[df_new.Sequence != "NA"]

df_final = dfnoNa.reset_index(drop=True)

df_final = df_final.astype({"Gene End": "int", "Gene Start": "int"})

df_final.to_csv(csvfile, index=False, header=False, sep="\t")
