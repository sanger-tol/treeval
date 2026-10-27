#!/usr/bin/env python3
import pandas as pd
import optparse

# Script originally developed by Yumi Sims (yy5@sanger.ac.uk)
# -------------------
# Update for BUSCO 5.5.0 - by we3 (Will Eagles)
# Reorder start and end so smallest always second column. Also, trim range from scaffold name in first column.
# -------------------

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

df_final["Sequence"] = df_final["Sequence"].str.replace(r":.*", "", regex=True)

df_final[["Gene Start", "Gene End"]] = df_final.apply(
    lambda row: (
        (row["Gene Start"], row["Gene End"])
        if row["Gene Start"] < row["Gene End"]
        else (row["Gene End"], row["Gene Start"])
    ),
    axis=1,
    result_type="expand",
)
df_final.to_csv(csvfile, index=False, header=False, sep="\t")
