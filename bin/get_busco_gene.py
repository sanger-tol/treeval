#!/usr/bin/env python
import pandas as pd
import os
import optparse


# Script originally developed by Yumi Sims (yy5@sanger.ac.uk)

parser = optparse.OptionParser(version="%prog 1.0")
parser.add_option('-i', '--input_fulltable', 
                  dest="input_fulltable", 
                  default="default.input",
                  )

parser.add_option('-o', '--output', 
                  dest="output", 
                  default="default.output",
                  )



options, remainder = parser.parse_args()

fulltable=options.input_fulltable
genefile=options.output



full_table = pd.read_csv(
    fulltable,
    sep="\t", 
    comment='#'
)

full_table.rename(columns = {'Busco id':'buscoID'}, inplace = True)

print(full_table)

column_names = list(full_table.columns)

df_a = full_table.loc[:,'Sequence':'Gene End']

print(df_a)



df_new=df_a.join(full_table[["buscoID"]]).join(full_table[["Score"]]).join(full_table[["Strand"]]).join(full_table[["OrthoDB url"]])


df_new.fillna('NA',inplace = True)

dfnoNa=df_new[df_new.Sequence != 'NA']

df_final = dfnoNa.reset_index(drop=True)

df_final = df_final.astype({'Gene End':'int','Gene Start':'int'})

df_final.to_csv(genefile, 
               index=False, 
               header=False,
               sep="\t")
