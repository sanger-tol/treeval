"""
---- Gene Alignment CSV Generator ----
      By Damon-Lee Pointon (dp24)

This script generates the csv files
 required by TreeVal (the gene alignment
 sub-workflows).
 
Script generates a csv per organism.assembly 
 in their respective classT folder

USAGE:
    python3 GA_data_sorting.py
    
TODO: ADD argparse + version data
"""

import argparse
import sys
import os
from os import path

def list_dir(dir_loc: str):
    return [os.path.join(dir_loc, i) for i in os.listdir(dir_loc) if path.isdir(os.path.join(dir_loc, i))]

def get_file_list(root: str):
    return [os.path.join(path, name) for path, subdirs, files in os.walk(root) for name in files]

def list_2_dict (file_list: list):
    file_dict = {}
    path_list = []
    for i in file_list:
        path_list = i.split('/')
        if path_list[-1].lower() in ["readme.txt", "readme"]:
            pass
        else:
            file_dict[path_list[-1]] = [path_list[-3], path_list[-2], i]
    return file_dict, path_list[-3]

def save_data(dict_of_data: dict, save_loc: str, org_accession: str):
    save_path = f'{save_loc}/csv_data/{org_accession}-data.csv'
    if os.path.exists(save_path):
        os.remove(save_path)
    else:
        pass
    print(f'Generating CSV for:\t{org_accession}\nSave Path:\t\t{save_path}')
    with open(save_path, 'w+') as new_csv:
        new_csv.write('org,type,data_file')
        for x, y in dict_of_data.items():
            new_csv.write(f'\n{y[0]},{y[1]},{y[2]}')

def main():
    gene_alignment_dir = sys.argv[1]
    clade_list = list_dir(gene_alignment_dir)
    master_list = []
    
    for i in clade_list:
        org_list = list_dir(i)
        for ii in org_list:
            accession_list = list_dir (ii)
            for iii in accession_list:
                print(f'============> {iii.split("/")[-1]} -- {i.split("/")[-1]}')
                data_list = list_dir (iii)
                master_list = []
                master_list += [get_file_list(j) for j in data_list]
                file_dict, org = list_2_dict([item for sublist in master_list for item in sublist])
                save_data(file_dict, i, org)
            print('GA_csv_gen: Complete')

if __name__ == '__main__':
    main()