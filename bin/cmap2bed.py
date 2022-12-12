#!/usr/bin/env python 

import optparse
import sys
import io
import string

parser = optparse.OptionParser()
parser.add_option('-t', '--inputfile', 
                  dest="inputfile", 
                  default="default.input",
                  )

parser.add_option('-z', '--enzyme', 
                  dest="enzyme", 
                  default="default.enzyme",
                  )

options, remainder = parser.parse_args()

enzyme=options.enzyme

def join2lines(previous_line, current_line):
    return ((previous_line.strip()+"\t"+current_line.strip()).split("\t"))

def get_fields(line):
    mylist=line.split("\t")
    return mylist[0]+"\t"+mylist[5]

def reformat_cmap(cmap,enzyme):

    my_file = open(cmap, "r")
    my_new_file=[]
    if my_file:
        current_line = my_file.readline()

    for line in my_file:
        previous_line = current_line
        current_line = line
        mylinelist=join2lines(get_fields(previous_line), get_fields(current_line))
        if mylinelist[0] == mylinelist[2]:
            my_new_file.append(mylinelist[0]+"\t"+str(int(float(mylinelist[1])))
            +"\t"+str(int(float(mylinelist[3])))+"\t"+enzyme
            +"\t"+str(int(float(mylinelist[3]))-int(float(mylinelist[1]))))
        else:
            my_new_file.append(mylinelist[2]+"\t"+"0"
            +"\t"+str(int(float(mylinelist[3])))+"\t"+enzyme+"\t"
            +str(int(float(mylinelist[3]))))

    firtLineList = my_new_file[0].split("\t")
    firtline=firtLineList[0]+"\t"+"0"+"\t"+firtLineList[1]+"\t"+enzyme+"\t"+firtLineList[1]
    my_new_file.insert(0, firtline)
    return my_new_file

mymap=reformat_cmap(options.inputfile, enzyme)

for line in mymap:
    print (line)