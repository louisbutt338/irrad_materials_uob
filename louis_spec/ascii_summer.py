from dataclasses import dataclass
import json
from math import pi, sqrt, log, exp
import numpy as np


####################################
############ USER INPUTS ###########
####################################


ascii_start_filename = ['fe','au']
ascii_end_filename = '_ubb_280324.Spe'

folder_path = '../short-lived/'

# writes a summed ASCII file from the input ASCII files in the array. Takes the header and footer parameters 
# (i.e. livetime, timings, calibration) from the first file in the array, so this will need to be edited 
# work on this is in progress

####################################
####################################

def parse_ascii(material):

    filename = f"{material}{ascii_end_filename}"
    with open(filename,'r') as ascii_data_file:
        #ascii_contents = ascii_data_file.read().strip().split('\n')
        ascii_contents = ascii_data_file.readlines()
        ascii_header = ascii_contents[:12]
        ascii_footer = ascii_contents[8204:]
    with open(filename,'r') as ascii_data_file:
        ascii_data_strings = ascii_data_file.read().replace(" ", "").strip().split('\n')[12:8204]
        ascii_data = [int(x) for x in ascii_data_strings]
    return ascii_header,ascii_data,ascii_footer

all_ascii_data = []
for m in ascii_start_filename:
    all_ascii_data.append(parse_ascii(m)[1])
ascii_histogram = [sum(x) for x in zip(*all_ascii_data)]

def write_ascii():
    print('writing summed ASCII...')
    filename = f"summed{ascii_end_filename}"
    with open(filename,'w') as ascii_histogram_file:
        for line in parse_ascii(ascii_start_filename[0])[0]:
            ascii_histogram_file.write(line)
        for line in ascii_histogram:
            ascii_histogram_file.write(f"{line}\n")
        for line in parse_ascii(ascii_start_filename[0])[2]:
            ascii_histogram_file.write(line)
        #ascii_histogram_file.write(ascii_histogram)
        #ascii_histogram_file.write(parse_ascii(ascii_start_filename[0])[2])
    
write_ascii()