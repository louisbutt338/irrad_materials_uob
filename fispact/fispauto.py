import os
import shutil
import numpy as np
import sys
import matplotlib.font_manager
import matplotlib
import matplotlib.pyplot as plt
#from sklearn.linear_model import LinearRegression
from scipy.stats import gaussian_kde
import seaborn
#plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = "Helvetica"
plt.rcParams["font.size"] = 22
plt.rcParams["font.weight"] = "normal"
from collections import Counter
from itertools import count
from random import randint
from pathlib import Path
import math

######################
######################
######################
######################

#materials = ['au', 'al', 'fe', 'in', 'nb', 'ni', 'rh', 'sc', 'y', 'dy', 'cd', 'cu'] # all mat ['li','ss','v44','v','au', 'al', 'fe', 'in', 'nb', 'ni', 'rh', 'sc', 'y', 'dy', 'cd','cu'] # other materials: ['ss', 'v44','v','li','lipb','macor']
materials = ['au', 'al', 'fe', 'in', 'nb', 'ni', 'rh', 'sc', 'y', 'dy', 'cd', 'cu']
ubb_or_other = 'ubb'
projectile = 'neutron' # select: proton neutron
irrad_time = [120] #[20,36.3,67,10.3,41.5] #[60*24] # set irrad times in minutes (if just one then do array of one)
dose_input = '2 0.3' # 1 or 2 0.3

# note: Set the stoich/mass/density in an original FISPACT {material}.i file. Need input files for each material in /{materials} folder

# for microbreeder flux:
fluxes = [67] # select uBB flux number #67
approach = '1' # select uBB approach
current_ua = [20] #[5.5,0,9,0,10] set the proton current in uA. (again if only one irradiation do array of one)

surface_or_cell = 'cell' #state surface or cell or volume

input_filetag = 'uBB' # or uBB
experiment_or_approx = 'approximation'

#flux_folder_path = '/Users/ljb841@student.bham.ac.uk/fispact/WORKSHOP/fluxes/ubb_design2/approach1/results/results_shielded/{}_flux'.format(projectile)
#flux_folder_path = '/Users/ljb841@student.bham.ac.uk/fispact/WORKSHOP/fluxes/experiment_shielded_2024'
flux_folder_path = '/Users/ljb841@student.bham.ac.uk/fispact/WORKSHOP/fluxes/proton_deuteron_comparison/deuteron'

# set up to use 

library = 'tendl21' # 'tendl21' or 'endfb8' or 'irdff2

flux_benchmark_factor = 0.02756

######################
######################
######################
######################

if library == 'irdff2':
    group = 725
if library == 'endfb8':
    group = 709
if library == 'tendl21':
    group = 1102


def fispact_setup(volume,material):

    # Sets up the FISPACT file to run - changing the flux normalisation, filename and fluxes file
    input_file_path = shutil.copyfile('/Users/ljb841@student.bham.ac.uk/fispact/WORKSHOP/uBB/materials/{}/uBB_{}.i'.format(experiment_or_approx,material), '{}_{}_{}{}.i'.format(input_filetag,material,surface_or_cell,volume))
    files_file_path = '/Users/ljb841@student.bham.ac.uk/fispact/WORKSHOP/uBB/files'
    arb_flux_file_path = '/Users/ljb841@student.bham.ac.uk/fispact/WORKSHOP/uBB/arb_flux_{}'.format(group)

    #if group == 1102:
    #    fluxes_file_path = '{}/{}_{}_FLUXES.dat'.format(flux_folder_path,surface_or_cell,volume)
    #else:
    #    fluxes_file_path = '{}/{}_{}_FLUXES_{}.dat'.format(flux_folder_path,surface_or_cell,volume,group)
    #fluxes_file_path = '{}/fluxes_{}/{}_foil.dat'.format(flux_folder_path,group,material)
    fluxes_file_path = '{}/fluxes_{}/fe_foil.dat'.format(flux_folder_path,group)


    with open(fluxes_file_path, 'r') as filename3:
        fluxes_file = filename3.readlines()

    if library == 'irdff2':
        initial_flux_value = float(fluxes_file[183].split()[2])
    if library == 'endfb8':
        initial_flux_value = float(fluxes_file[120].split()[2])
    if library == 'tendl21':
        initial_flux_value = float(fluxes_file[185].split()[2])
    fluence_factor = flux_benchmark_factor*6.24151e+12 

    final_flux_values = []
    for i in current_ua:
        proton_flux = (i)*fluence_factor
        final_flux_values.append(initial_flux_value*proton_flux)
    print('flux values are:',final_flux_values, 'n/cm2/s')
    if projectile == 'proton':
        with open(arb_flux_file_path,'r') as filename4:
            arb_flux_file = filename4.readlines()
            arb_flux_file[group+1:] = fluxes_file
        with open(arb_flux_file_path,'w') as filename4:
            filename4.writelines(arb_flux_file)

    with open(input_file_path,'r') as filename1:
        input_file = filename1.readlines()
    flux_keyword_line = 25+int(input_file[17-1].split()[2])
    input_file[10-1] = '<< ALLDISPEN 40 >> \n'
    input_file[13-1] = '* Approach {} - {} - vol{} - {} flux - first irrad current {}uA \n'.format(approach,material,volume,projectile,current_ua[0])
    if j == 'li':
        input_file[14-1] = '<< NUCGRAPH 1 0.01 1 1 >> \n'
    input_file[flux_keyword_line-6] = 'DOSE {} \n'.format(dose_input)
    input_file[flux_keyword_line-1] = 'FLUX {} \n'.format(final_flux_values[0])
    input_file[flux_keyword_line+1] = 'TIME {} MINS ATOMS \n'.format(irrad_time[0])
    if experiment_or_approx == 'experiment':   
        input_file[flux_keyword_line+2] = 'FLUX {} \n'.format(final_flux_values[1])
        input_file[flux_keyword_line+4] = 'TIME {} MINS ATOMS \n'.format(irrad_time[1])

        input_file[flux_keyword_line+5] = 'FLUX {} \n'.format(final_flux_values[2])
        input_file[flux_keyword_line+7] = 'TIME {} MINS ATOMS \n'.format(irrad_time[2])

        input_file[flux_keyword_line+8] = 'FLUX {} \n'.format(final_flux_values[3])
        input_file[flux_keyword_line+10] = 'TIME {} MINS ATOMS \n'.format(irrad_time[3])

        input_file[flux_keyword_line+11] = 'FLUX {} \n'.format(final_flux_values[4])
        input_file[flux_keyword_line+13] = 'TIME {} MINS ATOMS \n'.format(irrad_time[4])

    else:
        pass

    if projectile == 'neutron':
        input_file[5-1] = 'PROJ 1 \n'
        input_file[6-1] = '<< GRPCONVERT 1102 162 >>\n'
        input_file[8-1] = 'GETXS 1 {} \n'.format(group)
    if projectile == 'proton':
        input_file[5-1] = 'PROJ 3 \n'
        input_file[6-1] = 'GRPCONVERT 1102 162 \n'
        input_file[8-1] = 'GETXS 1 162 \n' 
    with open(input_file_path, 'w') as filename1:
        filename1.writelines(input_file)

    with open(files_file_path,'r') as filename2:
        files_file = filename2.readlines()
    if projectile == 'neutron':
        files_file[11-1] = 'fluxes {} \n'.format(fluxes_file_path)
        files_file[12-1] = '# arb_flux arb_flux \n'
        files_file[8-1] = 'prob_tab  ../../nuclear_data/tendl21data/tp-1102-294 \n'
        files_file[18-1] = '# enbins ../../nuclear_data/IRDFF-II/ebins_725 \n' 
        files_file[15-1] = 'dk_endf ../../nuclear_data/decay2020/decay_2020 \n'
        files_file[2-1] = 'ind_nuc  ../../nuclear_data/decay2020/decay_2020_index.txt \n'
        files_file[5-1] = 'xs_endf ../../nuclear_data/tendl21data/gendf-1102 \n' 
        if library == 'irdff2':
            files_file[5-1] = 'xs_endf ../../nuclear_data/IRDFF-II/irdff-II_725-n \n'
            files_file[18-1] = 'enbins ../../nuclear_data/IRDFF-II/ebins_725 \n' 
            files_file[15-1] = 'dk_endf ../../nuclear_data/IRDFF-II/decay_irdff-II \n'
            files_file[2-1] = 'ind_nuc  ../../nuclear_data/IRDFF-II/IRDFF_index \n'
        if library == 'endfb8':
            files_file[5-1] = 'xs_endf ../../nuclear_data/ENDFB80data/endfb80-n/gxs-709 \n' 
        if library == 'tend21':
            files_file[5-1] = 'xs_endf ../../nuclear_data/tendl21data/gendf-1102 \n' 


    if projectile == 'proton':
        files_file[5-1] = 'xs_endf ../../nuclear_data/p-tendl2019/gxs-162 \n'
        files_file[8-1] = 'prob_tab  ../../nuclear_data/tendl19data/tp-1102-294 \n'
        files_file[11-1] = 'fluxes {}/{}_{}_fluxes_162.dat \n'.format(flux_folder_path,surface_or_cell,volume)
        files_file[12-1] = 'arb_flux arb_flux \n'
    with open(files_file_path, 'w') as filename2:
        filename2.writelines(files_file)

# Runs FISPACT for all requested volumes    
for i in fluxes:
    for j in materials:
        print('Running FISPACT for {}_{}_{}{}.i'.format(input_filetag,j,surface_or_cell,str(i)))
        fispact_setup(str(i),j)
        os.system('$fispact {}_{}_{}{}.i'.format(input_filetag,j,surface_or_cell,str(i)))