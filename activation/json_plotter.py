import os
import shutil
import numpy as np
import sys
import matplotlib.font_manager
import matplotlib
import matplotlib.pyplot as plt
#plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = "Helvetica"
plt.rcParams["font.size"] = 22
plt.rcParams["font.weight"] = "normal"
import json

######################
######################
# USER INPUTS
######################
#####################

# input directory with the json file in (which is also the directory where the plots will be saved) + input json filename (not extension)
directory = '/Users/ljb841@student.bham.ac.uk/fispact/WORKSHOP/demo_hcll/1hr1g_fullflux_fe'
json_name = 'hcll_fw_fe'

######################
######################
######################
######################

dictionary = json.load(open('{}/{}.json'.format(directory,json_name), 'r'))
inventory_data = dictionary['inventory_data']
#print(inventory_data[0]['nuclides'][0])

timestep_array = [] # in days
total_activity_array = [] # in Bq
total_activity_normalised_array = [] # in Bq/g
total_dose_array = [] # in Sv/hr

for timestep in range(0,len(inventory_data)):
    if inventory_data[timestep]['irradiation_time'] != 0:
        timestep_array.append(inventory_data[timestep]['cooling_time']/(3600*24))
        total_activity_array.append(inventory_data[timestep]['total_activity'])
        total_activity_normalised_array.append(inventory_data[timestep]['total_activity']/(1e3*inventory_data[timestep]['total_mass']))
        total_dose_array.append(inventory_data[timestep]['dose_rate']['dose'])

print('********************','plotting activity in Bq/g and dose in Sv/hr', '********************')
print('note: FISPACT {} dose calculation was used'.format(inventory_data[timestep]['dose_rate']['type']))

#PLOT activity bq/g
fig, ax1 = plt.subplots()
ax1.set_xlabel('Decay time (days)') 
ax1.set_ylabel('Activity (Bq/g)')
ax1.tick_params(axis='y')
ax1.set_xlim(1e-3,1e3)
ax1.set_xscale("log")
#ax1.set_ylim(1e0,1e6)
ax1.set_yscale("log")
ax1.plot(timestep_array, total_activity_normalised_array , 'k-' ,     linewidth=1.5)
ax1.axhline(y=1e5, ls='-', c='green', lw=1.5)
ax1.axvline(x=0.04167, ls='--', c='grey')
ax1.axvline(x=1e0, ls='--', c='grey')
ax1.axvline(x=30, ls='--', c='grey')
ax1.text(1.02, 0.9, 'Approx. background', transform = ax1.transAxes, fontsize=12, c='green')
ax1.text(0.2, 1.02, '1 hour', transform = ax1.transAxes, fontsize=12, c='grey')
ax1.text(0.44, 1.02, '1 day', transform = ax1.transAxes, fontsize=12, c='grey')
ax1.text(0.654, 1.02, '1 month', transform = ax1.transAxes, fontsize=12, c='grey')
ax1.grid(which='major')
#ax1.legend(loc="upper left", bbox_to_anchor=(1, 1), borderaxespad=0, frameon=False, fontsize=11)
fig.set_size_inches((8, 6))
fig.savefig(os.path.join(directory, 'total_activity_{}.png'.format(json_name)), transparent=False, bbox_inches='tight')

#PLOT ISOTOPIC DOSE
fig, ax1 = plt.subplots()
ax1.set_xlabel('Decay time (days)') 
ax1.set_ylabel('Dose (Sv/h)')
ax1.tick_params(axis='y')
ax1.set_xlim(1e-3,1e3)
ax1.set_xscale("log")
#ax1.set_ylim(1e-20,1e-5)
ax1.set_yscale("log")
ax1.plot(timestep_array, total_dose_array,'k-' ,    linewidth=1.5)
ax1.axhline(y=1e-6, ls='-', c='green', lw=1.5)
ax1.axvline(x=0.04167, ls='--', c='grey')
ax1.axvline(x=1e0, ls='--', c='grey')
ax1.axvline(x=30, ls='--', c='grey')
ax1.text(1.02, 0.9, 'Approx. background', transform = ax1.transAxes, fontsize=12, c='green')
ax1.text(0.2, 1.02, '1 hour', transform = ax1.transAxes, fontsize=12, c='grey')
ax1.text(0.44, 1.02, '1 day', transform = ax1.transAxes, fontsize=12, c='grey')
ax1.text(0.654, 1.02, '1 month', transform = ax1.transAxes, fontsize=12, c='grey')
ax1.grid(which='major')
#ax1.legend(loc="upper left", bbox_to_anchor=(1, 1), borderaxespad=0, frameon=False, fontsize=11)
fig.set_size_inches((8, 6))
fig.savefig(os.path.join(directory, 'total_dose_{}.png'.format(json_name)), transparent=False, bbox_inches='tight')
