import matplotlib.pyplot as plt
import numpy as np
import json
import sys
import os

library = 'tendl21'

fispact_results_folder = '/Users/ljb841@student.bham.ac.uk/fispact/WORKSHOP/uBB/040924_foils_fe_flux_analysis/{}'.format(library)
collected_results_folder = '/Users/ljb841@student.bham.ac.uk/fispact/WORKSHOP/analysis/calculated_activities'

#materials = ['au', 'al', 'fe', 'in', 'nb', 'ni', 'rh', 'sc', 'y', 'dy', 'cd','cu']
materials = 'fe'

for m in materials:
    calculated_result_json = json.load(open('{}/uBB_{}_cell12.json'.format(fispact_results_folder,m)))
    activity = calculated_result_json['neutron_cell_flux'][3]['value']
    activity_uncert = calculated_result_json['neutron_cell_flux'][3]['value']

    print(activity,activity_uncert)

