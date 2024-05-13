from dataclasses import dataclass
import json
from math import pi, sqrt, log, exp
import numpy as np
import actigamma as ag
from scipy.integrate import quad
from datetime import datetime

data_dictionary = {}

################################################################################
########### USER INPUTS ########################################################
################################################################################
################################################################################

# insert measurement information and the count rates and uncertainties for the top 5 high-intensity peaks in lists (in descending intensity order):
# use example format: data_dictionary['Isotope'] = [ live time(s), datetime(year,month,day,hour,minute,second), [count rates], [uncertainties] ]

#EXAMPLE:
data_dictionary['Mn56'] = [703,datetime(2024,3,28, 18,42,32),[12530,1663,744,0,0],[114,44,30,0,0]]
data_dictionary['Au196'] = [418,datetime(2024,3,28, 22,43,46), [107,0,0,0,0], [23,0,0,0,0]]
data_dictionary['Au198'] = [418,datetime(2024,3,28, 22,43,46), [11584,56,0,0,0], [109,11,0,0,0]] 
data_dictionary['Na24'] = [8644,datetime(2024,3,28, 20,0,42),[8514,4319,0,0,0],[98,68,0,0,0]]
data_dictionary['Ni65'] = [929,datetime(2024,3,28, 18,58,33),[83,101,0,0,0],[11,18,0,0,0]]
data_dictionary['Cu64'] = [929,datetime(2024,3,28, 18,58,33),[66],[9]]
                           
# create a folder in working directory to save the results
folder_path = '/Users/ljb841@student.bham.ac.uk/fispact/WORKSHOP/uBB/analysis/activities_test' 

# gamma spec distance (cm). Script will run for all isotopes in your library for this distance. Current options are 1, 6, 10, 15 and 34 cm
# ADVISE DOING OWN CALIBRATION FOR DISTANCE YOU ARE USING. 
measurement_distance = 1

# input total irradiation time (s)
irrad_time = (20+67+41)*60 + 30

# input timestamp for the end of the irradiation in the format datetime(year,month,day,hour,minute,second)
irradiation_end = datetime(2024,3,28, 18,17,32)

###########################################################################################
###########################################################################################
###########################################################################################
##########################################################################################

# function for calculating the decay time from input irradiation end and measurement start
def decay_time(ts):
    return (ts-irradiation_end).total_seconds()

# read the pypact actigamma decay2012 database and get halflife, intensities,peaks for specified isotope
def get_decay_database(isotope_name):
    SPECTYPE = "gamma"
    db = ag.Decay2012Database()
    half_life = db.gethalflife(isotope_name) 
    intensity = db.getintensities(isotope_name,spectype=SPECTYPE)
    peak_energy_kev = (db.getenergies(isotope_name, spectype=SPECTYPE))*1e-3
    sorted_lists = sorted(zip(peak_energy_kev, intensity),key=lambda x: x[1])
    e, i = zip(*reversed(sorted_lists))
    return(i,e,half_life)

# for convenience
@dataclass
class FispactOutput:
    name: str
    activity: float

# equation for the efficiency curves used below
def efficiency_eqn(energy:float,n1:float,n2:float) -> float:
     return n1 * (energy)** (n2)

# use the measurement distance and efficiency curves to calculate activity over the live time
def activity_livetime(c,i,e) :  
    if measurement_distance == 1:
        selected_efficiency = efficiency_eqn(e,8.7655,-0.934) # louis fit with ba133,cs137,co60 (omitting ba133 81keV peak)
    if measurement_distance == 6:
        selected_efficiency = efficiency_eqn(e,0.6264,-0.749) # Kyle fit
    if measurement_distance == 10:
        selected_efficiency = efficiency_eqn(e,0.3755,-0.765) # Kyle NEW fit
    if measurement_distance == 15:
        selected_efficiency = efficiency_eqn(e,0.1966,-0.763) # Kyle fit 
    if measurement_distance == 34:
        selected_efficiency = efficiency_eqn(e,0.074,-0.868) # louis fit with ba133,cs137,co60 (omitting ba133 81keV peak)
    activity = (c) / ((i)
        * selected_efficiency
    )
    return activity

# integral to calculate initial activity at end of irradiation (activity_0) from the measured activity over a live time
def activity_integrand(t,half_life):
    return exp(- log(2) * (t/half_life))
def activity_0(c,i,e) :
    activity = activity_livetime(c,i,e) / (quad(activity_integrand, decay_time(data_dictionary[isotope_name][1]), (decay_time(data_dictionary[isotope_name][1])+data_dictionary[isotope_name][0]), args=(get_decay_database(isotope_name)[2])))
    return activity

# calculate number of atoms generated over the irradiation (using the a0 calculated from the highest-intensity peak)
def no_of_atoms(a):
    no_atoms = activity_0(a,get_decay_database(isotope_name)[0][0],get_decay_database(isotope_name)[1][0])[0] / ( log(2) / get_decay_database(isotope_name)[2] )
    return no_atoms

# calculate reaction rate over the entire irradiation (using the end-of-irradiation activity a0 calculated from the highest intensity peak)
def reaction_rates(a, irrad_time):
    rr_ave = activity_0(a,get_decay_database(isotope_name)[0][0],get_decay_database(isotope_name)[1][0])[0] / (1 - activity_integrand(irrad_time,get_decay_database(isotope_name)[2]))
    return rr_ave


# print and save results 
for isotope_name in list(data_dictionary.keys())[:]:
    print('************ Calculating activities and reaction rates for', isotope_name, '************')
    with open('{}/{}_activities.txt'. format(folder_path,isotope_name), 'w') as output_file:
        for n in range(len(get_decay_database(isotope_name)[0][:5])):
            if data_dictionary[isotope_name][2][n] != 0:
                print('(e=',get_decay_database(isotope_name)[1][n],'keV, i=',get_decay_database(isotope_name)[0][n],') activity at end of irradiation is', activity_0(data_dictionary[isotope_name][2][n],get_decay_database(isotope_name)[0][n],get_decay_database(isotope_name)[1][n])[0],'+-',activity_0(data_dictionary[isotope_name][3][n],get_decay_database(isotope_name)[0][n],get_decay_database(isotope_name)[1][n])[0], 'Bq')
                output_file.writelines(f" e={get_decay_database(isotope_name)[1][n]} keV, i={get_decay_database(isotope_name)[0][n]}: activity at end of irradiation is {activity_0(data_dictionary[isotope_name][2][n],get_decay_database(isotope_name)[0][n],get_decay_database(isotope_name)[1][n])[0]} +- {activity_0(data_dictionary[isotope_name][3][n],get_decay_database(isotope_name)[0][n],get_decay_database(isotope_name)[1][n])[0]} Bq \n")    
        
    # output average reaction rates for each isotope in a text file for unfolding purposes
    print('Average reaction rate over irradiation:', '{:.3e}'.format(reaction_rates(data_dictionary[isotope_name][2][0],irrad_time)), '+-', '{:.3e}'.format(reaction_rates(data_dictionary[isotope_name][3][0],irrad_time)))
    with open('{}/reaction_rates.txt'. format(folder_path), 'a') as output_file:        
        output_file.writelines(f"{reaction_rates(data_dictionary[isotope_name][2][0],irrad_time)} \n")
    with open('{}/reaction_rate_uncertainties.txt'. format(folder_path), 'a') as output_file:    
        output_file.writelines(f"{reaction_rates(data_dictionary[isotope_name][3][0],irrad_time)} \n") 