from dataclasses import dataclass
import json
from math import pi, sqrt, log, exp
import numpy as np
import actigamma as ag
from scipy.integrate import quad
from datetime import datetime

################################################################################
########### USER INPUTS ########################################################
################################################################################
################################################################################

# choose whether you want reaction rates or cross-sections along with the activities
reaction_rate_calculator = True
cross_section_calculator = False

# choose whether to run all FOILS isotopes ('foils'), TARGET isotopes ('target') 
# or a specific isotope ('isotope'):
automation = 'foils'

# choose peak analysis library 
peak_library = 'maestro'

# folder to save the results
folder_path = f"/Users/ljb841@student.bham.ac.uk/fispact/WORKSHOP/uBB/analysis/experimental_activities/{peak_library}" 

# total irradiation time secs
irrad_time = (20+67+41)*60 + 30

# datetime timestamp for the end of your irradiation
irradiation_end = datetime(2024,3,28, 18,17,32)

# insert live time, datetime of gamma spec measurement:
data_dictionary = {}
data_dictionary['Be7'   ] = [9091.4, datetime(2024,4,15, 12,39,1)]
data_dictionary['Zn65'  ] = [9091.4, datetime(2024,4,15, 12,39,1)]

data_dictionary['Mn56'  ] = [703,datetime(2024,3,28, 18,42,32)]
data_dictionary['Au196' ] = [418,datetime(2024,3,28, 22,43,46)]
data_dictionary['Au198' ] = [418,datetime(2024,3,28, 22,43,46)] 
data_dictionary['Na24'  ] = [8644,datetime(2024,3,28, 20,0,42)]
data_dictionary['Ni65'  ] = [929,datetime(2024,3,28, 18,58,33)]
data_dictionary['Cu64'  ] = [929,datetime(2024,3,28, 18,58,33)]
data_dictionary['Cd111m'] = [1536,datetime(2024,3,28, 19,24,23)]
data_dictionary['In117' ] = [1536,datetime(2024,3,28, 19,24,23)] # do not use - isotope clash
data_dictionary['Cd117' ] = [1536,datetime(2024,3,28, 19,24,23)] # do not use - isotope clash
data_dictionary['Cd115' ] = [1536,datetime(2024,3,28, 19,24,23)]
data_dictionary['In115m'] = [459,datetime(2024,3,28, 19,51,24)]
data_dictionary['In116m'] = [459,datetime(2024,3,28, 19,51,24)]
data_dictionary['Ni57'  ] = [274848.9,datetime(2024,3,29, 11,29,26)]
data_dictionary['Co58'  ] = [274848.9,datetime(2024,3,29, 11,29,26)] # do not use - isotope clash
data_dictionary['Dy165' ] = [44220,datetime(2024,3,28, 23,8,18)]
data_dictionary['Dy157' ] = [44220,datetime(2024,3,28, 23,8,18)]
data_dictionary['Nb92m' ] = [24133,datetime(2024,4,3, 10,16,14)]
#    data_dictionary['Sc44m'] = [56363,datetime(2024,4,4, 6,33,49)]
#    data_dictionary['Rh102m'] = [185374,datetime(2024,4,5, 13,50,32)]
#    data_dictionary['Y88'] = [189824.6,datetime(2024,4,12, 11,57,12)]

# insert count rates and uncertainties for the top 5 high-intensity peaks in lists (in descending intensity order):

# MAESTRO dictionary
if peak_library == 'maestro':
    data_dictionary.setdefault('Be7'   ,[]).extend([[52479],[302]])
    data_dictionary.setdefault('Zn65'  ,[]).extend([[353229],[615]])
    data_dictionary.setdefault('Mn56'  ,[]).extend([[12530,1663,744,0,0],[114,44,30,0,0]])
    data_dictionary.setdefault('Au196' ,[]).extend([[107,0,0,0,0],[23,0,0,0,0]])
    data_dictionary.setdefault('Au198' ,[]).extend([[11584,56,0,0,0],[109,11,0,0,0]]) 
    data_dictionary.setdefault('Na24'  ,[]).extend([[8514,4319,0,0,0],[98,68,0,0,0]])
    data_dictionary.setdefault('Ni65'  ,[]).extend([[83,101,0,0,0],[11,18,0,0,0]])
    data_dictionary.setdefault('Cu64'  ,[]).extend([[66],[9]])
    data_dictionary.setdefault('Cd111m',[]).extend([[13537,5794],[124,94]])
    data_dictionary.setdefault('In117' ,[]).extend([[280,814,0],[25,59,0]]) 
    data_dictionary.setdefault('Cd117' ,[]).extend([[219,0,0,0,31],[40,0,0,0,22]])
    data_dictionary.setdefault('Cd115' ,[]).extend([[228,48,0,0,0],[26,21,0,0,0]])
    data_dictionary.setdefault('In115m',[]).extend([[14256,0,0,0,0],[243,0,0,0,0]])
    data_dictionary.setdefault('In116m',[]).extend([[50630,37609,34604,7151,9558],[372,351,333,141,233]])
    data_dictionary.setdefault('Ni57'  ,[]).extend([[7959,0,897,0,0],[251,0,84,0,0]])
    data_dictionary.setdefault('Co58'  ,[]).extend([[212850,0,655],[512,0,85]]) 
    data_dictionary.setdefault('Dy165' ,[]).extend([[19022,2148,913,667,1430],[324,136,91,84,179]])
    data_dictionary.setdefault('Dy157' ,[]).extend([[2463,0,0,0,0],[154,0,0,0,0]])
    data_dictionary.setdefault('Nb92m' ,[]).extend([[3935,0,0,0,0],[82,0,0,0,0]])

# interspec library
if peak_library == 'interspec':
    data_dictionary.setdefault('Be7'   ,[]).extend([[9091.4*5.824],[9091.4*0.027]])
    data_dictionary.setdefault('Zn65'  ,[]).extend([[9091.4*48.87],[9091.4*0.13]])
    data_dictionary.setdefault('Mn56'  ,[]).extend([[12549,1641,739,0,0],[703*0.16,703*0.059,703*0.039,0,0]])
    data_dictionary.setdefault('Au196' ,[]).extend([[0.2961*418,0,0,0,0], [0.0324*418,0,0,0,0]])
    data_dictionary.setdefault('Au198' ,[]).extend([[12777,47,0,0,0], [418*0.44,418*0.0181,0,0,0]]) 
    data_dictionary.setdefault('Na24'  ,[]).extend([[8644*0.9896,8644*0.492,0,0,0],[8644*0.0109,8644*0.0076,0,0,0]])
    data_dictionary.setdefault('Ni65'  ,[]).extend([[929*0.08525,929*0.1381,0,0,0],[929*0.01006,929*0.0131,0,0,0]])
    data_dictionary.setdefault('Cu64'  ,[]).extend([[929*0.04651],[929*0.00753]])
    data_dictionary.setdefault('Cd111m',[]).extend([[1536*8.821,1536*3.78],[1536*0.077,1536*0.054]])
    data_dictionary.setdefault('In117' ,[]).extend([[297,820,0],[25,59,0]]) 
    data_dictionary.setdefault('Cd117' ,[]).extend([[219,0,0,0,31],[40,0,0,0,22]]) 
    data_dictionary.setdefault('Cd115' ,[]).extend([[1536*0.163,1536*0.04644,0,0,0],[1536*0.0116,1536*0.00778,0,0,0]])
    data_dictionary.setdefault('In115m',[]).extend([[459*53.49,0,0,0,0],[459*0.48,0,0,0,0]])
    data_dictionary.setdefault('In116m',[]).extend([[602.4*459,497*459,416.8*459,71.46*459,71.7*459],[5.2*459,5.9*459,5.9*459,2.12*459,4.41*459]])
    data_dictionary.setdefault('Ni57'  ,[]).extend([[2.748e5*0.027,0,2.748e5*0.0033,0,0],[2.748e5*0.00038,0,2.748e5*0.000164,0,0]])
    data_dictionary.setdefault('Co58'  ,[]).extend([[2.748e5*0.7757,0,2.748e5*0.002439],[2.748e5*0.0017,0,2.748e5*0.000174]]) 
    data_dictionary.setdefault('Dy165' ,[]).extend([[44220*0.4192,44220*0.04384,44220*0.01968,44220*0.01465,44220*0.03389],[44220*0.0045,44220*0.000169,442200*0.0011,44220*0.00099,44220*0.00227]])
    data_dictionary.setdefault('Dy157' ,[]).extend([[44220*0.0533,0,0,0,0],[44220*0.00199,0,0,0,0]])
    data_dictionary.setdefault('Nb92m' ,[]).extend([[24133*0.1614,0,0,0,0],[24133*0.0028,0,0,0,0]])

# ROOT library
if peak_library == 'root':
    data_dictionary.setdefault('Be7'   ,[]).extend([[52870],[253.3]])
    data_dictionary.setdefault('Zn65'  ,[]).extend([[9091.4*48.87],[9091.4*0.13]])
    data_dictionary.setdefault('Mn56'  ,[]).extend([[12520,1641,0,0,0],[113.2,41.6,0,0,0]]) 
    data_dictionary.setdefault('Au196' ,[]).extend([[125.8,0,0,0,0], [14.8,0,0,0,0]])
    data_dictionary.setdefault('Au198' ,[]).extend([[11590,0,0,0,0], [107.4,0,0,0,0]]) 
    data_dictionary.setdefault('Na24'  ,[]).extend([[8539.92,4245.37,0,0,0],[93.807,65.6465,0,0,0]])
    data_dictionary.setdefault('Ni65'  ,[]).extend([[75,118.6,0,0,0],[10.8,12.3,0,0,0]])
    data_dictionary.setdefault('Cu64'  ,[]).extend([[41.92],[8.33]])
    data_dictionary.setdefault('Cd111m',[]).extend([[13510,0],[119,0]])
    data_dictionary.setdefault('Cd115' ,[]).extend([[253.8,0,0,0,0],[18.6,0,0,0,0]])
    data_dictionary.setdefault('Dy165' ,[]).extend([[0,1811,0,0,0],[0,77.4,0,0,0]])
    data_dictionary.setdefault('Dy157' ,[]).extend([[2153,0,0,0,0],[90.4,0,0,0,0]])
    data_dictionary.setdefault('Nb92m' ,[]).extend([[3861,0,0,0,0],[67.5,0,0,0,0]])


###########################################################################################
###########################################################################################
###########################################################################################
##########################################################################################

# setting library isotopes to analyse and the distance they were analysed at
if automation == 'target':
    isotope_run_list = ['Be7','Zn65']
if automation == 'foils':
    isotope_run_list = list(data_dictionary.keys())[2:]
else:
    isotope_run_list = list(automation.split(" "))

if automation in ('Be7','Zn65','target'):
    measurement_distance = 34
else:
    measurement_distance = 1

# function for calculating the decay time from input irradiation end and measurement start
def decay_time(ts):
    return (ts-irradiation_end).total_seconds()

# read the pypact actigamma decay2012 database and get halflife, intensities,peaks for specified isotope
def get_decay_database(isotope_name):
    SPECTYPE = "gamma"
    db = ag.Decay2012Database()
    if isotope_name == 'Co58':
        half_life = db.gethalflife('Co58m') 
    else:
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

# solid angle approximation (knoll, p120)
def solid_angle(crystal_radius: float, distance: float) -> float:
    return 2 * pi * (1 - distance / sqrt(distance**2 + crystal_radius**2))

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

# integral to calculate initial activity from the measured activity over a live time
def activity_integrand(t,half_life):
    return exp(- log(2) * (t/half_life))
def activity_0(c,i,e) :
    activity = activity_livetime(c,i,e) / (quad(activity_integrand, decay_time(data_dictionary[isotope_name][1]), (decay_time(data_dictionary[isotope_name][1])+data_dictionary[isotope_name][0]), args=(get_decay_database(isotope_name)[2])))
    #activity = (activity_livetime(c,i,e)/data_dictionary[isotope_name][0]) * exp(get_decay_database(isotope_name)[2] / (log(2) * decay_time(data_dictionary[isotope_name][1])))
    return activity

# calculate reaction rate from the activity a0 at time t0 under irradiation 
def reaction_rates(a, irrad_time):
    rr_ave = activity_0(a,get_decay_database(isotope_name)[0][0],get_decay_database(isotope_name)[1][0])[0] / (1 - activity_integrand(irrad_time,get_decay_database(isotope_name)[2]))
    #rr_ave = activity_0(a,get_decay_database(isotope_name)[0][0],get_decay_database(isotope_name)[1][0])[0] / irrad_time
    rr_max = activity_0(a,get_decay_database(isotope_name)[0][0],get_decay_database(isotope_name)[1][0])[0] * (10/(8.778*128.5))
    return rr_ave,rr_max


# print and save results for individual isotope activities
for isotope_name in isotope_run_list:
    print(f"************ activities for {isotope_name} ************")
    with open(f"{folder_path}/{isotope_name}_activities.txt", 'w') as output_file:
        for n in range(len(get_decay_database(isotope_name)[0][:5])):
            if data_dictionary[isotope_name][2][n] != 0:
                print(f"(e={get_decay_database(isotope_name)[1][n]}keV, i={get_decay_database(isotope_name)[0][n]}) activity at end of irradiation is {activity_0(data_dictionary[isotope_name][2][n],get_decay_database(isotope_name)[0][n],get_decay_database(isotope_name)[1][n])[0]} +- {activity_0(data_dictionary[isotope_name][3][n],get_decay_database(isotope_name)[0][n],get_decay_database(isotope_name)[1][n])[0]} Bq")
                output_file.writelines(f" e={get_decay_database(isotope_name)[1][n]} keV, i={get_decay_database(isotope_name)[0][n]}: activity at end of irradiation is {activity_0(data_dictionary[isotope_name][2][n],get_decay_database(isotope_name)[0][n],get_decay_database(isotope_name)[1][n])[0]} +- {activity_0(data_dictionary[isotope_name][3][n],get_decay_database(isotope_name)[0][n],get_decay_database(isotope_name)[1][n])[0]} Bq \n")    

# print and save activities and uncertainties for all analysed isotopes as one nice txt
    with open(f"{folder_path}/exp_activities.txt", 'a') as output_file:
        for n in range(len(get_decay_database(isotope_name)[0][:1])):
            if data_dictionary[isotope_name][2][n] != 0:
                output_file.writelines(f"{activity_0(data_dictionary[isotope_name][2][n],get_decay_database(isotope_name)[0][n],get_decay_database(isotope_name)[1][n])[0]}\n")
    with open(f"{folder_path}/exp_uncertainties.txt", 'a') as output_file:
        for n in range(len(get_decay_database(isotope_name)[0][:1])):
            if data_dictionary[isotope_name][2][n] != 0:
                output_file.writelines(f"{activity_0(data_dictionary[isotope_name][3][n],get_decay_database(isotope_name)[0][n],get_decay_database(isotope_name)[1][n])[0]}\n")

# print and save average reaction rates for each isotope in one nice txt file for unfolding
    if reaction_rate_calculator == True:
        print(f"-------- reaction rates for {isotope_name} ----------")
        pathway_prob = [1]
        if isotope_name == 'Cu64':
            pathway_prob = [0.342,0.658]
        if isotope_name == 'Cd111m':
            pathway_prob = [0.02776,0.80484,0.16737]
        if isotope_name == 'Cd115':
            pathway_prob = [0.72795,0.27205]
        if isotope_name == 'Dy157':
            pathway_prob = [0.73696,0.26304]
        for p in pathway_prob:
            print(f"Average (fraction={p}) reaction rate over irradiation:{p*reaction_rates(data_dictionary[isotope_name][2][0],irrad_time)[0]:.3e} +- {p*reaction_rates(data_dictionary[isotope_name][3][0],irrad_time)[0]:.3e}")
            if automation == 'foils':
                with open(f"{folder_path}/reaction_rates.txt", 'a') as output_file:        
                    output_file.writelines(f"{p*reaction_rates(data_dictionary[isotope_name][2][0],irrad_time)[0]} \n")
                with open(f"{folder_path}/reaction_rate_uncertainties.txt", 'a') as output_file:    
                    output_file.writelines(f"{p*reaction_rates(data_dictionary[isotope_name][3][0],irrad_time)[0]} \n") 
    
# output cross-sections 
#    if cross_section_calculator == True:
