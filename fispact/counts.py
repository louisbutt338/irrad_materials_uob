from dataclasses import dataclass
import json
from math import pi, sqrt, log, exp
import numpy as np


####################################
############ USER INPUTS ###########
####################################

timestep = 4 # 7 is approx 1.1 days, 4 is approx 1.1 hours, 10 is approx 46 days, 8 a week

flux_number = 67
materials = ['fe','au','al','cu','in','nb','ni','rh','sc','y','dy','cd']

total_counts_minimum = 1e4
no_of_timesteps = 14
surface_or_cell = 'cell'

folder_path = '/Users/ljb841@student.bham.ac.uk/fispact/WORKSHOP/uBB/model_results/approach_1_d/neutrons/foils_hpge_scoping_121124'
results_folder_name = 'timestep{}_analysis'.format(timestep)

####################################
######### TONY FUNCTIONS ###########
####################################

# for convenience
@dataclass
class FispactOutput:
    name: str
    activity: float


# load a json file to dict
def load_json(path: str):
    """Load pre-processed decay data JSON library"""
    with open(path, "r") as f:
        json_dict = json.load(f)
    return json_dict


# Convert nuclide names into a common format
def common_name(name: str) -> str:
    """Convert nuclide names into a common format

    Designed to accommodate various styles of isotope definitions for
    convenience. For example Cobalt 60 can be 'Co-60', 'co60', and 'co -60-'
    will all become 'co60'.
    """
    # remove all spaces, get rid of dashes etc...
    common = name.lower()
    common = common.replace("-", "")
    common = common.replace("_", "")
    common = common.replace(" ", "")
    return common


# check if a nuclide even has photon emissions
def is_nuclide_relevant(library: dict, nuclide: str) -> bool:
    if nuclide in library.keys():
        radtypes = list(library[nuclide].keys())
        # just get any photon emissions
        if any(x in ["gamma", "x"] for x in radtypes):
            return True
    return False


# collect all energy/intensity values from a
def collect_energy_intensity(nuclide_data: dict) -> tuple:

    energy = []
    intensity = []

    for radtype, data in nuclide_data.items():
        # skip over anything not relevant like alpha or beta
        if not radtype in ["gamma", "x"]:
            continue

        for e, i, n in zip(data["energy"], data["intensity"], data["norm"]):
            energy.append(e)
            intensity.append(i * n)

    # Sort both lists based on the energy
    sorted_lists = sorted(zip(energy, intensity))
    e, i = zip(*sorted_lists)
    return e, i


# solid angle approximation (knoll, p120)
def solid_angle(crystal_radius: float, distance: float) -> float:
    return 2 * pi * (1 - distance / sqrt(distance**2 + crystal_radius**2))


# taken from calibration curve of standard germanium at ukaea
def knee_energy() -> float:
    return 233.00


# taken from calibration curve of standard germanium at ukaea
def above_knee(energy: float) -> float:
    ln_efficiency = -4.4729 - 0.152145 * log(energy) - 0.0408071 * (log(energy)) ** 2
    return exp(ln_efficiency)


# taken from calibration curve of standard germanium at ukaea
def below_knee(energy: float) -> float:
    ln_efficiency = -19.3723 + 5.584864 * log(energy) - 0.591785 * (log(energy)) ** 2
    return exp(ln_efficiency)


# rough germanium efficiency for a given peak energy
def efficiency(energy: float) -> float:
    if energy > knee_energy():
        return max(above_knee(energy), 0.0)
    else:
        return max(below_knee(energy), 0.0)


def relevant_nuclides_only(nuclides: list, activity: list, library: dict) -> list:
    # sort out the name to match this horrible method
    names = [common_name(n) for n in fispact_nuclides]

    return [
        FispactOutput(name, act)
        for name, act in zip(names, fispact_activity)
        if is_nuclide_relevant(library, name)
    ]

####################################
###### LOUIS FUNCTIONS #############
####################################
####################################

def grn_isotope_finder(foil_grn_data, isotopes_line):
    isotopes_array = []
    for x in foil_grn_data[isotopes_line-1:isotopes_line]:
        isotopes_array.append(x.split()[5::1])
        isotopes_array2=isotopes_array[0]
        isotopes_array3=[]
        for y in np.arange(0,len(isotopes_array2),1):
            if len(isotopes_array2[y]) > 4:
                isotopes_array3.append(isotopes_array2[y])
            if isotopes_array2[y].isdigit() == True:
                continue
            if isotopes_array2[y].isalpha() == True:
                isotopes = ''.join(isotopes_array2[y:y+2])
                isotopes_array3.append(isotopes)
            else:
                continue 
    return isotopes_array3

def grn_data_finder(foil_grn_data, isotopes_line):
    activities_array = []
    for i in foil_grn_data[isotopes_line+timestep-1:isotopes_line+timestep]:
        activities_array.append(i.split()[4::1])
        activities_array2 = np.array(activities_array[0])
        activities_array3 = [j*1e-3 for j in activities_array2.astype(float)]
        dose_array = [j for j in activities_array2.astype(float)]
    return activities_array3,dose_array     

def collect_fispact_info(foil_grn_filepath):
    with open(foil_grn_filepath,'r') as foil_grn_file:
        foil_grn_data = foil_grn_file.readlines()
        act_isotopes_array = grn_isotope_finder(foil_grn_data,8)
        dose_isotopes_array = grn_isotope_finder(foil_grn_data,2*no_of_timesteps+22)
        activities_array = grn_data_finder(foil_grn_data,8)[0]
        dose_array = grn_data_finder(foil_grn_data,2*no_of_timesteps+22)[1]

    return(act_isotopes_array,activities_array,dose_isotopes_array,dose_array)


#################################### 
for m in materials:
    counts_results_filepath = '{}/{}/{}_counts_ts{}_counts.txt'.format(folder_path,results_folder_name,m,str(timestep))
    dose_results_filepath =   '{}/{}/{}_counts_ts{}_doses.txt'. format(folder_path,results_folder_name,m,str(timestep))
    grn_filepath = '{}/uBB_{}_{}{}.grn'.format(folder_path,m,surface_or_cell,str(flux_number))

     # gimme the basic info from fispact
    fispact_nuclides = collect_fispact_info(grn_filepath)[0]
    #fispact_nuclides = ["h3", "au196m", "au196n", "au198", "au196","al28","na24","mg27"]
    fispact_activity = collect_fispact_info(grn_filepath)[1]

    fispact_dose_nuclides = collect_fispact_info(grn_filepath)[2]
    fispact_doses = collect_fispact_info(grn_filepath)[3]

      # load library of gamma lines takes from printlib5 output
    library = load_json("./decay2020.json")

     # just collect any that are relevant by checking for photon emission data
    nuclides = relevant_nuclides_only(fispact_nuclides, fispact_activity, library)

     # for convenience
    @dataclass
    class TableRow:
        nuclide_name: str
        energy: float
        intensity: float
        activity: float
        count_rate: float

     # grab every peak
    table_rows=[]
    total_count_rates = []
    for n in nuclides:
        print(f"Processing nuclide: {n.name}")
        energy, intensity = collect_energy_intensity(library[n.name])

        for e, i in zip(energy, intensity):
            counts = (
                (n.activity * i)
                * efficiency(e * 1e-3)
                * (solid_angle(crystal_radius=3.5, distance=1.0) / (4 * pi)) # change detector size and distance here
            )
            table_rows.append(
                TableRow(
                    nuclide_name=n.name,
                    energy=e,
                    intensity=i,
                    activity=n.activity,
                    count_rate=counts,
                )
            )
            
    # sort by the number of expected counts in the spectrum
    table_rows = sorted(table_rows, key=lambda row: row.count_rate, reverse=True)
    total_count_rate = sum([r.count_rate for r in table_rows])
  
    fispact_nuclides_sorted = [fispact_nuclides for _, fispact_nuclides in sorted(zip(fispact_activity,fispact_nuclides),reverse=True)]
    fispact_activity_sorted = sorted(fispact_activity, reverse=True)
    fispact_dose_nuclides_sorted = [fispact_dose_nuclides for _, fispact_dose_nuclides in sorted(zip(fispact_doses,fispact_dose_nuclides),reverse=True)]
    fispact_doses_sorted = sorted(fispact_doses, reverse=True)
    total_activity_sum = sum(fispact_activity)
    total_dose_sum = sum(fispact_doses)

    print('\nCount rate info:')
    for r in table_rows[:10]:
        per_count = 100.0 * (r.count_rate / total_count_rate)
        measurement_times = total_counts_minimum*((r.count_rate)**-1)
        print(
            f"{r.nuclide_name:<7} {r.energy*1e-3:<10.3f} keV   {r.intensity:<10.8f}   {r.count_rate:.3e} counts/s  ({per_count:.2f}%)"
        )
        with open(counts_results_filepath, 'a') as output_file:
            output_file.writelines(f"{r.nuclide_name:<7} {r.energy*1e-3:<10.3f} keV   {r.intensity:<10.8f}   {r.count_rate:.3e} counts/s  ({per_count:.2f}%)   time to 1e4 counts = {measurement_times:.3e} s \n")    
    
    print('\nActivity info:')
    with open(dose_results_filepath, 'a') as output_file:
        output_file.writelines("Activity info: \n")
    for r in np.arange(0,min(len(fispact_nuclides),10),1):
        per_activity = 100.0 * (fispact_activity_sorted[r] / total_activity_sum)
        print(
            f"{fispact_nuclides_sorted[r]:<7}   {fispact_activity_sorted[r]:.3e} Bq/g  ({per_activity:.2f}%)"
        )
        with open(dose_results_filepath, 'a') as output_file:
            output_file.writelines(f"{fispact_nuclides_sorted[r]:<7}   {fispact_activity_sorted[r]:.3e} Bq/g  ({per_activity:.2f}%) \n")
    print('\nDose info:')   
    with open(dose_results_filepath, 'w') as output_file:
        output_file.writelines("Dose info: \n")           
    for r in np.arange(0,min(len(fispact_dose_nuclides),10),1):
        per_dose = 100.0 * (fispact_doses_sorted[r] / total_dose_sum)
        print(
            f"{fispact_nuclides_sorted[r]:<7}   {fispact_doses_sorted[r]:.3e} Sv/hr  ({per_dose:.2f}%)"
        )
        with open(dose_results_filepath, 'a') as output_file:
            output_file.writelines(f"{fispact_nuclides_sorted[r]:<7}   {fispact_doses_sorted[r]:.3e} Sv/hr  ({per_dose:.2f}%) \n")    
            
    
     #        print(
    #            f"{e*1e-3:<10.3f} keV    {i:<10.8f}   {n.activity:.3e} Bq    {counts:.3e} counts/s"
    #        )

