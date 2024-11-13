
import numpy as np
from f4enix.output.mctal import Mctal
import matplotlib.pyplot as plt
from pathlib import Path

##################################
# user inputs
##################################

group_structure = 1102
tally_number = 154
cells_tallied = 1
current_ua = 20 # current in uA
cell_number = [67] # enter  cell number (starting with 1) that you would like to plot. if wish to plot all cells enter 0

folder = "proton_deuteron_comparison/deuteron"

filename="run0.m"

##################################
##################################

cells_array = [58,81,82,67,68,69,70,71,72,73,74,75,76,77,78]
def which_cells(cell_number_array):
    for cell1 in cell_number_array:
        for count,cell2 in enumerate(cells_array):
            if cell1 == cell2:
                return count+1
            else:
                continue


def flux_calculator():
    proton_flux_scaling = 0.02756*6.24151e12 # protons/s for 1uA
    if which_cells(cell_number) == None: 
        print('Cell number is not listed')
    else:
        n = which_cells(cell_number) - 1 
        flux = np.array(value[(n*group_structure)+n:((n+1)*group_structure)+n])
        normalised_flux = flux*proton_flux_scaling*current_ua
    return flux,normalised_flux,n

def get_tally(mctal: Mctal, tally: int) -> tuple:
    # grab the relevant data from the pandas dataframe and make a few lists
    energy = mctal.tallydata.get(tally).get("Energy").tolist()
    value = mctal.tallydata.get(tally).get("Value").tolist()
    relative_error = mctal.tallydata.get(tally).get("Error").tolist()
 
    return energy, value, relative_error
 
 
mctal = Mctal('{}/{}'.format(folder,filename))
energy, value, rel_err = get_tally(mctal, tally_number)
 
#print(energy[:group_structure+2])
#print(value[:group_structure+2]) 
#print(energy[1*group_structure+3:2*group_structure+5]) # n =1 
#value[n*(group_structure+3):(n+1)*(group_structure)+2+(n*3)]
#print(rel_err)


def dump_fluxes_file(output: Path):

    total_flux = sum(flux_calculator()[0])

    with open(output, "w") as f:
        count = 1
        # legacy fispact likes reverse order
        for flx in reversed(flux_calculator()[0]):
            f.write(f"{flx:.5e} ")
            count += 1
            if count > 6:
                f.write("\n")
                count = 1

        f.write(f"\n1.000")
        f.write(f"\nTotal = {total_flux:.5e} [n/cm2/src proton]")

def plot_flux(path,cell,groups):

    if groups > 500:
        energy_binning = energy[:group_structure]
    else:
        energy_binning = energy[:group_structure+2]

    fig, ax1 = plt.subplots(figsize=(8,6))
    ax1.step(energy_binning, flux_calculator()[1], where='post', label='Cell {}'.format(flux_calculator()[2]+1))
    #ax1.set_xscale('log')
    ax1.set_yscale('log')
    #ax1.set_ylim(1e6,1e9)
    ax1.set_xlim(0,40)
    ax1.grid()
    ax1.set_ylabel('Flux (n cm$^{-2}$ s$^{-1}$)')
    ax1.set_xlabel('Neutron energy (MeV)')
    fig.tight_layout()
    #fig.legend(loc="upper left", bbox_to_anchor=(0.8, 0.8), borderaxespad=0, frameon=True, fontsize=12,fancybox=False,facecolor='white',framealpha=1)
    plt.subplots_adjust(wspace=0.04, hspace=0.1)
    #plt.show()
    plt.savefig(f"{path}/cell{cell}_{groups}_plot.png")



for i, groups in enumerate([group_structure]):
    path_fluxes = Path(f"{folder}/fluxes_{groups}/")
    path_fluxes.mkdir(exist_ok=True)

    for cell in cell_number:
        flux_filename = f"cell{cell}_{groups}.dat"
        output = path_fluxes.joinpath(flux_filename)
        dump_fluxes_file(output)

        plot_flux(path_fluxes,cell,groups)