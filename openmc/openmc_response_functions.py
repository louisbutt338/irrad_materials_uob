import matplotlib.backend_tools
import openmc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rc("font", family="sans-serif",weight='normal')
import numpy as np
import pandas as pd
import os
import os.path
from matplotlib.pyplot import cm
import actigamma as ag

# choose data library: endfb8 or tendl21 or irdff2
xs_library = 'tendl21'

####################

xs_folder_path = '/Users/ljb841@student.bham.ac.uk/MCNP/MCNP_DATA'
energy_group = openmc.mgxs.EnergyGroups('VITAMIN-J-175')


response_matrix_file = 'response_matrix_{}.csv'.format(xs_library)
if os.path.isfile(response_matrix_file) == True:
    os.remove(response_matrix_file)

def irdff2_xs_extraction(irdff_ace_filepath,mt_number,energy_bins):
    ace_table = openmc.data.ace.get_table(irdff_ace_filepath)
    nxs = ace_table.nxs
    jxs = ace_table.jxs
    xss = ace_table.xss
    lmt = jxs[3]
    nmt = nxs[4]
    lxs = jxs[6]
    mts = xss[lmt : lmt+nmt].astype(int)
    print(mts)
    locators = xss[lxs : lxs+nmt].astype(int)

    cross_sections = {}
    for mt, loca in zip(mts, locators):
    # Determine starting index on energy grid
        nr = int(xss[jxs[7] + loca - 1])
        if nr == 0:
            breakpoints = None
            interpolation = None
        else:
            breakpoints = xss[jxs[7] + loca : jxs[7] + loca + nr].astype(int)
            interpolation = xss[jxs[7] + loca + nr : jxs[7] + loca + 2*nr].astype(int)

        # Determine number of energies in reaction
        ne = int(xss[jxs[7] + loca + 2*nr])

        # Read reaction cross section
        start = jxs[7] + loca + 1 + 2*nr
        energy = xss[start : start + ne] * 1e6
        xs = xss[start + ne : start + 2*ne]
        cross_sections[mt] = openmc.data.Tabulated1D(energy, xs, breakpoints, interpolation)

    return cross_sections[mt_number](energy_bins)

def tendl_extraction(isotope):
    isotope_tendl_format = isotope
    if len(isotope) == 4:
        isotope_tendl_format = isotope[:2] + '0' + isotope[2:]
    return isotope_tendl_format

def reaction_info(isotope, foil, ace_filename, mt_number, density, mass):

    #energy_bins = isotope_data.energy['294K']
    energy_bins = energy_group.group_edges
    material = openmc.Material()
    material.set_density('g/cm3', density)
    material.add_element(foil, 1) 
    isotope_atom_density = material.get_nuclide_atom_densities()[isotope]
    foil_volume = mass / material.density

    if xs_library == 'endfb8':
        PATH_ACE = '{}/Lib80x/{}/{}'.format(xs_folder_path,foil,ace_filename)
        isotope_data = openmc.data.IncidentNeutron.from_ace(PATH_ACE)
        cross_section = (isotope_data[mt_number].xs['294K'](energy_bins))
    if xs_library == 'tendl21':
        PATH_ACE = '{}/tendl21c/tendl21c/{}'.format(xs_folder_path,tendl_extraction(isotope))
        isotope_data = openmc.data.IncidentNeutron.from_ace(PATH_ACE)
        cross_section = (isotope_data[mt_number].xs['294K'](energy_bins))
    if xs_library == 'irdff2':
        if isotope not in ['Cu65', 'Dy164']: #in ['In115','Au197','Fe56','Al27','Cu65','Ni58']:
            filename_irdff_format = ace_filename[:6] + '34y'
            PATH_ACE = '{}/IRDFF-II/endf_format/{}'.format(xs_folder_path,filename_irdff_format)
            cross_section = (irdff2_xs_extraction(PATH_ACE,mt_number,energy_bins))
        else:
            PATH_ACE = '{}/tendl21c/tendl21c/{}'.format(xs_folder_path,tendl_extraction(isotope))
            isotope_data = openmc.data.IncidentNeutron.from_ace(PATH_ACE)
            cross_section = (isotope_data[mt_number].xs['294K'](energy_bins))

#    if xs_library == 'irdff2':
#        if isotope not in ['Cu65', 'Dy164']: #in ['In115','Au197','Fe56','Al27','Cu65','Ni58']:
#            cross_section = (irdff2_xs_extraction(PATH_ACE,mt_number,energy_bins))
#        else:
#            xs_library == 'tendl21'
#            isotope_data = openmc.data.IncidentNeutron.from_ace(PATH_ACE)
#            cross_section = (isotope_data[mt_number].xs['294K'](energy_bins))
#    else:
#        isotope_data = openmc.data.IncidentNeutron.from_ace(PATH_ACE)
#        cross_section = (isotope_data[mt_number].xs['294K'](energy_bins))

    response_function = isotope_atom_density*foil_volume*cross_section 
    return energy_bins,isotope_atom_density,foil_volume,cross_section,response_function


data_dictionary = {}
data_dictionary['${}^{115}$In(n,$\gamma$)'] = reaction_info('In115','In','49115.800nc',102,7.29 ,0.6452)
data_dictionary['${}^{164}$Dy(n,$\gamma$) *'] = reaction_info('Dy164','Dy','66164.800nc',102,8.551,0.0435)
#data_dictionary['${}^{156}$Dy(n,$\gamma$)'] = reaction_info('Dy156','Dy','66156.800nc',102,8.551,0.0435)
#data_dictionary['${}^{63}$Cu(n,$\gamma$)']  = reaction_info('Cu63' ,'Cu','29063.800nc',102,8.83 ,0.7877)
## data_dictionary['${}^{116}$Cd(n,$\gamma$)'] = reaction_info('Cd116','Cd','48116.800nc',102,8.65 ,0.1328)
#data_dictionary['${}^{114}$Cd(n,$\gamma$)'] = reaction_info('Cd114','Cd','48114.800nc',102,8.65 ,0.1328)
#data_dictionary['${}^{110}$Cd(n,$\gamma$)'] = reaction_info('Cd110','Cd','48110.800nc',102,8.65 ,0.1328)
data_dictionary['${}^{197}$Au(n,$\gamma$)'] = reaction_info('Au197','Au','79197.800nc',102,19.3 ,0.1498)
#data_dictionary['${}^{111}$Cd(n,n\')']      = reaction_info('Cd111','Cd','48111.800nc',4  ,8.65 ,0.1328)
data_dictionary['${}^{115}$In(n,n\')']      = reaction_info('In115','In','49115.800nc',11004  ,7.29 ,0.6452)
## data_dictionary['${}^{58}$Ni(n,p) ']        = reaction_info('Ni58', 'Ni','28058.800nc',103,8.90 ,1.5845)
data_dictionary['${}^{65}$Cu(n,p) *']        = reaction_info('Cu65' ,'Cu','29065.800nc',103,8.83 ,0.7877)
data_dictionary['${}^{56}$Fe(n,p)']         = reaction_info('Fe56' ,'Fe','26056.800nc',103,7.874,1.4440)
data_dictionary['${}^{27}$Al(n,$\\alpha$)'] = reaction_info('Al27' ,'Al','13027.800nc',107,2.7  ,0.2334)
data_dictionary['${}^{197}$Au(n,2n) ']      = reaction_info('Au197','Au','79197.800nc',16 ,19.3 ,0.1498)
data_dictionary['${}^{93}$Nb(n,2n)']        = reaction_info('Nb93', 'Nb','41093.800nc',11016 ,8.57 ,0.7317)
#data_dictionary['${}^{116}$Cd(n,2n) ']      = reaction_info('Cd116','Cd','48116.800nc',16 ,8.65 ,0.1328)
#data_dictionary['${}^{112}$Cd(n,2n)']       = reaction_info('Cd112','Cd','48112.800nc',16 ,8.65 ,0.1328)
#data_dictionary['${}^{158}$Dy(n,2n)']       = reaction_info('Dy158','Dy','66158.800nc',16 ,8.551,0.0435)
#data_dictionary['${}^{65}$Cu(n,2n)  ']      = reaction_info('Cu65' ,'Cu','29065.800nc',16 ,8.83 ,0.7877)
data_dictionary['${}^{58}$Ni(n,2n) ']       = reaction_info('Ni58', 'Ni','28058.800nc',16 ,8.90 ,1.5845)

def reaction_rate_calc(reaction,isotope_name):
    flux_175 = [4318.042765398374, 14814.46317739122, 5424.852974161929, 6754.3096536459925, 6346.05242475784, 5301.239926105321, 8895.653569952787, 6523.479350686184, 8044.727643484646, 8823.973092601569, 12612.04961694252, 10265.727825025502, 10633.756929839874, 15304.037454044434, 11670.104453440292, 14086.982969321552, 17663.533946366624, 14696.413065256687, 18073.180149736407, 12945.152371755808, 18364.95991078865, 19128.094401195816, 18909.833129096736, 16694.48548106331, 21263.03902981031, 19647.419851724255, 20733.861636918817, 24521.65068334278, 21458.57732311212, 17935.685667634418, 19035.256007201464, 23457.01003276177, 20994.601272393094, 20989.360389449375, 27157.778506265982, 32774.798146072266, 11234.826131315629, 10833.760567397634, 5337.641125012645, 3485.0321188072217, 11118.952476613209, 13320.488661253605, 15400.955596406351, 22062.011663123467, 36113.760138317906, 44164.733517385466, 28952.417939833264, 25932.2083808752, 21040.349859699443, 54833.8554522872, 45311.06366176897, 22695.27135221706, 16514.805909369174, 7477.59875850594, 10416.521536114557, 23957.22407282016, 14516.966516166065, 11705.859560951701, 24004.699641443884, 20071.716106869644, 48824.74868646387, 48445.214605223155, 46763.699909132025, 24684.970350031177, 74519.32215918563, 44415.90012315226, 59181.22462580889, 36365.90778817834, 25246.470519814462, 94485.20080027204, 85628.42934213844, 37564.32586974248, 47034.50920096144, 53905.68566136057, 66169.85232899502, 57118.15694864307, 47159.155634978815, 56792.07538679848, 66997.58427948358, 67516.8462777175, 84638.64831161182, 78486.90281200547, 73674.9265131542, 85241.4005848564, 90963.25789746315, 94270.58313279644, 95969.15896991812, 217898.28437580162, 99008.26767223435, 54838.208159306996, 20729.347620622462, 9863.778793456177, 26817.14859470878, 282533.0410186167, 310133.5429241854, 168747.60931439567, 166608.00663289684, 338653.9198041544, 383936.41621479555, 200053.5497398328, 208276.88681730337, 223281.39642896515, 237222.36347443878, 264729.33655443246, 289344.8568829838, 307691.3215264174, 314640.5491916379, 316046.53676740185, 330289.61805773084, 329058.0581188052, 345117.96936477354, 409680.73590194574, 321057.6723774492, 750685.3889990158, 414449.461630013, 420846.24951817794, 410665.5537802989, 430433.4449913736, 442971.4739120111, 429648.47753512295, 423919.1204986436, 400166.0671688591, 427210.74693818204, 415253.8910555242, 395595.4265657715, 387987.85696875764, 391128.1303745563, 393014.0289898152, 255928.83536634577, 128727.11236621159, 65027.0979994264, 65690.63981574701, 246803.15538173972, 353996.4242326478, 344760.54129587876, 338080.7554083941, 324402.7371664211, 308329.0184122137, 286060.61247931165, 539422.916151527, 494539.9362459795, 441856.2175021492, 204011.13140385912, 192698.19524281044, 178653.5798099801, 170875.52653844195, 162645.51849320525, 150505.5617290813, 144576.16795118924, 96976.54040034405, 44824.38058668863, 133176.00812006954, 132270.21736616406, 123244.05804903673, 125092.29169735612, 128455.0692537148, 125783.54252846059, 110978.41064445073, 100431.04532500621, 80194.58635974528, 40478.22342049091, 14153.934294073422, 9579.679418773727, 7020.0838589020095, 15597.93030573221, 467751.789563248, 528756.0481820796, 365355.3067896321, 103478.480713344, 1905.6536264222939, 0.0, 0.0, 0.0, 0.0, 0.0]
    flux_175_conc = np.concatenate((flux_175, flux_175[-1:]))
    cross_section_175 = data_dictionary[reaction][3]
    sigma_phi =  np.dot(flux_175_conc,cross_section_175)
    no_of_atoms = data_dictionary[reaction][1]*data_dictionary[reaction][2]
    decay_fraction = 1 - np.exp(-(128.5*60)*( np.log(2) / ag.Decay2012Database().gethalflife(isotope_name)))
    print(no_of_atoms*sigma_phi*decay_fraction)

   # a_0 = N sigma phi (1-exp(-t (irrad) lambda))
#reaction_rate_calc(list(data_dictionary.keys())[8],'Nb92m')

color = iter(cm.rainbow(np.linspace(0, 1, len(data_dictionary.keys()))))

fig, (ax1,ax2) = plt.subplots(1,2,figsize=(14,7),gridspec_kw={'width_ratios': [2, 3.5]})
for reaction in list(data_dictionary.keys())[0:]:
    c=next(color)
    ax1.step(data_dictionary[reaction]      [0]/1e6, data_dictionary[reaction]    [4], label= reaction,c=c)    
    ax2.step(data_dictionary[reaction]      [0]/1e6, data_dictionary[reaction]    [4], label= reaction,c=c)
    df_response_matrices = pd.DataFrame(data_dictionary[reaction]   [4][1:]).T
    df_response_matrices.to_csv('response_matrix_{}.csv'.format(xs_library),index=False,header=False,mode='a')
df_reactions = pd.DataFrame(data_dictionary.keys())
df_reactions.to_csv('reaction_rate_labels_{}.csv'.format(xs_library),index=False,header=False,mode='w')
df_energygroup = pd.DataFrame(energy_group.group_edges)
df_energygroup.to_csv('group_structure.csv',index=False,header=False,mode='w')


ax1.set_xlim(1e-8,1e0)
ax1.set_ylim(1e-11,1e3)
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.grid()
ax2.set_xlim(1e0,16)
ax2.set_ylim(1e-11,1e3)
ax2.tick_params(axis='y',left=False,labelleft=False)
ax2.set_yscale('log')
ax2.grid()  
ax2.legend(loc="upper left", bbox_to_anchor=(0.025, 0.955), borderaxespad=0, frameon=True, fontsize=10,fancybox=False,facecolor='white',framealpha=1,ncol=4)
fig.supylabel('Response function Rn(E)',y=0.5)
fig.supxlabel('Neutron energy (MeV)',x=0.4,y=0.03)
fig.tight_layout()
plt.savefig('response_functions_{}_TEST.png'.format(xs_library))
