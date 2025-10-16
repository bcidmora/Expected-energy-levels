import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
import energy_functions as enf
import ensembles_info as ens

ensChoices = enf.MENU_ENSEMBLES()

the_markers_list = ["s","o","^","*","+",">","<"]
the_colors_list=["#b90f22", "#5d83d5","#ffa500","#008000","#c44601","#f57600","#5ba300","#e6308a" ]

the_min_y=100000000
the_max_y=0

the_x_axis_labels = []
the_y_axis_vals = []
the_ensemble_names = []

item = 0 
hadronsInfo_0 = [['pi', '%s'%str(ensChoices[item]["pion_mass"]), '%s'%str(enf.FM_TO_MEV(ensChoices[item]["pion_mass"],ensChoices[item]["latt_spacing"]))], 
    ['KB','%s'%str(ensChoices[item]["kaon_mass"]), '%s'%str(enf.FM_TO_MEV(ensChoices[item]["kaon_mass"],ensChoices[item]["latt_spacing"]))], 
    ['eta','%s'%str(ensChoices[item]["eta_mass"]), '%s'%str(enf.FM_TO_MEV(ensChoices[item]["eta_mass"],ensChoices[item]["latt_spacing"]))], 
    ['N', '%s'%str(ensChoices[item]["nucleon_mass"]), '%s'%str(enf.FM_TO_MEV(ensChoices[item]["nucleon_mass"],ensChoices[item]["latt_spacing"]))], 
    ['etaprime-958','%s'%str(ensChoices[item]["etaprime_mass"]), '%s'%str(enf.FM_TO_MEV(ensChoices[item]["etaprime_mass"],ensChoices[item]["latt_spacing"]))], 
    ['Lambda', '%s '%str(ensChoices[item]["lambda_mass"]), '%s'%str(enf.FM_TO_MEV(ensChoices[item]["lambda_mass"],ensChoices[item]["latt_spacing"]))], 
    ['Sigma','%s'%str(ensChoices[item]["sigma_mass"]), '%s'%str(enf.FM_TO_MEV(ensChoices[item]["sigma_mass"],ensChoices[item]["latt_spacing"]))], 
    ['Delta-1232', '%s'%str(ensChoices[item]["delta_mass"]), '%s'%str(enf.FM_TO_MEV(ensChoices[item]["delta_mass"],ensChoices[item]["latt_spacing"]))], 
    ['Xi', '%s'%str(ensChoices[item]["xi_mass"]), '%s'%str(enf.FM_TO_MEV(ensChoices[item]["xi_mass"],ensChoices[item]["latt_spacing"]))], 
    ['Sigma-1385', '%s'%str(ensChoices[item]["sigmastar_mass"]), '%s'%str(enf.FM_TO_MEV(ensChoices[item]["sigmastar_mass"],ensChoices[item]["latt_spacing"]))], 
    ['Xi-1530', '%s'%str(ensChoices[item]["xistar_mass"]), '%s'%str(enf.FM_TO_MEV(ensChoices[item]["xistar_mass"],ensChoices[item]["latt_spacing"]))], 
    ['Omega', '%s'%str(ensChoices[item]["omega_mass"]), '%s'%str(enf.FM_TO_MEV(ensChoices[item]["omega_mass"],ensChoices[item]["latt_spacing"]))]]

for had in hadronsInfo_0:
    the_x_axis_labels.append(enf.PLOT_SINGLE_HADRON_NAMES(had[0]))

the_x_axis = np.arange(0,len(the_x_axis_labels)*12,12)

for item in range(len(ensChoices)):
    hadronsInfo = [['pi', '%s'%str(ensChoices[item]["pion_mass"]), '%s'%str(enf.FM_TO_MEV(ensChoices[item]["pion_mass"],ensChoices[item]["latt_spacing"]))], 
    ['KB','%s'%str(ensChoices[item]["kaon_mass"]), '%s'%str(enf.FM_TO_MEV(ensChoices[item]["kaon_mass"],ensChoices[item]["latt_spacing"]))], 
    ['eta','%s'%str(ensChoices[item]["eta_mass"]), '%s'%str(enf.FM_TO_MEV(ensChoices[item]["eta_mass"],ensChoices[item]["latt_spacing"]))], 
    ['N', '%s'%str(ensChoices[item]["nucleon_mass"]), '%s'%str(enf.FM_TO_MEV(ensChoices[item]["nucleon_mass"],ensChoices[item]["latt_spacing"]))], 
    ['etaprime-958','%s'%str(ensChoices[item]["etaprime_mass"]), '%s'%str(enf.FM_TO_MEV(ensChoices[item]["etaprime_mass"],ensChoices[item]["latt_spacing"]))], 
    ['Lambda', '%s '%str(ensChoices[item]["lambda_mass"]), '%s'%str(enf.FM_TO_MEV(ensChoices[item]["lambda_mass"],ensChoices[item]["latt_spacing"]))], 
    ['Sigma','%s'%str(ensChoices[item]["sigma_mass"]), '%s'%str(enf.FM_TO_MEV(ensChoices[item]["sigma_mass"],ensChoices[item]["latt_spacing"]))], 
    ['Delta-1232', '%s'%str(ensChoices[item]["delta_mass"]), '%s'%str(enf.FM_TO_MEV(ensChoices[item]["delta_mass"],ensChoices[item]["latt_spacing"]))], 
    ['Xi', '%s'%str(ensChoices[item]["xi_mass"]), '%s'%str(enf.FM_TO_MEV(ensChoices[item]["xi_mass"],ensChoices[item]["latt_spacing"]))], 
    ['Sigma-1385', '%s'%str(ensChoices[item]["sigmastar_mass"]), '%s'%str(enf.FM_TO_MEV(ensChoices[item]["sigmastar_mass"],ensChoices[item]["latt_spacing"]))], 
    ['Xi-1530', '%s'%str(ensChoices[item]["xistar_mass"]), '%s'%str(enf.FM_TO_MEV(ensChoices[item]["xistar_mass"],ensChoices[item]["latt_spacing"]))], 
    ['Omega', '%s'%str(ensChoices[item]["omega_mass"]), '%s'%str(enf.FM_TO_MEV(ensChoices[item]["omega_mass"],ensChoices[item]["latt_spacing"]))]]
    
    the_ensemble_names.append(ensChoices[item]["ens_name"])
    
    the_y_axis = []
    for had in hadronsInfo:
        the_y_axis.append(float(had[2])/float(hadronsInfo[3][2]))
    the_y_axis_vals.append(the_y_axis)

offset = .5
the_plot = plt.figure()

for ii in range(len(the_y_axis_vals)):
    if the_min_y > min(the_y_axis_vals[ii]): 
        the_min_y = min(the_y_axis_vals[ii])
    if the_max_y < max(the_y_axis_vals[ii]): 
        the_max_y = max(the_y_axis_vals[ii])
    plt.plot(the_x_axis+offset*ii, the_y_axis_vals[ii], marker = the_markers_list[ii], label = the_ensemble_names[ii], color = the_colors_list[ii], ms=7.5,ls="None")
the_new_offset = (offset*(len(the_y_axis_vals)-1))/2
the_x_axis = np.arange(the_x_axis[0]+the_new_offset, the_x_axis[-1]+the_new_offset+1,12)
plt.xticks(the_x_axis, the_x_axis_labels, fontsize=20)
plt.tick_params(axis='y', labelsize=11)
plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
plt.ylim([the_min_y*.5, the_max_y*1.1])
plt.ylabel(r'$m_{\mathrm{H}}/m_{\mathrm{N}}$ [MeV]',fontsize=16)
plt.tight_layout()
fig=plt.gcf()
axes=fig.axes
top_y = max(ax.get_position().y1 for ax in axes)
plt.legend(fontsize=13, loc='upper left',handletextpad=0.3)
plt.show()
the_plot.savefig("single_hadrons_ensemble_comparison.pdf", bbox_inches='tight')
