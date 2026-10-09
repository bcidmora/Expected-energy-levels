from pylatex import Document, LongTable, MultiColumn, Math, Table, Tabularx, Tabular, Section, Center, Alignat, Subsection, Subsubsection
from pylatex.utils import NoEscape

import matplotlib.pyplot as plt
import numpy as np

from ensembles_info import ensE250
import ensembles_info as ens
import list_of_lists as lol


### reads a txt file with the info of the expected energy levels.
def READ_TXT_FILE(the_name_file):
    with open(the_name_file) as file: 
        read_file = file.readlines()
    return read_file


#--------------  MODIFICATIONS TO THE HADRON NAME SCHEME  ---------------- 
#-------------------  OR SOMETHING LIKE THAT ---------------- 

#########################################################################
#                                                                       #
#                                                                       #
#       ** Hadrons info and momentum:                                   #
#        Extracts the momentum for each hadron in the string and        #
#        it includes the corresponding momentum to those without        #
#        explicit momentum.                                             #
#                                                                       #
#       ** Splitting the operators:                                     #
#        It separates the line in Colin's file that has the energy,     #
#        the multiplicity and the hadrons in that level.                #
#                                                                       #
#                                                                       #
#                                                                       #
#########################################################################




### adds the momentum to the operator when there is no momentum in the operator name.
# the_operator: is a string with all the operators in it such as pi_PSQ0_Sigma_PSQ3_isosinglet_0
# the_sqred_lattice_momentum: is the squared momentum in the center of mass frame
def HADRONS_INFO_AND_MOMENTUM(the_operator, the_sqred_lattice_momentum):
    the_decomposed_operators = []
    if 'PSQ' not in the_operator:
        the_decomposed_operators.append([the_operator, int(the_sqred_lattice_momentum)])
    else:
        the_list_the_operators = the_operator.split('_')
        for i in range(0, len(the_list_the_operators)-2, 2):
            da_mom = the_list_the_operators[i+1][3:]
            if len(da_mom)==1:
                da_mom = the_list_the_operators[i+1][3]
            elif len(da_mom)==2:
                if da_mom[-1]=='A' or da_mom[-1]=='B':
                    da_mom = the_list_the_operators[i+1][3]
                else:
                    da_mom = the_list_the_operators[i+1][3:]
            elif len(da_mom)==3:
                if da_mom[-1]=='A' or da_mom[-1]=='B':
                    da_mom = the_list_the_operators[i+1][3:-1]
                else:
                    da_mom = the_list_the_operators[i+1][3:]
            the_decomposed_operators.append([the_list_the_operators[i], int(da_mom)])
    return the_decomposed_operators



### Comments:
# This function returns the names of the irreps in tex format for plotting
def PLOT_HADRON_LABELINGS(the_irrep_name):
    the_irrep_name = the_irrep_name.replace(" ","")
    return lol.the_list_irreps[the_irrep_name].texName




###  --------------  PLOTTING FUNCTIONS  ----------------


####################################################
#                                                  #
#                                                  #
#  This routine plots the thresholds using all     #
#  the inputs of the ensembles, only for the       #
#  Lambda, it looks at the KN and PS thresholds.   #
#                                                  #
#                                                  #
####################################################


def PLOT_ENERGY_LEVELS(list_of_energies, the_ref_levels, the_y_axis_label, the_name_plot, the_nr_levels):
    line_styles=["-","--","-.",":"]
    line_colors=["#b90f22", "#5d83d5","#ffa500","#008000","#c44601","#f57600","#5ba300","#e6308a" ]
    the_plot = plt.figure()
    the_x_axis_labels = []
    for ii in range(len(list_of_energies)):
        the_name_irrep = PLOT_HADRON_LABELINGS(list_of_energies[ii][0])
        x_axis, y_axis = f'{the_name_irrep}({list_of_energies[ii][1]})', list_of_energies[ii][2]
        
        if x_axis not in the_x_axis_labels: the_x_axis_labels.append(x_axis)
        
        plt.plot(x_axis, y_axis, marker='_', ls='None', ms=14, markeredgewidth=3, lw=0.5, zorder=3, color = line_colors[1])
    
    for ii in range(len(the_ref_levels)):
        plt.axhline(the_ref_levels[ii][0], ls=line_styles[ii], lw=1.5, color=line_colors[0], label=the_ref_levels[ii][1])
        
    plt.xlabel('Irreducible representations',fontsize=16)
    the_x_axis = np.arange(0, len(the_x_axis_labels))
    plt.xticks(the_x_axis, the_x_axis_labels, rotation=45, fontsize=12)
    plt.tick_params(axis='y', labelsize=12)
    plt.ylabel(the_y_axis_label, fontsize=16)
    plt.tight_layout()
    # the_plot.legend(loc='upper center', bbox_to_anchor=(0.5, top_y+0.22), ncol = len(the_ref_levels),fontsize=16, columnspacing=1.5, handletextpad=0.3)
    # the_plot.legend(ncol = len(the_ref_levels),fontsize=16, columnspacing=1.75, handletextpad=0.3, loc='upper center')
    # plt.show()
    plt.margins(x=0.1, y=0.1)
    plt.tight_layout()
    the_plot.savefig(the_name_plot, bbox_inches='tight')
    


### This function returns the names of the hadrons in tex format for plotting
def PLOT_SINGLE_HADRON_NAMES(the_hadron_name):
    return lol.the_list_hadrons[the_hadron_name].texName
    
    

###############################################################
#                                                             #
#       Extra function to reorganize, to prepare things       #
#       for plots, and some menu options                      # 
#                                                             #
###############################################################



### It reorganizes the list of selected levels to plot them in increasing sqred momentum.
# the_irrep: The info of the irreducible representation
# the_sqred_mom: the squared momentum in this irrep
# the_energy_list: the energy list to be included in the tables
# the_energy_list_plot: the list for the plots.
def LIST_FOR_PLOT(the_irrep, the_sqred_mom, the_energy_list, the_energy_list_plot):
    for zz in range(len(the_energy_list)):
        the_energy_list_plot.append([the_irrep,the_sqred_mom,the_energy_list[zz][0], the_energy_list[zz][2]])
    return the_energy_list_plot




### It counts how many different irreps are there to plot.
# an_energy_list: the final list, it counts how many levels to plot in the end.
def HOW_MANY_LEVELS_PLOT(an_energy_list):
    the_first_irrep=an_energy_list[0][0]
    how_many_levels=1
    for ii in range(len(an_energy_list)):
        if an_energy_list[ii][0]!=the_first_irrep: 
            how_many_levels+=1
            the_first_irrep=an_energy_list[ii][0]
        else: continue
    return how_many_levels


### It shows the list of quantum numbers to choose from.
def MENU_QUANTUM_NUMBERS():
    print("Choose: (I: isospin, S: strangeness, B: Baryon nr.)")    
    for ii, the_qn in enumerate(lol.the_quantum_numbers):
        the_temp_s=(f'-{the_qn.theStrangeness[2:]}' if 'm' in the_qn.theStrangeness else the_qn.theStrangeness[1:])
        print(f"[{ii+1}] {the_qn.theParticleType} {the_qn.theIsospinStr} S={the_temp_s} B={the_qn.theBaryonNumber[-1]}")
    print(f"[{len(lol.the_quantum_numbers)+1}] All")
    
    the_quantum_number=input("The choice can be an integer or several separated by ',', e.g. 1,8...\n")
    the_choices=[int(item) for item in list(the_quantum_number.split(","))]
    
    if (len(lol.the_quantum_numbers)+1) in the_choices:
        return lol.the_quantum_numbers.copy()
    else:
        return [lol.the_quantum_numbers[item-1] for item in the_choices]            
            

### It gives the selection of ensembles that are available to be computed.
def MENU_ENSEMBLES():
    print("For which ensemble(s) do you want to obtain the expected energy levels?")
    for ii in range(len(ens.ensembleList)):
        print(f'    [{ii+1}] {ens.ensembleList[ii]['ens_name']}')
    print(f'    [{len(ens.ensembleList)+1}] All')
    whichEnsemble=str(input('The choice can be an integer or op1,op2,...\n'))
    ensChoices=[]
    if "," in whichEnsemble:
        ensChoices_pre=whichEnsemble.split(",")
        for num in ensChoices_pre:
            ensChoices.append(ens.ensembleList[int(num)-1])
    elif whichEnsemble=="All" or whichEnsemble==f'{len(ens.ensembleList)+1}':
        ensChoices=ens.ensembleList
    elif int(whichEnsemble)<=len(ens.ensembleList)+1:
        ensChoices.append(ens.ensembleList[int(whichEnsemble)-1])
    else:
        print("Incorrect choice.")
    return ensChoices




### It gives the masses of the hadrons that are available in ALL ensembles, it excludes the ones that are not in every ensemble.
def MENU_HADRONS():
    print("Choose the hadron(s):")
    for ii, (key, hadron) in enumerate(lol.the_list_hadrons.items()):
        print(f'[{ii+1}] {hadron.longName}')
    
    choice=str(input("You can choose one or more hadrons, eg. 1,1,6: "))
    if len(choice)>1:
        choice=choice.split(",")
    choice=list(choice)
    
    the_hadron=[]
    hadrons = list(lol.the_list_hadrons.values())
    
    for item in choice:
        item=int(item)
        hadron=hadrons[item-1]
        the_hadron.append([hadron.mass, hadron.texName])
    return the_hadron



### It constructs a list of thresholds, which is also a list of hadrons to obtain the threhsolds over the nucleon mass.
def CHOICE_THRESHOLDS_PLOT():
    how_many_thresholds=int(input("How many thresholds do you want to plot?\n"))
    the_hadrons=[]
    for jj in range(how_many_thresholds):
        print(f"For threshold nr. {jj+1}: ")
        the_hadrons.append(MENU_HADRONS())
    return the_hadrons
