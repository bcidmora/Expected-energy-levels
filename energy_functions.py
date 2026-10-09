from pylatex import Document, LongTable, MultiColumn, Math, Table, Tabularx, Tabular, Section, Center, Alignat, Subsection, Subsubsection
from pylatex.utils import NoEscape

import numpy as np

from ensembles_info import ensE250
import ensembles_info as ens
import labels_and_plots as lap


fmToMeV = 197.327 # This is to convert to MeV


#--------------  MATH ROUTINES ---------------- 

########################################################
#                                                      #
#                                                      #
#    Routines transform from fm to MeV and viceversa   #
#                                                      #
#   Computes the momentum in the Lattice:              #
#           p = 2 * pi * d / L                         #
#                                                      #
#   Relativistic Relation:                             #
#           E^2 = m^2 + p^2                            #
#                                                      #
#   Read a txt file: gets the info from E250           # 
#                                                      #
#   Shifts the energy from Elab to Ecm                 #
#                                                      #
#                                                      #
########################################################


### transforms the momentum in lattice units to MeV units, where the_latt_spacing=a.
# the_quantity: the variable to transform
# the_latt_spacing: lattice spacing a
def FM_TO_MEV(the_quantity, the_latt_spacing):
    return the_quantity*fmToMeV/the_latt_spacing


### transforms from MeV units to Lattice units, where the_latt_spacing=a.
# the_quatity: the variable to transform
# the_latt_spacing: the lattice spacing a
def MEV_TO_FM(the_quantity, the_latt_spacing):
    return the_quantity*the_latt_spacing/fmToMeV


### transforms the units of momentum in a volume (in the lattice) to units  of momentum in energy
# the_sqred_mom_units: The units of squared momentum, this is an integer
# the_latt_size: the lattice extent L
def MOMENTUM_COMPUTATION(the_sqred_mom_units, the_latt_size):
    return 2.* np.pi * (np.sqrt(float(the_sqred_mom_units))) / float(the_latt_size)


### returns the energy squared where the momentum and mass are in the same units
# the_had_mass: the hadron mass
# the_had_mom: the hadron momentum 
def RELATIVISTIC_RELATION(the_had_mass, the_had_mom):
    return the_had_mass**2 + the_had_mom**2


### shifts the lab frame-energy to the CM frame
# the_lab_energy: energy in the laboratory frame
# the_momentum: the momentum in the center of mass frame
def SHIFT_TO_ECM(the_lab_energy, the_momentum):
    return np.sqrt(the_lab_energy**2 - the_momentum**2)


### This function returns the possible momentum combinations along x,y,z axis that give that total momentum squared
# the_hadron: this is the hadron written as pi(0) for example:
def POSSIBLE_MOMENTUM_HADRONS(the_hadron):
    the_mom = int(the_hadron[the_hadron.index('(')+1:the_hadron.index(')')])
    if the_mom==0:
        return 'P = (0, 0, 0)'
    elif the_mom==1:
        return 'P = (k, 0, 0); k = 1/-1'
    elif the_mom==2:
        return 'P = (k, k, 0); k = 1/-1'
    elif the_mom==3:
        return 'P = (k, k, k); k = 1/-1'
    elif the_mom==4:
        return 'P = (k, 0, 0); k = 2/-2'
    elif the_mom==5:
        return 'P = (2, 1, 0)'
    elif the_mom==6:
        return 'P = (2, 1, 1)'



### separates the full operator string into separated hadrons with their momenta
## the_operator_row: It containes the energy levels with its multiplicity and the energy value and it splits them into a list (new_operator_row)
def SPLITTING_THE_OPERATORS_ROW(the_operator_row):
    the_pos_1=the_operator_row.index("(")
    the_pos_2=the_operator_row.index(")")
    the_pos_3=the_operator_row.index("\n")
    the_new_operator_row=[]
    the_energy_val=float(the_operator_row[0:the_pos_1])
    the_multiplicity_val=int(the_operator_row[the_pos_1+1:the_pos_2])
    the_new_operator_row.append(the_energy_val)
    the_new_operator_row.append(the_multiplicity_val)
    if "*" in the_operator_row:
        the_pos_4=the_operator_row.index("*")
        the_operators_string=the_operator_row[the_pos_2+1:the_pos_4]
        while " " in the_operators_string: the_operators_string=the_operators_string.replace(" ", "")
        the_new_operator_row.append(the_operators_string)
        the_new_operator_row.append(the_operator_row[the_pos_4:the_pos_4+3])
    elif "#" in the_operator_row:
        the_pos_4=the_operator_row.index("#")
        the_operators_string=the_operator_row[the_pos_2+1:the_pos_4]
        while " " in the_operators_string: the_operators_string=the_operators_string.replace(" ", "")
        the_new_operator_row.append(the_operators_string)
        the_new_operator_row.append(the_operator_row[the_pos_4:the_pos_4+4])
    else:
        the_operators_string=the_operator_row[the_pos_2+1:the_pos_3]
        while " " in the_operators_string: the_operators_string=the_operators_string.replace(" ", "")
        the_new_operator_row.append(the_operators_string)
    return the_new_operator_row





#--------------  ENERGY AND STUFF  ----------------

####################################################
#                                                  #
#                                                  #
#  Calculates the energy in fm Units,              #
#                                                  #
#                                                  #
####################################################



### This function calculates the final energy of a certain non-interacting level in lattice units.
### If there is no information about the hadron in that ensemble, then it uses the values of E250
# the_operators: [hadron name, squared momentum]
# the_list_of_masses: Is the list of masses known for this ensemble
# the_latt_extent: the lattice extent L
# the_sqred_lattice_momentum: is the squared momentum in the center of mass frame
# the_unknown: these are the unknown masses for a specific ensemble.
# the_ecm: If True, it moves the energy calculation from lab frame to the center of mass frame
def CALCULATING_FINAL_ENERGY(the_operators, the_list_of_masses, the_latt_extent, the_sqred_lattice_momentum, the_unknown, the_ecm):
    the_reordered_masses=list(zip(*the_list_of_masses))
    the_total_energy = 0.
    for item in the_operators:
        the_had_mom = MOMENTUM_COMPUTATION(float(item[1]), the_latt_extent)
        if item[0] not in the_reordered_masses[0]:
            if the_unknown: the_had_mass = ensE250[item[0]]
            else: the_total_energy=0.;break
        else: the_had_mass = float(the_reordered_masses[1][the_reordered_masses[0].index(item[0])])
        the_energy_for_one_hadron = np.sqrt(RELATIVISTIC_RELATION(the_had_mass, the_had_mom))
        the_total_energy += float(the_energy_for_one_hadron)
    if the_ecm:
        if the_unknown:
            the_mom=MOMENTUM_COMPUTATION(float(the_sqred_lattice_momentum),ensE250["s_extent"])
            the_energy_cm=SHIFT_TO_ECM(the_total_energy,the_mom)
        else:
            if the_total_energy==0.:the_energy_cm=0.
            else:
                the_mom=MOMENTUM_COMPUTATION(float(the_sqred_lattice_momentum),the_latt_extent)
                the_energy_cm=SHIFT_TO_ECM(the_total_energy,the_mom)
        return the_energy_cm/float(the_reordered_masses[1][the_reordered_masses[0].index("N")])
    else: return the_total_energy/float(the_reordered_masses[1][the_reordered_masses[0].index("N")])





### It calculates all the energy levels and returns them in a list: {energy, multiplicity, operators}
# the_flavor_sector: this is the line in the txt file that contains the energy, multiplicity and operators
# the_list_of_masses: the list of known masses for a certain ensemble
# the_latt_extent: spatial extent
# the_sqred_lattice_momentum: squared of momentum in the center of mass frame
# the_unknown: the unknown masses for this ensemble are taken from the E250
# the_ecm: It calculates the energy in the center of mass frame
def ENERGY_LIST_RAW(the_flavor_sector, the_list_of_masses, the_latt_extent, the_sqred_lattice_momentum, the_unknown, the_ecm):
    the_energy_multi_hads_list=[]
    for ii in range(len(the_flavor_sector)):
        the_row = SPLITTING_THE_OPERATORS_ROW(the_flavor_sector[ii])
        the_hadrons_row = the_row[2]
        hadrons_and_momenta = lap.HADRONS_INFO_AND_MOMENTUM(the_hadrons_row, the_sqred_lattice_momentum)
        the_energy = CALCULATING_FINAL_ENERGY(hadrons_and_momenta, the_list_of_masses, the_latt_extent, the_sqred_lattice_momentum, the_unknown, the_ecm) 
        if the_energy==0.: ii=ii+1
        else:
            if len(the_row)>3: the_energy_multi_hads_list.append([the_energy, the_row[1], the_row[2] + " " + the_row[3]])
            else: the_energy_multi_hads_list.append([the_energy, the_row[1], the_row[2]])
    new_multi_hads_list= sorted(the_energy_multi_hads_list, key=lambda k:[k[0], k[1], k[2]])
    return new_multi_hads_list




### It prepares everything to be in a filtered list like: {energy, multiplicity, operators} (The cutoff happens here)
# the_hadrons_energy_list: This is the list to put in the tables
# the_threshold: This is the point where the full list of energies is cut
# the_threeparticle: is a boolean that tells to cut one level after the three particle threshold
def ENERGY_LIST_TABLES(the_hadrons_energy_list, the_threshold, the_threeparticle):
    the_final_energy_list=[]
    for ii in range(len(the_hadrons_energy_list)):
        if the_threeparticle==True:
            if float(the_hadrons_energy_list[ii][0])>the_threshold:
                if float(the_hadrons_energy_list[ii][0])<=the_threshold*1.05:
                    the_final_energy_list.append(the_hadrons_energy_list[ii])
                    the_final_energy_list.append(the_hadrons_energy_list[ii+1])
                    break
                else:break
            else: the_final_energy_list.append(the_hadrons_energy_list[ii])
        else:
            if float(the_hadrons_energy_list[ii][0])>the_threshold:break
            else: the_final_energy_list.append(the_hadrons_energy_list[ii])
    return the_final_energy_list



### Checks all irreps in the tables
def FINAL_LIST_OF_OPERATORS(the_hadrons_list, the_sqred_lattice_momentum):
    the_hadrons_and_momenta=[]
    for ii in range(len(the_hadrons_list)):
        the_hads_mom = lap.HADRONS_INFO_AND_MOMENTUM(the_hadrons_list[ii][2],the_sqred_lattice_momentum)
        for item in the_hads_mom:
            the_new_item = f'{item[0]}({item[1]})'
            if the_new_item in the_hadrons_and_momenta: continue
            else:the_hadrons_and_momenta.append(the_new_item)
    return the_hadrons_and_momenta



### it sums up all the hadron masses given (over the nucleon mass). Usually for a threshold
def SUMMING_HADRON_MASSES(the_threshold,the_ensemble):
    result=0.
    result_name=''
    for hadron in the_threshold:
        result+=the_ensemble[hadron[0]]/the_ensemble["nucleon_mass"]
        result_name+=hadron[1]
    return [result, result_name]


### It actually goes over all thresholds and sums up each hadron to get each threshold value.
def GETTING_THRESHOLDS(the_ensemble,the_thresholds_list):
    the_final_thresholds=[]
    for threshold in the_thresholds_list:
        the_final_thresholds.append(SUMMING_HADRON_MASSES(threshold,the_ensemble))
    return the_final_thresholds
    
