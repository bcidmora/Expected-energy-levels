import numpy as np
from pylatex import Document, LongTable, MultiColumn, Math, Table, Tabularx, Tabular, Section, Center, Alignat, Subsection, Subsubsection
from pylatex.utils import NoEscape
import matplotlib.pyplot as plt
from ensembles_info import ensE250
import ensembles_info as ens


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

    # transforms the momentum in lattice units to MeV units, where latt_size=a.
def FM_TO_MEV(the_quantity, latt_size):
    return the_quantity*fmToMeV/latt_size


    # transforms from MeV units to Lattice units, where latt_size=a.
def MEV_TO_FM(the_quantity, latt_size):
    return the_quantity*latt_size/fmToMeV


    # transforms the units of momentum in a volume (in the lattice) to units  of momentum in energy
def MOMENTUM_COMPUTATION(the_sqred_mom_units, latt_extent):
    return 2.* np.pi * (np.sqrt(float(the_sqred_mom_units))) / float(latt_extent)


    # returns the energy squared
def RELATIVISTIC_RELATION(the_had_mass, the_had_mom):
    return the_had_mass**2 + the_had_mom**2


    # reads a txt file with the info of the expected energy levels.
def READ_TXT_FILE(name_file):
    with open(name_file) as file: 
        read_file = file.readlines()
    return read_file


    # shifts the lab frame-energy to the CM frame
def SHIFT_TO_ECM(the_lab_energy, the_momentum):
    return np.sqrt(the_lab_energy**2 - the_momentum**2)



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




    # adds the momentum to the operator when there is no momentum in the operator name.
def HADRONS_INFO_AND_MOMENTUM(the_operator, the_sqred_lattice_momentum):
    decomposed_operators = []
    if 'PSQ' not in the_operator:
        decomposed_operators.append([the_operator, int(the_sqred_lattice_momentum)])
    else:
        list_the_operators = the_operator.split('_')
        for i in range(0, len(list_the_operators)-2, 2):
            da_mom = list_the_operators[i+1][3:]
            if len(da_mom)==1:
                da_mom = list_the_operators[i+1][3]
            elif len(da_mom)==2:
                if da_mom[-1]=='A' or da_mom[-1]=='B':
                    da_mom = list_the_operators[i+1][3]
                else:
                    da_mom = list_the_operators[i+1][3:]
            elif len(da_mom)==3:
                if da_mom[-1]=='A' or da_mom[-1]=='B':
                    da_mom = list_the_operators[i+1][3:-1]
                else:
                    da_mom = list_the_operators[i+1][3:]
            decomposed_operators.append([list_the_operators[i], int(da_mom)])
    return decomposed_operators

    


    # separates the full operator string into separated hadrons with their momenta
def SPLITTING_THE_OPERATORS_ROW(the_operator_row):
    pos_1=the_operator_row.index("(")
    pos_2=the_operator_row.index(")")
    pos_3=the_operator_row.index("\n")
    new_operator_row=[]
    energy_val=float(the_operator_row[0:pos_1])
    multiplicity_val=int(the_operator_row[pos_1+1:pos_2])
    new_operator_row.append(energy_val)
    new_operator_row.append(multiplicity_val)
    if "*" in the_operator_row:
        pos_4=the_operator_row.index("*")
        operators_string=the_operator_row[pos_2+1:pos_4]
        while " " in operators_string: operators_string=operators_string.replace(" ", "")
        new_operator_row.append(operators_string)
        new_operator_row.append(the_operator_row[pos_4:pos_4+3])
    elif "#" in the_operator_row:
        pos_4=the_operator_row.index("#")
        operators_string=the_operator_row[pos_2+1:pos_4]
        while " " in operators_string: operators_string=operators_string.replace(" ", "")
        new_operator_row.append(operators_string)
        new_operator_row.append(the_operator_row[pos_4:pos_4+4])
    else:
        operators_string=the_operator_row[pos_2+1:pos_3]
        while " " in operators_string: operators_string=operators_string.replace(" ", "")
        new_operator_row.append(operators_string)
    return new_operator_row








#--------------  ENERGY AND STUFF  ----------------

####################################################
#                                                  #
#                                                  #
#  Calculates the energy in fm Units,              #
#                                                  #
#                                                  #
####################################################



    # This function calculates the final energy of a certain non-interacting level in lattice units.
    # If there is no information about the hadron in that ensemble, then it uses the values of E250
def CALCULATING_FINAL_ENERGY(the_operators, the_list_of_masses, latt_extent, the_sqred_lattice_momentum,the_unknown,the_ecm):
    reordered_masses=list(zip(*the_list_of_masses))
    the_total_energy = 0.
    for item in the_operators:
        the_had_mom = MOMENTUM_COMPUTATION(float(item[1]), latt_extent)
        if item[0] not in reordered_masses[0]:
            if the_unknown: the_had_mass = ensE250[item[0]]
            else: the_total_energy=0.;break
        else: the_had_mass = float(reordered_masses[1][reordered_masses[0].index(item[0])])
        the_energy_for_one_hadron = np.sqrt(RELATIVISTIC_RELATION(the_had_mass, the_had_mom))
        the_total_energy += float(the_energy_for_one_hadron)
    if the_ecm:
        if the_unknown:
            the_mom=MOMENTUM_COMPUTATION(float(the_sqred_lattice_momentum),ensE250["s_extent"])
            the_energy_cm=SHIFT_TO_ECM(the_total_energy,the_mom)
        else:
            if the_total_energy==0.:the_energy_cm=0.
            else:
                the_mom=MOMENTUM_COMPUTATION(float(the_sqred_lattice_momentum),latt_extent)
                the_energy_cm=SHIFT_TO_ECM(the_total_energy,the_mom)
        return the_energy_cm/float(reordered_masses[1][reordered_masses[0].index("N")])
    else: return the_total_energy/float(reordered_masses[1][reordered_masses[0].index("N")])





    # It calculates all the energy levels and returns them in a list: {energy, multiplicity, operators}
def ENERGY_LIST_RAW(the_flavor_sector, the_list_of_masses, latt_extent, the_sqred_lattice_momentum,the_unknown,the_ecm):
    energy_multi_hads_list=[]
    for i in range(len(the_flavor_sector)):
        the_row = SPLITTING_THE_OPERATORS_ROW(the_flavor_sector[i])
        the_hadrons_row = the_row[2]
        hadrons_and_momenta = HADRONS_INFO_AND_MOMENTUM(the_hadrons_row,the_sqred_lattice_momentum)
        the_energy = CALCULATING_FINAL_ENERGY(hadrons_and_momenta, the_list_of_masses, latt_extent, the_sqred_lattice_momentum,the_unknown,the_ecm) 
        if the_energy==0.: i=i+1
        else:
            if len(the_row)>3: energy_multi_hads_list.append([the_energy, the_row[1], the_row[2] + " " + the_row[3]])
            else: energy_multi_hads_list.append([the_energy, the_row[1], the_row[2]])
    new_multi_hads_list= sorted(energy_multi_hads_list, key=lambda k:[k[0], k[1], k[2]])
    return new_multi_hads_list




    # It prepares everything to be in a filtered list like: {energy, multiplicity, operators} (The cutoff happens here)
def ENERGY_LIST_TABLES(the_hadrons_energy_list,the_threshold,threeparticle):
    the_final_energy_list=[]
    for i in range(len(the_hadrons_energy_list)):
        if threeparticle==True:
            if float(the_hadrons_energy_list[i][0])>the_threshold:
                if float(the_hadrons_energy_list[i][0])<=the_threshold*1.05:
                    the_final_energy_list.append(the_hadrons_energy_list[i])
                    the_final_energy_list.append(the_hadrons_energy_list[i+1])
                    break
                else:break
            else: the_final_energy_list.append(the_hadrons_energy_list[i])
        else:
            if float(the_hadrons_energy_list[i][0])>the_threshold:break
            else: the_final_energy_list.append(the_hadrons_energy_list[i])
    return the_final_energy_list




###############################################################
#                                                             #
#       Extra function to reorganize, to prepare things       #
#       for plots, and some menu options                      # 
#                                                             #
###############################################################




     # It reorganizes the list of selected levels to plot them in increasing sqred momentum.
def LIST_FOR_PLOT(the_irrep,the_sqred_mom,the_energy_list, the_energy_list_plot):
    for zz in range(len(the_energy_list)):
        the_energy_list_plot.append([the_irrep,the_sqred_mom,the_energy_list[zz][0], the_energy_list[zz][2]])
    return the_energy_list_plot




    # It counts how many different irreps are there to plot.
def HOW_MANY_LEVELS_PLOT(an_energy_list):
    the_first_irrep=an_energy_list[0][0]
    how_many_levels=1
    for ii in range(len(an_energy_list)):
        if an_energy_list[ii][0]!=the_first_irrep: 
            how_many_levels+=1
            the_first_irrep=an_energy_list[ii][0]
        else: continue
    return how_many_levels



def MENU_QUANTUM_NUMBERS():
    print("Choose: (I: isospin, S: strangeness, B: Baryon nr.)")
    print("   [1] fermionic I=1/2 S=0 B=1 \n   [2] fermionic I=1/2 S=-2 B=1 \n   [3] fermionic I=3/2 S=0 B=1 \n   [4] fermionic I=3/2 S=-2 B=1 \n [5] fermionic I=5/2 S=0 B=1 \n   [6] fermionic I=0 S=-1 B=1 \n   [7]  fermionic I=0 S=-3 B=1 \n   [8] fermionic I=1 S=-1 B=1 \n   [9] fermionic I=1 S=-3 B=1 \n   [10] fermionic I=2 S=-1 B=1 \n    [11] bosonic I=0 S=0 B=0 \n    [12] bosonic I=0 S=2 B=0'] [13] bosonic I=1/2 S=1 B=0 \n   [14] bosonic I=3/2 S=1 B=0 \n   [15] bosonic I=5/2 S=1 B=0 \n   [16] bosonic I=1 S=0 B=0 \n   [17] bosonic I=1 S=2 B=0 \n   [18] bosonic I=2 S=0 B=0 \n   [19] bosonic I=2 S=2 B=0 \n    [20] bosonic I=3 S=0 B=0 \n    [21] All")
    the_quantum_number=int(input("The choice can be an integer or several separated by ',', e.g. 1,8..."))
    
    the_meson_baryon_list=[]
    if len(the_quantum_number)>1:
        the_quantum_number=the_quantum_number.split(",")
    for item in the_quantum_number:
        if item==1: the_meson_baryon_list.append(['fermionic', '2I1', 'S0', 'B1'])
        elif item==2: the_meson_baryon_list.append(['fermionic', '2I1', 'Sm2', 'B1'])
        elif item==3: the_meson_baryon_list.append(['fermionic', '2I3', 'S0', 'B1'])
        elif item==4: the_meson_baryon_list.append(['fermionic', '2I3', 'Sm2', 'B1'])
        elif item==5: the_meson_baryon_list.append(['fermionic', '2I5', 'S0', 'B1'])
        elif item==6: the_meson_baryon_list.append(['fermionic', 'I0', 'Sm1', 'B1'])
        elif item==7: the_meson_baryon_list.append(['fermionic', 'I0', 'Sm3', 'B1'])
        elif item==8: the_meson_baryon_list.append(['fermionic', 'I1', 'Sm1', 'B1'])
        elif item==9: the_meson_baryon_list.append(['fermionic', 'I1', 'Sm3', 'B1'])
        elif item==10: the_meson_baryon_list.append(['fermionic', 'I2', 'Sm1', 'B1']) 
        elif item==11: the_meson_baryon_list.append(['bosonic', 'I0', 'S0', 'B0'])
        elif item==12: the_meson_baryon_list.append(['bosonic', 'I0', 'S2', 'B0'])
        elif item==13: the_meson_baryon_list.append(['bosonic', '2I1', 'S1', 'B0'])
        elif item==14: the_meson_baryon_list.append(['bosonic', '2I3', 'S1', 'B0'])
        elif item==15: the_meson_baryon_list.append(['bosonic', '2I5', 'S1', 'B0'])
        elif item==16: the_meson_baryon_list.append(['bosonic', 'I1', 'S0', 'B0'])
        elif item==17: the_meson_baryon_list.append(['bosonic', 'I1', 'S2', 'B0'])
        elif item==18: the_meson_baryon_list.append(['bosonic', 'I2', 'S0', 'B0'])
        elif item==19: the_meson_baryon_list.append(['bosonic', 'I2', 'S2', 'B0'])
        elif item==20: the_meson_baryon_list.append(['bosonic', 'I3', 'S0', 'B0'])
        elif item==21: 
            the_meson_baryon_list=[
                ['fermionic', '2I1', 'S0', 'B1'],
                ['fermionic', '2I1', 'Sm2', 'B1'],
                ['fermionic', '2I3', 'S0', 'B1'],
                ['fermionic', '2I3', 'Sm2', 'B1'], 
                ['fermionic', '2I5', 'S0', 'B1'],
                ['fermionic', 'I0', 'Sm1', 'B1'],
                ['fermionic', 'I0', 'Sm3', 'B1'],
                ['fermionic', 'I1', 'Sm1', 'B1'],
                ['fermionic', 'I1', 'Sm3', 'B1'],
                ['fermionic', 'I2', 'Sm1', 'B1'],
                ['bosonic', 'I0', 'S0', 'B0'],
                ['bosonic', 'I0', 'S2', 'B0'],
                ['bosonic', '2I1', 'S1', 'B0'],
                ['bosonic', '2I3', 'S1', 'B0'],
                ['bosonic', '2I5', 'S1', 'B0'],
                ['bosonic', 'I1', 'S0', 'B0'],
                ['bosonic', 'I1', 'S2', 'B0'],
                ['bosonic', 'I2', 'S0', 'B0'],
                ['bosonic', 'I2', 'S2', 'B0'],
                ['bosonic', 'I3', 'S0', 'B0']];break
    return the_meson_baryon_list
            
            


    # It gives the selection of ensembles that are available to be computed.
def MENU_ENSEMBLES():
    print("For which ensemble(s) do you want to obtain the expected energy levels?")
    print("   [1] E250 \n   [2] B451 \n   [3] B452 \n   [4] B450 \n   [5] N452 \n   [6] N451 \n   [7] N200 \n   [8] N201 \n   [9] X252 \n   [10] X253 \n   [11] X451 \n   [12] D251 \n   [13] All")
    whichEnsemble=str(input('The choice can be an integer or op1,op2,...\n'))
    ensChoices=[]
    if "," in whichEnsemble:
        ensChoices_pre=whichEnsemble.split(",")
        for num in ensChoices_pre:
            ensChoices.append(ens.ensembleList[int(num)-1])
    elif whichEnsemble=="All" or "13" in whichEnsemble:
        ensChoices=ens.ensembleList
    elif int(whichEnsemble)<14:
        ensChoices.append(ens.ensembleList[int(whichEnsemble)-1])
    else:
        print("Incorrect choice.")
    return ensChoices



    # It gives the masses of the hadrons that are available in ALL ensembles, it excludes the ones that are not in every ensemble.
def MENU_HADRONS():
    print("Choose the hadron(s): \n   [1] Pion\n   [2] Kaon\n   [3] eta \n   [4] Nucleon\n   [5] eta-prime\n   [6] Lambda \n   [7] Sigma \n   [8] Xi \n   [9] Delta \n   [10] Sigma-star \n   [11] Xi-star\n   [12] Omega\n")
    choice=str(input("You can choose one or more hadrons, eg. 1,1,6: "))
    if len(choice)>1:
        choice=choice.split(",")
    choice=list(choice)
    the_hadron=[]
    for item in choice:
        item=int(item)
        if item==1:
            the_hadron.append(["pion_mass", r'$\pi$'])
        elif item==2:
            the_hadron.append(["kaon_mass", r'$K$'])
        elif item==3:
            the_hadron.append(["eta_mass",r'$\eta$'])
        elif item==4:
            the_hadron.append(["nucleon_mass", r'$N$'])
        elif item==5:
            the_hadron.append(["etaprime_mass",r"$\eta'$"])
        elif item==6:
            the_hadron.append(["lambda_mass", r'$\Lambda$'])
        elif item==7:
            the_hadron.append(["sigma_mass",r'$\Sigma$'])
        elif item==8:
            the_hadron.append(["xi_mass",r'$\Xi$'])
        elif item==9:
            the_hadron.append(["delta_mass",r'$\Delta$'])
        elif item==10:
            the_hadron.append(["sigmastar_mass",r'$\Sigma*$'])
        elif item==11:
            the_hadron.append(["xistar_mass",r'$\Xi*$'])
        elif item==12:
            the_hadron.append(["omega_mass",r'$\Omega$'])
    return the_hadron
           


    # It constructs a list of thresholds, which is also a list of hadrons to obtain the threhsolds over the nucleon mass.
def CHOICE_THRESHOLDS_PLOT():
    how_many_thresholds=int(input("How many thresholds do you want to plot?\n"))
    the_hadrons=[]
    for j in range(how_many_thresholds):
        print("For threshold nr. %s: "%(j+1))
        the_hadrons.append(MENU_HADRONS())
    return the_hadrons


    # it sums up all the hadron masses given (over the nucleon mass)
def SUMMING_HADRON_MASSES(the_threshold,the_ensemble):
    result=0.
    result_name=''
    for hadron in the_threshold:
        result+=the_ensemble[hadron[0]]/the_ensemble["nucleon_mass"]
        result_name+=hadron[1]
    return [result,result_name]


    # It actually sums up each hadron to get the threshold value.
def GETTING_THRESHOLDS(the_ensemble,the_thresholds_list):
    the_final_thresholds=[]
    for threshold in the_thresholds_list:
        the_final_thresholds.append(SUMMING_HADRON_MASSES(threshold,the_ensemble))
    return the_final_thresholds
    


def FINAL_LIST_OF_OPERATORS(the_hadrons_list):
    the_hadrons_and_momenta=[]
    for ii in range(len(the_hadrons_list)):
        the_hads_mom = HADRONS_INFO_AND_MOMENTUM(the_hadrons_list[ii][3],the_hadrons_list[ii][1])
        for item in the_hads_mom:
            if item in the_hadrons_and_momenta: continue
            else:the_hadrons_and_momenta.append(item)
    return the_hadrons_and_momenta


