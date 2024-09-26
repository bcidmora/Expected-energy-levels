from pylatex import Document, LongTable, MultiColumn, Math, Table, Tabularx, Tabular, Section, Center, Alignat, Subsection, Subsubsection, NewPage, Command, Figure, Hyperref, Package
from pylatex.utils import italic, NoEscape
import os
import numpy as np
import energy_functions as enf
import tables_latex as tl
import ensembles_info as ens


    # Directory where the ExpectedEnergyLevels produced by Colins are
energyLevelsLocation = os.path.expanduser('~')+'/ExpectedEnergyLevels/96x96x96_phys/'


unkwnownHads=False
noEnergyCut=False
threeParticleThreshold=False
cmEnergy=False
plotEnergyLevels=True

ensChoices=enf.MENU_ENSEMBLES()

    # Info about the final state of interest
f_2I1_S0_B1 = ['fermionic', '2I1', 'S0', 'B1'] 
f_2I1_Sm2_B1 = ['fermionic', '2I1', 'Sm2', 'B1'] 
f_2I3_S0_B1 = ['fermionic', '2I3', 'S0', 'B1'] 
f_2I3_Sm2_B1 = ['fermionic', '2I3', 'Sm2', 'B1'] 
f_2I5_S0_B1 = ['fermionic', '2I5', 'S0', 'B1'] 
f_I0_Sm1_B1 = ['fermionic', 'I0', 'Sm1', 'B1'] 
f_I0_Sm3_B1 = ['fermionic', 'I0', 'Sm3', 'B1'] 
f_I1_Sm1_B1 = ['fermionic', 'I1', 'Sm1', 'B1'] 
f_I1_Sm3_B1 = ['fermionic', 'I1', 'Sm3', 'B1'] 
f_I2_Sm1_B1 = ['fermionic', 'I2', 'Sm1', 'B1'] 

b_I0_S0_B0 = ['bosonic', 'I0', 'S0', 'B0'] 
b_I0_S2_B0 = ['bosonic', 'I0', 'S2', 'B0'] 
b_2I1_S1_B0 = ['bosonic', '2I1', 'S1', 'B0'] 
b_2I3_S1_B0 = ['bosonic', '2I3', 'S1', 'B0'] 
b_2I5_S1_B0 = ['bosonic', '2I5', 'S1', 'B0']  
b_I1_S0_B0 = ['bosonic', 'I1', 'S0', 'B0'] 
b_I1_S2_B0 = ['bosonic', 'I1', 'S2', 'B0']  
b_I2_S0_B0 = ['bosonic', 'I2', 'S0', 'B0'] 
b_I2_S2_B0 = ['bosonic', 'I2', 'S2', 'B0'] 
b_I3_S0_B0 = ['bosonic', 'I3', 'S0', 'B0'] 


#allMesonBaryonLevels = [f_2I1_S0, f_2I1_Sm2, f_2I3_S0, f_2I3_Sm2, f_2I5_S0, f_I0_Sm1, f_I0_Sm3, f_I1_Sm1, f_I1_Sm3, f_I2_Sm1]

#allBosonicLevels=[b_I0_S0_B0, b_I0_S2_B0 ,b_2I1_S1_B0, b_2I3_S1_B0, b_2I5_S1_B0, b_I1_S0_B0, b_I1_S2_B0, b_I2_S0_B0, b_I2_S2_B0, b_I3_S0_B0]

#allMesonBaryonLevels=MENU_QUANTUM_NUMBERS()

allMesonBaryonLevels=[f_I0_Sm1_B1]

#momCombinations = ['000', '001', '002', '003', '011', '012', '022', '111', '112', '122'] # All momentum combinations
momCombinations = ['000', '001', '002', '003', '011', '022', '111'] 


cutChoice=int(input("[1] Three Particle Threshold \n[2] Energy value \n[3] No cut (includes all levels)\n"))


if cutChoice==1:
    threeParticleThreshold=True
    print("Now choosing the three particles...")
    chosenThresholdList=enf.MENU_HADRONS()
elif cutChoice==2:
    threeParticleThreshold=False
    chosenThreshold=1.68
    newChosenThreshold=input("Default cutoff is: %s \nEnter a new value if you want to change it: "%str(chosenThreshold))
    try:
        chosenThreshold=float(newChosenThreshold)
    except ValueError:
        chosenThreshold=float(chosenThreshold)
    
elif cutChoice==3:
    noEnergyCut=True
    threeParticleThreshold=False
    chosenThreshold=10


if cmEnergy:
    energyHeader=r'$E_{cm}/m_{N}$'
else:
    energyHeader=r'$E/m_{N}$'
    

if plotEnergyLevels:
    print("Choosing thresholds to be in the plot...")
    refHadrons=enf.CHOICE_THRESHOLDS_PLOT()

## --------------------------------------------------------------------




# This is the beginning of the script
if __name__=='__main__':
    
    for each_hadron in allMesonBaryonLevels: 
        hadronType = each_hadron[0] # Final state: bosonic or fermionic
        totalIsospin = each_hadron[1]  # Isospin value
        totalStrangeness =  each_hadron[2] # Strangeness 
        numberBaryons =  each_hadron[3] # Number of Baryons
        
        levelsPlot=[]

        # This is about the document
        geometry_options = {"tmargin": "1cm", "lmargin": "1.8cm", "rmargin": "1.8cm"}
        doc = Document(page_numbers=True, geometry_options=geometry_options) 
        
        # This is the setup for the table of contents
        doc.append(NoEscape(r'\setcounter{tocdepth}{1}'))
        

        # Setting up the title page
        doc.preamble.append(Command('title', 'Expected Energy Levels: Colin Morningstar Group Theory input'))
        doc.preamble.append(Command('author', 'B. Cid-Mora'))
        doc.preamble.append(Command('date', NoEscape(r'\today')))
        
        doc.append(NoEscape(r'\maketitle'))
        doc.append(NoEscape(r'{\Large{Hadron: %s}}'%(hadronType + ' ' + totalIsospin + ' '+ totalStrangeness + ' ' +numberBaryons)))
        
        doc.append(NoEscape(r'\tableofcontents'))
        
        doc.append(NewPage())
        with doc.create(Section('Some introduction')):   
            doc.append(NoEscape(r'\label{sec:intro}'))
            
            doc.append(NoEscape('This document is intended to obtain the expected energy levels for a few ensembles. Below some tables show the momentum combinations to study the $\Lambda(1405)$. The following relations are used:\n'))
            doc.append('\n1. Relativistic relation:\n')
            with doc.create(Alignat(numbering=True, escape=False)) as eqn:
                eqn.append(r'E^{2} (\mathbf{d})=  m_{H}^{2} + \left(\frac{2 \pi |\Vec{\mathbf{d}}|}{L}\right)^{2},')
            
            doc.append(NoEscape('where $L:$ lattice extent, and $\mathbf{d}:$ units of momentum in the lattice.\n'))
            
            doc.append('\n2. Then the energy is transformed into MeV as follows:\n')
            with doc.create(Alignat(numbering=True, escape=False)) as eqn:
                eqn.append(r'm_{\rm had} [\rm MeV] = m_{\rm latt} \cdot \frac{%s}{a_{\rm latt}}'%str(enf.fmToMeV))
            
            doc.append('Then one obtains the expected energy levels over the nucleon mass \n')
                
            with doc.create(Alignat(numbering=True, escape=False)) as eqn:
                eqn.append(r'\frac{E_{\rm expected}( \mathbf{d})}{m_{N}} = \frac{1}{m_{N}}\cdot \sum_{i}^{\rm nr. hads} \sqrt{ m_{i}^{2} + \left( \frac{2\pi \mathbf{d}_{i}}{L} \right)^{2} }')
                
            with doc.create(Alignat(numbering=False, escape=False)) as eqn:
                eqn.append(r'E_{\rm cm} =  \sqrt{ E_{\rm lab}^{2} - \left( \frac{2\pi \mathbf{d}_{\rm tot}}{L}  \right)^{2} }')
        
            
        
        for item in range(len(ensChoices)):
            
            hadronsPlot=[]
            
            doc.append(NewPage())
            with doc.create(Section('Ensemble %s'%ensChoices[item]["ens_name"])):        
                doc.append('LATTICE PROPERTIES AND OTHERS\n')

                hadronsInfo = [['pi', '%s'%str(ensChoices[item]["pion_mass"]), '%s'%str(enf.FM_TO_MEV(ensChoices[item]["pion_mass"], ensChoices[item]["latt_spacing"])), NoEscape('$0^{-}$')],
                            
                            ['K','%s'%str(ensChoices[item]["kaon_mass"]), '%s'%str(enf.FM_TO_MEV(ensChoices[item]["kaon_mass"], ensChoices[item]["latt_spacing"])), NoEscape('$0^{-}$') ],
                            
                            ['KB','%s'%str(ensChoices[item]["kaon_mass"]), '%s'%str(enf.FM_TO_MEV(ensChoices[item]["kaon_mass"], ensChoices[item]["latt_spacing"])), NoEscape('$0^{-}$') ],
                            
                            ['eta','%s'%str(ensChoices[item]["eta_mass"]), '%s'%str(enf.FM_TO_MEV(ensChoices[item]["eta_mass"], ensChoices[item]["latt_spacing"])), NoEscape('$0^{-}$') ],
                            
                            ['N', '%s'%str(ensChoices[item]["nucleon_mass"]), '%s'%str(enf.FM_TO_MEV(ensChoices[item]["nucleon_mass"], ensChoices[item]["latt_spacing"])), NoEscape('$1/2^{+}$') ],
                            
                            ['etaprime-958','%s'%str(ensChoices[item]["etaprime_mass"]), '%s'%str(enf.FM_TO_MEV(ensChoices[item]["etaprime_mass"], ensChoices[item]["latt_spacing"])), NoEscape('$0^{-}$') ],
                            
                            ['Lambda', '%s '%str(ensChoices[item]["lambda_mass"]), '%s'%str(enf.FM_TO_MEV(ensChoices[item]["lambda_mass"], ensChoices[item]["latt_spacing"])), NoEscape('$1/2^{+}$') ],
                            
                            ['Sigma','%s'%str(ensChoices[item]["sigma_mass"]), '%s'%str(enf.FM_TO_MEV(ensChoices[item]["sigma_mass"], ensChoices[item]["latt_spacing"])), NoEscape('$1/2^{+}$')],
                            
                            ['Delta-1232', '%s'%str(ensChoices[item]["delta_mass"]), '%s'%str(enf.FM_TO_MEV(ensChoices[item]["delta_mass"], ensChoices[item]["latt_spacing"])), NoEscape('$3/2^{+}$')],
                            
                            ['Xi', '%s'%str(ensChoices[item]["xi_mass"]), '%s'%str(enf.FM_TO_MEV(ensChoices[item]["xi_mass"], ensChoices[item]["latt_spacing"])), NoEscape('$1/2^{+}$')],
                            
                            ['Sigma-1385', '%s'%str(ensChoices[item]["sigmastar_mass"]), '%s'%str(enf.FM_TO_MEV(ensChoices[item]["sigmastar_mass"], ensChoices[item]["latt_spacing"])), NoEscape('$3/2^{+}$')],
                            
                            ['Xi-1530', '%s'%str(ensChoices[item]["xistar_mass"]), '%s'%str(enf.FM_TO_MEV(ensChoices[item]["xistar_mass"], ensChoices[item]["latt_spacing"])), NoEscape('$3/2^{+}$')],
                            
                            ['Omega', '%s'%str(ensChoices[item]["omega_mass"]), '%s'%str(enf.FM_TO_MEV(ensChoices[item]["omega_mass"], ensChoices[item]["latt_spacing"])), NoEscape('$3/2^{-}$')]
                            ]
                try:
                    hadronsInfo.append(['rho-770', '%s'%str(ensChoices[item]["rho_mass"]), '%s'%str(enf.FM_TO_MEV(ensChoices[item]["rho_mass"], ensChoices[item]["latt_spacing"])), NoEscape('$1^{-}$') ])
                except KeyError:
                    print('rho-770 is not included in %s ensemble'%str(ensChoices[item]["ens_name"]))
                    
                try:
                    hadronsInfo.append(['Lambda-1405', '%s'%str(ensChoices[item]["lambda1405_mass"]), '%s'%str(enf.FM_TO_MEV(ensChoices[item]["lambda1405_mass"], ensChoices[item]["latt_spacing"])), NoEscape('$1/2^{-}$')])
                except KeyError:
                    print('Lambda(1405) is not included in %s ensemble'%str(ensChoices[item]["ens_name"]))
                
                newHadronsInfo=sorted(hadronsInfo, key=lambda k:[k[1], k[0], k[2],k[3]])
                
                tl.TABLE_FOUR_COLS_JPC_ENSEMBLES(doc, ensChoices[item]["cls_name"], ensChoices[item]["latt_spacing"], ensChoices[item]["s_extent"], ensChoices[item]["t_extent"], ensChoices[item]["beta"], ensChoices[item]["ud_kappa"], ensChoices[item]["s_kappa"], newHadronsInfo, ensChoices[item]["cnfgs"], 'Hadrons info and lattice size (a) and extent (L).')
                
                doc.append(NewPage())
                
                for mom in momCombinations:
                    momentumFileDir = energyLevelsLocation + 'mom_%s/'%mom + hadronType + '_' + totalIsospin + '_' + totalStrangeness + '_' + numberBaryons + '_levels.txt'
                    
                    if not os.path.isfile(momentumFileDir): 
                        continue
                    else:
                        energyLevelsFile = enf.READ_TXT_FILE(momentumFileDir)
                        doc.append(NewPage())
                        with doc.create(Subsection(str(energyLevelsFile[5]))):
                            
                            psqMom=int(mom[0])**2 + int(mom[1])**2 + int(mom[2])**2
                            
                            if threeParticleThreshold:
                                chosenThreshold=enf.SUMMING_HADRON_MASSES(chosenThresholdList,ensChoices[item])[0]
                            
                            doc.append("Threshold: ")
                            doc.append(NoEscape(r'$E = %s$'%str(chosenThreshold)))
                            
                            indexFlavorInData = [i for i, x in enumerate(energyLevelsFile) if "Flavor" in x]
                            
                            tabHeaders = [NoEscape(energyHeader), 'Degeneracy', 'Operators']
                            
                            for jj in range(len(indexFlavorInData)):
                                tabCaption = '(%s) '%str(ensChoices[item]["ens_name"]) + str(energyLevelsFile[indexFlavorInData[jj]]) + ('. [d %s]')%str(energyLevelsFile[5][23:])
                                
                                posIrrepStr=int(str(energyLevelsFile[indexFlavorInData[jj]]).index("Irrep"))+7
                                disIrrep = str(energyLevelsFile[indexFlavorInData[jj]])[posIrrepStr: -1]
                                
                                if jj!=(len(indexFlavorInData)-1):  
                                     final_hadrons_list=enf.ENERGY_LIST_RAW(energyLevelsFile[int(indexFlavorInData[jj])+2:int(indexFlavorInData[jj+1])-2], newHadronsInfo,ensChoices[item]["s_extent"], psqMom, unkwnownHads,cmEnergy)
                                     
                                else:
                                    final_hadrons_list=enf.ENERGY_LIST_RAW(energyLevelsFile[int(indexFlavorInData[jj])+2:], newHadronsInfo,ensChoices[item]["s_extent"], psqMom, unkwnownHads,cmEnergy)
                                
                                if not noEnergyCut:
                                    final_energy_list=enf.ENERGY_LIST_TABLES(final_hadrons_list,chosenThreshold,threeParticleThreshold)
                                    if len(final_energy_list)==0:break
                                else: 
                                    final_energy_list=final_hadrons_list 
                                tl.CONSTRUCTING_TABLES(doc,final_energy_list,tabHeaders, tabCaption)
                                
                                levelsPlot=enf.LIST_FOR_PLOT(disIrrep,str(psqMom),final_energy_list,hadronsPlot)

                    doc.append(NewPage())
            
            
            if plotEnergyLevels: 
                namePlot="EnergyLevels_plot_"+hadronType+"_"+totalIsospin+"_"+totalStrangeness+"_"+numberBaryons+"_%s.pdf"%ensChoices[item]["ens_name"]
                
                newLevelsPlot=list(sorted(levelsPlot, key=lambda k: [ k[1], k[0], k[2], k[3]]))
                
                howManyLevels=enf.HOW_MANY_LEVELS_PLOT(newLevelsPlot)

                refHadronLevels=enf.GETTING_THRESHOLDS(ensChoices[item],refHadrons)
                
                tl.PLOT_ENERGY_LEVELS(newLevelsPlot, refHadronLevels,energyHeader,namePlot,howManyLevels)
                doc.append(NewPage())
                
                with doc.create(Subsection('Summary Energy Levels')):
                    plot = os.path.join(os.path.dirname(__file__), namePlot)
                    with doc.create(Figure(position='h!')) as energies_plot:
                        energies_plot.add_image(plot, width='400px')
                        energies_plot.add_caption(NoEscape('Summary of expected energy levels (%s). Orange lines are relevant thresholds. Irreps are in the x-axis, and momentum in the form $P_{tot}^{2}=P_{x}^{2} + P_{y}^{2} + P_{z}^{2}$ is in parenthesis. '%ensChoices[item]["ens_name"]))
                
        doc.generate_pdf("Operators_"+hadronType+"_"+totalIsospin+"_"+totalStrangeness+"_"+numberBaryons, clean_tex=True)
