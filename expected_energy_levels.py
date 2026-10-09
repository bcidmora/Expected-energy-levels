from pylatex import Document, LongTable, MultiColumn, Math, Table, Tabularx, Tabular, Section, Center, Alignat, Subsection, Subsubsection, NewPage, Command, Figure, Hyperref, Package
from pylatex.utils import italic, NoEscape

import os
import numpy as np

import energy_functions as enf
import tables_latex as tl
import ensembles_info as ens
import labels_and_plots as lap
import list_of_lists as lol


    ### Directory where the ExpectedEnergyLevels produced by Colins are
energyLevelsLocation = os.path.expanduser('~')+'${YOUR_PATH}/96x96x96_phys/'

unkwnownHads=False
noEnergyCut=False
threeParticleThreshold=False
cmEnergy=True
plotEnergyLevels=True

ensChoices=lap.MENU_ENSEMBLES()

    ### Info about the final state of interest
allMesonBaryonLevels=lap.MENU_QUANTUM_NUMBERS()

momCombinations = ['000', '001', '002', '003', '011', '012', '022', '111', '112', '122'] # All momentum combinations

cutChoice=int(input("[1] Three Particle Threshold \n[2] Energy value \n[3] No cut (includes all levels)\n"))

if cutChoice==1:
    threeParticleThreshold=True
    print("Now choosing the three particles...")
    chosenThresholdList=lap.MENU_HADRONS()
elif cutChoice==2:
    threeParticleThreshold=False
    chosenThreshold=1.71
    newChosenThreshold=input(f"Default cutoff is: {chosenThreshold} \nEnter a new value if you want to change it: ")
    try:
        chosenThreshold=float(newChosenThreshold)
    except ValueError:
        chosenThreshold=float(chosenThreshold)
elif cutChoice==3:
    noEnergyCut=True
    threeParticleThreshold=False
    chosenThreshold=10
else:
    print('Invalid option. Exiting.')
    exit()


if cmEnergy:
    energyHeader=r'E$_{\mathrm{cm}}/m_{N}$'
else:
    energyHeader=r'E$/m_{N}$'
    

if plotEnergyLevels:
    print("Choosing thresholds to be in the plot...")
    refHadrons=lap.CHOICE_THRESHOLDS_PLOT()


### --------------------------------------------------------------------

### This is the beginning of the script
if __name__=='__main__':
    
    for each_hadron in allMesonBaryonLevels: 
        hadronType = each_hadron.theParticleType # Final state: bosonic or fermionic
        totalIsospin = each_hadron.theIsospin  # Isospin value
        totalStrangeness =  each_hadron.theStrangeness # Strangeness 
        numberBaryons =  each_hadron.theBaryonNumber # Number of Baryons

        levelsPlot=[]
        
        for item in range(len(ensChoices)):
            hadronsPlot=[]
            finalListHadrons=[]
        
            ### This is about the document
            geometry_options = {"tmargin": "1cm", "lmargin": "1.8cm", "rmargin": "1.8cm"}
            doc = Document(page_numbers=True, geometry_options=geometry_options) 
            
            ### This is the setup for the table of contents
            doc.append(NoEscape(r'\setcounter{tocdepth}{1}'))
            

            ### Setting up the title page
            doc.preamble.append(Command('title', 'Expected Energy Levels: Colin Morningstar Group Theory input'))
            doc.preamble.append(Command('author', 'B. Cid-Mora'))
            doc.preamble.append(Command('date', NoEscape(r'\today')))
            
            doc.append(NoEscape(r'\maketitle'))
            doc.append(NoEscape(f'Hadron: {hadronType} {totalIsospin} {totalStrangeness} {numberBaryons}'))
            
            doc.append(NoEscape(r'\tableofcontents'))
            
            doc.append(NewPage())
            with doc.create(Section('Some introduction')):   
                doc.append(NoEscape(r'\label{sec:intro}'))
                
                doc.append(NoEscape(rf'This document is intended to obtain the expected energy levels for a few ensembles. Below some tables show the momentum combinations to study the $\Lambda(1405)$. The following relations are used:'))
                doc.append('\n\n1. Relativistic relation:')
                with doc.create(Alignat(numbering=True, escape=False)) as eqn:
                    eqn.append(r'E^{2} (\mathbf{d})=  m_{H}^{2} + \left(\frac{2 \pi |\Vec{\mathbf{d}}|}{L}\right)^{2},')
                
                doc.append(NoEscape(r'where $L:$ lattice extent, and $\mathbf{d}:$ units of momentum in the lattice.'))
                
                doc.append('\n\n2. Then the energy is transformed into MeV as follows:')
                with doc.create(Alignat(numbering=True, escape=False)) as eqn:
                    eqn.append(r'm_{\rm had} [\rm MeV] = m_{\rm latt} \cdot \frac{%s}{a_{\rm latt}}'%str(enf.fmToMeV))
                    
                doc.append('Then one obtains the expected energy levels over the nucleon mass \n')
                    
                with doc.create(Alignat(numbering=True, escape=False)) as eqn:
                    eqn.append(r'\frac{E_{\rm expected}( \mathbf{d})}{m_{N}} = \frac{1}{m_{N}}\cdot \sum_{i}^{\rm nr. hads} \sqrt{ m_{i}^{2} + \left( \frac{2\pi \mathbf{d}_{i}}{L} \right)^{2} }')
                
                doc.append('\nThe expected energy level in this case matches the Elab, which can be shifted as')
                    
                with doc.create(Alignat(numbering=False, escape=False)) as eqn:
                    eqn.append(r'E_{\rm cm} =  \sqrt{ E_{\rm lab}^{2} - \left( \frac{2\pi \mathbf{d}_{\rm tot}}{L}  \right)^{2} }')
            
            doc.append(NewPage())
            with doc.create(Section(f'Ensemble {ensChoices[item]["ens_name"]}')):        
                doc.append('LATTICE PROPERTIES AND OTHERS\n')
                
                hadronsInfo = []
                
                for key, hadron in lol.the_list_hadrons.items():
                    
                    try:
                        the_mass_lattice = ensChoices[item][hadron.mass]
                        the_mass_MeV = enf.FM_TO_MEV(the_mass_lattice, ensChoices[item]['latt_spacing'])
                    except KeyError:
                        print(f'{key} is not included in {ensChoices[item]["ens_name"]} ensemble')
                    
                    hadronsInfo.append([key, str(the_mass_lattice), str(the_mass_MeV),  NoEscape((hadron.spinParity))])
                
                newHadronsInfo=sorted(hadronsInfo, key=lambda k:[k[1], k[0], k[2],k[3]])
                
                tl.TABLE_FOUR_COLS_JPC_ENSEMBLES(doc, ensChoices[item]["cls_name"], ensChoices[item]["latt_spacing"], ensChoices[item]["s_extent"], ensChoices[item]["t_extent"], ensChoices[item]["beta"], ensChoices[item]["ud_kappa"], ensChoices[item]["s_kappa"], newHadronsInfo, ensChoices[item]["cnfgs"], 'Hadrons info and lattice size (a) and extent (L).')
                
                doc.append(NewPage())
                
                for mom in momCombinations:
                    momentumFileDir = f'{energyLevelsLocation}mom_{mom}/{hadronType}_{totalIsospin}_{totalStrangeness}_{numberBaryons}_levels.txt'
                    
                    if not os.path.isfile(momentumFileDir): 
                        continue
                    else:
                        energyLevelsFile = lap.READ_TXT_FILE(momentumFileDir)
                        doc.append(NewPage())
                        with doc.create(Subsection(str(energyLevelsFile[5]))):
                            
                            psqMom=int(mom[0])**2 + int(mom[1])**2 + int(mom[2])**2
                            
                            if threeParticleThreshold:
                                chosenThreshold=enf.SUMMING_HADRON_MASSES(chosenThresholdList, ensChoices[item])[0]
                            
                            doc.append("Threshold: ")
                            doc.append(NoEscape(rf'$E = {chosenThreshold}$'))
                            
                            indexFlavorInData = [i for i, x in enumerate(energyLevelsFile) if "Flavor" in x]
                            
                            tabHeaders = [NoEscape(energyHeader), 'Degeneracy', 'Operators']
                            
                            for jj in range(len(indexFlavorInData)):
                                
                                tabCaption = f'({ensChoices[item]["ens_name"]}) {energyLevelsFile[indexFlavorInData[jj]]}. [d {energyLevelsFile[5][23:]}]'
                                
                                posIrrepStr=int(str(energyLevelsFile[indexFlavorInData[jj]]).index("Irrep"))+7
                                disIrrep = str(energyLevelsFile[indexFlavorInData[jj]])[posIrrepStr: -1]
                                
                                if jj!=(len(indexFlavorInData)-1):  
                                     final_hadrons_list=enf.ENERGY_LIST_RAW(energyLevelsFile[int(indexFlavorInData[jj])+2:int(indexFlavorInData[jj+1])-2], newHadronsInfo, ensChoices[item]["s_extent"], psqMom, unkwnownHads, cmEnergy)
                                     
                                else:
                                    final_hadrons_list=enf.ENERGY_LIST_RAW(energyLevelsFile[int(indexFlavorInData[jj])+2:], newHadronsInfo, ensChoices[item]["s_extent"], psqMom, unkwnownHads,cmEnergy)
                                
                                if not noEnergyCut:
                                    final_energy_list=enf.ENERGY_LIST_TABLES(final_hadrons_list,chosenThreshold,threeParticleThreshold)
                                else: 
                                    final_energy_list=final_hadrons_list 
                                tl.CONSTRUCTING_TABLES(doc,final_energy_list,tabHeaders, tabCaption)
                                
                                levelsPlot=lap.LIST_FOR_PLOT(disIrrep,str(psqMom),final_energy_list,hadronsPlot)
                                
                                irrepHadronsList = enf.FINAL_LIST_OF_OPERATORS(final_energy_list,str(psqMom))
                                for ll in range(len(irrepHadronsList)):
                                    if irrepHadronsList[ll] not in finalListHadrons: finalListHadrons.append(irrepHadronsList[ll])
                    
                    doc.append(NewPage())
            
            
            for xx in range(len(finalListHadrons)):
                finalListHadrons[xx] = [finalListHadrons[xx], enf.POSSIBLE_MOMENTUM_HADRONS(finalListHadrons[xx])]
                
            if plotEnergyLevels: 
                
                namePlot=f"EnergyLevels_plot_{hadronType}_{totalIsospin}_{totalStrangeness}_{numberBaryons}_{ensChoices[item]["ens_name"]}.pdf"
                
                newLevelsPlot=list(sorted(levelsPlot, key=lambda k: [ k[1], k[0], k[2], k[3]]))

                howManyLevels=lap.HOW_MANY_LEVELS_PLOT(newLevelsPlot)

                refHadronLevels=enf.GETTING_THRESHOLDS(ensChoices[item],refHadrons)
                
                the_plot_size=10*howManyLevels+260
                
                lap.PLOT_ENERGY_LEVELS(newLevelsPlot, refHadronLevels, energyHeader, namePlot, howManyLevels)
                doc.append(NewPage())
                
                with doc.create(Subsection('Summary Energy Levels')):
                    plot = os.path.join(os.path.dirname(__file__), namePlot)
                    with doc.create(Figure(position='h!')) as energies_plot:
                        energies_plot.add_image(plot, width=f'{the_plot_size}px')
                        energies_plot.add_caption(NoEscape(rf'Summary of expected energy levels ({ensChoices[item]["ens_name"]}). Orange lines are relevant thresholds. Irreps are in the $x$-axis, and momentum in the form $P_{{tot}}^{{2}}=P_{{x}}^{{2}} + P_{{y}}^{{2}} + P_{{z}}^{{2}}$ is in parenthesis. '))
            
            tl.TABLE_TWOCOLS(doc, list(sorted(list(sorted(finalListHadrons, key=lambda k: [k[0][k[0].index('(')+1: k[0].index(')')], k[1] ])), key=lambda k: [k[0], k[1]])), ['Hadron', 'Momentum'], 'All needed hadrons with their corresponding momentum.')
        
            doc.generate_pdf(f"Operators_{hadronType}_{totalIsospin}_{totalStrangeness}_{numberBaryons}_{ensChoices[item]["ens_name"]}", clean_tex=True)
