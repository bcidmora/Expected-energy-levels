# Expected energy levels

To run this script, the following libraries are needed:

    + pylatex
    + numpy
    + os
    + matplotlib
    + dataclasses
    + Colin's expected energy levels for the E250
    
These scripts are organized as:

** expected_energy_levels.py ** 
    This is the main script to run. The directory at the top must be changed. It asks for which ensembles you want to get the expected energy levels, the quantum numbers of interest, and nearby thresholds of interest for plotting. It returns a pdf file per ensemble per quantum number with the info of the ensemble, the hadron masses in lattice and in physical units and at the end, a plot with the energy levels and a summary of all the single hadrons required to compute. 
    
    - unkwnownHads: if True, the unknown masses for a certain ensemble will be include in the calculation of energy levels using the E250 hadron masses. If False, then the energy levels that include these hadrons are skipped and not included in the tables.
    
    - noEnergyCut: If False, the full list of energy levels is reduced to those that are within the chosen cutoff, otherwise (True) it includes all possible levels.
    
    - threeParticleThreshold: if True, three particles must be chosen and then it cuts the energy two levels after this threshold. 
    
    - cmEnergy: if True the final energy levels are in terms of the E_{cm}. If False, all the levels are given in the Elab frame.
    
    - plotEnergyLevels: If True, all the selected energy levels are plotted at the end of the section. Else, no plot is generated.
    
    - ensChoices: Ensembles to be studied. More ensembles can be added to the dictionary in "ensembles_info.py"
    
    - chosenThreshold: The cutoff value (three particles, a specific value or no cut)
    
    - hadronsInfo: This list contains all the hadrons that are known in all the ensembles in "ensembles_info.py". More hadron values can be added, but their mass values must be included for the ensemble in "ensembles_info.py". 

    - momCombinations: These are all the momentum combinations of interest. One can focus only on certain frames. 

Later, the process is the following: 

        * All energy values are computed using Colin's group theory (ENERGY_LIST_RAW). 
        * Then, a selection of levels is made according to the chosen cutoff (ENERGY_LIST_TABLES). 
        * Finally, all tables are constructed according to this selection, if the tables are too long, they are separated in several pages (CONSTRUCTING_TABLES). 
        * When plotEnergyLevels==True, the plots can include thresholds.    

** tables_latex.py **
    Builds all the tables that appear in the pdf file. These tables are written in a TeX file and compiled. 

** labels_and_plots.py **
    It has all the routines displaying the Menu for hadrons, ensembles, etc. It also has the plot routines.

** list_of_lists.py **
    It contains all the lists of hadrons, irreps and quantum numbers as lists. 

** energy_functions.py **
    It has all the routines that calculate stuff. It gets the energies for the new ensemble following Colin's group theory.  
    
** ensembles_info.py **
    This file contains all the information about the ensembles one wants to get the energy levels from. Please add more hadron masses in this file if you have more values, and share them with the group.
    
