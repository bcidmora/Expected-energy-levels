from pylatex import Document, LongTable, MultiColumn, Math, Table, Tabularx, Tabular, Section, Center, Alignat, Subsection, Subsubsection, NewPage
from pylatex.utils import NoEscape

import numpy as np

from energy_functions import *



#--------------  TEX TABLES HADRONS  ----------------

### It creates the TeX tables with 2 columns
def TABLE_TWOCOLS(the_doc, the_hadrons_list, the_headers, the_caption):
    with the_doc.create(Table(position = 'h!')) as main_table:
        with the_doc.create(Tabularx('l l', col_space='.75cm')) as table:            
            table.add_hline()
            table.add_row((the_headers[0], the_headers[1]))
            table.add_hline()
            table.add_hline()
            for i in range(len(the_hadrons_list)):
                table.add_row(str(the_hadrons_list[i][0]), the_hadrons_list[i][1])
            table.add_hline()
        main_table.add_caption(the_caption) 
        
### It creates the TeX tables with 3 columns
def TABLE_THREECOLS_ENERGIES(the_doc, the_hadrons_list, the_headers, the_caption):
    with the_doc.create(Table(position = 'h!')) as main_table:
        with the_doc.create(Tabularx('c c l', col_space='.75cm')) as table:            
            table.add_hline()
            table.add_row((the_headers[0], the_headers[1], the_headers[2]))
            table.add_hline()
            table.add_hline()
            for i in range(len(the_hadrons_list)):
                table.add_row(str(the_hadrons_list[i][0]), the_hadrons_list[i][1], the_hadrons_list[i][2])
            table.add_hline()
        main_table.add_caption(the_caption) 
        

### It puts the tables in the document. 
def CONSTRUCTING_TABLES(the_doc, the_hadrons_list, the_headers, the_caption):
    if len(the_hadrons_list)<=45:
        TABLE_THREECOLS_ENERGIES(the_doc,the_hadrons_list,the_headers, the_caption)
        the_doc.append(NewPage())
    elif len(the_hadrons_list)>45:
        for ii in range(0,len(the_hadrons_list),45):
            if ii+45<len(the_hadrons_list):
                TABLE_THREECOLS_ENERGIES(the_doc,the_hadrons_list[ii:ii+45],the_headers, the_caption)
                the_doc.append(NewPage())
            else:
                TABLE_THREECOLS_ENERGIES(the_doc,the_hadrons_list[ii:],the_headers, the_caption)
                the_doc.append(NewPage())





#--------------  TEX TABLES HADRONS AND ENSEMBLES INFO  ----------------


####################################################
#                                                  #
#                                                  #
#  These are tables summarizing all the info of a  #
#  ensemble, the hadron masses and other consts.   #
#                                                  #
#                                                  #
####################################################

# This routine is pretty much the same than the ones above, but includes the JPC numbers (the must be included in the the_had_data as the last entry for each item.)

# the_a: lattice spacing in fm.
# the_N_lattice: lattice extent L.
# the_T_lattice: lattice extent T.
# the_had_data: list of lists, [name of hadron, mass of hadron ]
# the_caption: whatever one wants to put in the caption of this table.
def TABLE_FOUR_COLS_JPC_ENSEMBLES(the_doc, the_cls_name, the_a, the_N_lattice, the_T_lattice, the_beta_vals, the_kappa_u, the_kappa_s, the_had_data, the_nr_configs, the_caption):
     with the_doc.create(Table(position = 'h!')) as main_table:
        with the_doc.create(Tabular('|lccc|', col_space='1.1cm', booktabs=True)) as table:
            table.add_row(('Properties of the Lattice', '','',  'Values'))
            table.add_hline()
            table.add_row(('Lattice size', '', '', the_a))
            table.add_row((NoEscape(r'Lattice extent $L^{3}$'),'',  '',  NoEscape(f'${the_N_lattice}^{3}$')))
            table.add_row((NoEscape(r'Lattice extent $T$'), '', '', NoEscape(f'${the_T_lattice}$' )))
            table.add_row((NoEscape(r'$\beta$'), '', '',  the_beta_vals))
            table.add_row((NoEscape(r'$\kappa_{u}$'), '', '', the_kappa_u))
            table.add_row((NoEscape(r'$\kappa_{s}$'), '', '', the_kappa_s))
            table.add_row(('CLS name', '', '', the_cls_name))
            table.add_row(('Nr. Gauge Configs', '', '', the_nr_configs))
            table.add_hline()
            table.add_hline()
            table.add_row(('Hadrons', NoEscape(r'$J^{P}$'), 'Masses [am]', NoEscape('Masses [MeV]')))
            table.add_hline()
            for i in range(len(the_had_data)):
                table.add_row(the_had_data[i][0], the_had_data[i][3],  str(np.round(float(the_had_data[i][1]), 3)), str(np.round(float(the_had_data[i][2]), 1)))
            table.add_hline()
        main_table.add_caption(the_caption)
        
            

def TABLE_TWOCOLS_OPERATORS(the_doc, the_hadrons_list, the_headers, the_caption):
    with the_doc.create(Table( )) as main_table:
        with the_doc.create(Tabularx('c c')) as table:            
            table.add_hline()
            table.add_row((the_headers[0], the_headers[1]))
            table.add_hline()
            table.add_hline()
            for i in range(len(the_hadrons_list)):
                table.add_row(str(the_hadrons_list[i][0]), the_hadrons_list[i][1])
            table.add_hline()
        main_table.add_caption(the_caption) 
