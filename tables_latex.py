import numpy as np
from pylatex import Document, LongTable, MultiColumn, Math, Table, Tabularx, Tabular, Section, Center, Alignat, Subsection, Subsubsection, NewPage
from pylatex.utils import NoEscape
import matplotlib.pyplot as plt
from energy_functions import *



#--------------  TEX TABLES HADRONS  ----------------

    # It creates the TeX tables
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
        
    # It creates the TeX tables
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

# This routine is pretty much the same than the one above, but includes the JPC numbers (the must be included in the the_had_data as the last entry for each item.)

# a: Lattice size in fm.
# N_lattice: lattice extent L.
# t_lattice: lattice extent T.
# the_had_data: list of lists, [name of hadron, mass of hadron ]
# the_caption: whatever one wants to put in the caption of this table.

def TABLE_FOUR_COLS_JPC_ENSEMBLES(the_doc, the_cls_name, a, N_lattice, t_lattice, the_beta_vals, the_kappa_u, the_kappa_s, the_had_data, the_nr_configs, the_caption):
     with the_doc.create(Table(position = 'h!')) as main_table:
        with the_doc.create(Tabular('|lccc|', col_space='1.1cm', booktabs=True)) as table:
            table.add_row(('Properties of the Lattice', '','',  'Values'))
            table.add_hline()
            table.add_row(('Lattice size', '', '', a))
            table.add_row((NoEscape('Lattice extent $L^{3}$'),'',  '',  NoEscape('$%s^{3}$'%str(N_lattice))))
            table.add_row((NoEscape('Lattice extent $T$'), '', '', NoEscape('$%s$'%t_lattice )))
            table.add_row((NoEscape(r'$\beta$'), '', '',  the_beta_vals))
            table.add_row((NoEscape('$\kappa_{u}$'), '', '', the_kappa_u))
            table.add_row((NoEscape('$\kappa_{s}$'), '', '', the_kappa_s))
            table.add_row(('CLS name', '', '', the_cls_name))
            table.add_row(('Nr. Gauge Configs', '', '', the_nr_configs))
            table.add_hline()
            table.add_hline()
            table.add_row(('Hadrons', NoEscape('$J^{P}$'), 'Masses [am]', NoEscape('Masses [MeV]')))
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


    
    


def PLOT_ENERGY_LEVELS(list_of_energies,the_ref_levels,the_y_axis_label,the_name_plot,the_nr_levels):
    line_styles=["--","-","-.",":"]
    line_colors=["#b90f22", "#5d83d5","#ffa500","#008000","#c44601","#f57600","#5ba300","#e6308a" ]
    the_plot = plt.figure()
    for ii in range(len(list_of_energies)):
        the_name_irrep = PLOT_HADRON_LABELINGS(list_of_energies[ii][0])
        x_axis, y_axis = the_name_irrep+"(%s)"%str(list_of_energies[ii][1]), list_of_energies[ii][2]
        plt.rc('axes', labelsize=10) 
        plt.plot(x_axis, y_axis, marker='_', ls='None', ms=20, markeredgewidth=2.5, lw=0.95, zorder=3, color = line_colors[1])
    for ii in range(len(the_ref_levels)):
        plt.hlines(the_ref_levels[ii][0], xmin=-2, xmax=the_nr_levels+1, ls=line_styles[ii],lw=1.25,color=line_colors[0],label=the_ref_levels[ii][1])
    plt.xlim([-.5, the_nr_levels+.05])
    plt.xlabel('Irreducible representations',fontsize=12)
    plt.xticks(rotation=45,fontsize=11)
    plt.ylabel(the_y_axis_label,fontsize=11)

    fig=plt.gcf()
    axes=fig.axes
    top_y = max(ax.get_position().y1 for ax in axes)
    the_plot.legend(loc='upper center', bbox_to_anchor=(0.5, top_y+0.08), ncol = len(the_ref_levels))
    the_plot.savefig(the_name_plot, bbox_inches='tight')
    
    
    
