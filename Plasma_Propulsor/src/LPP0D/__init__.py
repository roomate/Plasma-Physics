"""Petite description"""

__author__ = "Benjamin Esteves"
__mail__ = "benjamin.esteves@lpp.polytechnique.fr"
__version__ = '0.2.0'

from .objects_generator import generate_all_objects, Reaction
#from .diff_generator import gen_diff
#from .diff_generator_gas_heating import gen_diff_gas_heating
#from .diff_generator_GRONDEIN import gen_diff_GRONDEIN
#from .diff_generator_cell_test_GRONDEIN import gen_diff_cell_test_GRONDEIN
#from .diff_generator_cell_test import gen_diff_cell_test
from .diff_generator_all_gases import gen_diff_all_gases
#from .diff_generator_all_gases_manuscrit_2Temps import gen_diff_all_gases_manuscrit_2Temps
from .solve_and_post_proc import solve, solve_gas_heating, solve_GRONDEIN,solve_all_gases,solve_cell_test_GRONDEIN,solve_cell_test, solve_all_gases_manuscrit_2Temps
from .toolbox_functions import v_th, uB, mass, sccm_to_persecond, h0, Vplasma, Vsheath, alpha_s, nu_m, R_ind, nu_m_GRONDEIN, R_ind_GRONDEIN, uB_GRONDEIN, hL_GRONDEIN, hR_GRONDEIN, R_ind_noble_gas,R_ind_cell_test_GRONDEIN,hl,hr,Kappa_I,R_ind_cell_test,Kappa_noble_gas,R_ind_Istar_included,nu_m_Istar_included,list_of_species,Omega_I_I2,Twall_star_cell,R_ind_cell_test_Istar_included,Kappa_I2

__ALL__ = ['generate_all_objects', 'Reaction', 'gen_diff', 'gen_diff_gas_heating', 'gen_diff_GRONDEIN','gen_diff_all_gases','gen_diff_cell_test_GRONDEIN','solve','gen_diff_all_gases_manuscrit_2Temps', 'solve_gas_heating','solve_GRONDEIN','solve_cell_test_GRONDEIN', 'v_th', 'uB', 'mass', 'h0','nu_m', 'R_ind','nu_m_GRONDEIN', 'R_ind_GRONDEIN', 'Vplasma', 'Vsheath', 'alpha_s','uB_GRONDEIN', 'hL_GRONDEIN', 'hR_GRONDEIN','R_ind_Ar','R_ind_Xe','Twall_star_cell','R_ind_Kr','R_ind_noble_gas','R_ind_cell_test_GRONDEIN','R_ind_cell_test','hl','hr','Kappa_I','Kappa_noble_gas','R_ind_Istar_included','nu_m_Istar_included',
'list_of_species','Omega_I_I2','R_ind_cell_test_Istar_included','Kappa_I2']
