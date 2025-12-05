"""Petite description"""

__author__ = "Florian Marmuse, Romain Lucken, Pascaline Grondein, Pascal Chabert"
__mail__ = ""
__version__ = '0.2.0'

from .objects_generator import generate_all_objects, Reaction
from .diff_generator import gen_diff
from .solve_and_post_proc import solve
from .toolbox_functions import v_th, uB, mass, sccm_to_persecond

__ALL__ = ['generate_all_objects', 'Reaction', 'gen_diff', 'solve', 'v_th', 'uB', 'mass']
