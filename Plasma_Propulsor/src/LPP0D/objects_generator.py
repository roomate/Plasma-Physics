########################################################################
#                                                                      #
#    #######  ########        ## ########  ######  ########  ######    #
#   ##     ## ##     ##       ## ##       ##    ##    ##    ##    ##   #
#   ##     ## ##     ##       ## ##       ##          ##    ##         #
#   ##     ## ########        ## ######   ##          ##     ######    #
#   ##     ## ##     ## ##    ## ##       ##          ##          ##   #
#   ##     ## ##     ## ##    ## ##       ##    ##    ##    ##    ##   #
#    #######  ########   ######  ########  ######     ##     ######    #
#                                                                      #
########################################################################

# This file encompasses the classes objects, thruster, and propellant.

# Written by Florian Marmuse, January 2018

# Requirements for the input text file :
# # is to comment
# $$ shall only appear in section titles
# finish with $$ ENDFILE $$

# Best practices for the input text file :
# Write e- first in the species.

# TBD: remove the absolute path to files.
# raise a warning is a species in the reaction set is not in species

import numpy as np
from scipy.constants import e, pi, m_e, k
import pandas as pd

from .toolbox_functions import persecond_to_sccm, line_number, mass, reaction_elements, lines_in_section, get, v_Te

import os
dir_path = os.path.dirname(os.path.realpath(__file__))

###############################################################################
#                                                                             #
# ######## ##     ## ########  ##     ##  ######  ######## ######## ########  #
#    ##    ##     ## ##     ## ##     ## ##    ##    ##    ##       ##     ## #
#    ##    ##     ## ##     ## ##     ## ##          ##    ##       ##     ## #
#    ##    ######### ########  ##     ##  ######     ##    ######   ########  #
#    ##    ##     ## ##   ##   ##     ##       ##    ##    ##       ##   ##   #
#    ##    ##     ## ##    ##  ##     ## ##    ##    ##    ##       ##    ##  #
#    ##    ##     ## ##     ##  #######   ######     ##    ######## ##     ## #

###############################################################################


class NewThruster:
    def __init__(self, input_txt):
        """"Constructor"""
        self.R               = get('R',             'THRUSTER', input_txt)
        self.L               = get('L',             'THRUSTER', input_txt)
        self.beta_neutrals   = get('beta_neutrals', 'THRUSTER', input_txt)
        self.beta_ions       = get('beta_ions',     'THRUSTER', input_txt)
        self.R_coil          = get('R_coil',        'THRUSTER', input_txt)
        self.w               = get('w',             'THRUSTER', input_txt)
        self.N               = get('N',             'THRUSTER', input_txt)
        self.Vgrids          = get('Vgrids',        'THRUSTER', input_txt)
        self.grid_dist       = get('grid_dist',     'THRUSTER', input_txt)
        self.r_coil          = get('r_coil',     'THRUSTER', input_txt)
        self.l_coil          = get('l_coil',     'THRUSTER', input_txt)
        self.UQ_coef = 1  # for UQ

    def volume(self):
        """Volume of the thruster."""
        return pi * self.R**2 * self.L

    def total_area(self):
        """Total area of the thruster."""
        return 2 * pi * self.R**2 + 2 * pi * self.R * self.L

    def open_area_neutrals(self):
        """Open area of the thruster for the gas."""
        return self.beta_neutrals * pi * self.R**2

    def open_area_ions(self):
        """Open area of the thruster for the ions."""
        return self.beta_ions * pi * self.R**2

    def heat_diff_length(self):
        """Heat diffusion length, from Chabert, 2012."""
        return np.sqrt((self.R / 2.405)**2 + (self.L / pi)**2)

    def __str__(self):
        """Readable description of the object, to be called with print."""
        return "The thruster\'s parameters are:\n" \
            "    R = {} m\n" \
            "    L = {} m\n" \
            "    R_coil = {} Ohm\n" \
            "    w = {} Hz\n" \
            "    N = {} \n" \
            "    Vgrids = {} V\n" \
            "    grid_dist = {} m\n" \
            "    r_coil = {} m\n" \
            "    l_coil = {} m\n" \
            "    beta_neutrals = {} \n" \
            "    beta_ions = {} \n\n".format(self.R,
                                             self.L,self.R_coil,self.w,self.N,self.Vgrids,self.grid_dist,
                                             self.r_coil,self.l_coil,self.beta_neutrals,self.beta_ions)


###################################################################################################
#
# ########     ###    ########     ###    ##     ## ######## ######## ######## ########   ######
# ##     ##   ## ##   ##     ##   ## ##   ###   ### ##          ##    ##       ##     ## ##    ##
# ##     ##  ##   ##  ##     ##  ##   ##  #### #### ##          ##    ##       ##     ## ##
# ########  ##     ## ########  ##     ## ## ### ## ######      ##    ######   ########   ######
# ##        ######### ##   ##   ######### ##     ## ##          ##    ##       ##   ##         ##
# ##        ##     ## ##    ##  ##     ## ##     ## ##          ##    ##       ##    ##  ##    ##
# ##        ##     ## ##     ## ##     ## ##     ## ########    ##    ######## ##     ##  ######
#
###################################################################################################


class NewParameters:
    def __init__(self, input_txt):
        """Constructor."""
        # self.B0                 = get('B0',               'PARAMETERS', input_txt)
        self.Q0                 = get('Q0',               'PARAMETERS', input_txt)
        self.wall_temperature   = get('wall_temperature', 'PARAMETERS', input_txt)
        self.Pabs               = get('Pabs',             'PARAMETERS', input_txt)
        # self.stoch_heat         = get('stoch_heat',       'PARAMETERS', input_txt)
        self.Tg                 = get('Tg',               'PARAMETERS', input_txt)
        # self.pressure           = get('pressure',         'PARAMETERS', input_txt)
        # self.EPS                = get('EPS',              'PARAMETERS', input_txt)
        self.gamma              = get('gamma',            'PARAMETERS', input_txt)
        self.I_coil             = get('I_coil',            'PARAMETERS', input_txt)
        self.sigma              = get('sigma',             'PARAMETERS', input_txt)


    def __str__(self):
        """Readable description of the object, to be called with print."""
        return "The parameters are: \n" \
               "    Q0 = {:.2e} persecond\n" \
               "    Q0 = {:.2e} sccm\n" \
               "    Pabs = {:.2e} W/m3\n" \
               "    I_coil = {:.2e} A\n" \
               "    gamma = {}\n" \
               "    T_wall = {:.1f} K\n"\
               "    sigma = {:.1e} m2".format(self.Q0, persecond_to_sccm(self.Q0), self.Pabs,self.I_coil, self.gamma, self.wall_temperature*e/k,self.sigma)


###############################################################
#
#  ######  ########  ########  ######  #### ########  ######
# ##    ## ##     ## ##       ##    ##  ##  ##       ##    ##
# ##       ##     ## ##       ##        ##  ##       ##
#  ######  ########  ######   ##        ##  ######    ######
#       ## ##        ##       ##        ##  ##             ##
# ##    ## ##        ##       ##    ##  ##  ##       ##    ##
#  ######  ##        ########  ######  #### ########  ######
#
###############################################################


class Species:
    """Contains the species. Contains class attributes, class methods and static methods."""

    number_of_species = 0

    def __init__(self, species_name):
        self.name = species_name
        self.mass = mass(species_name)

        Species.number_of_species += 1

    def __str__(self):
        return "Species {} with mass {}.".format(self.name, self.mass)


##################################################################################
#
# ########  ########    ###     ######  ######## ####  #######  ##    ##  ######
# ##     ## ##         ## ##   ##    ##    ##     ##  ##     ## ###   ## ##    ##
# ##     ## ##        ##   ##  ##          ##     ##  ##     ## ####  ## ##
# ########  ######   ##     ## ##          ##     ##  ##     ## ## ## ##  ######
# ##   ##   ##       ######### ##          ##     ##  ##     ## ##  ####       ##
# ##    ##  ##       ##     ## ##    ##    ##     ##  ##     ## ##   ### ##    ##
# ##     ## ######## ##     ##  ######     ##    ####  #######  ##    ##  ######
#
##################################################################################

class Chemistry:
    """Carries R in the parametric analyses."""
    def __init__(self, input_txt):
        self.R = dict()
        self.init_vector = self.make_init_vec(input_txt)

        for line in lines_in_section('REACTIONS', input_txt):
            tmp = Reaction(line)

            try:
                self.R[tmp.reactants[0]][tmp.reactants[1]][tmp.reactype] = tmp
            except KeyError:
                try:
                    self.R[tmp.reactants[0]][tmp.reactants[1]] = dict()
                    self.R[tmp.reactants[0]][tmp.reactants[1]][tmp.reactype] = tmp
                except KeyError:
                    try:
                        self.R[tmp.reactants[0]] = dict()
                        self.R[tmp.reactants[0]][tmp.reactants[1]] = dict()
                        self.R[tmp.reactants[0]][tmp.reactants[1]][tmp.reactype] = tmp
                    except KeyError:
                        print("Init of reaction ", tmp.__str__, "failed.")

    def make_init_vec(self, input_txt):
        init_vec = []

        species_and_init = lines_in_section('SPECIES', input_txt)
        for a in species_and_init:
            init_vec.append(float(a.split('=')[-1]))

        temp_and_init = lines_in_section('TEMPERATURES', input_txt)
        for b in temp_and_init:
            init_vec.append(float(b.split('=')[-1]))

        return init_vec


class Reaction:
    """Contains a reaction."""
    def __init__(self, line_in_the_reaction_section):

        self.UQ_coef = 1

        self.reaction_elements = reaction_elements(line_in_the_reaction_section)

        self.reactants = self.reaction_elements[0]      # example ['Xe', 'e^-']
        self.products = self.reaction_elements[1]       # example ['Xe^+', 'e^-', 'e^-']

        self.reactype = self.reaction_elements[2]       # example 'IONIZATION'

        cwd = dir_path  # os.getcwd()  # current working directory
        # print(cwd)

        self.ratefile = cwd+'/../lxcat/{}_{}_{}_rates_and_energy.csv'.format(self.reactants[0], self.reactants[1], self.reactype).replace('e^-', 'e')

        try:
            self.df_data = pd.read_csv(self.ratefile, skiprows=1, names=['Te', 'K', 'energy', 'Kdruy'])  # forth column with druy rate.
        except:
            print("FileNotFoundError while building rate table for {} with {}.".format(self.reactype, self.reactants))

    def __str__(self):
        """Readable description of the object."""
        return "Reaction description:\n" \
               "    {}\n" \
               "    {}\n" \
               "    UQ coef: {}" \
               "    K(3ev): {}".format(self.reaction_elements, self.ratefile, self.UQ_coef, self.K(3))

    def energy_loss(self, T):
        """Energy loss (eV) of a given reaction, to calculate Ploss. T in eV.
        Corresponds to the energy in line 2 of lxcat columns.
        """
        return self.UQ_coef * np.interp(T, self.df_data.Te, self.df_data.energy)

    def K(self, T):  # TODO: sort K(Te), K(Tg), K(T-)...
        """reaction rate, for T in eV."""
        return self.UQ_coef * np.interp(T, self.df_data.Te, self.df_data.K)
        #return self.UQ_coef * np.interp(T, self.df_data.Te, self.df_data.Kdruy)


########################################################################################
#
#  ######   ######## ##    ## ######## ########     ###    ########  #######  ########
# ##    ##  ##       ###   ## ##       ##     ##   ## ##      ##    ##     ## ##     ##
# ##        ##       ####  ## ##       ##     ##  ##   ##     ##    ##     ## ##     ##
# ##   #### ######   ## ## ## ######   ########  ##     ##    ##    ##     ## ########
# ##    ##  ##       ##  #### ##       ##   ##   #########    ##    ##     ## ##   ##
# ##    ##  ##       ##   ### ##       ##    ##  ##     ##    ##    ##     ## ##    ##
#  ######   ######## ##    ## ######## ##     ## ##     ##    ##     #######  ##     ##
#
#########################################################################################


def generate_all_objects(text_file):
    """Generates thruster, params, and all species and reactions."""
    my_thruster = NewThruster(text_file)
    my_params = NewParameters(text_file)
    my_chem = Chemistry(text_file)

    return my_thruster, my_params, my_chem