"""Aux fun used in global models.

Benjamin Esteves
04/2022 dernière MAJ
"""

import numpy as np
from scipy.constants import pi, e, m_e, epsilon_0, k, c, m_u
# import scipy.special
# import pandas as pd
from scipy.special import erf, jv
from .molar_masses import masses_dictionary
import os
dir_path = os.path.dirname(os.path.realpath(__file__))

##############
# CONVERSION #
##############


def sccm_to_persecond(sccm):
    """Convert sccm to particle per second."""
    return sccm * 4.477962E17


def persecond_to_sccm(Q):
    """Convert particle per second to sccm."""
    return Q / 4.477962E17


def power_to_volumic_power(power, volume):
    """Convert power to volumic power."""
    return power / volume


def volumic_power_to_power(volumic_power, volume):
    """Convert volumic power to power."""
    return volumic_power * volume


def find_nearest(array, value):
    """Find nearest element."""
    idx = (np.abs(array - value)).argmin()
    return array[idx]


##########
# PLASMA #
##########

def remove_charge(element):
    """Return the neutral element corresponding to an ion.
    Example: remove_charge('Xe')
    remove_charge('Xe^+') = 'Xe'
    """

    return element.split('^')[0]  # the whole part before the ^, or the whole part if there is not ^.


def relevant_lines(text_file):
    """Return an ordered list of the non-empty lines, without the comments.
    Example: relevant_lines('input_xenon').
    """

    with open(text_file, 'r') as input_txt:
        return [
               line.split('#')[0].strip()  # the part before # of each line
               for line in input_txt.readlines()
               if (line.rstrip() is not '' and line.startswith('#') is False)
               ]  # if line not empty and not starting with #


def get(variable, section, text_file):
    """Return the value of a variable in a given section of a text file.

    Quite specific to the input file for LPP0D.
    Case sensitive.

    Example: get('R', 'THRUSTER', 'input_xenon.txt')
    """

    for line in lines_in_section(section, text_file):
        if line.startswith(variable):
            try:
                return eval(line.split('=')[1].strip())  # main return here
            except NameError:
                return line.split('=')[1].strip()

    # if nothing found:
    return 'This variable is missing.'


def list_of_reactions(text_file):
    """Return a list of the reactions, i.e. a list containing the lines of the "REACTIONS" section."""
    return lines_in_section('REACTIONS', text_file)


def reaction_elements(reaction):
    """Return a list [[reactants], [products], rate] for each reaction line you give it."""

    lhs = reaction.split('->')[0]
    reactants = [item.strip() for item in lhs.split(' + ')]

    rhs = reaction.split('->')[1].split('!')[0]
    products = [item.strip() for item in rhs.split(' + ')]
    # print(products)

    # rate = reaction.rsplit('!')[-1].strip()
    reaction_type = reaction.rsplit('!')[-1].split('&')[0].strip()

    #rate_given = float(reaction.rsplit('!')[-1].split('&')[1].strip()) if '&' in reaction.rsplit('!')[-1] else False

    # return [reactants, products, rate]
    return [reactants, products, reaction_type]


def format_rate(reaction):
    [reactants, products, type] = reaction_elements(reaction)
    return 'R[\'{}\'][\'{}\'][\'{}\'].K(Te)'.format(reactants[0], reactants[1], type).replace('\'e\'', '\'e^-\'')


def format_energy(reaction):
    [reactants, products, type] = reaction_elements(reaction)
    return 'R[\'{}\'][\'{}\'][\'{}\'].energy_loss'.format(reactants[0], reactants[1], type).replace('\'e\'', '\'e^-\'')


def all_reactions(text_file):
    """Return a list of the reactions under the format [reac1_reac2_type]"""
    return [l.split('!')[-1].strip().split('(')[0][2:] for l in list_of_reactions(text_file)]


def mass(part):
    """ Returns the mass of a particle, in kg.

    mass('I2p') = 126.905 * 2 * m_u
    """

    while part[-1] == 'p' or part[-1] == 'm' or part[-1] == '+' or part[-1] == '-' or part[-1] == '^' or part[-1] == '*':  # for I2p or Im I guess.
        part = part[:-1]  # remove last character

    if part == 'e':  # TODO nomenclature précise des espèces, conventions de notations.
        return m_e

    else:

        if str.isdigit(part[-1]):
            number_of_elements = eval(part[-1])
            part = part[:-1]

        else:
            number_of_elements = 1

        return masses_dictionary[part] * number_of_elements * m_u

#####


def line_number(reaction_name, lxfile):
    """Return the line number from which to read the data for a given reaction type."""

    with open(lxfile, 'r') as f:

        reaction_number = eval(reaction_name[-1]) if reaction_name[-1].isdigit() else 1

        lines = f.readlines()
        good_line_numbers = []

        for (a, b) in enumerate(lines):
            if b.startswith(reaction_name):
                good_line_numbers.append(a)
        try:
            return good_line_numbers[reaction_number-1]
        except IndexError:
            print('reaction {} not found'.format(reaction_name))

#######5te

def uB(ne, nIm, Te, Tg, part):
    """Bohm velocity.
    INPUT : ne, nIm, Te, Tg, part
    """
    alpha = nIm / ne
    gamma = Te/Tg
    ub_pos = np.sqrt(e * Te / mass(part))
    return ub_pos * np.sqrt((1+alpha_s(alpha, gamma))/(1 + gamma * alpha_s(alpha, gamma)))
    #return ub_pos

def uB_GRONDEIN(Te,part):
    """Bohm velocity.
    INPUT : Te, part
    """
    ub_pos = np.sqrt(e * Te / mass(part))
    return ub_pos 

def hL_GRONDEIN(L,ng):
    """Chabert 2012 """
    sigma_i = 1.e-18
    li = (ng * sigma_i) ** (-1)
    l = L/2
    return 0.86 * (3 + l/li) ** (-1 / 2) 

def hR_GRONDEIN(R,ng):
    """Chabert 2012 """
    sigma_i = 1.e-18
    li = (ng * sigma_i) ** (-1)
    return 0.8 * (4 + R/li) ** (-1 / 2)




def hl(L, ng, ne, nIm, Te, Tg, sigma):
    """Chabert 2016 """
    alpha = float(nIm/ne)
    li = (ng * sigma) ** (-1)
    l = L/2
    gamma = Te/Tg
    term1 = 0.86 * (3 + l/li + (1+alpha)**(1/2) * 0.2 * (0.1 / Te) * (2*l/li)**2) ** (-1 / 2)
    term2 = (((gamma-1)/(gamma*(1+alpha)**2))+(1/gamma))**(1/2)
    return term1 * term2



# def hr(R, ng, ne, nIm, Te, Tg):
#     """Edge-to-center plasma density ratio.
#     Answer: Yes, we should have different sigma_i for each ion-atom couple."""
#     alpha = nIm / ne
#     sigma_i = 1.e-18
#     li = (ng * sigma_i) ** (-1)
#     return (1 + alpha) ** (-1) * 0.8 * (4 + R / li + (Tg / Te) * (R/li) ** 2) ** (-1 / 2)

def hr(R, ng, ne, nIm, Te, Tg,sigma):
    """inspired from Chabert 2016 """
    alpha = nIm / ne
    li = (ng * sigma) ** (-1)
    gamma = Te/Tg
    term1 = 0.8 * (4 + R/li + (1+alpha)**(1/2) * 0.41 * (0.1/Te) * (R/li)**2) ** (-1 / 2)
    term2 = (((gamma-1)/(gamma*(1+alpha)**2))+(1/gamma))**(1/2)
    return term1 * term2

hl=np.vectorize(hl)
hr=np.vectorize(hr)

def h0(L, R, ng, ne, nIm, Te, Tg):
    """Global h factor
    INPUT : L, R, ng, ne, nIm, Te, Tg
    """
    total_area = 2 * pi * R * R + 2 * pi * R * L
    effective_area = hr(R, ng, ne, nIm, Te, Tg,sigma) * 2 * pi * L * R + hl(L, ng, ne, nIm, Te, Tg,sigma) * 2 * pi * R ** 2
    return effective_area / total_area


def alpha_s(alpha, gamma):  # Thorsteinsson 2010
    a1 = 0.607
    a2 = 5.555
    a3 = -11.16
    a4 = 1.634
    a5 = 12e-3
    a6 = -107e-3

    rho = np.abs(alpha + a5 * (np.exp(a6 * (gamma - 50))-1))

    num = alpha * a1 * erf(a2 * rho + a3) * np.exp(-a4 / rho**1.35)
    denom = np.exp((gamma-1)/(2*gamma) - 0.49)

    return np.maximum(num/denom, 0)

def nu_m(nI,nI2,Te,R):
    """Collision frequency.
    INPUT : nI,nI2,Te
    OUPUT : Elastic collision frequency
    """
    for reaction in R['I']['e^-'].keys():
        if reaction[:5] == 'ELAST':
            K_I = float(R['I']['e^-'][reaction].K(Te))
    K_I2 = float(R['I2']['e^-']['ELASTIC'].K(Te))
    #return (nI * K_I + nI2 * K_I2)
    return (nI * 2.5e-13 + nI2 * K_I2)

def nu_m_Istar_included(nI,nIstar,nI2,Te,R):
    """Collision frequency.
    INPUT : nI,nIstar,nI2,Te
    OUPUT : Elastic collision frequency
    """
    for reaction in R['I']['e^-'].keys():
        if reaction[:5] == 'ELAST':
            K_I = float(R['I']['e^-'][reaction].K(Te))
    for reaction in R['I']['e^-'].keys():
        if reaction[:5] == 'ELAST':
            K_Istar = float(R['I*']['e^-'][reaction].K(Te))
    K_I2 = float(R['I2']['e^-']['ELASTIC'].K(Te))
    return (nI * K_I + nIstar * K_Istar + nI2 * K_I2)
    #return (nI * 2.5e-13 + nIstar * K_Istar + nI2 * K_I2)

nu_m_Istar_included=np.vectorize(nu_m_Istar_included)

def R_ind_cell_test_GRONDEIN(ne,nI,nI2,Te,Radius,L,R):
    """Resistive part of plasma impedance, Ohm.
    INPUT : ne,nI,nI2,Te,Tg,Radius R, length L, omega frequency w, N turns
    OUPUT : Inductive resistance from Grondein et al.
    """
    
    #wpe = np.sqrt(ne * e**2 / (m_e * epsilon_0))
    nu_eff = nu_m_GRONDEIN(nI,nI2,Te,R) #+ nu_stoch(u) * simu.stoch_heat
    V = np.pi * Radius**2 * L
    A = np.pi * Radius**2

    return m_e*nu_eff*V/(A**2*ne*e**2)

def R_ind_cell_test(ne,nI,nI2,Te,Radius,L,R):
    """Resistive part of plasma impedance, Ohm.
    INPUT : ne,nI,nI2,Te,Tg,Radius R, length L, omega frequency w, N turns
    OUPUT : Inductive resistance from Grondein et al.
    """
    
    #wpe = np.sqrt(ne * e**2 / (m_e * epsilon_0))
    nu_eff = nu_m(nI,nI2,Te,R) #+ nu_stoch(u) * simu.stoch_heat
    V = np.pi * Radius**2 * L
    A = np.pi * Radius**2

def R_ind_cell_test_Istar_included(ne,nI,nIstar,nI2,Te,Radius,L,R,Tg):
    """Resistive part of plasma impedance, Ohm.
    INPUT : ne,nI,nIstar,nI2,Te,Radius R, length L, chem R
    OUPUT : Inductive resistance from Grondein et al.
    """
    
    #wpe = np.sqrt(ne * e**2 / (m_e * epsilon_0))
    #nu_eff = nu_m(nI,nI2,Te,R) #+ nu_stoch(u) * simu.stoch_heat
    nu_eff = nu_m_Istar_included(nI,nIstar,nI2,Te,R)
    V = np.pi * Radius**2 * L
    A = np.pi * Radius**2

    return m_e*nu_eff*V/(A**2*ne*e**2)

def R_ind(ne,nI,nI2,Te,Tg,coil_Radius,coil_Length,w,N,R):
    """Resistive part of plasma impedance, Ohm.
    INPUT : ne,nI,nI2,Te,Tg,coil_Radius,coil_Length,w frequency,N turns,R chemistry
    OUPUT : Inductive resistance from Grondein et al.
    """
    
    wpe = np.sqrt(ne * e**2 / (m_e * epsilon_0))
    nu_eff = nu_m(nI,nI2,Te,R) #+ nu_stoch(u) * simu.stoch_heat

    eps_p = 1 - (wpe**2 / (w * (w - 1j*nu_eff)))
    k2 = (w/c) * np.sqrt(eps_p)

    Rind_1 = (2 * pi * N**2)/(coil_Length * w * epsilon_0)
    Rind_2 = 1j * k2 * coil_Radius * jv(1, k2 * coil_Radius)
    Rind_3 = eps_p * jv(0, k2 * coil_Radius)
    # verifier pour nu_m = 0
    return Rind_1 * np.real(Rind_2 / Rind_3)

def R_ind_Istar_included(ne,nI,nIstar,nI2,Te,Tg,coil_Radius,coil_Length,w,N,R):
    """Resistive part of plasma impedance, Ohm.
    INPUT : ne,nI,nIstar,nI2,Te,Tg,coil_Radius,coil_Length,w frequency,N turns,R chemistry
    OUPUT : Inductive resistance from Grondein et al.
    """
    
    wpe = np.sqrt(ne * e**2 / (m_e * epsilon_0))
    nu_eff = nu_m_Istar_included(nI,nIstar,nI2,Te,R) #+ nu_stoch(u) * simu.stoch_heat

    eps_p = 1 - (wpe**2 / (w * (w - 1j*nu_eff)))
    k2 = (w/c) * np.sqrt(eps_p)

    Rind_1 = (2 * pi * N**2)/(coil_Length * w * epsilon_0)
    Rind_2 = 1j * k2 * coil_Radius * jv(1, k2 * coil_Radius)
    Rind_3 = eps_p * jv(0, k2 * coil_Radius)
    # verifier pour nu_m = 0
    return Rind_1 * np.real(Rind_2 / Rind_3)


def Kappa_I(Tg):
    return 6.32e-5 * (Tg*e/k)**(0.768)

def Kappa_I2(T_I2):
    return 2.42e-5 * (T_I2*e/k)**(0.849)

def Kappa_noble_gas(Tg,main_element):
    if main_element == 'Xe':
        return 7.98e-5 * (Tg*e/k)**(0.745)
    if main_element == 'Ar':
        return 3.59e-4 * (Tg*e/k)**(0.681)
    if main_element == 'Kr':
        return (3.59e-4 * (Tg*e/k)**(0.681)+7.98e-5 * (Tg*e/k)**(0.745))/2


def R_ind_noble_gas(ne,ng,Te,Tg,Radius,L,w,N,R,main_element):
    """Resistive part of plasma impedance, Ohm.
    INPUT : ne,ng,Te,Tg,coil_Radius,coil_Length,w frequency,N turns,R chemistry,main_element
    OUPUT : Inductive resistance from Grondein et al.
    """
    for reaction in R[main_element]['e^-'].keys():
        if reaction[:5] == 'ELAST':
            nu_eff = float(R[main_element]['e^-'][reaction].K(Te))*ng 

    wpe = np.sqrt(ne * e**2 / (m_e * epsilon_0))
    eps_p = 1 - (wpe**2 / (w * (w - 1j*nu_eff)))
    k2 = (w/c) * np.sqrt(eps_p)

    Rind_1 = (2 * pi * N**2)/(L * w * epsilon_0)
    Rind_2 = 1j * k2 * Radius * jv(1, k2 * Radius)
    Rind_3 = eps_p * jv(0, k2 * Radius)
    # verifier pour nu_m = 0
    return Rind_1 * np.real(Rind_2 / Rind_3)

def nu_m_GRONDEIN(nI,nI2,Te,R):
    """Collision frequency.
    INPUT : nI,nI2,Te
    OUPUT : Elastic collision frequency
    """
    K_I = float(R['I']['e^-']['ELASTIC_GRONDEIN'].K(Te))
    K_I2 = float(R['I2']['e^-']['ELASTIC_GRONDEIN'].K(Te))
    #return (nI + nI2) * K_I                     # in Grondein paper
    return (nI*K_I + nI2*K_I2)

def R_ind_GRONDEIN(ne,nI,nI2,Te,Tg,Radius,L,w,N,R):
    """Resistive part of plasma impedance, Ohm.
    INPUT : ne,nI,nI2,Te,Tg,Radius R, length L, omega frequency w, N turns
    OUPUT : Inductive resistance from Grondein et al.
    """
    
    wpe = np.sqrt(ne * e**2 / (m_e * epsilon_0))
    nu_eff = nu_m_GRONDEIN(nI,nI2,Te,R) #+ nu_stoch(u) * simu.stoch_heat

    eps_p = 1 - (wpe**2 / (w * (w - 1j*nu_eff)))
    k2 = (w/c) * np.sqrt(eps_p)

    Rind_1 = (2 * pi * N**2)/(L * w * epsilon_0)
    Rind_2 = 1j * k2 * Radius * jv(1, k2 * Radius)
    Rind_3 = eps_p * jv(0, k2 * Radius)
    # verifier pour nu_m = 0
    return Rind_1 * np.real(Rind_2 / Rind_3)


def Omega_I_I2(T_I,T_I2):
    T_I2_I = 1/3.*(T_I2+2*T_I)
    m_I2_I = mass('I')*mass('I2')/(mass('I')+mass('I2'))
    d_I=280e-12
    d_I2=2*d_I
    d_I2_I=0.5*(d_I+d_I2)
    #print(np.sqrt(e*T_I2_I/(2*np.pi*m_I2_I))*np.pi*d_I2_I**2)
    return np.sqrt(e*T_I2_I/(2*np.pi*m_I2_I))*np.pi*d_I2_I**2

def Twall_star_cell(Pabs):
    return 293+Pabs*1.4

# potentiel plasma qui varie avec les densité, comment je gère ça ? pas aujourd'hui en tous cas.
# def mean_ub(text_file):
#     species = list_of_species(text_file)
#
#     for specie in spe

def Vplasma(Te, Tg, ne, nIm):  # plasma potential, always positive
    alpha = nIm / ne
    gamma = Te/Tg
    alphas = alpha_s(alpha, gamma)

    return 0.5 * (1 + alphas)/(1 + gamma*alphas) * Te


def Vsheath(Te, Tg, ne, nI, nI2, nIp, nI2p, nIm):  # sheath potential, always positive (abs value)

    meanub = np.average([uB(ne, nIm, Te, Tg, 'I'), uB(ne, nIm, Te, Tg, 'I2')], weights=[nIp, nI2p], axis=0)
    meanv = np.average([v_th(Tg, 'I'), v_th(Tg, 'I2')], weights=[nI, nI2], axis=0)

    alpha = nIm / ne
    gamma = Te / Tg

    alphas = alpha_s(alpha, gamma)

    ve = np.sqrt(8 * e * Te / (pi * m_e))

    result = Te * np.log(4 * (meanub / ve) * (1 + alphas) / (1 + alphas * (meanv / ve) ** 2))

    return np.maximum(-result, 0)


def v_th(Tg, part):
    """Thermal velocity."""
    return np.sqrt((8 * e * Tg) / (pi * mass(part)))


def wpe(n):
    """plasma frequency [rad/s]"""
    return np.sqrt(e**2*n/(epsilon_0*m_e)) 


def l_Debye(n, T):
    """Debye length [m]"""
    return np.sqrt(epsilon_0*T/(e*n))   


def v_Te(T):
    """ Electron 'thermal' velocity """
    return (e*T/m_e)**0.5


def init(text_file):
    """Return a vector of initial densities and temperature"""
    lines = relevant_lines(text_file)
    section_line_numbers = [i for i, x in enumerate(lines) if '$$' in x]

    index_species = lines.index('$$ SPECIES $$')
    index_end_species = section_line_numbers[section_line_numbers.index(index_species)+1]

    index_temperatures = lines.index('$$ TEMPERATURES $$')
    index_end_temperatures = section_line_numbers[section_line_numbers.index(index_temperatures)+1]

    irange = np.append(np.arange(index_species+1, index_end_species),
                       np.arange(index_temperatures+1, index_end_temperatures))
    return [eval(lines[i].split('=')[1].strip()) for i in irange]


###################


def lines_in_section(section, text_file):
    """return the list of the lines in a given section"""
    lines = relevant_lines(text_file)
    section_line_numbers = [i for i, x in enumerate(lines) if '$$' in x]

    index_section = lines.index('$$ {} $$'.format(section))  # number of the line where $$ SECTION $$ is.
    index_end_section = section_line_numbers[
                                        section_line_numbers.index(index_section) + 1]  # number of the next section.
    return [lines[i] for i in np.arange(index_section + 1, index_end_section)]


def list_of_species(text_file):
    """Return a list of the species declared, for example ['e^-','Xe','Xe^+'].
    Example: list_of_species('input_xenon.txt')
    """

    lines = lines_in_section('SPECIES', text_file)
    return [line.split('=')[0].strip() for line in lines]
