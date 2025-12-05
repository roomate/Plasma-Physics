"""This is the solve function. And where the post-processing used to be."""


from scipy.constants import m_e, e, epsilon_0, c, pi
import numpy as np
import datetime as dt
from scipy.special import jv
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# all power are power densities.
# all temperatures are in eV.


def solve(thruster, params, chem, init, temporal=False):
    """Solver for the scikits.odes cvode method."""
    # extra_options = {'old_api': False, 'max_steps': 1E9}  # 5000 is enough for Xenon.

    # remove "_iodine_thesis" if your are not me !

    #from generated_diff_iodine_thesis import differentials
    from generated_diff import differentials

    def rhs(t, y):
            return differentials(t, y, [thruster, params, chem])

    #     ode_solver = ode('cvode', rhs, **extra_options)

    if temporal:

        tf = 1
        time_span = np.arange(0, tf, 1e-3)
        output = solve_ivp(rhs, [0, tf], init, method='BDF', dense_output=True, t_eval=time_span)

        return output

    else:

        output = solve_ivp(rhs, [0, 10], init, method='BDF', dense_output=True)

        return output



def solve_gas_heating(thruster, params, chem, init, temporal=False):
    """Solver for the scikits.odes cvode method."""
    # extra_options = {'old_api': False, 'max_steps': 1E9}  # 5000 is enough for Xenon.

    # remove "_iodine_thesis" if your are not me !

    #from generated_diff_iodine_thesis import differentials
    from generated_diff_gas_heating import differentials

    def rhs(t, y):
            return differentials(t, y, [thruster, params, chem])

    #     ode_solver = ode('cvode', rhs, **extra_options)

    if temporal:

        tf = 0.5
        time_span = np.arange(0, tf, 1e-3)
        output = solve_ivp(rhs, [0, tf], init, method='BDF', dense_output=True, t_eval=time_span)

        return output

    else:

        output = solve_ivp(rhs, [0, 1], init, method='BDF', dense_output=True)

        return output


def solve_GRONDEIN(thruster, params, chem, init, temporal=False):
    """Solver for the scikits.odes cvode method."""
    # extra_options = {'old_api': False, 'max_steps': 1E9}  # 5000 is enough for Xenon.

    from generated_diff_GRONDEIN import differentials

    def rhs(t, y):
            return differentials(t, y, [thruster, params, chem])

    #     ode_solver = ode('cvode', rhs, **extra_options)

    if temporal:

        tf = 0.5
        time_span = np.arange(0, tf, 1e-3)
        output = solve_ivp(rhs, [0, tf], init, method='BDF', dense_output=True, t_eval=time_span)

        return output

    else:

        output = solve_ivp(rhs, [0, 1], init, method='BDF', dense_output=True)

        return output

def solve_cell_test_GRONDEIN(thruster, params, chem, init, temporal=False):
    """Solver for the scikits.odes cvode method."""
    # extra_options = {'old_api': False, 'max_steps': 1E9}  # 5000 is enough for Xenon.

    from generated_diff_cell_test_GRONDEIN import differentials

    def rhs(t, y):
            return differentials(t, y, [thruster, params, chem])

    #     ode_solver = ode('cvode', rhs, **extra_options)

    if temporal:

        tf = 0.5
        time_span = np.arange(0, tf, 1e-3)
        output = solve_ivp(rhs, [0, tf], init, method='BDF', dense_output=True, t_eval=time_span)

        return output

    else:

        output = solve_ivp(rhs, [0, 1], init, method='BDF', dense_output=True)

        return output

def solve_cell_test(thruster, params, chem, init, temporal=False):
    """Solver for the scikits.odes cvode method."""
    # extra_options = {'old_api': False, 'max_steps': 1E9}  # 5000 is enough for Xenon.

    from generated_diff_cell_test import differentials

    def rhs(t, y):
            return differentials(t, y, [thruster, params, chem])

    #     ode_solver = ode('cvode', rhs, **extra_options)

    if temporal:

        tf = 1
        time_span = np.arange(0, tf, 1e-3)
        output = solve_ivp(rhs, [0, tf], init, method='BDF', dense_output=True, t_eval=time_span)

        return output

    else:

        output = solve_ivp(rhs, [0, 1], init, method='BDF', dense_output=True)

        return output

def solve_all_gases(thruster, params, chem, init, temporal=False):
    """Solver for the scikits.odes cvode method."""
    # extra_options = {'old_api': False, 'max_steps': 1E9}  # 5000 is enough for Xenon.

    from generated_diff_all_gases import differentials

    def rhs(t, y):
            return differentials(t, y, [thruster, params, chem])

    # def rhs(t, y):
    #      return [0, 0, 0, 0, 0]

    #     ode_solver = ode('cvode', rhs, **extra_options)

    if temporal:

        tf = 1.0
        time_span = np.arange(0, tf, 1e-3)
        output = solve_ivp(rhs, [0, tf], init, method='BDF', dense_output=True, t_eval=time_span)

        return output

    else:

        output = solve_ivp(rhs, [0, 1], init, method='BDF', dense_output=True)

        return output

def solve_all_gases_manuscrit_2Temps(thruster, params, chem, init, temporal=False):
    """Solver for the scikits.odes cvode method."""
    # extra_options = {'old_api': False, 'max_steps': 1E9}  # 5000 is enough for Xenon.

    from generated_diff_all_gases_manuscrit_2Temps import differentials

    def rhs(t, y):
            return differentials(t, y, [thruster, params, chem])

    #     ode_solver = ode('cvode', rhs, **extra_options)

    if temporal:

        tf = 1
        time_span = np.arange(0, tf, 1e-3)
        output = solve_ivp(rhs, [0., tf], init, method='BDF', dense_output=True, t_eval=time_span)

        return output

    else:

        output = solve_ivp(rhs, [0, 1], init, method='BDF', dense_output=True)

        return output
# FIXME: Post-processing to clean-up!


    #########################
    # ## Post processing ## #
    #########################







    # return [convergeFlag,
            # simu.Pabs,
            # simu.Pabs * thruster.volume(),
            # simu.Q0,
            # thruster.R,
            # thruster.L,
            # thruster.volume(),
            # deff(out),
            # simu.EPS,
            # thruster.Vgrids,
            # Y_steady[0],
            # Y_steady[1],
            # Y_steady[2],
            # Y_steady[3],
            # pressure(out),
            # uB(out),
            # RFpower(out),
            # thrustpower(out),
            # totalpower(out),
            # J_i(out),
            # Jcl(out),
            # thrust(out),
            # thrust_ions(out),
            # thrust_neutrals(out),
            # zeta(out),
            # gamma(out),
            # ksi(out),
            # eta(out),
            # thruster_eff(out),
            # discharge_loss(out),
            # R_ind(out),
            # nu_m(out),
            # nu_stoch(out),
            # lambda_d(out),
            # sheath_thickness(out),
            # hl(out),
            # hl_new(out),
            # hr(out)
            # ]


# uB, lD, sheath lengths



def pressure(u):
    """Neutral partial pressure considered total, mTorr."""
    # was wrong before..
    [ne, ng, Te, Tg] = unnormalize(u)
    return (ng * e * Tg) * 1e-2 * 1e3 * (1.33)**(-1)


def lambda_d(u):
    """Debye length."""  # not checked yet, caution Te in eV
    [ne, ng, Te, Tg] = unnormalize(u)
    return np.sqrt((epsilon_0 * Te)/(ne * e))


def sheath_thickness(u):
    """Child Langmuir Sheath thickness, from Chabert 3.16."""
    [ne, ng, Te, Tg] = unnormalize(u)
    return (2**(1/2)/3) * lambda_d(u) * ((2*thruster.Vgrids)/Te)**(3/4)


def nu_m(u):
    """"Collision frequency."""
    [ne, ng, Te, Tg] = unnormalize(u)
    return ng * Kel(Te)


def nu_stoch(u):    # stochastic heating :
    """Stochastic collision frequency."""
    [ne, ng, Te, Tg] = unnormalize(u)
    v_mean_e = ((8 * e * Te)/(pi * m_e))**(0.5)
    wpe = np.sqrt(ne * e**2 / (m_e * epsilon_0))
    # delta = c / wpe
    delta_anomalous = ((c**2 * v_mean_e)/(thruster.w * wpe**2))**(1/3)
    # nouvelle formule ?
    return v_mean_e / (4 * delta_anomalous)  # dans le Chabert.


def R_ind(u):
    """Resistive part of plasma impedance, Ohm."""
    [ne, ng, Te, Tg] = unnormalize(u)

    R = thruster.R
    w = thruster.w

    wpe = np.sqrt(ne * e**2 / (m_e * epsilon_0))
    nu_eff = nu_m(u) + nu_stoch(u) * simu.stoch_heat

    eps_p = 1 - (wpe**2 / (w * (w - 1j*nu_eff)))
    k2 = (w/c) * np.sqrt(eps_p)

    Rind_1 = (2 * pi * thruster.N**2)/(thruster.L * w * epsilon_0)
    Rind_2 = 1j * k2 * R * jv(1, k2 * R)
    Rind_3 = eps_p * jv(0, k2 * R)
    # verifier pour nu_m = 0
    return Rind_1 * np.real(Rind_2 / Rind_3)


def I_coil(u):
    """Current in the coil (A)."""
    return np.sqrt(2 * simu.Pabs * thruster.volume() / R_ind(u))


def RFpower(u):
    """RF power."""
    return 0.5 * (R_ind(u) + thruster.R_coil) * I_coil(u)**2


def J_i(u):
    """Ion current, en A/m2."""
    [ne, ng, Te, Tg] = unnormalize(u)
    return e * ne * hl(u) * uB(u)


def Jcl(u):
    """Child langmuir ion current limit, en A/m2."""  # uB changed for v_beam
    [ne, ng, Te, Tg] = unnormalize(u)
    M = prop.m1
    Vgrids = thruster.Vgrids
    s = thruster.grid_dist
    return (4/9) * epsilon_0 * np.sqrt(2 * e / M) * Vgrids**(3/2) * s**(-2)

def thrustpower(u):
    """Power used to accelerate: ion current times acceleration voltage."""
    [ne, ng, Te, Tg] = unnormalize(u)
    return J_i(u) * thruster.open_area_ions() * thruster.Vgrids


def totalpower(u):
    """Sum of RF power + acceleration power."""
    [ne, ng, Te, Tg] = unnormalize(u)
    return RFpower(u) + thrustpower(u)


def thrust_neutrals(u):
    """Thrust from the neutrals, mN."""
    [ne, ng, Te, Tg] = unnormalize(u)
    M = prop.m1
    v_mean_ion = np.sqrt((8 * e * Tg)/(pi * prop.m1))
    v_mean_gas = v_mean_ion
    GammaG = 0.25 * ng * v_mean_gas  # neutral flux through grids
    Ag = thruster.open_area_neutrals()
    return GammaG * M * v_mean_gas * Ag * 1E3


def thrust_ions(u):
    """Thrust from the ions, mN."""
    [ne, ng, Te, Tg] = unnormalize(u)
    M = prop.m1
    v_beam = np.sqrt(2 * e * thruster.Vgrids / M)
    Gammai = ne * hl(u) * uB(u)
    Ai = thruster.open_area_ions()
    return Gammai * M * v_beam * Ai * 1E3


def thrust(u):
    """Total thrust, mN."""
    return thrust_ions(u) + thrust_neutrals(u)


def zeta(u):
    """Inductive coupling transfer efficiency."""
    return R_ind(u) / (R_ind(u) + thruster.R_coil)


def gamma(u):
    """Thrust power efficiency."""
    [ne, ng, Te, Tg] = unnormalize(u)
    v_beam = np.sqrt(2 * e * thruster.Vgrids / prop.m1)
    ion_thrust_power = 0.5 * v_beam * (thrust_ions(u)*1E-3)  # mN to N
    v_mean_gas = np.sqrt((8 * e * Tg)/(pi * prop.m1))
    neutral_thrust_power = 0.5 * v_mean_gas * (thrust_neutrals(u)*1E-3)  # mN to N

    num = ion_thrust_power + neutral_thrust_power
    return (num)/(num + RFpower(u))


def ksi(u):
    """Thrust efficiency, in mN per kW."""
    return (thrust_ions(u)+thrust_neutrals(u))/(RFpower(u)*1e-3)


def eta(u):
    """Mass utilization efficiency."""
    [ne, ng, Te, Tg] = unnormalize(u)
    Gammai = ne * hl(u) * uB(u)
    return (Gammai * thruster.open_area_ions())/simu.Q0
    # return (Gammai * thruster.open_area_ions())/((ne+ng)*thruster.volume())  # neutral depletion, not working.


def thruster_eff(u):
    """"Thruster efficiency (Goebel 2.5-4).

    Real thrust means we consider zero beam divergence.
    """
    mfr = simu.Q0 * prop.m1
    return (thrust(u)*1e-3)**2 / (2 * mfr * totalpower(u))


def discharge_loss(u):
    """"Discharge loss (Goebel 2.5-2)."""
    return (simu.Pabs*thruster.volume()) / (J_i(u) * thruster.open_area_ions())


# ############################### #
# ## Formula for simplified GM ## #
# ############################### #

def deff(u):
    """Effective distance of the thruster."""
    [ne, ng, Te, Tg] = unnormalize(u)
    Aeff = (2 * hr(u) * pi * thruster.L * thruster.R
            + 2 * hl(u) * pi * thruster.R**2)
    return thruster.volume() / Aeff


# def power_density(u):
#     """Power density."""
#     [ne, ng, Te, Tg] = unnormalize(u)
#     return simu.Peasy / thruster.volume()


def part_balance(u):
    """Particle balance. Must equal 0. Positive if particle creation."""
    [ne, ng, Te, Tg] = unnormalize(u)
    Aeff = (2 * hr(u) * pi * thruster.L * thruster.R
            + 2 * hl(u) * pi * thruster.R**2)
    deff = thruster.volume() / Aeff
    return K_iz(Te) - uB(u) / (ng * deff)


def power_balance(u):
    """Power gain minus losses."""
    [ne, ng, Te, Tg] = unnormalize(u)

    Vs = Te * 0.5 * np.log(prop.m1 / (2 * pi * m_e))  # sheat potential at a floating wall
    # print(np.log(prop.m1 / (2 * pi * m_e)))
    Eki = 0.5 * Te + Vs  # ions: pre-sheath + sheath drop

    Eke = 2 * Te  # Liebermann 2.4.11

    Eiz = prop.E_iz
    Eexc = prop.E_exc

    Ec_exc = (K_exc(Te) / K_iz(Te)) * Eexc

    Ecoll = (3 * m_e * Te * Kel(Te)) / (prop.m1 * K_iz(Te))

    Ec = Ecoll + Ec_exc + Eiz

    Et = Eki + Eke + Ec  # energy total_area

    Aeff = (2 * hr(u) * pi * thruster.L * thruster.R
            + 2 * hl(u) * pi * thruster.R**2)

    # print(type(Aeff))
    return (simu.Pabs * thruster.volume()) - e * ne * uB(u) * Aeff * Et

# ############################### #
# ## Parametric studies        ## #
# ############################### #

def power_range(Pmax, res, params, thruster, propel):
    Pmin = params.Pabs
    
    # output file header
    date = str(dt.datetime.now())+"\n"
    thr = "Thruster: R[m]="+str(thruster.R)+' ; L[m]='+str(thruster.L)+'\n'
    param = "Parameters Q0[s-1]="+str(params.Q0)+'\n'
    prange = "Power range [W/m3]="+str(Pmin)+' - '+str(Pmax)+'\n'    
    ##################

    Nparam = 4
    dPow   = (Pmax-Pmin)/Nparam
    paramtab = []    
    for i in range(Nparam):
        paramtab.append(np.append([params.Pabs], res.values.y[-1]))
        print("P[W/m3] = ", params.Pabs)
        params.Pabs += dPow
        res = solve_cvode(propel, thruster, params, res.values.y[-1])
    np.savetxt('data/paramP.txt', paramtab, header = date+thr+param+prange)     
    return paramtab
               



#######################################
# ## Plots                         ## #
#######################################

def multi_plot(results):
    fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(3, 2, sharex='col')
    ax1.plot(results.values.t, results.values.y[:,0])
    ax1.set_yscale('Log')
    ax1.set_ylabel('n_e')
    
    ax2.plot(results.values.t, results.values.y[:,1])
    ax2.set_yscale('Log')
    ax2.set_ylabel('n_g')
    
    ax3.plot(results.values.t, results.values.y[:,2])
    ax3.set_yscale('Log')
    ax3.set_ylabel('n_Xe+')
    
    ax4.plot(results.values.t, results.values.y[:,3])
    ax4.set_ylabel('T_e')
    
    ax5.plot(results.values.t, results.values.y[:,4])
    ax5.set_ylabel('Tg')
    plt.tight_layout()
    plt.show() 

def power_plot(fil, thruster):
    fil = 'data/'+fil
    data = np.genfromtxt(fil)

    data[:,0] *= thruster.volume()

    fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(3, 2, sharex='col')
    ax1.plot(data[:,0], data[:,1])
    #ax1.set_yscale('Log')
    ax1.set_ylabel('n_e')

    ax2.plot(data[:,0], data[:,2])
    #ax2.set_yscale('Log')
    ax2.set_ylabel('n_g')

    ax3.plot(data[:,0], data[:,3])
    #ax3.set_yscale('Log')
    ax3.set_ylabel('n_Xe+')

    ax4.plot(data[:,0], data[:,4])
    ax4.set_ylabel('T_e')

    ax5.plot(data[:,0], data[:,5])
    ax5.set_ylabel('Tg')
    ax5.set_xlabel('Power [W]')

    ax6.set_xlabel('Power [W]')  

    plt.tight_layout()
    plt.show()


