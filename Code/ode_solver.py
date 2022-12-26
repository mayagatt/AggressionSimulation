import numpy as np
import pandas as pd
from scipy.integrate import odeint
from scipy.integrate import trapezoid


def odefun_creator(params):
    """
    DEfining the 4 diff equation for the in-batch development of the nutrient c, toxin x, aggressor rho_A, victim rho_V
    :param t: time
    :param y: vector (c, x, rho_A, rho_V)
    :param params: must have: K (for g[c]), K_L (for g[x]), E (max consumption), DeltaE (resources allocated for toxin
    production), deltaL (not sure todo??) , gtype (str, 'hill' or 'linear'), h (hill func param)
    :return: dydt [= f(y)]
    """
    gfun = get_g_fun(params)

    def odefun(y, t):
        dydt = np.zeros_like(y)

        # c = y[0]
        # x = y[1]
        # rho_A = y[2]
        # rho_V = y[3]

        gc = gfun(y[0], params['K'])  # g[c]

        # dc/dt = -E*g[c](rho_A + rho_B)
        dydt[0] = -params['E'] * gc * y[2] - params['E'] * gc * y[3]  # assuming params['h']=1 => gfun = monod
        # dx/dt = rho_A * deltaE * g[c]
        dydt[1] = y[2] * params['DeltaE'] * gc
        # d rho_A/dt = rho_A * deltaE * g[c]
        dydt[2] = y[2] * (params['E'] - params['DeltaE']) * gc
        # d rho_V/dt = rho_V E g[c] - rho_V deltaL k[x]; k[x]=x/(x+K_L)
        dydt[3] = y[3] * params['E'] * gc - y[3] * params['deltaL'] * gfun(y[1], params['K_L'])

        return dydt

    return odefun


def get_g_fun(params):
    '''
    :param params:
    :return: function, g(c)
    '''
    if params['gtype'] == 'hill':
        def gfun(x, K):
            return (x ** params['h']) / ((K ** params['h']) + (x ** params['h']))
    elif params['gtype'] == 'linear':
        def gfun(x, K):
            return x / K
    else:
        def gfun(x, K):
            return None
    return gfun


def two_species_batches(params):
    '''

    :param params: all params necessary foe a run, including initial conditions
    :return:
    '''
    t = np.arange(0, params['tf'], 0.01)  # todo: change to adaptive mesh, save t in sol
    odefun = odefun_creator(params)
    g_fun = get_g_fun(params)

    rhoAs_t0 = []  # saving rhoA(0) for batch b=1,..., until steady state
    rhoVs_t0 = []
    xs_tf = []
    rhoAs_tf = []  # saving rhoA(0) for batch b=0,1,..., until steady state
    rhoVs_tf = []
    rhoV_deads_ft = []
    dilution_factors = []
    nutrient_integrals = []
    toxin_integrals = []

    rhoA_ts = []  # array of arrays, including time evolution within batches [...[rhoA(t=0)(b),...,rhoA(tf)(b)],...]
    rhoV_ts = []
    x_ts = []

    eps = 1e-8
    cond = True
    y0 = {'log10c0': params['log10c0'], 'x': params['x'], 'rhoA': params['rhoA'], 'rhoV': params['rhoV']}
    i = 0
    while cond:  # cond
        i += 1
        sol = run_solver(odefun, t, y0)
        dilution_factor = (params['rhoA'] + params['rhoV']) / (
                params['rhoA'] + params['rhoV'] + 10 ** params['log10c0'])

        # continue simulation while both species are still changing by more than eps between batches
        cond = max(abs(y0['rhoV'] - (sol['rhoV'][-1] * dilution_factor)),
                   abs(y0['rhoA'] - (sol['rhoA'][-1] * dilution_factor))) > eps

        # save y0 at t=0
        rhoAs_t0.append(y0['rhoA'])
        rhoVs_t0.append(y0['rhoV'])
        dilution_factors.append(dilution_factor)

        # save y0 at t=tf
        rhoAs_tf.append(sol['rhoV'][-1])
        rhoVs_tf.append(sol['rhoV'][-1])
        xs_tf.append(sol['x'][-1])

        # update y0
        y0['x'] = 0 if params['is_restart_antibiotics'] else sol['x'][-1]
        y0['rhoV'] = sol['rhoV'][-1] * dilution_factor  # rhoV(tf)
        y0['rhoA'] = sol['rhoA'][-1] * dilution_factor

        # integrate g(c) per batch
        temp_integral = trapezoid(y=g_fun(10 ** sol['log10c0'], K=params['K']), x=t)
        nutrient_integrals.append(temp_integral)

        # integrate k(x) per batch
        k_x = g_fun(sol['x'], K=params['K_L'])
        toxin_integral = trapezoid(y=k_x, x=t)
        toxin_integrals.append(toxin_integral)

        # integrate rhoV dead per batch - deltaL * rhoV * k[x]
        dead_integral = trapezoid(y=params['deltaL'] * k_x * sol['rhoV'], x=t)
        rhoV_deads_ft.append(dead_integral)

        # save all values of the batch
        rhoA_ts.append(sol['rhoA'])
        rhoV_ts.append(sol['rhoV'])
        x_ts.append(sol['x'])

    tf_data_base = pd.DataFrame(
        {'rhoA_t0': np.array(rhoAs_t0), 'rhoV_t0': np.array(rhoVs_t0), 'rhoA_tf': np.array(rhoAs_tf),
         'rhoV_tf': np.array(rhoVs_tf), 'x_tf': xs_tf, 'dilution_factor': np.array(dilution_factors),
         'nutrient_integral': nutrient_integrals, 'toxin_integral': toxin_integrals})
    full_time_series = pd.DataFrame({'rhoA_ts': rhoA_ts, 'rhoV_ts': rhoV_ts, 'x_ts': x_ts})
    return tf_data_base, full_time_series


def run_solver(odefun, t, y0):
    sol = odeint(odefun, [10 ** y0['log10c0'], y0['x'], y0['rhoA'], y0['rhoV']], t)
    c_t = np.array(sol[:, 0])
    x_t = sol[:, 1]
    rhoA_t = sol[:, 2]
    rhoV_t = sol[:, 3]
    log10c0 = np.log10(np.where(c_t > 1e-10, c_t, 1e-10))
    return {'log10c0': log10c0, 'x': x_t, 'rhoA': rhoA_t, 'rhoV': rhoV_t}


if __name__ == "__main__":
    import argparse
    import json

    PARSER = argparse.ArgumentParser()
    PARSER.add_argument('-o', type=str, required=True, help='Output file')
    PARSER.add_argument('-d', type=str, required=True, help='Default params')
    PARSER.add_argument('-s', type=str, required=True, help='Setup file')  # every row is a single run
    PARSER.add_argument('-n', type=int, required=True, help='row in setup file to run')
    parser = PARSER.parse_args()
    out_file = parser.o
    default_params_file = parser.d
    setup_file = parser.s
    n = parser.n
    with open(default_params_file, "r") as in_file:
        default_params = json.load(in_file)
    params_df = pd.read_csv(setup_file, index_col=0)[n:n + 1]
    params = params_df.to_dict('records')[0]

    out_df, full_time_series = two_species_batches(params)
    params_df = (pd.DataFrame(np.repeat(params_df.values, out_df.shape[0], axis=0), columns=params_df.columns))
    out_df = pd.concat([out_df, params_df], axis=1)  # add params to out file, in every row
    out_df.to_pickle(out_file + '.pkl')
