import numpy as np
import pandas as pd
from scipy.integrate import odeint
from scipy.integrate import trapezoid
import matplotlib.pyplot as plt
import copy


def odefun_creator(params):
    '''
    DEfining the 4 diff equation for the in-batch development of the nutrient c, toxin x, aggressor rho_A, victim rho_V
    :param t: time
    :param y: vector (c, x, rho_A, rho_V)
    :param params: must have: K (for g[c]), K_L (for g[x]), E (max consumption), DeltaE (resources allocated for toxin
    production), deltaL (not sure todo??) , gtype (str, 'hill' or 'linear'), h (hill func param)
    :return: dydt [= f(y)]
    '''
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
    match params['gtype']:  # set gfun as the growth function
        case 'hill':
            def gfun(x, K):
                return (x ** params['h']) / ((K ** params['h']) + (x ** params['h']))
        case 'linear':
            def gfun(x, K):
                return x / K
        case _:
            print('invalid growth function type')

            def gfun(x, K):
                return None
    return gfun


def two_species_batches(y_initial, params):
    '''

    :param params:
    :param y_initialt:
    :return:
    '''
    t = np.arange(0, params['tf'], 0.01)  # todo: change to adaptive mesh, save t in sol
    odefun = odefun_creator(params)
    g_fun = get_g_fun(params)

    rhoAs = []  # saving rhoA(0) for batch b=0,1,..., until steady state
    rhoVs = []
    xs = []
    dilution_factors = []
    nutrient_integrals = []
    toxin_integrals = []

    rhoA_ts = []  # array of arrays, including time evolution within batches [...[rhoA(t=0)(b),...,rhoA(tf)(b)],...]
    rhoV_ts = []
    x_ts = []

    eps = 1e-8
    cond = True
    y0 = copy.deepcopy(y_initial)
    i = 0
    while cond:  # cond
        i += 1
        sol = run_solver(odefun, t, y0)
        dilution_factor = (y_initial['rhoA'] + y_initial['rhoV']) / (
                y_initial['rhoA'] + y_initial['rhoV'] + y_initial['c'])

        cond = max(y0['rhoV'] - (sol['rhoV'][-1] * dilution_factor),
                   y0['rhoA'] - (sol['rhoA'][-1] * dilution_factor)) > eps

        # update y0
        y0['x'] = 0 if params['is_restart_antibiotics'] else sol['x'][-1]
        y0['rhoV'] = sol['rhoV'][-1] * dilution_factor  # rhoV(tf)
        y0['rhoA'] = sol['rhoA'][-1] * dilution_factor

        # integrate g(c) per batch
        temp_integral = trapezoid(y=g_fun(sol['c'], K=params['K']), x=t)
        nutrient_integrals.append(temp_integral)

        # integrate k(x) per batch
        toxin_integral = trapezoid(y=g_fun(sol['x'], K=params['K_L']), x=t)
        toxin_integrals.append(toxin_integral)

        # save values at tf
        rhoAs.append(y0['rhoA'])
        rhoVs.append(y0['rhoV'])
        xs.append(sol['x'][-1])  # note! this is x in tf but rho_A,V  are in t=0 of batch b+1
        dilution_factors.append(dilution_factor)

        # save all values of the batch
        rhoA_ts.append(sol['rhoA'])
        rhoV_ts.append(sol['rhoV'])
        x_ts.append(sol['x'])
    # plot_batch(params, sol, t, y0)
    tf_data_base = pd.DataFrame(
        {'rhoA': np.array(rhoAs), 'rhoV': np.array(rhoVs), 'dilution_factor': np.array(dilution_factors), 'x': xs,
         'nutrient_integral': nutrient_integrals, 'toxin_integral': toxin_integrals})
    full_time_series = pd.DataFrame({'rhoA_ts': rhoA_ts, 'rhoV_ts': rhoV_ts, 'x_ts': x_ts})
    return tf_data_base, full_time_series


def run_solver(odefun, t, y0):
    sol = odeint(odefun, [y0['c'], y0['x'], y0['rhoA'], y0['rhoV']], t)
    c_t = sol[:, 0]
    x_t = sol[:, 1]
    rhoA_t = sol[:, 2]
    rhoV_t = sol[:, 3]
    return {'c': c_t, 'x': x_t, 'rhoA': rhoA_t, 'rhoV': rhoV_t}


def plot_batch(i, out_df):
    tf = len(out_df['rhoV_ts'][0])
    x = np.arange(i * tf, (i + 1) * tf)
    plt.plot(x, out_df['rhoV_ts'][i], label='rho_V', c='b', linewidth=0.8)
    plt.plot(x, out_df['rhoA_ts'][i], label='rho_A', c='r', linewidth=0.8)
    plt.xlabel('t')
    plt.ylabel('Species Density')
    plt.legend()
    plt.savefig('density_evolution_over_80_batches')
    plt.show()


def plot_species_density_over_batch_number(out_df, params):
    plt.plot(out_df['rhoA'], label='rhoA', color='r')
    plt.plot(out_df['rhoV'], label='rhoV', color='b')
    plt.plot(y0['c'] - out_df['rhoV'] - out_df['rhoA'], label='wasted biomass', color='grey')
    # plt.plot(out_df['dilution_factor'], label='dilut')
    plt.xlabel('Batch number')
    plt.ylabel('Species Density')
    plt.title(f'delta_L={params["deltaL"]}')
    plt.legend()
    plt.show()


def plot_species_density_over_time(out_df):
    plt.figure(dpi=300)
    # # print(type(out_df['rhoA_ts']))
    # rhoA_ts_bs = np.concatenate(out_df['rhoA_ts'].values)
    # print(rhoA_ts_bs)
    # plt.scatter(np.arange(len(rhoA_ts_bs)), rhoA_ts_bs, label='rhoA_ts')
    # plt.scatter(np.arange(len(rhoA_ts_bs)), np.arange(len(rhoA_ts_bs)), s=0.05)
    tf = len(out_df['rhoV_ts'][0])
    for i in range(len(out_df['rhoV_ts'])):
        x = np.arange(i * tf, (i + 1) * tf)
        plt.plot(x, out_df['rhoV_ts'][i], c='b', linewidth=0.8)
        plt.plot(x, out_df['rhoA_ts'][i], c='r', linewidth=0.8)
    plt.xlabel('t')
    plt.ylabel('Species Density')
    # plt.legend()
    plt.grid()
    plt.savefig('density_evolution_over_80_batches')
    plt.show()


# def search_critical_delta_l(delta_min, delta_max):
#     y0 = {'c': 1, 'x': 0, 'rhoA': 0.5, 'rhoV': 0.5}  # initial cond: c0=1, x0=0, rho_A=1/2, rho_V=1/2
#     params = {'K': 1, 'K_L': 1, 'E': 1, 'DeltaE': 0.2, 'deltaL': 0.025, 'gtype': 'hill', 'h': 1, 'tf': 10,
#               'is_restart_antibiotics': True}
#     # delta_Ls = np.arange(delta_min, delta_max, d_delta)
#     params['deltaL'] = delta_L
#     out_df, full_time_series = two_species_batches(y0, params)
#     if out_df['rhoA'] - out_df['rhoV'] > 0:  # rhoA wins
#
#     for delta_L in delta_Ls:
#         plot_species_density_over_time(full_time_series)


if __name__ == "__main__":
    import argparse
    import json

    PARSER = argparse.ArgumentParser()
    PARSER.add_argument('-o', type=str, required=True, help='Output file')
    PARSER.add_argument('-d', type=str, required=True, help='Default setup file')
    PARSER.add_argument('-i', type=str, required=True, help='Setup file')
    s = PARSER.parse_args()
    out_file = s.o
    default_setup_file = s.d
    setup_file = s.i
    print(setup_file)
    with open(default_setup_file, "r") as in_file:
        default_setup_data = json.load(in_file)

    with open(setup_file, "r") as in_file:
        setup_data = json.load(in_file)

    # merge params with priority to setup_data
    params = default_setup_data["y0"].copy()
    y0 = default_setup_data["params"].copy()
    for key, val in setup_data["params"].items():
        params[key] = val
    for key, val in setup_data["y0"].items():
        y0[key] = val

    print(out_file)

    out_df, full_time_series = two_species_batches(default_setup_data["y0"], default_setup_data["params"])
    out_df.to_csv(out_file)

# if __name__ == "__main__":
#     # later: read params from file, run solver, save to data file
#
#     # delta_L = 0.221660525  # pretty critical value
#     integral_ratios = []
#     dE_dL_ratios = []
#
#     # delta_L = 0.142857
#     # delta_L = 0.2696325
#     delta_L = 0.25
#     y0 = {'c': 1, 'x': 0, 'rhoA': 0.5, 'rhoV': 0.5}  # initial cond: c0=1, x0=0, rho_A=1/2, rho_V=1/2
#     params = {'K': 1, 'K_L': 1, 'E': 1, 'DeltaE': 0.1, 'deltaL': delta_L, 'gtype': 'hill', 'h': 1, 'tf': 10,
#               'is_restart_antibiotics': True}
#     out_df, full_time_series = two_species_batches(y0, params)
#     # plot_species_density_over_time(full_time_series)
#     plot_species_density_over_batch_number(out_df, params)
#     print(out_df['rhoA'][58])
#
#     # shifts = np.linspace(-1e-06, 5*1e-06, 9)
#     # print(shifts)
#     # for ii, temp in enumerate(shifts):
#     #     print(ii)
#     #     delta_L = 0.2696325 + temp  # pretty critical value for
#     #     # delta_L = 0.25
#     #     y0 = {'c': 0.1, 'x': 0, 'rhoA': 0.5, 'rhoV': 0.5}  # initial cond: c0=1, x0=0, rho_A=1/2, rho_V=1/2
#     #     params = {'K': 1, 'K_L': 1, 'E': 1, 'DeltaE': 0.2, 'deltaL': delta_L, 'gtype': 'hill', 'h': 1, 'tf': 10,
#     #               'is_restart_antibiotics': True}
#     #     out_df, full_time_series = two_species_batches(y0, params)
#     #
#     #     toxin_integral_ss = np.array(out_df['toxin_integral'])[-1]
#     #     nut_integral_ss = np.array(out_df['nutrient_integral'])[-1]
#     #     integral_ratios.append(toxin_integral_ss / nut_integral_ss)
#     #     dE_dL_ratios.append(params['DeltaE'] / params['deltaL'])
#     #     plot_species_density_over_batch_number(out_df,params)
#
#     # plot_species_density_over_time(full_time_series)
#     # plot_species_density_over_batch_number(out_df)
