import numpy as np
import matplotlib.pyplot as plt


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
    plt.plot(10 ** params['log10c0'] - out_df['rhoV'] - out_df['rhoA'], label='wasted biomass', color='grey')
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
