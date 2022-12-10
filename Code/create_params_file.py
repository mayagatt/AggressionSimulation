import numpy as np
import pandas as pd


def change_two_params(setup_file_path, params, param1, param1_arr, param2, param2_arr):
    """
    writes to a csv file a DataFrame containing the same value in each row for all the params, accepts param1 and param2.
    For every value of param1 we  match all the values of param2, e.g.:
    param1,param2,all_other_params
    val_1,val_1,val
    val_1,val_2,val
    ...
    val_1,val_n,val
    val_2,val_1,val
    val_2,val_2,val
    ...
    val_2,val_n,val
    ...
    :param setup_file_path:
    :param params: all parameters necessary for a single run
    :param param1: 1st parameters to be changed.
    :param param1_arr: the different values param1 should have
    :param param2: 2nd parameters to be changed.
    :param param2_arr: the different values param2 should have
    """
    df = pd.DataFrame([params])
    tot_len = len(param1_arr) * len(param2_arr)
    df = (pd.DataFrame(np.repeat(df.values, tot_len, axis=0), columns=df.columns))  # repeat same row in all df
    df[param1] = np.repeat(param1_arr, len(param2_arr))  # replace param1
    df[param2] = np.round(np.tile(param2_arr, len(param1_arr)), 1)  # replace param2, todo: beware of the rounding
    df.to_csv(setup_file_path)


if __name__ == "__main__":
    import argparse
    PARSER = argparse.ArgumentParser()
    PARSER.add_argument('-o', type=str, required=True, help='Output file')
    parser = PARSER.parse_args()
    out_file = parser.o

    params = {'log10c0': 1, 'x': 0, 'rhoA': 0.5, 'rhoV': 0.5, 'K': 1, 'K_L': 1, 'E': 1, 'DeltaE': 0.1, 'deltaL': 1,
              'gtype': 'hill', 'h': 1, 'tf': 10,
              'is_restart_antibiotics': True}
    min_deltaL = 0
    max_deltaL = 1
    tot_deltaL = 101
    min_log10c0 = -3
    max_log10c0 = 3
    tot_log10c0 = 31
    # path = ".\Data\Raw\\v_a_vs_deltaL_and_c0\\to_run.tsv"
    deltaLs = np.linspace(min_deltaL, max_deltaL, tot_deltaL)
    log10c0s = np.linspace(min_log10c0, max_log10c0, tot_log10c0)
    change_two_params(out_file, params, param1='deltaL', param2='log10c0', param1_arr=deltaLs, param2_arr=log10c0s)
