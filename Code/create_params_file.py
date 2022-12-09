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
    df = pd.DataFrame()
    for key, val in params.items():
        if key == param1:
            array = np.repeat(param1_arr, len(param2_arr))
        elif key == param2:
            array = np.tile(param2_arr, len(param1_arr))
        else:
            array = np.tile(val, len(param1_arr) * len(param2_arr))
        df[key] = array
    df.to_csv(setup_file_path)


if __name__ == "__main__":
    params = {'log10c0': 1, 'x': 0, 'rhoA': 0.5, 'rhoV': 0.5, 'K': 1, 'K_L': 1, 'E': 1, 'DeltaE': 0.1, 'deltaL': 1, 'gtype': 'hill', 'h': 1, 'tf': 10,
              'is_restart_antibiotics': True}
    min_deltaL = 0
    max_deltaL = 1
    tot_deltaL = 2
    min_log10c0 = 0
    max_log10c0 = 3
    tot_log10c0 = 4
    path = "..\Data\Raw\\v_a_vs_deltaL_and_c0\\to_run.tsv"
    # change_deltaL_c0(path, params, min_deltaL, max_deltaL, tot_deltaL, min_log10c0, max_log10c0, tot_log10c0)
    deltaLs = np.linspace(min_deltaL, max_deltaL, tot_deltaL)
    log10c0s = np.linspace(min_log10c0, max_log10c0, tot_log10c0)
    change_two_params(path, params, param1='deltaL', param2='log10c0', param1_arr=deltaLs, param2_arr=log10c0s)

