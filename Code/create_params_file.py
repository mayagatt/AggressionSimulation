import numpy as np
import json


def different_delta_L(min, max, step_size, params):
    """
    :param min: min value of delta_L
    :param max: max value of delta_L
    :param step_size: step size of delta_L
    :param params: all parameters needed for a successful run. Only the delta_L entry will be changed
    :return: str, the name of a directory containing all the generated files. Each file contains the same params except
    for delta_L.
    """
    delta_Ls = np.arange(min, max, step_size)
    dir_name = 'setup_files_delta_Ls\\'
    for delta_L in delta_Ls:
        params[
            'deltaL'] = delta_L  # beware, it is changing the value of the original dict. it works only when changing one param at a time
        file_name = dir_name + params_to_file_name(params) + '.json'
        with open(file_name, 'w') as fp:
            json.dump({'params': params}, fp)

    return dir_name


def params_to_file_name(params):
    file_name = ''
    for key, val in params.items():
        file_name = file_name + key + '_' + str(val) + '__'
    return file_name[:-2]


def change_deltaL_c0(setup_file_path, params, min_deltaL, max_deltaL, tot_deltaL, min_log10c0, max_log10c0,
                     tot_log10c0):
    deltaLs = np.linspace(min_deltaL, max_deltaL, tot_deltaL)
    log10c0s = np.linspace(min_log10c0, max_log10c0, tot_log10c0)
    for deltaL in deltaLs:
        for log10c0 in log10c0s:


if __name__ == "__main__":
    params = {'K': 1, 'K_L': 1, 'E': 1, 'DeltaE': 0.1, 'deltaL': 1, 'gtype': 'hill', 'h': 1, 'tf': 10,
              'is_restart_antibiotics': True}
    min_deltaL = 0
    max_deltaL = 1
    d_deltaL = 0.5
