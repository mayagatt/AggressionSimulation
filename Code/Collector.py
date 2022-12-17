import pandas as pd
import os

if __name__ == "__main__":
    import argparse

    PARSER = argparse.ArgumentParser()
    PARSER.add_argument('-d', type=str, required=True, help='Raw data dir')
    PARSER.add_argument('-o', type=str, required=True, help='Out file')  # without a suffix!
    parser = PARSER.parse_args()
    raw_data_dir = parser.d
    out_file = parser.o

    collected_df = pd.DataFrame()
    i = 0
    for file in os.listdir(raw_data_dir):
        filename = raw_data_dir + "/" + os.fsdecode(file)
        if filename.endswith(".pkl"):
            try:
                ss_data_df = pd.read_pickle(filename)[-2:-1]
                if i == 0:
                    collected_df = ss_data_df
                    i += 1
                else:
                    collected_df = pd.concat([collected_df, ss_data_df], ignore_index=True)
            except:
                print('no results for run')

    collected_df.to_pickle(out_file + '.pkl')
