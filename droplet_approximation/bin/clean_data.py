import os

import sys

import pandas as pd
import numpy as np

from droplet_approximation.data import *

number_epochs = 7

def main( argv ):
    if len( argv ) != 3:
        print("Usage: <input_path> <output_path>")
        sys.exit(1)

    be_dump_path = argv[1]
    data_save_path = argv[2]

    df = read_NTLP_file(be_dump_path)

    temperature_range = [df["input temperature"].quantile(0.00005),df["input temperature"].quantile(0.99995)]
    radius_range = [1.99e-7, df["input radius"].quantile(0.99995)]

    data = df.loc[(df["output temperature"] > temperature_range[0]) & (df["output temperature"] < temperature_range[1]) & (df["input temperature"] > temperature_range[0]) & (df["input temperature"] < temperature_range[1]) & (df["output radius"] < radius_range[1]) & (df["air density"] == 1.0) & (df["relative humidity"] < 5.0)][["input radius", "input temperature", "salinity", "air temperature", "relative humidity", "air density", "integration time", "output radius", "output temperature"]].to_numpy(dtype=np.float32)

    data.tofile(data_save_path)

if __name__ == "__main__":
    sys.exit( main( sys.argv ) )
