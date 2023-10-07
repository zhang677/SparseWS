import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os
import sys
import argparse
import re

def plot_memory(csv_file):
    """
    csv_file is in format time,memory,stage
    Plot memory usage. X-axis is time, Y-axis is memory usage. Each record is a point.
    """
    print("plotting memory usage")
    df = pd.read_csv(csv_file, dtype = {'time': float, 'memory': int, 'stage': str})
    sns.set(style="darkgrid")
    sns.lineplot(x="time", y="memory", data=df)
    # save to pdf
    plt.savefig(os.path.splitext(csv_file)[0] + '.pdf')

def analyse_stages(csv_file):
    df = pd.read_csv(csv_file)
    # df is in format time,memory,stage
    # There are 10 types of stages:
    # 1. 'Init'
    # 2. r'(T(\d+)\+A(\d+))'
    # 3. r'T(\d+)'
    # 4. r'FC_\d+_\d+_\d+/\d+'
    # 5. 'FM'
    # 6. 'MA(\d+)'
    # 7. 'MR(\d+)'
    # 8. 'M(\d+)'
    # 9. 'FF'
    # 10. 'OR'

def extract_main_stages(csv_file):
    print("extracting main stages")
    df = pd.read_csv(csv_file, dtype = {'time': float, 'memory': int, 'stage': str})
    # Extract stages S1, FC_0_\d+_\d+/\d+, FM, MA, MR, M, FF, OR from df and store them in a new dataframe df2.
    # df2 is in format time,memory,stage
    df2 = pd.DataFrame(columns=['time', 'memory', 'stage'])
    lst = []
    for index, row in df.iterrows():
        if row['stage'] == 'Init':
            lst.append(row)
        elif re.match(r'FC_0_\d+_\d+/\d+', row['stage']):
            lst.append(row)
        elif row['stage'] == 'FM':
            lst.append(row)
        elif re.match(r'MA(\d+)', row['stage']):
            lst.append(row)
        elif re.match(r'MR(\d+)', row['stage']):
            lst.append(row)
        elif re.match(r'M(\d+)', row['stage']):
            lst.append(row)
        elif row['stage'] == 'FF':
            lst.append(row)
        elif row['stage'] == 'OR':
            lst.append(row)
        else:
            continue
    df_extended = pd.DataFrame(lst, columns=['time', 'memory', 'stage'])
    df2 = pd.concat([df2, df_extended])
    df2.to_csv(os.path.splitext(csv_file)[0] + '_succinct.csv', index=False)

def count_firstmode_collisions(firstmode, csv_file):
    df = pd.read_csv(csv_file, dtype = {'time': float, 'memory': int, 'stage': str})
    # [std, ...] of each cycle to measure the collision of first mode. The higher the std, the more collisions.
    std_list = []
    in_cycle = False
    cycle_start = r'FC_0_\d+_\d+/\d+'
    cycle_dur = r'FC_\d+_\d+_\d+/\d+'
    cycle_end = {"FM", "MA"}
    for index, row in df.iterrows():
        match_start = re.match(cycle_start, row['stage'])
        if match_start:
            in_cycle = True
            lst = np.zeros(firstmode)
            lst[int(match.group(1))] = int(match.group(3))
        elif row['stage'] in cycle_end:
            in_cycle = False
            std_list.append(np.std(lst))
        elif in_cycle:
            match = re.match(cycle_dur, row['stage'])
            lst[int(match.group(1))] = int(match.group(3))
        else:
            continue
    # store std_list to a csv file
    with open(os.path.splitext(csv_file)[0] + '_occupancy_std.txt', 'w') as f:
        for item in std_list:
            f.write("%s\n" % item)



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--file', default="~/code/SparseWS/data/results/CSC_CSR_T/profile-coord-bucket_0.csv", type=str, help='csv file')
    parser.add_argument('--firstmode', default=1916, type=int, help='Number of first mode') # first mode for CAG_mat1916 is 1916
    args = parser.parse_args()
    plot_memory(args.file)
    #extract_main_stages(args.file)
    # count_firstmode_collisions(args.firstmode, args.file)



