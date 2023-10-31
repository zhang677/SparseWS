import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy

# read xlsx
df = pd.read_excel('/home/zgh23/code/SparseWS/data/ablation/CSC_CSR_T/memory/speedup-outputnnz.xlsx', header=None)
# x equals df['memory']
x = np.array(df[1].values.tolist()[1:])
# print(x)
# y equals df['speedup']
y = np.array(df[2].values.tolist()[1:])
# Split x into intervals of 0.1
# intervals = np.array([0.8,0.9,0.91,0.92,0.93,0.94,0.95,1])
# intervals = np.array([0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0001])
# intervals = np.array([0.0,0.2,0.4,0.6,0.8,1.0001])
intervals = np.array([0, 1e4, 1e5, 1e6, 1e7, 1e8])
y_values = []
x_values = []
# print(intervals)
# Calculate average, minimum, and maximum for each interval of x
for i in range(len(intervals)-1):
    #indices = np.where((np.array(x) >= interval - 0.05) & (np.array(x) < interval + 0.05))
    indices = np.where((x >= intervals[i]) & (x < intervals[i+1]))
    print(indices)
    interval_y = y[indices]
    # interval_x = x[indices]
    # x_values.append(scipy.stats.mstats.gmean(interval_x))
    # x_values.append(intervals[i] + 0.1)
    x_values.append(0.5*(intervals[i+1]-intervals[i]))
    y_values.append((scipy.stats.mstats.gmean(interval_y), np.min(interval_y), np.max(interval_y)))

# Extract x and y values for the plot
# x_plot = [(interval + 0.05) for interval in intervals]
# x_plot = [0.9,0.91,0.92,0.93,0.94,0.95,0.96]
# x_plot = x_values
x_plot = [4, 5, 6, 7, 8]
y_avg = [value[0] for value in y_values]
y_min = [value[1] for value in y_values]
y_max = [value[2] for value in y_values]

# Plotting the figure
# plt.bar(x_plot, y_max, width=0.1, align='center', alpha=0.5, color='r', label='Max')
# plt.bar(x_plot, y_min, width=0.1, align='center', alpha=0.5, color='b', label='Min')
plt.plot(x_plot, y_max, marker='o', linestyle='-', color='r', label='Max')
plt.plot(x_plot, y_min, marker='o', linestyle='-', color='b', label='Min')
plt.plot(x_plot, y_avg, marker='o', linestyle='-', color='g', label='Average')
# plt.xlabel('memory saving')
# plt.ylabel('speedup')
# plt.title('Speedup vs. Memory Saving')
# plt.xlabel('Collision Ratio')
plt.xlabel('Output NNZ')
plt.ylabel('Fusion SpeedUp')
# plt.title('Speedup vs. Collision Ratio')
plt.title('Speedup vs. Output NNZ (log10)')
plt.legend()
plt.grid(True)
plt.savefig('figure-outputnnz-2.png')