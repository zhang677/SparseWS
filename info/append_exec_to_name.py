import argparse

parser = argparse.ArgumentParser(description='Output matrix names from csv to a txt file.')
parser.add_argument('--input_file', default="/home/zgh23/code/SparseWS/info/band_names_1001.txt", type=str, help='Read path')
parser.add_argument('--output_file', default="/home/zgh23/code/SparseWS/info/band_names_1001_exec.txt", type=str, help='Read path')
args = parser.parse_args()

exec_list = ['/home/zgh23/code/SparseWS/bench-hash','/home/zgh23/code/SparseWS/bench-hash-flex','/home/zgh23/code/SparseWS/bench-hash-flex-mt','/home/zgh23/code/SparseWS/bench-hash-mt','/home/zgh23/code/SparseWS/bench-coord-c','/home/zgh23/code/SparseWS/bench-coord-cf','/home/zgh23/code/SparseWS/bench-bucket']
# Input file is a txt file with each line being a matrix name
# Example: regular_1000_20_0.02
# Output file is a txt file with each line being a mtrix name and a command connected by a comma
# Example: regular_1000_20_0.02,/home/zgh23/code/SparseWS/bench-hash-flex-mt
# The commands come from exec_list
with open(args.input_file, "r") as input_file, open(args.output_file, "w") as output_file:
    for line in input_file:
        line = line.strip()
        for exec in exec_list:
            output_file.write(line + "," + exec + "\n")