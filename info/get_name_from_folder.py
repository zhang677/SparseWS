import glob
import argparse

parser = argparse.ArgumentParser(description='Output matrix names from csv to a txt file.')
parser.add_argument('--folder_path', default="/home/zgh23/code/SparseWS/data/origin/band", type=str, help='Read path')
parser.add_argument('--output_file', default="/home/zgh23/code/SparseWS/info/band_names.txt", type=str, help='Read path')
args = parser.parse_args()
folder_path = args.folder_path  # Replace with the path to your folder

# Find all files with .mtx extension in the folder
file_names = glob.glob(folder_path + "/*.mtx")

# Write the file names to a.txt
with open(args.output_file, "w") as file:
    for file_name in file_names:
        # file_name is like "/home/zgh23/code/SparseWS/data/origin/band/band_125000_1250_0.019892.mtx"
        # I want to extract "band_125000_1250_0.019892"
        file_name = file_name.split("/")[-1]
        # Delete .mtx
        file_name = file_name[:-4]
        file.write(file_name + "\n")

