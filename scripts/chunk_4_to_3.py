import argparse
import random

def modify_each_line(input_file, output_file):
    with open(input_file, 'r') as in_file, open(output_file, 'w') as out_file:
        for line in in_file:
            line = line.split()
            line = line[:-2]
            # val is a random float between -0.5 and 0.5
            val = random.uniform(-0.5, 0.5)
            # convert val to a string, keep 4 digits and append it to the line
            line.append(str(round(val, 4)))
            line = ' '.join(line)
            out_file.write(line + '\n')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=str, default='/scratch/zgh23/sparse_ten/nips.tns')
    parser.add_argument('--output', type=str, default='/scratch/zgh23/sparse_ten/nips3.tns')
    args = parser.parse_args()

    modify_each_line(args.input, args.output)