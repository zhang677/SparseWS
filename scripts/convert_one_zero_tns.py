import argparse

# Input is a file, each line is composed of several ints and one float
# E.g. 
# 1 2 1 0.5
# 2 3 1 0.5
# Output is a file, where each line is composed of ints with the value minus 1 and the same float
# E.g.
# 0 1 0 0.5
# 1 2 0 0.5

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=str, default='/scratch/zgh23/sparse_ten/nell-2.tns')
    parser.add_argument('--output', type=str, default='/scratch/zgh23/sparse_ten/nell-2-zero.tns')
    args = parser.parse_args()
    with open(args.input, 'r') as f:
        lines = f.readlines()
    with open(args.output, 'w') as f:
        for line in lines:
            line = line.split()
            val = line[-1]
            line = [str(int(x) - 1) for x in line[:-1]]
            line.append(val)
            line = ' '.join(line)
            f.write(line + '\n')
