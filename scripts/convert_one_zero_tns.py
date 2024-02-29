import argparse

# Input is a file, each line is composed of several ints and one float
# E.g. 
# 1 2 1 0.5
# 2 3 1 0.5
# Output is a file, where each line is composed of ints with the value minus 1 and the same float
# E.g.
# 0 1 0 0.5
# 1 2 0 0.5

def read_all_lines(input_file, output_file):
    with open(input_file, 'r') as f:
        lines = f.readlines()
    with open(output_file, 'w') as f:
        for line in lines:
            line = line.split()
            val = line[-1]
            line = [str(int(x) - 1) for x in line[:-1]]
            line.append(val)
            line = ' '.join(line)
            f.write(line + '\n')

def read_each_line(input_file, output_file):
    with open(input_file, 'r') as f:
        with open(output_file, 'w') as g:
            for line in f:
                line = line.split()
                val = line[-1]
                line = [str(int(x) + 1) for x in line[:-1]]
                line.append(val)
                line = ' '.join(line)
                g.write(line + '\n')

def count_each_mode(input_file):
    order = 3
    dim = [0,0,0]
    nnz = 0
    with open(input_file, 'r') as f:
        for line in f:
            line = line.split()
            val = line[-1]
            line = [int(x) for x in line[:-1]]
            for i in range(order):
                dim[i] = max(dim[i], line[i])
            nnz += 1
    print(input_file)
    print("Max elements: ")
    print(dim)
    print("Number of non-zero elements: ")
    print(nnz)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=str, default='/scratch/zgh23/sparse_ten/nell-2.tns')
    parser.add_argument('--output', type=str)
    args = parser.parse_args()

    # read_all_lines(args.input, args.output)
    read_each_line(args.input, args.output)
    # count_each_mode(args.input)

