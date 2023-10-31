import argparse

def count_each_mode(input_file):
    order = 3
    dim = [0,0,0]
    with open(input_file, 'r') as f:
        for line in f:
            line = line.split()
            val = line[-1]
            line = [int(x) for x in line[:-1]]
            for i in range(order):
                dim[i] = max(dim[i], line[i])
    print(input_file)
    print("Max elements: ")
    print(dim)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=str, default='/scratch/zgh23/sparse_ten/nell-2.tns')
    args = parser.parse_args()

    count_each_mode(args.input)

