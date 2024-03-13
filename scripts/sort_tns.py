import argparse
def sort_all_lines(input_file, output_file):
    tuple_list = []
    with open(input_file, 'r') as f:
        lines = f.readlines()
    
        for line in lines:
            line = line.split()
            cur = (line[0], line[1], line[2], line[3])
            tuple_list.append(cur)
    sorted_list = sorted(tuple_list, key=lambda x: (x[0], x[1], x[2]))
    with open(output_file, 'w') as f:
        for line in sorted_list:
            line = ' '.join(line)
            f.write(line + '\n')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=str, default='/scratch/zgh23/sparse_ten/nell-2.tns')
    parser.add_argument('--output', type=str)
    args = parser.parse_args()  
    sort_all_lines(args.input, args.output)