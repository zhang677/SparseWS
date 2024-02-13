import argparse

def clean_zeros_ctf(input_file, output_file):
    with open(input_file, 'r') as f:
        with open(output_file, 'w') as g:
            for line in f:
                line_sp = line.split()
                if line_sp[2] == '0.000000':
                    continue
                g.write(line)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=str, default='/home/zgh23/code/ctf/A-facebook-4.txt')
    parser.add_argument('--output', type=str, default='/home/zgh23/code/ctf/A-facebook-4-clean.txt')

    args = parser.parse_args()
    clean_zeros_ctf(args.input, args.output)