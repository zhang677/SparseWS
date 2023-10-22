import re
import argparse

def extract_memory_data(input_string, functions=None):
    number_pattern = r'(\d{1,3}(?:,\d{3})*)'
    string_pattern = r'(\/[\w\/.-]+:\w+)'
    pattern = number_pattern+r'\s+'+number_pattern+r'\s+'+number_pattern+r'\s+'+number_pattern+r'\s+'+number_pattern+r'\s+'+number_pattern+r'\s+'+string_pattern
    matches = re.findall(pattern, input_string, re.MULTILINE)
    
    extracted_data = []
    for match in matches:
        # print(match)
        # # get the third number_pattern from match

        totB = match[2].replace(',', '')
        totFdB = match[4].replace(',', '')
        full_function_name = match[6]
        function_name = full_function_name.split(':')[-1]
        # print((totB, totFdB, function_name))
        if function_name in functions:
            extracted_data.append((function_name, totB, totFdB))
            functions.remove(function_name)
        if functions == []:
            break
    
    return extracted_data

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Get memory consumption from xtree memory report')
    parser.add_argument('--input', default="/home/zgh23/code/SparseWS/data/profile/memory/xtmem_hash_nofuse_gre_512.txt", type=str, help='Store path')
    args = parser.parse_args()
    with open(args.intput, 'r') as f:
        input_string = f.read()
    functions = {'compute_nofuse','transpose','CSC_CSR_T_hash_nonfuse'}
    extracted_data = extract_memory_data(input_string, functions)
    
    for data in extracted_data:
        print(data)