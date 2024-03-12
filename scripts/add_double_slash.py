def add_double_slash():
    file = "/home/zgh23/code/SparseWS/info/test_matrix_info_merged_latex.txt"
    outfile = "/home/zgh23/code/SparseWS/info/test_matrix_info_merged_latex_double_slash.txt"
    with open(file, "r") as f:
        with open(outfile, "w") as out:
            lines = f.readlines()
            for line in lines:
                # Add \\ to the end of each line
                line = line.strip() + " \\\\" + "\n"
                out.write(line)
    
def add_kind_name():
    file = "/home/zgh23/code/SparseWS/info/test_matrix_info_merged.txt"
    outfile = "/home/zgh23/code/SparseWS/info/test_matrix_info_merged_latex_double_slash_kind_name.txt"
    with open(file, "r") as f:
        with open(outfile, "w") as out:
            lines = f.readlines()
            for line in lines:
                # Add \\ to the end of each line
                name = line.split(",")[0]
                print(name)
                with open("/scratch/zgh23/sparse_mat/" + name + ".mtx", "r") as mtx:
                    for in_line in mtx:
                        if in_line.startswith("% kind: "):
                            kind_name = in_line.split(": ")[1].strip()
                            print(kind_name)
                            break
                line = line.strip() + " & " + kind_name + " \\\\" + "\n"
                out.write(line)

if __name__ == "__main__":
    # add_double_slash()
    add_kind_name()