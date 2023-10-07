import networkx as nx
import scipy as sp
import scipy.io
import io
import os
import argparse

# Generate random graphs
def gen_random_regular_graph(args):
    fh = io.BytesIO()
    G = nx.random_regular_graph(args.d, args.n, args.seed)
    m = nx.to_scipy_sparse_array(G)
    sp.io.mmwrite(fh, m, symmetry='general')
    density = float("{:.6f}".format(args.d / args.n)) # not the same as nx.density(G)
    filename = os.path.join(args.output, "regular", f"regular_{args.n}_{args.d}_{density}.mtx")
    file = open(filename, 'w')
    print(fh.getvalue().decode('utf-8'), file = file)
    file.close()

def gen_random_powerlaw_graph(args):
    fh = io.BytesIO()
    G = nx.powerlaw_cluster_graph(args.n, args.d, args.p, args.seed)
    m = nx.to_scipy_sparse_array(G)
    sp.io.mmwrite(fh, m, symmetry='general')
    density = float("{:.6f}".format(m.size / G.number_of_nodes() / G.number_of_nodes())) # not the same as nx.density(G)
    filename = os.path.join(args.output, "powerlaw", f"powerlaw_{args.n}_{args.d}_{density}.mtx")
    file = open(filename, 'w')
    print(fh.getvalue().decode('utf-8'), file = file)
    file.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate random regular graph')
    parser.add_argument('--m', default="regular", type=str, help='method')
    parser.add_argument('--n', default=6, type=int, help='# of nodes')
    parser.add_argument('--d', default=2, type=int, help='degree')
    parser.add_argument('--p', default=0.1, type=float, help='probability')
    parser.add_argument('--seed', default=42, type=int, help='Random seed')
    parser.add_argument('--output', default="/home/zgh23/code/SparseWS/data/origin", type=str, help='Store path')
    args = parser.parse_args()
    if args.m == "regular":
        gen_random_regular_graph(args)
    elif args.m == "powerlaw":
        gen_random_powerlaw_graph(args)
    else:
        print("Please enter a valid method")
        exit()