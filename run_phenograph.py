import argparse
import pandas as pd
import phenograph

# This script runs Phenograph 

# Get input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--data', help='Infile (rows = cells, columns = features)')
parser.add_argument('--out', help='Outfile (clusters)')
parser.add_argument('-k', help='Number of neighbors for kNN graph', type=int, default=50)
parser.add_argument('--metric', help='Distance metric to use', choices=['manhattan', 'euclidean', 'cosine', 'correlation'], default='cosine')
parser.add_argument('--ncores', help='Number of cores to use', type=int, default=-1)
args = parser.parse_args()

# Run phenograph
data = pd.read_table(args.data)
communities, graph, Q = phenograph.cluster(data, k=args.k, primary_metric=args.metric, n_jobs=args.ncores)

# Write output
out = open(args.out, 'w')
for xi in communities:
    out.write('%s\n' %(xi))
out.close()
