import BioSimSpace as BSS
from BioSimSpace import _Exceptions
import sys
import csv
import os
import numpy as np

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import colors
import matplotlib.pyplot as plt
import seaborn as sns

print ("%s %s %s %s" % (sys.argv[0], sys.argv[1], sys.argv[2], sys.argv[3]))
results_file_path = "./outputs/summary.csv"

def plotOverlapMatrix(overlap_matrix, savepath):
	"""
	Given a SOMD-style overlap numpy matrix, plot heatmap in best practices style.

	--args
	overlap matrix (array): numpy 2D array of (n,n)

	--returns
        None; saves plot instead.
	"""


	ovlp_mtx = overlap_matrix

	cmap = colors.ListedColormap(['#FBE8EB','#88CCEE','#78C592', '#117733'])
	bounds=[0.0, 0.025, 0.1, 0.3,0.8]
	norm = colors.BoundaryNorm(bounds, cmap.N, clip=False)
	cbar_kws=dict(ticks=[.025, .1, .3,0.8], label='Phase space overlap')
	ax = sns.heatmap(ovlp_mtx,annot=False, fmt='.2f', linewidths=.3,
	                 annot_kws={"size": 14},square=True,robust=True,cmap=cmap,
	                 norm=norm,cbar_kws=cbar_kws)
	ax.xaxis.tick_top()

	plt.savefig(savepath, dpi=200)

# simply load the FEP directory of the corresponding ligand using BSS.
# this function computes the binding free energy as well.

#print (sys.argv)
engine = sys.argv[3].rstrip()
#print ("@%s@" % sys.argv[1])
#print ("@%s@" % sys.argv[2])
#print ("@%s@" % engine)

path_to_dir = f"./outputs/{engine}/{sys.argv[1]}~{sys.argv[2]}"
print (path_to_dir)
#path_to_dir = "./outputs/%s/%s~%s" % (engine, sys.argv[1], sys.argv[2])
#print (path_to_dir)
freenrg_val = "NaN"
freenrg_err = "NaN"
try:
    pmf_bound, overlap_matrix_bound = BSS.FreeEnergy.Relative.analyse(path_to_dir+"/bound")
    pmf_free, overlap_matrix_free = BSS.FreeEnergy.Relative.analyse(path_to_dir+"/free")

    freenrg = BSS.FreeEnergy.Relative.difference(pmf_bound, pmf_free)
    freenrg_val = round(freenrg[0].value(), 4)
    freenrg_err = round(freenrg[1].value(), 4)

except _Exceptions.AnalysisError:
    freenrg_val = freenrg_err = "NaN"
    overlap_matrix_bound = overlap_matrix_free = None


data_point = [sys.argv[1], sys.argv[2], str(freenrg_val), str(freenrg_err), engine]

####### WRITING DATA

# use csv to open the results file.
with open(results_file_path, "a") as freenrg_writefile:
    writer = csv.writer(freenrg_writefile)
    
    # first, write a header if the file is created for the first time.
    if os.path.getsize(results_file_path) == 0:
        print(f"Starting {results_file_path} file.")
        writer.writerow(["lig_1", "lig_2", "freenrg", "error", "engine"])

    
with open(results_file_path, "r") as freenrg_readfile:
    # then, grab all of the data that is already in the file.
    reader = csv.reader(freenrg_readfile)
    data_entries = [ row for row in reader ]

# check if our data entry is not already in the results file. Raise an error if is.
if data_point in data_entries:
    raise Exception(f"Results for this run are already in {results_file_path}. Exiting.")

# at this point we know that we are writing a new entry in the results file. Append the line to the file.
# use csv to open the results file.
with open(results_file_path, "a") as freenrg_writefile:
    writer = csv.writer(freenrg_writefile)

    print("Writing MBAR results. Free energy of binding and error are (rsp.):")
    print(freenrg)
    writer.writerow(data_point)


# in case of SOMD, we will also have overlap matrices for both legs. These are helpful for troubleshooting, so store 
# them in ./logs/
if overlap_matrix_bound:
    np.save(f"logs/overlap_bound_{sys.argv[1]}~{sys.argv[2]}", np.matrix(overlap_matrix_bound))
    plotOverlapMatrix(np.matrix(overlap_matrix_bound), f"logs/overlap_bound_{sys.argv[1]}~{sys.argv[2]}.png")
else:
    print("Failed to write overlap matrix for bound leg.")
if overlap_matrix_free:
    np.save(f"logs/overlap_free_{sys.argv[1]}~{sys.argv[2]}", np.matrix(overlap_matrix_free))
    plotOverlapMatrix(np.matrix(overlap_matrix_free), f"logs/overlap_free_{sys.argv[1]}~{sys.argv[2]}.png")
else:
    print("Failed to write overlap matrix for free leg.")
