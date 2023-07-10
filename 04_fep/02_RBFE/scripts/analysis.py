# analysis script

import warnings
import BioSimSpace as BSS
import sys
import csv
import os as _os

# someway that this is for the transformations listed in the csv file made at the start for the setup
# or alternatively do this in bash
trans = sys.argv[1].rstrip()
lig_1 = trans.split("~")[0]
lig_2 = trans.split("~")[1]
engine = sys.argv[2].rstrip()

main_dir = _os.environ["MAINDIRECTORY"]
final_results_file_path = f"{main_dir}/outputs/final_summary_{engine}.csv"

path_to_dir = f"{main_dir}/outputs/{engine}/{trans}"

if not _os.path.exists(path_to_dir):
    raise OSError(f"{path_to_dir} does not exist.")

print(f"analysing results for {path_to_dir}")
print(f"using {chosen_method} and {chosen_estimator} for analysis")

try:
    pmf_bound, overlap_matrix_bound = BSS.FreeEnergy.Relative.analyse(
        f"{path_to_dir}/bound"
    )
except Exception as e:
    print(e)
    print(f"Unable to analyse values for bound in {path_to_dir}.")

try:
    pmf_free, overlap_matrix_free = BSS.FreeEnergy.Relative.analyse(
        f"{path_to_dir}/free"
    )
except Exception as e:
    print(e)
    print(f"Unable to analyse values for free in {path_to_dir}.")

freenrg_rel = BSS.FreeEnergy.Relative.difference(pmf_bound, pmf_free)
freenrg_val = freenrg_rel[0].value()
freenrg_err = freenrg_rel[1].value()

# ####### WRITING DATA for the final result
# data point for average
data_point = [lig_1, lig_2, str(freenrg_val), str(freenrg_err), engine]
print(data_point)

# use csv to open the results file.
with open(final_results_file_path, "a") as freenrg_writefile:
    writer = csv.writer(freenrg_writefile)

    # first, write a header if the file is created for the first time.
    if _os.path.getsize(final_results_file_path) == 0:
        print(f"Starting {final_results_file_path} file.")
        writer.writerow(["lig_1", "lig_2", "freenrg", "error", "engine"])


with open(final_results_file_path, "r") as freenrg_readfile:
    # then, grab all of the data that is already in the file.
    reader = csv.reader(freenrg_readfile)
    data_entries = [row for row in reader]

# check if our data entry is not already in the results file. Raise an error if is.
if data_point in data_entries:
    warnings.warn(
        f"Results for in {trans}, {engine} are already in {final_results_file_path}."
    )

else:
    # at this point we know that we are writing a new entry in the results file. Append the line to the file.
    # use csv to open the results file.
    with open(final_results_file_path, "a") as freenrg_writefile:
        writer = csv.writer(freenrg_writefile)
        print(
            f"Writing results. Free energy of binding is {freenrg_rel[0]} and the error is {freenrg_rel[1]} for {trans}, {engine}."
        )
        writer.writerow(data_point)
