# Relative Binding Free Energy (RBFE) Network Tutorial - Analysis

This Jupyter notebook is a tutorial on an analysis workflow for RBFE calculations using a Free Energy Perturbation (FEP) network with BioSimSpace, as well as the Freenrgworkflows and Cinnabar packages.

This notebook includes core as well as <span style="color:teal">extra</span> options. These extra sections are there to include some more functionality for when the concepts of these tutorials are applied to your own work.      

**<span style="color:teal">Reading Time:</span>**
~ 30 mins

### Maintainers
- [Anna Herz -- @annamherz](https://github.com/annamherz)

See [README.md](https://github.com/michellab/BioSimSpaceTutorials/blob/main/04_fep/README.md) for complete list of authors.

### Prerequisites
 - Basic Python
 - Part 1 of this workshop (An Introduction to setting up alchemical free energy calculations)
    - this should include basic knowledge of the principles behind RBFE
 - The 01_setup_rbfe notebook

### Learning Objectives
 - Analyse and plot the results of an FEP pipeline setup using BSS.

### Table of Contents
1. [Analysing the edges of the Network](#intro)    
    1.1 [Analysing repeats](#reps)     
    1.2 [Experimental binding affinities](#exp)      
    1.3 [Plotting](#plot)      
2. [Analysing per Ligand](#lig)   
    2.1 [Freenrgworkflows](#fwf)     
    2.2 [Plotting](#plot2)    
    2.3 [Statistical analysis](#stats)     
    2.4 [Outliers](#outliers)  
3. [Cinnabar](#cinnabar)   

### Further reading for this topic
- [LiveComs Best Practices for Alchemical Free Energy Calculations](https://livecomsjournal.org/index.php/livecoms/article/view/v2i1e18378).

**<span style="color:black">Jupyter Cheat Sheet</span>**
- To run the currently highlighted cell and move focus to the next cell, hold <kbd>&#x21E7; Shift</kbd> and press <kbd>&#x23ce; Enter</kbd>;
- To run the currently highlighted cell and keep focus in the same cell, hold <kbd>&#x21E7; ctrl</kbd> and press <kbd>&#x23ce; Enter</kbd>;
- To get help for a specific function, place the cursor within the function's brackets, hold <kbd>&#x21E7; Shift</kbd>, and press <kbd>&#x21E5; Tab</kbd>;
- You can find the full documentation at [biosimspace.org](https://biosimspace.org).

You can use `!` to access terminal commands: e.g. `! head -n 20 myfilename.dat` will display the first 20 lines of a file. 


### Exercises
Exercises are announced using an alert alert-success box in this way:
<div class="alert alert-success">
<b>Exercise 1.1: Write a function that computes bond lengths:</b>
</div>
and followed by an incomplete cell. All exercises should be numbered. 
Missing parts are indicated by:

```python
#FIXME
```
These are included whilst running through the workshops and also in dedicated sections.   



```python
# import libraries
import BioSimSpace as BSS
import os
import csv
import math
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from scipy.stats import sem as sem

# define all the folder locations
main_folder = os.getcwd()
```


```python
from get_tutorial import download
download("02")
```

## 1. Analysing the edges of the Network
<a id="intro"></a>

Once we have obtained our results, we want to analyse them. The basics of this analysis using `BSS.FreeEnergy.Relative.analyse()` for both MBAR and TI, as well as plotting overlap matrices, have already been covered in the introduction to alchemistry part of this workshop.

In this part, we will look at how to carry out a large scale network analysis. As it would take some time to run the analysis for each perturbation, the [runs from this paper](https://chemrxiv.org/engage/chemrxiv/article-details/62ec4b0eadfd35eddd272954) have already been analysed using BSS to give the MBAR RBFE result and error in a csv file format. It is best practice to run repeats of the simulations, which is why there are multiple results files, one for each repeat. The files 'analysis/repeat_{r}_tyk2.csv' contain the results for this 'analysis/network_full.dat'. These runs are from a very large network, with some different ligands than we used for the first part of this RBFE tutorial. This is included incase you want to create a different network and analyse how the output changes based on the edges chosen. However, in this tutorial we will not be analysing this whole large network. We have a different network file to describe the inital network we want to consider, generated using LOMAP, in 'analysis/network.dat'. This uses the ligands from 'inputs/ligands/analysis_ligands'.

There is also a csv file with the experimental results for these ligands. In cases where experimental data is available, for instance when benchmarking a new protein-ligand set, we would like to compare how well FEP is predicting with respect to these data.    

First, we want to set the file paths to our results and experimental data. We also want to create a list of the perturbations and of ligands so we can use/edit these throughout the analysis.



```python
# selected network
# we will define the engine here for ease
engine = "SOMD"

# experimental values (e.g. ic50/ki) for all ligands in our set.
exp_filepath = "analysis/exp_data_tyk2.dat"

# We also want to create a list of the perturbations in our network.
# create a list of the perturbations
perturbations = []

# create a list of ligands
ligands = []

# use the network file to find the ligands and perturbations
for line in open("analysis/network_lomap.dat", "r"):
    lig_0 = line.split()[0]
    lig_1 = line.split()[1]
    pert = f"{lig_0}~{lig_1}"
    perturbations.append(pert)
    if lig_0 not in ligands:
        ligands.append(lig_0)
    elif lig_1 not in ligands:
        ligands.append(lig_1)
    else:
        pass

# make a directory for all our output files
if not os.path.exists("analysis/outputs"):
    os.mkdir("analysis/outputs")
```

Similarly to how we used NetworkX in the setup to visualise our adjusted network, we can also use it here to view the network we are analysing here.


```python
# Generate the graph.
graph = nx.Graph()

# Loop over the nligands and add as nodes to the graph.
for lig in ligands:
    graph.add_node(lig, label=lig, labelloc="t")

# Loop over the edges in the dictionary and add to the graph.
for edge in perturbations:
    lig_0 = edge.split("~")[0]
    lig_1 = edge.split("~")[1]
    graph.add_edge(lig_0, lig_1)

# Plot the networkX graph.
pos = nx.kamada_kawai_layout(graph)
plt.figure(figsize=(8, 8), dpi=150)
nx.draw(
    graph,
    pos,
    edge_color="black",
    width=1,
    linewidths=1,
    node_size=1800,
    node_color="skyblue",
    font_size=12,
    labels={node: node for node in graph.nodes()},
)

plt.savefig("analysis/outputs/analysis_network.png", dpi=300)
plt.show()
```

For the network that we are considering, we want to get results files with these perturbations from the overall file that contains the large network results. We do not usually have to do this, but in this case it is neccessary to get the files that we need for analysis later.


```python
results_all_files = []
no_repeats = 3
# create a list of the results files
for r in list(range(0, no_repeats)):
    file_name = f"analysis/freenrg_repeat_{r}_tyk2.csv"
    results_all_files.append(file_name)

results_files = []

for file in results_all_files:
    new_file_name = f"analysis/outputs/results_{results_all_files.index(file)}.csv"
    with open(new_file_name, "w") as result_file:
        writer = csv.writer(result_file, delimiter=",")
        writer.writerow(["lig_1", "lig_2", "freenrg(kcal/mol)", "error(kcal/mol)", "engine"])

        for row, index in pd.read_csv(file).iterrows():
            pert = f"{index['lig_1']}~{index['lig_2']}"
            if pert in perturbations:
                writer.writerow(
                    [
                        index["lig_1"],
                        index["lig_2"],
                        index["freenrg(kcal/mol)"],
                        index["error(kcal/mol)"],
                        index["engine"],
                    ]
                )

        results_files.append(new_file_name)
```

### 1.1 Analysing repeats
<a id="reps"></a>  

For the first analysis, we will look at the reproducibility between different repeats. We will calculate the average and SEM for the computed runs, based on their repeats. Here, we have six repeats for each. The SEM can be useful to see how reproducible the different runs are, and a larger SEM would indicate a perturbation that has poorer convergence between repeats.   

In general, for our analyses, we will create a dictionary of the values. These can then be easily converted into a pandas dataframe for plotting. It is also a good idea to save output csv files, so that these are easily accessible in the future and can be loaded into python/excel again for any other plotting.



```python
# make a dictionary with the results of the files
comp_dict_list = {}

# append for results file
for res_file in results_files:
    res_df = pd.read_csv(res_file)
    for index, row in res_df.iterrows():
        lig_0 = row[0]
        lig_1 = row[1]
        pert = f"{lig_0}~{lig_1}"
        # if not comp_dict_list[pert]:
        #     print("oop")
        ddG = row[2]

        if pert in comp_dict_list:
            # Key exist in dict, check if is a list
            if not isinstance(comp_dict_list[pert], list):
                # If type is not list then make it list
                comp_dict_list[pert] = [comp_dict_list[pert]]
            # Append the value in list
            comp_dict_list[pert].append(ddG)
        else:
            # As key is not in dict,
            # so, add key-value pair
            comp_dict_list[pert] = ddG

# now calculate all the avg and SEM for the network perturbations
# put these into a dictionary
comp_diff_dict = {}

# write these to a csv file
with open("analysis/outputs/computed_perturbations_average.csv", "w") as comp_pert_file:
    writer = csv.writer(comp_pert_file, delimiter=",")
    writer.writerow(["lig_1", "lig_2", "freenrg(kcal/mol)", "error(kcal/mol)", "engine"])
    for pert in perturbations:
        ddGs = comp_dict_list[pert]
        lig_0 = pert.split("~")[0]
        lig_1 = pert.split("~")[1]
        comp_ddG = np.average(ddGs)
        comp_err = sem(ddGs)

        # update the dictionary for plotting later
        comp_diff_dict.update({pert: (comp_ddG, comp_err)})

        writer.writerow([lig_0, lig_1, comp_ddG, comp_err, engine])
```

### 1.2 Experimental binding affinities

Next, we want to visualise our results whilst comparing them to experimental values. In this example here, TYK2 has binding affinities in Ki, and can be converted using ΔG = RTlnK . It is important at this stage to make sure that the units match, so they are consequently converted into kcal/mol. We can carry this conversion out like below:


```python
# create a dictionary for the experimental values
exper_val_dict = {}

# set the temperature in Kelvin
temp = 300

with open(exp_filepath, "r") as exp_file:
    for line in exp_file:
        lig = line.split(",")[0]
        if lig != "lig":
            exp_val = line.split(",")[1].strip()
        if exp_val != "value":
            exp_val = float(exp_val)
            # gas constant in kcal per Kelvin per mol, exp val converted into M
            exp_kcal = 0.0019872041 * temp * np.log(exp_val / (10 ^ 9))
        err = line.split(",")[2].strip()
        if err != "error":
            err = float(err)
            # get the upper and lower values for the error
            exp_upper = exp_val + err
            exp_lower = exp_val - err
            exp_upper_kcal = 0.0019872041 * temp * np.log(exp_upper / (10 ^ 9))
            exp_lower_kcal = 0.0019872041 * temp * np.log(exp_lower / (10 ^ 9))
            err_kcal = abs(exp_upper_kcal - exp_lower_kcal) / 2
            exper_val_dict.update({lig: (exp_kcal, err_kcal)})
```


```python
# now that we have our dictionary,
# we can also create a dictionary with all the experimental values for the perturbations
exper_diff_dict = {}

# calculate the experimental RBFEs
# write these to a csv file
with open("analysis/outputs/experimental_perturbations.csv", "w") as exp_pert_file:
    writer = csv.writer(exp_pert_file, delimiter=",")
    writer.writerow(["lig_1", "lig_2", "freenrg(kcal/mol)", "error(kcal/mol)", "engine"])

    for pert in perturbations:
        lig_0 = pert.split("~")[0]
        lig_1 = pert.split("~")[1]
        exper_ddG = exper_val_dict[lig_1][0] - exper_val_dict[lig_0][0]
        exper_err = math.sqrt(
            math.pow(exper_val_dict[lig_0][1], 2)
            + math.pow(exper_val_dict[lig_1][1], 2)
        )
        exper_diff_dict.update({pert: (exper_ddG, exper_err)})

        writer.writerow([lig_0, lig_1, exper_ddG, exper_err, "experimental"])
```

### 1.3 Plotting
<a id="plot"></a>

Now we have our computational and experimental in dicitonary format, we can turn this into a pandas dataframe. For plotting it is typically easier to work with the pandas library, which is why this next piece of code will reshape it to this.   

Note that if pandas returns value errors at this step, it is likely there are ligands missing from either your FEP outputs or your experimental input. If a single repeat of a perturbation is missing, the code should be able to deal with this, but if a whole perturbation from the network file is missing there should be an error.



```python
freenrg_pert_dict = {}

# construct dict with experimental freenrg and error and computed
for pert in perturbations:
    exp_ddG = exper_diff_dict[pert][0]
    exp_err = exper_diff_dict[pert][1]
    comp_ddG = comp_diff_dict[pert][0]
    comp_err = comp_diff_dict[pert][1]
    freenrg_pert_dict[pert] = [exp_ddG, exp_err, comp_ddG, comp_err]

# want to put these in a joint dictionary so can convert into pandas df
freenrg_df_pert = pd.DataFrame(
    freenrg_pert_dict, index=["freenrg_exp", "err_exp", "freenrg_fep", "err_fep"]
).transpose()

# save our results to a file that can be opened in e.g. Excel.
freenrg_df_pert.to_csv("analysis/outputs/fep_diff_results_table.csv")
```

Now, we can plot our results against the experimental data. This is best done using a scatter plot:


```python
# plot a scatter plot
plt.rc("font", size=12)
fig, ax = plt.subplots(figsize=(10, 10))

# get these based on which column the data is in.
x = freenrg_df_pert["freenrg_exp"]
y = freenrg_df_pert["freenrg_fep"]
x_er = freenrg_df_pert["err_exp"]
y_er = freenrg_df_pert["err_fep"]

# plotting the scatterplot
scatterplot = [plt.scatter(x, y, zorder=10)]

# plotting error bars
plt.errorbar(
    x,
    y,
    yerr=y_er,
    # xerr=x_er,   # comment this line to hide experimental error bars \
    # as this can sometimes overcrowd the plot.
    ls="none",
    lw=0.5,
    capsize=2,
    color="black",
    zorder=5,
)

# plot 1/2 kcal bounds:
plt.fill_between(
    x=[-100, 100],
    y2=[-100.25, 99.75],
    y1=[-99.75, 100.25],
    lw=0,
    zorder=-10,
    alpha=0.3,
    color="grey",
)
# upper bound:
plt.fill_between(
    x=[-100, 100],
    y2=[-99.5, 100.5],
    y1=[-99.75, 100.25],
    lw=0,
    zorder=-10,
    color="grey",
    alpha=0.2,
)
# lower bound:
plt.fill_between(
    x=[-100, 100],
    y2=[-100.25, 99.75],
    y1=[-100.5, 99.5],
    lw=0,
    zorder=-10,
    color="grey",
    alpha=0.2,
)

# get the bounds. This can be done with min/max or simply by hand.
all_freenrg_values = np.concatenate(
    [freenrg_df_pert["freenrg_exp"].values, freenrg_df_pert["freenrg_fep"].values]
)
min_lim = min(all_freenrg_values)
max_lim = max(all_freenrg_values)

# for a scatterplot we want the axis ranges to be the same.
plt.xlim(min_lim * 1.3, max_lim * 1.3)
plt.ylim(min_lim * 1.3, max_lim * 1.3)

plt.axhline(color="black", zorder=1)
plt.axvline(color="black", zorder=1)

plt.ylabel("Computed $\Delta\Delta$G$_{bind}$ / kcal$\cdot$mol$^{-1}$")
plt.xlabel("Experimental $\Delta\Delta$G$_{bind}$ / kcal$\cdot$mol$^{-1}$")

plt.savefig("analysis/outputs/fep_vs_exp_scatterplot.png", dpi=300)
plt.show()
```

This is a very simple way to initially visualise our results. Below, we will consider more rigorous methods for analysis.

## 2. Network Analysis using Cinnabar
<a id="cinnabar"></a>

Instead of just computing ΔΔG values for each transformation, we would like to estimate the ΔΔG value for each individual ligand. There are some involved algorithms needed for these steps, and in this section of the tutorial, we will be using [cinnabar](https://github.com/OpenFreeEnergy/cinnabar) to care of that for us. This already has all the analysis and plotting in functions.


```python
from cinnabar import wrangle, plotting
```


```python
# write to a csv file
with open(f"analysis/outputs/cinnabar_data.csv", "w") as cinnabar_data_file:
    writer = csv.writer(cinnabar_data_file, delimiter=",")

    # first, write the experimental data
    writer.writerow(["# Experimental block"])
    writer.writerow(["# Ligand", "expt_DDG", "expt_dDDG"])

    for lig in exper_val_dict.keys():
        writer.writerow([lig, f"{exper_val_dict[lig][0]}", f"{exper_val_dict[lig][1]}"])

    # second write the perturbation data
    writer.writerow([" "])
    writer.writerow(["# Calculated block"])
    writer.writerow(
        ["# Ligand1", "Ligand2", "calc_DDG", "calc_dDDG(MBAR)", "calc_dDDG(additional)"]
    )

    # need to write the average of the data, otherwise cinnabar just uses the last entry
    comp_diff_dict = {}

    # make a dictionary with the results of the files
    comp_dict_list = {}
    comp_err_dict_list = {}

    # append for results file
    for res_file in results_files:
        res_df = pd.read_csv(res_file)
        for index, row in res_df.iterrows():
            lig_0 = row[0]
            lig_1 = row[1]
            pert = f"{lig_0}~{lig_1}"
            ddG = float(row[2])
            ddG_err = float(row[3])

            if pert in comp_dict_list:
                # Key exist in dict, check if is a list
                if not isinstance(comp_dict_list[pert], list):
                    # If type is not list then make it list
                    comp_dict_list[pert] = [comp_dict_list[pert]]
                if not isinstance(comp_err_dict_list[pert], list):
                    # If type is not list then make it list
                    comp_err_dict_list[pert] = [comp_err_dict_list[pert]]
                # Append the value in list
                comp_dict_list[pert].append(ddG)
                comp_err_dict_list[pert].append(ddG_err)
            else:
                # As key is not in dict,
                # so, add key-value pair
                comp_dict_list[pert] = [ddG]
                comp_err_dict_list[pert] = [ddG_err]

    # now calculate all the avg and SEM for the network perturbations
    for pert in perturbations:
        lig_0 = pert.split("~")[0]
        lig_1 = pert.split("~")[1]

        # check if the perturbations calculated are also those in the network file and if any are missing
        try:
            # find the values in the dictionary
            ddGs = comp_dict_list[pert]
            ddGs_error = comp_err_dict_list[pert]
            # calculate the average and the error
            comp_ddG = np.average(ddGs)
            # comp_ddG = np.average([ddG.value() for ddG in ddGs])
            if len(ddGs) == 1:
                comp_err = ddGs_error[0]
                # comp_err = ddGs_error.value()
            else:
                comp_err = sem(ddGs)
                # comp_err = sem([ddG.value() for ddG in ddGs])

        # if unable to calculate one of the perturbations, this is a None value.
        except:
            comp_ddG = None
            comp_err = None

        # update the dictionary
        comp_diff_dict.update({pert: (comp_ddG, comp_err)})

    # write to file
    for key in comp_diff_dict:
        lig_0 = key.split("~")[0]
        lig_1 = key.split("~")[1]
        comp_ddG = comp_diff_dict[key][0]
        comp_err = comp_diff_dict[key][1]

        if not comp_ddG:
            pass
        else:
            pert = f"{lig_0}~{lig_1}"
            if pert in perturbations:
                writer.writerow([lig_0, lig_1, comp_ddG, comp_err, "0.0"])
```

Next, we can use this file to instatiate a FEMap object, that we can then use to plot both the perturbations and also the analysis per ligand. These plots also return the statistical analysis.


```python
network = wrangle.FEMap(f"analysis/outputs/cinnabar_data.csv")
# plot the perturbations
plotting.plot_DDGs(network.graph, title="DDGs", filename=f"analysis/outputs/DDGs.png", figsize=6)
# plot the ligands
plotting.plot_DGs(network.graph, title="DGs", filename=f"analysis/outputs/DGs.png", figsize=6)
```

### 2.1. Outliers
<a id="outliers"></a>
In a network analysis, a badly computed perturbation can impact the overall quality of our network. Below, we want to annotate any outliers for our perturbations (on the scatter plot) so we can exclude them, rerun  our network analysis, and improve this. We can do this in one of two ways - either we can remove the values that disagree more substantially with the experimental, as below, or, if no experimental values are available, we can consider the error of the perturbation.

First, let's look at our scatter plot with the outlier:


```python
number_outliers_to_annotate = 3

# get an array of the MUE values comparing experimental and FEP values. Take the absolute values.
mue_values = abs(freenrg_df_pert["freenrg_exp"] - freenrg_df_pert["freenrg_fep"])

# find the n ligand names that are outliers.
outlier_names = mue_values.nlargest(number_outliers_to_annotate).index.values.tolist()
print(outlier_names)

# construct a list of labels to annotate the scatterplot with.
annot_labels = []
colours = []
for ligand in freenrg_df_pert.index.values:
    # if the ligand is an outlier, append the name to the annotation labels list.
    if ligand in outlier_names:
        annot_labels.append(ligand)
        colours.append("red")
    else:
        # if the ligand is not an outlier, append an empty string to the annotation labels list.
        annot_labels.append("")
        colours.append("blue")

# Create the same scatterplot as above. Can include some more of the formatting if needed.
plt.figure(figsize=(7, 7))

plt.scatter(
    freenrg_df_pert["freenrg_exp"], freenrg_df_pert["freenrg_fep"], zorder=10, c=colours
)
plt.plot((-10, 10), (-10, 10))
plt.xlim(min_lim * 1.3, max_lim * 1.3)
plt.ylim(min_lim * 1.3, max_lim * 1.3)

# then, after generating the figure, we can annotate:
for i, txt in enumerate(annot_labels):
    plt.annotate(
        txt,
        (
            freenrg_df_pert["freenrg_exp"].values.tolist()[i] + 0.1,  # x coords
            freenrg_df_pert["freenrg_fep"].values.tolist()[i] + 0.1,
        ),  # y coords
        size=20,
        color="crimson",
    )

plt.savefig("analysis/outputs/fep_vs_exp_outlier_plot.png", dpi=300)
plt.show()
```


```python
error_table = pd.DataFrame(freenrg_df_pert["err_fep"]).sort_values(by="err_fep")

# we can write the table to a csv file that can be opened in e.g. Excel.
error_table.to_csv("analysis/outputs/error_table.csv")
error_table
```

<div class="alert alert-success">
<b>Exercise 2.1.1: Excluding intermediates from the Network analysis</b>
</div>

Try exculding this perturbation and rerunning the above Network analysis. First, we need to remove the perturbation. Then, we need to make sure that our new output image is being saved using a different file path. Adjust these in the cells below where the #FIXME is.   



```python
# remove the perturbation
# FIXME
```

<details><summary {style='color:green;font-weight:bold'}> Click here to see solution to Exercise. </summary>

```python
perturbations.remove("lig_ejm42~lig_ejm49")
perturbations.remove("lig_ejm44~lig_ejm45")
perturbations.remove("lig_ejm31~lig_ejm49")
```

</details>


```python
# rerun the above cells for network analysis with a new file path
# FIXME
```

<details><summary {style='color:green;font-weight:bold'}> Click here to see solution to Exercise. </summary>

```python
# write to a csv file
with open(f"analysis/outputs/cinnabar_data_outliers_removed.csv", "w") as cinnabar_data_file:
    writer = csv.writer(cinnabar_data_file, delimiter=",")

    # first, write the experimental data
    writer.writerow(["# Experimental block"])
    writer.writerow(["# Ligand","expt_DDG","expt_dDDG"])

    for lig in exper_val_dict.keys():
        writer.writerow([lig,f"{exper_val_dict[lig][0]}",f"{exper_val_dict[lig][1]}"])


    # second write the perturbation data
    writer.writerow([" "])
    writer.writerow(["# Calculated block"])
    writer.writerow(["# Ligand1","Ligand2","calc_DDG","calc_dDDG(MBAR)", "calc_dDDG(additional)"])

    # need to write the average of the data, otherwise cinnabar just uses the last entry
    comp_diff_dict = {}

    # make a dictionary with the results of the files
    comp_dict_list = {}
    comp_err_dict_list = {}

    # append for results file
    for res_file in results_files:
        res_df = pd.read_csv(res_file)
        for index,row in res_df.iterrows():

            lig_0 = row[0]
            lig_1 = row[1]
            pert = f"{lig_0}~{lig_1}"
            ddG = float(row[2])
            ddG_err = float(row[3])
                
            if pert in comp_dict_list:
                # Key exist in dict, check if is a list
                if not isinstance(comp_dict_list[pert], list):
                    # If type is not list then make it list
                    comp_dict_list[pert] = [comp_dict_list[pert]]
                if not isinstance(comp_err_dict_list[pert], list):
                    # If type is not list then make it list
                    comp_err_dict_list[pert] = [comp_err_dict_list[pert]]
                # Append the value in list
                comp_dict_list[pert].append(ddG)
                comp_err_dict_list[pert].append(ddG_err)
            else:
                # As key is not in dict,
                # so, add key-value pair
                comp_dict_list[pert] = [ddG]
                comp_err_dict_list[pert] = [ddG_err]

    # now calculate all the avg and SEM for the network perturbations  
    for pert in perturbations:
        lig_0 = pert.split("~")[0]
        lig_1 = pert.split("~")[1]
        
        # check if the perturbations calculated are also those in the network file and if any are missing
        try:
            # find the values in the dictionary
            ddGs = comp_dict_list[pert]
            ddGs_error = comp_err_dict_list[pert]
            # calculate the average and the error
            comp_ddG = np.average(ddGs)
            # comp_ddG = np.average([ddG.value() for ddG in ddGs])
            if len(ddGs) == 1:
                comp_err = ddGs_error[0]
                # comp_err = ddGs_error.value()
            else:
                comp_err = sem(ddGs)
                # comp_err = sem([ddG.value() for ddG in ddGs])
    
        # if unable to calculate one of the perturbations, this is a None value.
        except:
            comp_ddG = None
            comp_err = None

        #update the dictionary
        comp_diff_dict.update({pert:(comp_ddG, comp_err)})

    # write to file
    for key in comp_diff_dict:
        lig_0 = key.split("~")[0]
        lig_1 = key.split("~")[1]
        comp_ddG = comp_diff_dict[key][0]
        comp_err = comp_diff_dict[key][1]            

        if not comp_ddG:
            pass
        else:
            pert = f"{lig_0}~{lig_1}"
            if pert in perturbations:
                writer.writerow([lig_0, lig_1, comp_ddG, comp_err, "0.0"])

network = wrangle.FEMap(f"analysis/outputs/cinnabar_data_outliers_removed.csv")
# plot the perturbations
plotting.plot_DDGs(network.graph, title="DDGs", filename=f"analysis/outputs/DDGs_outliers_removed.png", figsize=6)
# plot the ligands
plotting.plot_DGs(network.graph, title="DGs", filename=f"analysis/outputs/DGs_outliers_removed.png", figsize=6)
```

</details>

However, removing perturbations is not always a good idea, especially if we do not have very large outliers as is the case here. As can be seen in the sample output folder, when we remove [three of the perturbations](example_output/example_output_analysis/DGs_three_outliers_removed.png), we loose one of the ligands in our network (N=16). Removing only [one perturbation](example_output/example_output_analysis/DGs_three_outliers_removed.png) does not improve our statistics either. Generally, removing perturbations should only be done after very careful consideration.

### 2.2. Visualising the ligand

<div class="alert alert-success">
<b>Exercise 2.2.1: Visualising the ligand</b>
</div>

Based on our results, we may also want to visualise which of our ligands had the highest binding affinity.

**Hint:** First, obtain the node information from the cinnabar graph. Sort this resulting data frame to find the ligand with the best binding affinity and then use BSS to visualise it like in the introductionary workshop.


```python
# FIXME
```

<details><summary {style='color:green;font-weight:bold'}> Click here to see solution to Exercise. </summary>

```python

lig_dict = {entry[1]['name']:entry[1]['calc_DG'] for entry in network.graph.nodes.data()}
lig_df = pd.DataFrame.from_dict(lig_dict, orient="index", columns=["calc_DG"])
lig_df.to_csv("analysis/outputs/ligand_calc_DG.csv")
lig_df.sort_values("calc_DG")

```

The best binder based on this data is 'lig_jmc27'. We can visualise this using BSS:

```python

BSS.Notebook.View("inputs/ligands/analysis_ligands/jmc_27.sdf").system()

```
</details>

## 3. Congratulations
You have made it to the end! This notebook has covered the basics of some analysis that can be done on our RBFE results. There are many ways to customise this, and it is also a good idea to check out the further reading for a good overview of AFE data reporting best practices: [Section 8.7, LiveComs Best Practices for Alchemical Free Energy Calculations](https://livecomsjournal.org/index.php/livecoms/article/view/v2i1e18378).

As for the setup workshop, all the generated graphs and tables outputs are available in the 'example_output' folder.



```python
# Uncomment the line below to download example outputs
# download("03")
```
