# Absolute Binding Free Energies and Experimental Features in BioSimSpace: Analysis

This is the second of two jupyter notebooks on Absolute Binding Free Energy (ABFE) calculations in BioSimSpace. This notebook describes the analysis of ABFE calculations. The first notebook, "setup_abfe.ipynb" discussed the integration of experimental features into BioSimSpace through the *sandpit*, and the setup of absolute binding free energy calculations using BioSimSpace.

### Maintainers
 - Finlay Clark -- finlay.clark@ed.ac.uk (@fjclark)

### Prerequisites

 - Completion of the FEP `01_intro_to_alchemy` tutorial.
 - Completion of the first notebook in this tutorial, `setup_abfe.ipynb`
    
### Learning Objectives

 - Understand the functionality supplied by BioSimSpace for the analysis of ABFE calculations
 - Understand some general principles applying to the analysis of ABFE calculations
 
### TOC

- [1. Basic Analysis in BioSimSpace](#basic_ana)    
    -  [1.1 Overlap Matrices](#overlap)    
    -  [1.2 Potentials of Mean Force](#pmf)
-  [2. Further Analysis](#further_analysis)   
    - [2.1 Are My Results Reproducible?](#repro)     
    - [2.2 Is My Equilibration Period Sufficient and are my Results Converged?](#conv)    
    - [2.3 Are there Alternative Binding Poses?](#poses)    
    - [2.4 Do I Need Symmetry Corrections?](#sym)    
    - [2.5 Final Comments](#comments)
 

### Further reading for this topic
- [LiveComs Best Practices for Alchemical Free Energy Calculations](https://livecomsjournal.org/index.php/livecoms/article/view/v2i1e18378).

**<span style="color:black">Jupyter Cheat Sheet</span>**
- To run the currently highlighted cell and move focus to the next cell, hold <kbd>&#x21E7; Shift</kbd> and press <kbd>&#x23ce; Enter</kbd>;
- To run the currently highlighted cell and keep focus in the same cell, hold <kbd>&#x21E7; Ctrl</kbd> and press <kbd>&#x23ce; Enter</kbd>;
- To get help for a specific function, place the cursor within the function's brackets, hold <kbd>&#x21E7; Shift</kbd>, and press <kbd>&#x21E5; Tab</kbd>;

### Link to documentation:
You can find the full documentation at [biosimspace.org](https://biosimspace.org).


```python
# Import Exscientia sandpit as before
import BioSimSpace.Sandpit.Exscientia as BSS

# Imports for plotting and analysis
import matplotlib.pyplot as plt
import numpy as np
```

Let's download the tutorial dataset if necessary


```python
from get_tutorial import download
download("03")
```

## 1. Basic Analysis in BioSimSpace
<a id="basic_ana"></a>

Running a full set of ABFE calculations for the MIF / MIF-180 complex takes days, so we've provided some example output for simulations with SOMD in `example_output`. This includes:

- All required input files
- Compressed output energies for all stages of run 1 
- Compressed output energies for the bound discharge stages for runs 1 - 5. 

The Boresch restraints were based on the triazole ring (anchor points shown as spheres):

<img src="images/boresch_restr_2.png" width="400"/>


### 1.1 Overlap Matrices
<a id="overlap"></a>

As discussed in the relative binding free energy workshop, you can analyse your data to obtain the potential of mean force along lambda and overlap matrix through BioSimSpace. Normally, this would be done using the `AlchemicalFreeEnergy` object directly. For example, having run `restrain_fe_calc`:

```Python
restrain_fe_calc.analyse()
```

However, `analyse` can also be used directly through the AlchemicalFreeEnergy class, allowing us to analyse our pre-computed data:


```python
pmf, overlap = BSS.FreeEnergy.AlchemicalFreeEnergy.analyse(
    "example_output/bound/run001/restrain/output", estimator="MBAR"
)
```

<div class="alert alert-success">
<b>Exercise 1: Inspect the Overlap Matrices</b>
</div>

As discussed in the previous notebook, sufficient overlap is required for the MBAR estimator to be reliable. Plot the overlap matrices for run001, for all stages for the bound and free legs, using `BSS.Notebook.plotOverlapMatrix`. Have a sufficient number of $\mathrm{\lambda}$ values been used for all simulations? Can you suggest improvements to the spacing or number of $\mathrm{\lambda}$ windows for any stage?


```python
?BSS.Notebook.plotOverlapMatrix
```


```python
# FIXME: Obtain the overlap data and plot the overlap matrices for all stages using plotOverlapMatrix
```

<details><summary {style='color:green;font-weight:bold'}> Click here to see solution to Exercise. </summary>

```python
# Repeat for all stages
_, overlap = BSS.FreeEnergy.AlchemicalFreeEnergy.analyse("example_output/bound/run001/vanish/output", estimator="MBAR")
BSS.Notebook.plotOverlapMatrix(overlap, continuous_cbar=True)
```
</details>

###  1.2 Potentials of Mean Force
<a id="pmfs"></a>

`BSS.FreeEnergy.AlchemicalFreeEnergy.Analyse()` also returns the potential of mean force (PMF) with respect to $\mathrm{\lambda}$, which shows how the free energy changes as an interaction (restraints, charges, or Lennard-Jones) is scaled with $\mathrm{\lambda}$:


```python
# Plot the PMFs for all stages

# List of tuples of stages of form (leg, stage)
stages = [
    ("free", "discharge"),
    ("free", "vanish"),
    ("bound", "restrain"),
    ("bound", "discharge"),
    ("bound", "vanish"),
]

fig, axs = plt.subplots(1, 5, figsize=(18, 4), dpi=500)
axs = axs.flatten()

for ax, (leg, stage) in zip(axs, stages):
    pmf, _ = BSS.FreeEnergy.AlchemicalFreeEnergy.analyse(
        f"example_output/{leg}/run001/{stage}/output", estimator="MBAR"
    )
    lam_vals = [x[0] for x in pmf]
    delta_g = [x[1].value() for x in pmf]  # Need .value() to remove energy units
    ax.plot(lam_vals, delta_g)
    ax.set_ylabel("$\Delta \it{G}$ / kcal.mol$^-$$^1$")
    ax.set_xlabel("$\mathrm{\lambda}$")
    ax.set_title(f"{leg} {stage}")

fig.tight_layout()
```

<div class="alert alert-success">
<b>Exercise 2: Calculate the Free Energy of Biniding </b>
</div>

Estimate $\Delta\textit{G}_{\textrm{Bind}}^o$ for run 1 by adding up all the terms around the free energy cycle:

<img src="images/abfe_cycle_details.png" width="500"/>

Hints: 
- Remember to change the signs to match the cycle
- Remember to include the correction term for releasing the restraint
- Remember to include an estimate of uncertainty (check `BSS.FreeEnergy.AlchemicalFreeEnergy.analyse?`). The uncertainty returned is generated by pyMBAR and reflects one standard deviation for the normal distribution of the estimate in the asymptotic limit
- When extracting the free energies from the PMF, remember to use .value() to remove the unit

First, obtain the free energies from the simulations:


```python
# Initialise dicts to store free energies and uncertainties
free_energies = {"free":{"discharge":0, "vanish":0}, "bound":{"restrain":0, "discharge":0, "vanish":0}}
conf_ints = {"free":{"discharge":0, "vanish":0}, "bound":{"restrain":0, "discharge":0, "vanish":0}}

# Store the total energy (minus the restraint correction) and 95 % confidence intervals
delta_g_tot = 0
conf_int_tot = 0

for leg in free_energies:
    for stage in free_energies[leg]:
        # FIXME: Using BSS.FreeEnergy.AlchemicalFreeEnergy.analyse, find the free energy changes corresponding to each 
        # stage in the ABFE cycle shown. Remember to include the associated uncertainties.
```

<details><summary {style='color:green;font-weight:bold'}> Click here to see solution to Exercise. </summary>

```python
# Initialise dicts to store free energies and uncertainties
free_energies = {"free":{"discharge":0, "vanish":0}, "bound":{"restrain":0, "discharge":0, "vanish":0}}
conf_ints = {"free":{"discharge":0, "vanish":0}, "bound":{"restrain":0, "discharge":0, "vanish":0}}

# Store the total energy (minus the restraint correction) and 95 % confidence intervals
delta_g_tot = 0
conf_int_tot = 0

for leg in free_energies:
    for stage in free_energies[leg]:
        pmf, _ = BSS.FreeEnergy.AlchemicalFreeEnergy.analyse(f"example_output/{leg}/run001/{stage}/output",
                                                              estimator="MBAR")
        if leg == "bound": # Reverse the sign of the free energy estimate
            delta_g = - pmf[-1][1].value()
        else: # Free leg
            delta_g = pmf[-1][1].value()
        std_dev = pmf[-1][2].value() # Also remove unit
        conf_int = 1.96 * std_dev # Convert to 95 % conf int
        free_energies[leg][stage] = delta_g
        conf_ints[leg][stage] = conf_int
        delta_g_tot += delta_g
        uncert_tot = np.sqrt(conf_int_tot**2 + conf_int**2)

```
Note that BioSimSpace (using `alchemlyb`) automatically attempts to obtain uncorrelated energies by subsampling according to statistical inefficiency. If this fails, the uncertainty will be (further) underestimated.
</details>

Now, obtain the correction term from the restraint:


```python
# Import required units
from BioSimSpace.Sandpit.Exscientia.Units.Length import angstrom
from BioSimSpace.Sandpit.Exscientia.Units.Angle import radian, degree
from BioSimSpace.Sandpit.Exscientia.Units.Energy import kcal_per_mol
from BioSimSpace.Sandpit.Exscientia.Units.Temperature import kelvin

# The restraint object requires a system
system = BSS.IO.readMolecules(
    ["input/complex/mif_mif180.prm7", "input/complex/mif_mif180.rst7"]
)
lig = BSS.Align.decouple(system[0])
system.updateMolecule(0, lig)

# Assign three atoms from the protein
r1 = system.getAtom(935)
r2 = system.getAtom(933)
r3 = system.getAtom(937)

# Assign three atoms from the ligand
l1 = system.getAtom(10)
l2 = system.getAtom(13)
l3 = system.getAtom(20)

restraint_dict = {
    "anchor_points": {"r1": r1, "r2": r2, "r3": r3, "l1": l1, "l2": l2, "l3": l3},
    "equilibrium_values": {
        "r0": 8.14 * angstrom,
        "thetaA0": 2.06 * radian,
        "thetaB0": 1.89 * radian,
        "phiA0": 1.68 * radian,
        "phiB0": 1.52 * radian,
        "phiC0": 0.20 * radian,
    },
    "force_constants": {
        "kr": 10.92 * kcal_per_mol / angstrom**2,
        "kthetaA": 126.83 * kcal_per_mol / (radian * radian),
        "kthetaB": 98.43 * kcal_per_mol / (radian * radian),
        "kphiA": 189.35 * kcal_per_mol / (radian * radian),
        "kphiB": 58.81 * kcal_per_mol / (radian * radian),
        "kphiC": 100.72 * kcal_per_mol / (radian * radian),
    },
}
```


```python
# FIXME: Using the system and restraint dictionary created above,
# create a BSS.FreeEnergy.Restraint object and obtain the correction term
```

<details><summary {style='color:green;font-weight:bold'}> Click here to see solution to Exercise. </summary>

```python
# Get the free energy of turning off the restraint. Note that there is no contribution to the 
# uncertainty, as this depends only on the restraint dictionary and is identical between runs
restraint = BSS.FreeEnergy.Restraint(system, restraint_dict, 298*BSS.Units.Temperature.kelvin)
dg_rest = restraint.getCorrection().value() # remove units of kcal / mol
print(f"Free energy of releasing restraint: {dg_rest:.2f} kcal / mol")
```
</details>


```python
# FIXME: Calculate the free energy of binding for run 1.
```

<details><summary {style='color:green;font-weight:bold'}> Click here to see solution to Exercise. </summary>

```python
# Reverse the sign when combining with the overall free energy to be consistent with the cycle
delta_g_tot -= dg_rest

print(f"Overall free energy of binding from run 1: {delta_g_tot:.2f} +/- {uncert_tot:.2f} kcal / mol (95 % C.I. based on single run)")
```
</details>

The experimental binding free energy of MIF-180 to MIF is [-8.98 $\pm$ 0.28 kcal mol $^{-1}$](https://pubs.acs.org/doi/pdf/10.1021/acs.jpcb.9b07588). Compare your answer. Consider:

- Can we have confidence in our binding free energy estimate? 
    - If yes, on what basis?
    - If no, what additional checks would you suggest? Have any key considerations been omitted?
- Do you think your uncertainty estimate is a realistic representation of the the true uncertainty? 
- In the limit of perfect sampling, would you expect the free energy estimate to be exactly the same as the experimental value?

 ```Python
 #FIXME: Answer the above question
 ```

<details><summary {style='color:green;font-weight:bold'}> Click here to see solution to Exercise. </summary>        

- No. We have not carried out even the most basic tests on convergence and reproducibility. See "Further Analysis".     
- No. Inter-run deviations can be substantially larger than within-run fluctuations due to the tendency of trajectories to become trapped in metastable regions of configuration space. To adequately capture this source of random error in our stated uncertainty, and to obtain a reliable free energy estimate, [we should perform multiple independent replicate runs](https://pubs.acs.org/doi/full/10.1021/acs.jctc.7b01143).       
- No. There will be systematic error as a result of the force field and other approximations, for example LJ interactions may be ignorned beyond some cutoff.       
</details>

## 2. Further Analysis
<a id="further_analysis"></a>

We should regard our free energy estimate above, no matter how similar to experiment, as an extremely rough estimate which serves only as a sanity check on our calculations. To have confidence in our result, we must first perform a number of key checks.

### 2.1 Are My Results Reproducible? </span>
<a id="repro"></a>

A minimal test of reproducibility is performing replicate runs where all atoms are assigned different initial velocities. This was done to produce the data provided for the bound discharge stages for runs 1 - 5 in `example_output`. Even better, we could have used different starting structures extracted from a simulation of the complex.

<div class="alert alert-success">
<b>Exercise 3: Deviation Between Replicate Runs </b>
</div>

Analyse the data for the 5 replicate runs for the discharge stage of the bound leg to provide a revised estimate of the free energy change for this stage. Calculate a new measure of the uncertainty for this stage, and compare it to the uncertainty estimate obtained for a single run. Note that to save space, only the results for the bound discharge stage have been supplied for all 5 runs.


```python
# FIXME: Find the free energy change for "vanishing" the bound ligand based on the
# 5 repeat runs. Provide a revised uncertainty estimate for this stage.
```

<details><summary {style='color:green;font-weight:bold'}> Click here to see solution to Exercise. </summary>

```python
# Import for calculating C.I.
import scipy.stats as st

dgs = []

# Loop over 5 replicates
for i in range(1,6):
    pmf, _ = BSS.FreeEnergy.AlchemicalFreeEnergy.analyse(f"example_output/bound/run00{i}/discharge/output", 
                                                         estimator="MBAR")
    dgs.append(pmf[-1][1].value())

# Convert list to numpy array for easier handling
dgs = np.array(dgs)

# Calculate 95 % C.I. assuming Gaussian distribution of free energy estimates
uncert = dgs.mean() - st.t.interval(0.95, len(dgs)-1, dgs.mean(), scale=st.sem(dgs))[0]

# Now calculate single run uncertainty for first run based on MBAR estimate
pmf, _ = BSS.FreeEnergy.AlchemicalFreeEnergy.analyse(f"example_output/bound/run00{i}/discharge/output", 
                                                     estimator="MBAR")
    
dg_single_run = pmf[-1][1].value()
conf_int_single_run = pmf[-1][2].value() * 1.96

print(f"dg_bound_discharge estimate over 5 replicates: {dgs.mean():.2f} +/- {uncert:.2f} kcal / mol (95 % C.I. based on 5 replicates)")
print(f"dg_bound_discharge estimate from single run: {dg_single_run:.2f} +/- {conf_int_single_run:.2f} kcal / mol (95 % C.I. based on single run)")

```
</details>

Here the inter-run and intra-run uncertainties do not differ dramatically, but much larger deviations can occur when there are sampling issues. Performing several replicate runs is essential for obtaining a meaningful uncertainty estimate, and therefore a meaningful result.

We should note we've assumed that we have enough samples so that the MBAR uncertainty represents one standard deviation for the normal distribution of the estimate, and that our replicate free energies are normally distributed. This is likely to be incorrect, and this should be recognised when dealing with the uncertainties.

### 2.2 Is My Equilibration Period Sufficient and are my Results Converged?
<a id="conv"></a>

If our free energy estimate is correct, then we should have sampled all relevant regions of configurations space at equilibrium. Performing more sampling by increasing the simulation time should not add any new information and our estimate should not change. This can be checked by repeatedly truncating the data at a given simulation time and analysing the free energy change with the remaining data. The results are plotted for the bound stages, where the blue line shows the average over all runs and the beige fill indicates the 95 % t-based confidence intervals. BioSimSpace does not provide the functionality to perform such analysis, but the code used to generate these plots is available [here](https://github.com/fjclark/abfe_analysis/blob/main/abfe_analysis/get_data/convergence_data.py) (provided purely as an example - reuse is not suggested).

<img src="images/convergence_no_equil.png" width="1000"/>

There seems to be some drift with increasing simulation time for the vanish stage especially. In addition, within all runs, all stages evolve in the same direction with increasing simulation time at early values of total simulation time. This is because we have completely failed to include any equilibration time for each $\mathrm{\lambda}$ window! As a result, the systems do not have enough time to relax from starting configurations which would otherwise have a very low probability of being sampled. To avoid this biasing our simulations, we should discard an initial portion of each trajectory as equilibration (or perform short equilibration simulations for each $\mathrm{\lambda}$ window).

We can also check how the PMF changes with increasing total simulation time for each run. For example, for run 1 (note that the cumulative sampling bar should show a maximum of 8 ns for the vanish stage, but is correct for the restrain and discharge stages):

<img src="images/pmf_evolution.png" width="1000"/>

For the bound vanish stage, the rate of divergence of the PMF is by far the greatest around $\lambda = 0.4$, which tells us that $\lambda$ windows in this region take far longer to equilibrate than those before or afterwards. Further analysis reveals this is the region at which the ligand becomes "soft" enough and water begins to enter the binding site. Slow equilibration likely occurs because the partially vanished ligand needs to be "pushed" to the side of the binding site to allow water in, which happens slowly. At higher values of $\lambda$, the ligand's intermolecular interactions are extremely weak and water is able to quickly enter the binding site.

The analysis was repeated after discarding the initial 1 ns (of 6 ns) for the restrain and discharge stages, and the initial 3 ns (of 8 ns) for the vanish stage. Greater equilibration time was required for the vanish stage in order to allow water to enter the binding site.

<img src="images/convergence_equil.png" width="1000"/>

This removes the initial correlation between the time evolution of the repeat runs, and results in average estimates (blue line) for the vanish and discharge stages which show no clear drift with time, suggesting that we have achieved convergence. There appears to be some drift for the restrain stage, although the gradient is very shallow. Starting with the full set of sampled energies and discarding data in the equilibration periods mentioned above, we found $\Delta\textit{G}_{\textrm{Bind}}^o=-9.03 \pm 1.06$ kcal mol $^{-1}$ (95 % t-based confidence interval). Note that because `AlchemicalFreeEnergy.analyse()` does not discard data from the beginning of the simulation, this analysis shows that our previous results will be subject to bias resulting from including the unequilibrated portion of the trajectory.

Repeating this analysis in reverse and comparing with the analysis in the forward direction can be helpful for diagnosing insufficient equilibration time. More rigorous approaches for detecting equilibration are also available, such as that of [Chodera](https://pubs.acs.org/doi/abs/10.1021/acs.jctc.5b00784) and [Yang, Bitetti-Putzer, and Karplus](https://pubs.aip.org/aip/jcp/article-abstract/120/6/2618/597166/Free-energy-simulations-Use-of-reverse-cumulative?redirectedFrom=fulltext).

While reassuring, it should be noted that a lack of drift of the free energy estimate with respect to increasing simulation time does not guarantee converged sampling - it may be that transitions to other important configurations are not sampled on the timescale of the simulations.

<div class="alert alert-success">
<b>Exercise 4: Discrepancies Between Runs </b>
</div>

Although the average results based on all runs seem fairly stable, the individual trajectories do not appear to be converging to a single estimate, as we would expect if all relevant regions of configuration space were sampled at equilibrium during every run. Suggest an explanation for this.

 ```Python
 #FIXME: Answer the above question
 ```

<details><summary {style='color:green;font-weight:bold'}> Click here to see solution to Exercise. </summary>

 - The inter-run deviations show that individual trajectories have become trapped in metastable regions of configuration space - see [Figure 1 here](https://www.sciencedirect.com/science/article/pii/S0959440X20301366?casa_token=h0F3n8Vkdm4AAAAA:ZFmnAqOkA3lCeWBOhA-ZLWakJtUUOrOdK8wGnf1Xu8tsm6eQUR7OYC40FPPDriMiaOps47IDeA). To attempt to recover the correct estimate which would be provided by sampling of all of these states in equilibrium, we must run replicate simulations.
</details>

### 2.3 Are there Alternative Binding Poses?
<a id="poses"></a>

We have selected restraints based on a short unrestrained simulation of the fully interacting complex. These restraints closely confine the ligand to the binding pose sampled, and therefore our ABFE calculations make the implicit assumption that this is the single dominant binding pose. If there are other binding poses with similar or more favourable free energies of binding, our estimate based on a single binding pose will be wrong. 

In fact, we found another binding pose for MIF-180. The two poses are shown below:

<img src="images/alt_poses.png" width="1000"/>

While they appear similar, we found that interconversion occurred very infrequently on the timescale of the simulations. The bound leg calculations were repeated for the alternative binding pose, yielding a binding free energy of $\Delta\textit{G}_{\textrm{Bind}}^o=-9.06 \pm 1.12$ kcal mol $^{-1}$ (95 % t-based confidence interval).

Where there are alternative binding poses, [the free energies should be combined according to](https://aip.scitation.org/doi/full/10.1063/1.2221683):

\begin{align*}
\Delta G_{\text {Bind}}^{o}=-k_{B}T\ln[\exp(-\beta\Delta\textit{G}_{\textrm{Bind, 1}}^o) + \exp(-\beta\Delta\textit{G}_{\textrm{Bind, 2}}^o)]
\end{align*}

Where the subscripts 1 and 2 indicate alternative binding poses. 

<div class="alert alert-success">
<b>Exercise 5: Combining Results for Alternative Binding Poses </b>
</div>

Calculate the overall binding free energy by combining the estimates from the two different poses. The estimates from the two different poses are $\Delta\textit{G}_{\textrm{Bind}}^o=-9.03 \pm 1.06$ kcal mol $^{-1}$ and $-9.06 \pm 1.12$ kcal mol $^{-1}$ (95 % t-based confidence intervals). Propagate the error.


```python
# FIXME: Calculate the overall free energy of binding based on both binding poses
```

<details><summary {style='color:green;font-weight:bold'}> Click here to see solution to Exercise. </summary>

```python
dg1, dg2 = (-9.03, 1.06), (-9.06, 1.12) # kcal / mol, 95 % t-based C.I.s
beta = 1/0.593 # mol / kcal

# Get the overall free energy of binding
dg_overall = - (1/beta) * np.log(np.exp(-beta*dg1[0]) + np.exp(-beta*dg2[0]))

# Propagate the error by taking derivatives with respect to each of the dgs
deriv_dg1 = np.exp(-beta*dg1[0]) / (np.exp(-beta*dg1[0]) + np.exp(-beta*dg2[0]))
deriv_dg2 = np.exp(-beta*dg2[0]) / (np.exp(-beta*dg1[0]) + np.exp(-beta*dg2[0]))
uncert_overall = np.sqrt((deriv_dg1**2 * dg1[1]**2) + (deriv_dg2**2 * dg2[1]**2))

print(f"Combined dg_bind estimate = {dg_overall:.2f} +/- {uncert_overall:.2f} kcal /mol (95 % t-based C.I. based on 5 replicate runs)")
```
</details>

### 2.4 Do I Need Symmetry Corrections? 
<a id="sym"></a>

When the receptor-ligand restraints are imposed, sometimes a symmetry correction is required. However, symmetry corrections can be confusing and they are a common source of errors in the literature. [This article](https://aip.scitation.org/doi/full/10.1063/5.0046853) is recommended.

Symmetry corrections are required when there is some symmetry in the energy of they system which is broken by the introduction of restraints, but only when the symmetric energy minima are not sampled at all as the restraints are introduced. If there is equilibrium sampling of all symmetry-related regions of configurations space as the restraints are introduced then a correction is not required, and if there is partial sampling then this represents a sampling issue. For example, imagine we have a ligand with two symmetry-related binding modes. If restraints are introduced which completely prevent sampling of one of the poses, and this pose is not sampled at all as the restraints are turned on, then a symmetry correction of $-k_BT\ln2$ is required. This recovers the free energy change which would have been observed during the restraining simulations had the other binding pose been sampled.

<div class="alert alert-success">
<b>Exercise 6: The Final Result</b>
</div>

Look at the structure of MIF, a homotrimer:

<img src="images/mif_mif180.png" width="400"/>

Is a symmetry correction required for the restraints used for these simulations? Based on this, calculate your final estimate of $\Delta\textit{G}_{\textrm{Bind}}^o$, making sure to account for the alternative binding poses. Compare with experiment and your rough initial estimate.


```Python
# FIXME: Answer the question above
```

<details><summary {style='color:green;font-weight:bold'}> Click here to see solution to Exercise. </summary>

 - Yes, because MIF has three symmetrical binding sites, but our restraints artificially confine the ligand to a single site. If we were to ignore the remaining binding sites, we would underestimate the probability of the ligand binding by a factor of 3, and therefore underestimate the free energy of binding by $-k_BT\ln3 = 0.65 $ kcal mol $^{-1}$. Hence, the final free energy of binding is:      

 \begin{equation*}
-9.46 \pm 0.77 - 0.65 = -10.11 \pm 0.77 \;\mathrm{kcal} \; \mathrm{mol}^{-1} 
 \end{equation*}
        
- Superficially, this appears to be a much poorer answer than our rough initial estimate of $-8.65 \pm 0.23$ kcal mol $^{-1}$, which is closer to the experimental result ($-8.98 \pm 0.28$ kcal mol $^{-1}$) and has a smaller confidence interval. However, as we have discussed above, the "correct" answer for our model in the limit of infinite sampling will be different to the experimental value, and our final answer is a far more rigorous estimate of this value with a far more realistic estimate of the random error arising from incomplete sampling. While not within "chemical accuracy" (1 kcal mol $^{-1}$) of the experimental result, our final answer shows reasonable agreement with experiment.      
</details>

<div class="alert alert-success">
<b>Exercise 7: Other Symmetry Corrections </b>
</div>

Can you think of another case for this system where symmetry corrections might be required?

```Python
# FIXME: Answer the question above
```

<details><summary {style='color:green;font-weight:bold'}> Click here to see solution to Exercise. </summary>        

- In the first tutorial, the restraints were based on the phenol group. In the case that the following criteria are satisfied:      
    - There is a symmetry in the energy of the system (flip of the phenol gives equivalent energies)        
    - The restraints make one symmetrical orientation of the phenol group so high in energy that its contribution to the partition function is negligible       
    - The orientation raised in energy by the restraints is not sampled at all as the restraints are turned on      

A symmetry correction of $-k_BT\ln2$ should be applied to account for failure to sample the orientation of the phenol which is penalised by the restraints as the restraints are turned on.         
</details>  

### 2.5 Final Comments
<a id="comments"></a>

Please note that only a small subset of potential issues and suggested analyses have been discussed - for a more comprehensive guide, check out the [best practices article](https://livecomsjournal.org/index.php/livecoms/article/view/v2i1e18378/963).
