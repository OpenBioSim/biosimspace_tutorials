Author: Lester Hedges<br>
Email:&nbsp;&nbsp; lester.hedges@bristol.ac.uk

# Molecular dynamics

The companion notebook for this section can be found [here](https://github.com/michellab/BioSimSpaceTutorials/blob/4844562e7d2cd0b269cead56562ec16a3dfaef7c/01_introduction/03_molecular_dynamics.ipynb)

## Introduction

In this section we will learn how to use BioSimSpace to configure and run some basic molecular dynamics simulations.


## Protocols

One of the key goals of BioSimSpace was to start a conversation regarding _best practice_ within the biomolecular simulation community and to facilitate the codifying of shareable, re-usable, and extensible simulation protocols.

The [BioSimSpace.Protocol](https://biosimspace.org/api/index_Protocol.html) package defines protocols for a range of common molecular dynamics simulations. We can query the package to ee what protocols are available:


```python
import BioSimSpace as BSS
BSS.Protocol.protocols()
```




    ['Equilibration',
     'FreeEnergy',
     'Metadynamics',
     'Minimisation',
     'Production',
     'Steering']



Since we require protocols to be _interoperable_, the classes listed above are simple objects that allow you to configure a _limited_ set of options that are handled by _all_ of the molecular dynamics engines that we support. This might seem quite restrictive, but we will see later how it is possible to fully customise a simulation for a particular molecular dynamics engine.

Each protocol comes with some default options. To see what those are we can instantiate an object using the default constructor. For example, let's explore the [Equilibration](https://biosimspace.org/api/generated/BioSimSpace.Protocol.Equilibration.html#BioSimSpace.Protocol.Equilibration) protocol.


```python
protocol = BSS.Protocol.Equilibration()
print(protocol)
```

    <BioSimSpace.Protocol.Equilibration: timestep=2.0000 fs, runtime=0.2000 ns, temperature_start=300.0000 K, temperature_end=300.0000 K, pressure=None, report_interval=100, restart_interval=500,restraint=None>


Here we can see that the default protocol performs an equlibration at fixed temperature (`temperature_start == temperature_end`) in the NVT ensemble (`pressure=None`) with no restraints (`restraint=None`). The total simulation time is 0.2 nanoseconds with an integration timestep of 2 femtoseconds. The `report_interval` and `restart_interval` govern how frequently information is written to log and restart (and/or trajectory) files respectively.

If any of these defaults are unsuitable, then you are free to change them by passing in appropriate values for each of the arguments when instantiating the object. In some cases it might be desirable to _override_ the default protocols and set specific values of the arguments that are suitable for a particular project or team. This can be achieved by defining a set of function wrappers that configure and return the protocols using your own defaults.

As an example, the following configuration could be used to provide alternative defaults for NVT and NPT equlibration protocols. (You could simply re-use the existing protocol name, but here we provide two separate protocols for convenience.)

```python
# myconfig/Protocol.py
import BioSimSpace as BSS

# Override the equilibration protocol with some custom defaults. Ideally all
# arguments to the BioSimSpace function would be mapped, but here we use a
# subset for simplicity.

# A custom equlibration in the NVT ensemble.
def EquilibrationNVT(runtime=5*BSS.Units.Time.nanosecond,
                     report_interval=2500,
                     restart_interval=250000,
                     restraint="backbone"):
    return BSS.Protocol.Equilibration(runtime=runtime,
                                      report_interval=report_interval,
                                      restart_interval=restart_interval,
                                      restraint=restraint)

# A custom equlibration in the NPT ensemble.
def EquilibrationNPT(runtime=5*BSS.Units.Time.nanosecond,
                     pressure=BSS.Units.Pressure.atm,
                     report_interval=2500,
                     restart_interval=25000,
                     restraint="backbone"):
    return BSS.Protocol.Equilibration(runtime=runtime,
                                      pressure=pressure,
                                      report_interval=report_interval,
                                      restart_interval=restart_interval,
                                      restraint=restraint)

```

We could then import the customised protocols from our local configuraton and use them instead, e.g.:


```python
from myconfig.Protocol import *

protocol = EquilibrationNVT()
print(protocol)
```

    <BioSimSpace.Protocol.Equilibration: timestep=2.0000 fs, runtime=5.0000 ns, temperature_start=300.0000 K, temperature_end=300.0000 K, pressure=None, report_interval=2500, restart_interval=250000,restraint='backbone'>


## Processes

Once you have created a molecular system and chosen a protocol, then it is time to create a simulation _process_. The [BioSimSpace.Process](https://biosimspace.org/api/index_Process.html) package provides functionality for configuring and running processes with several common molecular dynamics engines.

Let's query the package to see what engines are available:


```python
BSS.Process.engines()
```




    ['Amber', 'Gromacs', 'Namd', 'OpenMM', 'Somd']



Before creating a process let us once again load our example alanine-dipeptide system from file:


```python
system = BSS.IO.readMolecules("inputs/ala*")
```

As a simple example, let us use a short minimisation protocol:


```python
protocol = BSS.Protocol.Minimisation(steps=1000)
```

We'll now create a process to apply the `Protocol` to the `System` using the AMBER molecular dynamics engine:


```python
process = BSS.Process.Amber(system, protocol)
```

A lot of complexity is hidden in this line. BioSimSpace has automatically found an AMBER executable on the underlying operating system, has automatically written AMBER format molecular input files, generated an AMBER configuration file for the minimisation protocol, and configured any command-line arguments that are required.

By default, processes are run inside of a temporary working directory hidden away from the user. To see where this is, run:


```python
process.workDir()
```




    '/tmp/tmpjuzioj_o'



N.B. If you want to use a different temporary directory, e.g. one with a faster disk, then simply set the `TMPDIR` environment variable. Alternatively, you can pass the `work_dir` argument to the `Process` constructor to explicitly specify the path. This can be useful when you want named directories, or want to examine the intermediate files from the `Process` for debugging purposes.

To see what executable was found, run:


```python
process.exe()
```




    '/home/lester/sire.app/bin/sander'



To see the list of autogenerated input files:


```python
process.inputFiles()
```




    ['/tmp/tmpjuzioj_o/amber.cfg',
     '/tmp/tmpjuzioj_o/amber.rst7',
     '/tmp/tmpjuzioj_o/amber.prm7']



If you like, we could zip up the input files to use on another occasion. When working on a notebook server it's possible to return a file link so that we can download them:


```python
process.getInput(file_link=True)
```




<a href='amber_input.zip' target='_blank' download='amber_input.zip'>amber_input.zip</a><br>



We can query also query the list of configuration file options:


```python
process.getConfig()
```




    ['Minimisation',
     ' &cntrl',
     '  imin=1,',
     '  ntx=1,',
     '  ntxo=1,',
     '  ntpr=100,',
     '  irest=0,',
     '  maxcyc=1000,',
     '  ncyc=1000,',
     '  cut=8.0,',
     ' /']



And also get command-line argument string for the process:


```python
process.getArgString()
```




    '-O -i amber.cfg -p amber.prm7 -c amber.rst7 -o stdout -r amber.crd -inf amber.nrg'



If you're an expert in a particular package then BioSimSpace allows you to fully customise the process by tweaking the configuration options and command-line arguments. Read the help documentation for `process.setConfig` and `process.setArgs` if you are interested. Once again, it's possible to wrap the instantiation of `Process` objects in your own custom functions, allowing you to tweak the default configuration options for your own requirements. For example, if you always want to wrap coordinates to the minimum image when using AMBER, then this could be achieved as follows:

```python
# myconfig/Process.py
import BioSimSpace as BSS

# Wrap the instantiation of BSS.Process.Amber objects to configure them
# such that coordinates are always wrapped to the minimum image.
def Amber(system, protocol, exe=None, name="amber",
            work_dir=None, seed=None, property_map={}):

    # Create process using the passed parameters.
    process = BSS.Process.Amber(system,
                                protocol,
                                exe=exe,
                                name=name,
                                work_dir=work_dir,
                                seed=seed,
                                property_map=property_map)
    
    # Get the config.
    config = process.getConfig()
    
    # Add coordinate wrapping to the end of the config.
    config[-1] = "  iwrap=1,"
    config.append(" /")

    # Set the new config.
    process.setConfig(config)
    
    # Return the process.
    return process
```

Let us know create our custom AMBER process and check the configuration:


```python
from myconfig.Process import *

process = Amber(system, protocol)
process.getConfig()
```




    ['Minimisation',
     ' &cntrl',
     '  imin=1,',
     '  ntx=1,',
     '  ntxo=1,',
     '  ntpr=100,',
     '  irest=0,',
     '  maxcyc=1000,',
     '  ncyc=1000,',
     '  cut=8.0,',
     '  iwrap=1,',
     ' /']



N.B. You might want to add additional configuration details to your `Process` wrappers, e.g. to ensure that a specific executable is used.

Now that we have a process, let's go ahead and start it:


```python
process.start()
```




    BioSimSpace.Process.Amber(<BioSimSpace.System: nMolecules=631>, <BioSimSpace.Protocol.Custom>, exe='/home/lester/sire.app/bin/sander', name='amber', work_dir='/tmp/tmpm8kx_jms', seed=None)



BioSimSpace has now launched a minimisation process in the background! When in an interactive session you carry on working and periodically check in on the process to see how its doing.

To check whether the process is running:


```python
process.isRunning()
```




    True



We can see how many minutes it has been running for:


```python
process.runTime()
```




    0.1960 mins



Since this is a short minimisation it will likely finish pretty quickly. Let's print the final energy of the system and return the minimised molecular configuration.


```python
print(process.getTotalEnergy(block=True))
minimised = process.getSystem()
```

    -6954.7000 kcal/mol


When working interactively, any time we query a running process we get back the _latest_ information that has been written to disk. This means that we can get an update on how things are progressing, then immediately carry on with what we were doing in our notebook. By passing `block=True`, as we do when we call `getTotalEnergy` above, we request that the process finishes running before returning a result. This means we get the _final_ energy, and the minimised system that is returned afterwards represents the _final_ snapshot that was saved.

Let's now re-run the simulation, instead using GROMACS as the MD engine.


```python
process = BSS.Process.Gromacs(system, protocol)
```

When the process is instantiated, BioSimSpace takes the system that was read from AMBER format files and converts it to GROMACS format ready for simulation. Let's take a look at the list of input files that were autogenerated for us:


```python
process.inputFiles()
```




    ['/tmp/tmpsxf5hoej/gromacs.mdp',
     '/tmp/tmpsxf5hoej/gromacs.gro',
     '/tmp/tmpsxf5hoej/gromacs.top',
     '/tmp/tmpsxf5hoej/gromacs.tpr']



Let's start the process running and, once again, wait for it to finish before getting the minimised system.


```python
process.start()
minimised = process.getSystem(block=True)
```

## Interactive molecular dynamics

The example in the previous section was finished almost as soon as it began. Let's run a more complicated equilibration protocol so that we can learn more about how to monitor processes interactively using BioSimSpace.


```python
protocol = BSS.Protocol.Equilibration(runtime=20*BSS.Units.Time.picosecond,
                                      temperature_start=0*BSS.Units.Temperature.kelvin,
                                      temperature_end=300*BSS.Units.Temperature.kelvin,
                                      restraint="backbone")
```

This protocol will equlibrate a system for 20 picoseconds, while heating it from 0 to 300 Kelvin and restraining any atoms in the backbone of the molecule. Note that some of the parameters passed have units, e.g. the temperatures are in Kelvin. BioSimSpace has a built in type system for handling variables with units. The `BSS.Units` package provides a convenient way of declaring these, for example `10*BSS.Units.Temperature.kelvin` creates an object of type `BSS.Types.Temperature` with a magnitude of 10 and unit of Kelvin. This allows the user to pass parameters with whatever unit they like. BioSimSpace will simply convert it to the correct unit for the chosen MD engine internally.

One again, we now need a `Process` in order to run our simulation. Exectute the cell below to initialise an AMBER process and start it immediately. Note that we pass in the minimised system from the last example, along with our new protocol.


```python
process = BSS.Process.Amber(minimised, protocol).start()
```

We can monitor the time, temperature, and energy as the process runs. If you run this multiple times using "CTRL+Return" you'll see the temperature slowly increasing.


```python
print(process.getTime(), process.getTemperature(), process.getTotalEnergy())
```

    2.4000 ps 30.6900 K -6936.6583 kcal/mol


Since all of the values returned above are typed we can easily convert them to other units:


```python
print(process.getTime().nanoseconds(), process.getTemperature().celsius(), process.getTotalEnergy().kj_per_mol())
```

    0.0050 ns -204.7600 C -2.7995e+04 kJ/mol


It's possible to query many other thermodynamic records. What's available depends on type of protocol and the MD package that is used to run the protocol. To get more information, run:

N.B. Certain functionality is specific to the process in question, i.e. `BSS.Process.Amber` will have different options to `BSS.Process.Gromacs`, but, for the purposes of interoperability, there is a core set of functionality that is consistent across all `Process` classes, e.g. all classes implement a `getSystem` method.)

### Plotting time series data

As well as querying the most recent records we can also get a time series of results by passing the `time_series` keyword argument to any of the data record getter methods, e.g.

```python
# Get a time series of pressure records.
pressure = process.getPressure(time_series=True)
```

The `BSS.Notebook` package provides several useful tools that are available when working inside of a Jupyter notebook. One of these is the plot function, that allows us to create simple x/y plots of time-series data.

Let's grab the same record data as above and use it to make some graphs of the data.


```python
# Generate a plot of time vs temperature.
plot1 = BSS.Notebook.plot(process.getTime(time_series=True), process.getTemperature(time_series=True))

# Generate a plot of time vs energy.
plot2 = BSS.Notebook.plot(process.getTime(time_series=True), process.getTotalEnergy(time_series=True))
```

![Time-series plots](https://github.com/michellab/BioSimSpaceTutorials/blob/dd5a24e58778af21612ade7febe5ba7fd98f9885/01_introduction/assets/03_time_series.png)
    
(Note that, by default, the axis labels axis labels are automatically generated from the types and units of the x and y data that are passed to the function.)

Re-run the cell using "CTRL+Return" to see the graphs update as the simulation progesses. (Occasionally, you might see a warning that the x and y data sets are mismatched in length, this is because the data was extracted before all records were written to disk.)

Being able to query a process in real time is an incredibly useful tool. This could enable us to check for convergence, or spot errors in the simulation. If you ever need to kill a running process (perhaps it was configured incorrectly), run:

```python
process.kill()
```

### Visualising the molecular system

Another useful tool that is available when working inside of a notebook is the `View` class that can be used to visualise the molecular system while a process is running. To create a `View` object we must attach it to a process (or a molecular system), e.g.:


```python
view = BSS.Notebook.View(process)
```

We can now visualise the system:


```python
view.system()
```

![Visualise the system](https://github.com/michellab/BioSimSpaceTutorials/blob/dd5a24e58778af21612ade7febe5ba7fd98f9885/01_introduction/assets/03_view_system.png)

(If you see an empty view, try re-executing the cell.)

To only view a specific molecule:


```python
view.molecule(0)
```

![Visualise a molecule](https://github.com/michellab/BioSimSpaceTutorials/blob/86442df77e2ad33ae79f62e214a53af42cb320ec/01_introduction/assets/03_view_molecule.png)

To view a list of molecules:


```python
view.molecules([0, 5, 10])
```

![Visualise some molecules](https://github.com/michellab/BioSimSpaceTutorials/blob/86442df77e2ad33ae79f62e214a53af42cb320ec/01_introduction/assets/03_view_molecules.png)


If a particular view was of interest it can be reloaded as follows:


```python
# Reload the original view.
view.reload(0)
```

![Visualise the system](https://github.com/michellab/BioSimSpaceTutorials/blob/dd5a24e58778af21612ade7febe5ba7fd98f9885/01_introduction/assets/03_view_system.png)

To save a specific view as a PDB file:


```python
view.savePDB("my_view.pdb", index=0)
```

### Reading and analysing trajectory data¶ 

The `BSS.Trajectory` package comes with a set of tools for reading and analysis trajectory files. Files can be loaded directly, or if supported, can be read from a running process.

For example, to get the trajectory from the process, run:


```python
traj = process.getTrajectory()
```

(If you get an error, then the trajectory file may be in the process of being written. Simply try again.)

To get the current number of frames:


```python
traj.nFrames()
```




    20



To get all of the frames as a list of `System` objects:


```python
frames = traj.getFrames()
```

Specific frames can be extracted by passing a list of indices, e.g. the first and last:


```python
frames = traj.getFrames([0, -1])
```

Like most things in BioSimSpace, the `Trajectory` class is simply a wrapper around existing tools. Internally, trajectories are stored as an [MDTraj](http://mdtraj.org) object. This can be obtained, allowing the user direct access to the full power of MDTraj:


```python
mdtraj = traj.getTrajectory()
type(mdtraj)
```




    mdtraj.core.trajectory.Trajectory



Alternatively, a trajectory can be returned in [MDAnalysis](https://www.mdanalysis.org) format:


```python
mdanalysis = traj.getTrajectory(format="mdanalysis")
type(mdanalysis)
```




    MDAnalysis.core.universe.Universe



The `Trajectory` class also provides wrappers around some basic MDTraj analysis tools, allowing the user to compute quantities such as the root mean squared displacement (RMSD).

Let's measure the RMSD of the alanine-dipeptide molecule with a reference to its configuration in the first trajectory frame. To extract the alanine-dipeptide, we search the system for a residue named ALA. We'll also plot the RMSD for each frame of the trajectory.



```python
# Search the system for a residue named ALA. Since there is a single match,
# we take the first result.
molecule = system.search("mol with resname ALA")[0]

# Get the indices of the atoms in the molecule, relative to the original system.
indices = [system.getIndex(x) for x in molecule.getAtoms()]

# Compute the RMSD for each frame and plot the result.
BSS.Notebook.plot(y=process.getTrajectory().rmsd(frame=0, atoms=indices),
                  xlabel="Frame", ylabel="RMSD")
```

![RMSD vs frame index](https://github.com/michellab/BioSimSpaceTutorials/blob/dd5a24e58778af21612ade7febe5ba7fd98f9885/01_introduction/assets/03_rmsd.png)
