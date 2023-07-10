from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from metadynamics import *
from glob import glob
import os
import shutil

# Load the topology and coordinate files.
prmtop = AmberPrmtopFile("openmm.prm7")
inpcrd = AmberInpcrdFile("openmm.rst7")

# Initialise the molecular system.
system = prmtop.createSystem(
    nonbondedMethod=PME, nonbondedCutoff=1 * nanometer, constraints=HBonds
)

# Add a barostat to run at constant pressure.
barostat = MonteCarloBarostat(1.01325 * bar, 300.0 * kelvin)
system.addForce(barostat)

# Set funnel variables.
p1 = [1049, 1068, 1204, 1223, 1237, 1256, 2603, 2625, 2641]
p2 = [
    537,
    552,
    571,
    590,
    601,
    615,
    626,
    637,
    649,
    1272,
    1298,
    1305,
    1324,
    1331,
    1406,
    1472,
    1491,
    1498,
    1897,
    1904,
    2099,
    2115,
    2134,
    2582,
    2589,
    2603,
    2625,
]
lig = [x for x in range(0, 18)]

# Create the bias variable for the funnel projection.
projection = CustomCentroidBondForce(3, "distance(g1,g2)*cos(angle(g1,g2,g3))")
projection.addGroup(lig)
projection.addGroup(p1)
projection.addGroup(p2)
projection.addBond([0, 1, 2])
projection.setUsesPeriodicBoundaryConditions(True)
sigma_proj = 0.025
proj = BiasVariable(projection, 0.3, 3.7, 0.025, False, gridWidth=200)

# Create the bias variable for the funnel extent.
extent = CustomCentroidBondForce(3, "distance(g1,g2)*sin(angle(g1,g2,g3))")
extent.addGroup(lig)
extent.addGroup(p1)
extent.addGroup(p2)
extent.addBond([0, 1, 2])
extent.setUsesPeriodicBoundaryConditions(True)
sigma_ext = 0.05
ext = BiasVariable(extent, 0.0, 0.95, 0.05, False, gridWidth=200)

# Add restraints.
k1 = 10000 * kilojoules_per_mole
k2 = 1000 * kilojoules_per_mole
lower_wall = 0.5 * nanometer
upper_wall = 3.5 * nanometer

# Upper wall.
upper_wall_rest = CustomCentroidBondForce(
    3, "(k/2)*max(distance(g1,g2)*cos(angle(g1,g2,g3)) - upper_wall, 0)^2"
)
upper_wall_rest.addGroup(lig)
upper_wall_rest.addGroup(p1)
upper_wall_rest.addGroup(p2)
upper_wall_rest.addBond([0, 1, 2])
upper_wall_rest.addGlobalParameter("k", k1)
upper_wall_rest.addGlobalParameter("upper_wall", upper_wall)
upper_wall_rest.setUsesPeriodicBoundaryConditions(True)
system.addForce(upper_wall_rest)

# Sides of the funnel.
wall_width = 0.6 * nanometer
wall_buffer = 0.15 * nanometer
beta_cent = 1.5
s_cent = 2.0 * nanometer
dist_restraint = CustomCentroidBondForce(
    3,
    "(k/2)*max(distance(g1,g2)*sin(angle(g1,g2,g3)) - (a/(1+exp(b*(distance(g1,g2)*cos(angle(g1,g2,g3))-c)))+d), 0)^2",
)
dist_restraint.addGroup(lig)
dist_restraint.addGroup(p1)
dist_restraint.addGroup(p2)
dist_restraint.addBond([0, 1, 2])
dist_restraint.addGlobalParameter("k", k2)
dist_restraint.addGlobalParameter("a", wall_width)
dist_restraint.addGlobalParameter("b", beta_cent)
dist_restraint.addGlobalParameter("c", s_cent)
dist_restraint.addGlobalParameter("d", wall_buffer)
dist_restraint.setUsesPeriodicBoundaryConditions(True)
system.addForce(dist_restraint)

# Lower wall.
lower_wall_rest = CustomCentroidBondForce(
    3, "(k/2)*min(distance(g1,g2)*cos(angle(g1,g2,g3)) - lower_wall, 0)^2"
)
lower_wall_rest.addGroup(lig)
lower_wall_rest.addGroup(p1)
lower_wall_rest.addGroup(p2)
lower_wall_rest.addBond([0, 1, 2])
lower_wall_rest.addGlobalParameter("k", k1)
lower_wall_rest.addGlobalParameter("lower_wall", lower_wall)
lower_wall_rest.setUsesPeriodicBoundaryConditions(True)
system.addForce(lower_wall_rest)

# Initialise the metadynamics object.
bias = 10.0
meta = Metadynamics(
    system,
    [proj, ext],
    300.0 * kelvin,
    10.0,
    1.5 * kilojoules_per_mole,
    1000,
    biasDir=".",
    saveFrequency=1000,
)

# Define the integrator.
integrator = LangevinIntegrator(300.0 * kelvin, 1 / picosecond, 0.002 * picoseconds)

# Set the simulation platform.
platform = Platform.getPlatformByName("CUDA")
properties = {"CudaDeviceIndex": "0"}

# Initialise and configure the simulation object.
simulation = Simulation(prmtop.topology, system, integrator, platform, properties)
simulation.context.setPositions(inpcrd.positions)

# Setting intial system velocities.
simulation.context.setVelocitiesToTemperature(300.0)

# Look for a restart file.
if os.path.isfile("openmm.chk"):
    simulation.loadCheckpoint("openmm.chk")
    shutil.copy("openmm.out", "old_openmm.out")
    sim_log_file = [line[:-2] for line in open("openmm.out").readlines()]
    current_steps = int(sim_log_file[-1].split(",")[1])
    steps -= current_steps
    shutil.copy("COLVAR.npy", "old_COLVAR.npy")
    shutil.copy("HILLS", "old_HILLS")
    shutil.copy("openmm.dcd", "old_openmm.dcd")

# Add reporters.
simulation.reporters.append(DCDReporter("openmm.dcd", 50000))
simulation.reporters.append(
    StateDataReporter(
        "openmm.log",
        1000,
        step=True,
        time=True,
        potentialEnergy=True,
        kineticEnergy=True,
        totalEnergy=True,
        volume=True,
        temperature=True,
        totalSteps=True,
        separator=" ",
    )
)
simulation.reporters.append(CheckpointReporter("openmm.chk", 1000))

# Create PLUMED compatible HILLS file.
file = open("HILLS", "w")
file.write("#! FIELDS time pp.proj pp.ext sigma_pp.proj sigma_pp.ext height biasf\n")
file.write("#! SET multivariate false\n")
file.write("#! SET kerneltype gaussian\n")

# Initialise the collective variable array.
current_cvs = np.array(
    list(meta.getCollectiveVariables(simulation)) + [meta.getHillHeight(simulation)]
)
colvar_array = np.array([current_cvs])

# Write the inital collective variable record.
line = colvar_array[0]
time = 0
write_line = f"{time:15} {line[0]:20.16f} {line[1]:20.16f}          {sigma_proj}           {sigma_ext} {line[2]:20.16f}            {bias}\n"
file.write(write_line)

# Run the simulation.
steps = 50000000
cycles = 50000
steps_per_cycle = int(50000000 / cycles)
for x in range(0, cycles):
    meta.step(simulation, steps_per_cycle)
    current_cvs = np.array(
        list(meta.getCollectiveVariables(simulation)) + [meta.getHillHeight(simulation)]
    )
    colvar_array = np.append(colvar_array, [current_cvs], axis=0)
    np.save("COLVAR.npy", colvar_array)
    line = colvar_array[x + 1]
    time = int((x + 1) * 0.002 * steps_per_cycle)
    write_line = f"{time:15} {line[0]:20.16f} {line[1]:20.16f}          {sigma_proj}           {sigma_ext} {line[2]:20.16f}            {bias}\n"
    file.write(write_line)
