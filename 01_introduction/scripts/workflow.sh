#!/usr/bin/env bash

# Exit immediately on error.
set -e

# Remove existing output.
rm -f parameterised.*
rm -f solvated.*
rm -f minimised.*

echo "Parameterising..."
python nodes/parameterise.py --pdb inputs/methanol.pdb --forcefield gaff

echo "Solvating..."
python nodes/solvate.py --files parameterised.* --water_model tip3p

echo "Minimising..."
python nodes/minimise.py --files solvated.* --steps 1000

echo "Equilibrating..."
python nodes/equilibrate.py --files minimised.* --restraint heavy

echo "Done!"
