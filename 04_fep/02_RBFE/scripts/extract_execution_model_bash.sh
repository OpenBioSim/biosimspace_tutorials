#! /bin/bash 
# script to turn the execution model files into bash suitable arrays and variables.

# remove any ^M from end of file lines
dos2unix "$lig_file"
dos2unix "$net_file"
dos2unix "$prot_file"

# For all the ligands in ligands.dat, make an array
lig_array=()
while read lig; do 
lig_clean=$(sed 's/\r$//' <<< $lig); lig_array+=("${lig_clean}");
done < $lig_file

# Make a list of the transformations from the network.dat file
trans_array=()
eng_array=()
win_array=()
IFS=' '
while read trans; do
while read -a tra; do tran=${tra[0]}~${tra[1]}; eng=${tra[-1]}; win=${tra[2]};
trans_array+=("$tran"); eng_array+=("$eng"); win_array+=("$win"); done <<< $trans
done < $net_file

# check that the trans, eng, and win arrays are the same length.
