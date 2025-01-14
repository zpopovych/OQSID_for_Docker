#!/usr/bin/env bash
set -ex

# Check if the volume mount has the license file
if [ -f /results/mosek.lic ]; then
  echo "Copying MOSEK license to /root/mosek..."
  mkdir -p /root/mosek
  cp /results/mosek.lic /root/mosek/mosek.lic
else
  echo "MOSEK license file not found!" >&2
  exit 1
fi

# Define the directory where the Julia scripts are located
#SCRIPT_DIR="code"

# Define an array of Julia script names
JULIA_SCRIPTS=(
  "Kossak_CONSTR_POP.jl"
  "Kossak_SLSQP.jl"
  "Lindblad_CONSTR_POP.jl"
  "Lindblad_SLSQP.jl"
  "NonMarkovianity.jl"
  "Kossak_CONSTR_POP_MULT_DUR_TRAIN.jl"
  "Lindblad_CONSTR_POP_MULT_DUR_TRAIN.jl"
  "LSID_DMD_Bloch4.jl"
  "LSID_DMDvsERA.jl"
)

# Define the Python script to run after Julia scripts
PYTHON_SCRIPTS=(
"FIG_1_Plot_POP_Violins.py"
"FIG_2_Plot_by_train_duration.py"
"FIG_3_Plot_LSID_Violins.py"
)

# Run each Julia script
for J_SCRIPT in "${JULIA_SCRIPTS[@]}"; do
  echo "Running $J_SCRIPT..."
  julia "$J_SCRIPT"
  if [ $? -ne 0 ]; then
    echo "Error occurred while running $J_SCRIPT"
    exit 1
  fi
done

echo "Julia scripts executed successfully."

# Run each Python script
for P_SCRIPT in "${PYTHON_SCRIPTS[@]}"; do
  echo "Running $P_SCRIPT..."
  python "$P_SCRIPT"
  if [ $? -ne 0 ]; then
    echo "Error occurred while running $P_SCRIPT"
    exit 1
  fi
done


echo "All scripts executed successfully."