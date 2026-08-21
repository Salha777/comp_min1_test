#!/bin/bash
# =============================================================================
# SLURM job script — Stage 1 MD (CLASSICAL, CPU)
# LAMMPS stock CPU build: lammps/20230802.3-openmpi-5.0.3-gcc-12.2.0
#   -> has MOLECULE + KSPACE + RIGID (needed for atom_style full / pppm / fix shake)
# NOTE: the lab GPU build (lammps-mace-plumed-gpu) lacks MOLECULE, so it CANNOT run
#       this classical stage. Keep it for the later MACE metadynamics stage instead.
# =============================================================================
#
# Submit with:  sbatch submit_stage1.sh
# Monitor:      squeue -u $USER
# Cancel:       scancel <jobid>
# =============================================================================

#SBATCH -J mg_stage1
#SBATCH -o mg_stage1_%j.out
#SBATCH -e mg_stage1_%j.err
#SBATCH -n 32                      # MPI ranks; ~15k-atom box scales fine to 32-64 cores
#SBATCH --nodes=1                  # start on one node; raise only if you go multi-node
#SBATCH --mem-per-cpu=2G
#SBATCH -t 240:0:0                 # trim to your CPU partition's max walltime if needed
# CPU job: NO --gres / --constraint / GPU account. Apocrita CPU jobs use the default
# partition; add e.g. "#SBATCH -p short" ONLY if your cluster requires an explicit one.

mkdir -p logs

# --- Modules -----------------------------------------------------------------
module purge
module load lammps/20230802.3-openmpi-5.0.3-gcc-12.2.0

# Fail fast if lmp isn't on PATH after the module load
if ! command -v lmp &>/dev/null; then
    echo "ERROR: 'lmp' not found after 'module load'. Check the module name with 'module avail lammps'."
    exit 1
fi
echo "Using LAMMPS: $(which lmp)"

# --- Run ---------------------------------------------------------------------
echo "=== Starting Stage 1 MD (CPU): $(date) ==="
echo "Job $SLURM_JOB_ID on $SLURM_JOB_NODELIST  |  $SLURM_NTASKS MPI ranks"

mpirun -np "$SLURM_NTASKS" lmp -in stage1.in \
    2>&1 | tee "logs/stage1_${SLURM_JOB_ID}.log"

# Real lmp exit status — PIPESTATUS[0], NOT $? (which would be tee's status = always 0).
# This is the fix for the earlier "exited with code 0 / completed successfully" on a hard error.
EXIT_CODE=${PIPESTATUS[0]}
echo "=== LAMMPS exited with code $EXIT_CODE: $(date) ==="

# --- Post-run checks ---------------------------------------------------------
if [ "$EXIT_CODE" -eq 0 ]; then
    echo "Run completed successfully. Output files:"
    ls -lh stage1_prod.lammpstrj stage1_prod_final.data thermo_prod.dat 2>/dev/null
else
    echo "ERROR: LAMMPS exited non-zero. Check logs/stage1_${SLURM_JOB_ID}.log"
fi

# =============================================================================
# WORKFLOW SEQUENCE
# =============================================================================
# Before submitting:
#   1.  packmol < packmol.inp        # generates stage1_packed.pdb
#   2.  python build_topology.py     # generates stage1.data + mg_o_table.dat
#   3.  sbatch submit_stage1.sh      # this script
#
# To restart a timed-out run:
#   Change "read_data stage1.data" -> "read_restart stage1_equil.restart"
#   (or stage1_prod_final.restart) in stage1.in, then resubmit.
#
# Perf note: CPU is slower than the A100 GPU build would be, but the GPU build lacks
# MOLECULE so it isn't an option for the classical stage. If throughput is too low,
# raise -n (and --nodes) — pppm + lj/cut/coul/long parallelise well across cores.
# =============================================================================
