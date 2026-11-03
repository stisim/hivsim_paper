#!/bin/bash

#SBATCH --nodes=1

#SBATCH --ntasks-per-node=1

#SBATCH --cpus-per-task=1

#SBATCH --time=05-00:00:00
#SBATCH --partition=cpu_medium
#SBATCH --mem=8GB

#SBATCH --job-name=MWI_HIVControl
#SBATCH --chdir=/gpfs/data/bershteynlab/EMOD/citrod01/MWI_calib_20220422

#SBATCH --mail-type=END

#SBATCH --mail-user=daniel.citron@nyulangone.org

#SBATCH --output=/gpfs/data/bershteynlab/EMOD/citrod01/MWI_calib_20220422/slurm_%j.out

#python run_scenarios.py -c optim_script_transmission.py --resample-method provided --samples resampled_parameter_sets_20230424.csv --output-dir MWI_transmission --suite-name MWI_transmission --table scenarios_transmission.csv

python run_scenarios.py -c optim_script.py --resample-method provided --samples resampled_parameter_sets_20230424.csv --output-dir hivcontrol --suite-name hivcontrol --table scenarios_ai.csv