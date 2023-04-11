import argparse
import os
import tomli

parser = argparse.ArgumentParser(description='Runs a simulation')

parser.add_argument('--sim_name', type=str)  # Name of simulation
parser.add_argument('--run_local', action='store_true')  # Whether to run on cluster or locally
args = parser.parse_args()

for entry in os.scandir(f'configs/{args.sim_name}_configs'):

    with open(f'configs/{args.sim_name}_configs/{entry.name}', 'rb') as f:
        config = tomli.load(f)
    name_split = entry.name.split('_')
    config_index = name_split[-1].split('.')[0][4:]
    if args.run_local:
        os.system(f'bash _localbash {args.sim_name} {config["reals"]} {config_index}')
    else:
        os.system(f'sbatch --array=1-{config["reals"]} submit_to_queue {args.sim_name} {config_index}')
