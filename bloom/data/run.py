import os
import sys
import subprocess
from pathlib import Path

def run_with_nohup(script_path: str, operation: str):
    """
    Runs the specified Python script on a remote node via SSH.
    
    Parameters:
    -----------
    node : str
        Hostname of the node to run the script on (e.g., "compute02").
    script_path : str
        Absolute path of the Python script on the remote node.
    """

    # Paths
    script_path = Path(script_path).resolve()
    script_directory = script_path.parent

    # Existing cluster paths
    micromamba_path = "$(micromamba shell hook --shell bash)"
    bashrc_path = "/home/egusmao/.bashrc"
    conda_env = "bio"

    command = (
        f"nohup bash -c 'eval {micromamba_path} && "
        f"source {bashrc_path} && "
        f"cd {script_directory} && "
        f"conda activate {conda_env} && "
        f"python {script_path} {operation}' "
        f"> {operation}.log 2>&1 &"
    )

    print(f"Running: {command}")
    subprocess.run(command)

# Example of usage:
if __name__ == "__main__":

    script = "/storage2/egusmao/projects/Bloom/bloom/data/download_geo.py"
    operation = "metadata" # metadata download make run merge delete delete_all_sra_files

    run_with_nohup(compute_node, script, operation)

