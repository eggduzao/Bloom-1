import os
import sys
import subprocess
from pathlib import Path

def run_on_remote_node(node: str, script_path: str, operation: str):
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
    conda_path = "/sw/miniconda3/etc/profile.d/conda.sh"
    bashrc_path = "/home/egusmao/.bashrc"
    conda_env = "ml"

    command = (
        f"nohup bash -c 'source {conda_path} && "
        f"source {bashrc_path} && "
        f"cd {script_directory} && "
        f"conda activate {conda_env} && "
        f"python {script_path} {operation}' "
        f"> {script_directory}/{node}_download.txt 2>&1 &"
    )

    ssh_command = ["ssh", node, command]

    print(f"Running remotely on node {node}: {command}")
    subprocess.run(ssh_command)

# Example of usage:
if __name__ == "__main__":

    compute_node = "node1"
    operation = "make"
    script = "/storage2/egusmao/projects/Bloom/bloom/data/download_geo.py"

    run_on_remote_node(compute_node, script, operation)

