import os
import sys
import subprocess
from pathlib import Path

    def run_on_remote_node(node: str, script_path: str):
        """
        Runs the specified Python script on a remote node via SSH.
        
        Parameters:
        -----------
        node : str
            Hostname of the node to run the script on (e.g., "compute02").
        script_path : str
            Absolute path of the Python script on the remote node.
        """

        command = (
            f"nohup python {script_path} "
            f"> nodeX_download.txt 2>&1 &"
        )

        ssh_command = ["ssh", node, command]

        print(f"Running remotely on node {node}: {command}")
        subprocess.run(ssh_command)

# Example of usage:
if __name__ == "__main__":

    compute_node = "node3"
    script = "/storage2/egusmao/projects/Bloom/bloom/data/download_geo.py"

    run_on_remote_node(compute_node, script)

