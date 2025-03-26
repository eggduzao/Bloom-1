
# Import
import os
import sys
import subprocess
from pathlib import Path

# Run with nohup
def run_with_nohup(script_path: str, operation: str):
    """
    Runs the specified Python script on the head node.
    
    Parameters:
    -----------
    script_path : str
        Absolute path of the Python script on the remote node.
    operation : str
        The operation to perform: "metadata", "download", "make", "run",
        "merge", "delete" or "delete_all_sra_files".
    """

    # Resolve script path
    script_path = Path(script_path).resolve()

    # Construct command as a proper shell string
    command = f"nohup python {script_path} {operation} > {operation}.log 2>&1 &"

    # Log message
    print(f"Running: {command}")

    # Run the command inside a shell
    subprocess.run(command, shell=True, check=True)

# Entry Point
if __name__ == "__main__":

    script = "/storage2/egusmao/projects/Bloom/bloom/data/download_geo.py"
    operation = "make" # metadata download make run merge delete delete_all_sra_files

    run_with_nohup(script, operation)

