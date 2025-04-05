import os
import sys
import subprocess
from pathlib import Path

# Define the base path of your project
current_path = Path().resolve()
project_root = current_path.parent.parent # Moves 2 levels up from Bloom/bloom/data/ to Bloom

# Add project root to sys.path
sys.path.append(str(project_root))

# Now I can import the module"
from bloom.data.geo_sra_downloader import GEODataDownloader

class RunPBS:
    
    def __init__(self, geo_id, ncbi_path, download_path, temp_path, email=None, api_key=None):

        self.geo_id = geo_id
        self.ncbi_path = Path(ncbi_path).resolve()
        self.download_path = Path(download_path).resolve()
        self.temp_path = Path(temp_path).resolve()
        self.email = email
        self.api_key = api_key
        self.geo_downloader = GEODataDownloader(self.geo_id,
                                                self.download_path,
                                                self.email,
                                                self.api_key,
                                                self.ncbi_path)
        self.sra_list = self.geo_downloader.get_prefetch_sra_id_list(check_fastq_exists=False,
                                                                     check_prefetch_exists=True)
        
    def create_metadata(self):

        self.geo_downloader.create_metadata_table()

    def download(self):

        self.geo_downloader.download_raw_data()

    def create_files(self):

        # Current location
        current_dir = str(Path(__file__).resolve().parent)

        # Create PBS scripts
        for i, sra_id in enumerate(self.sra_list):

            # Parameters
            sra_file_name = self._chech_sra_file_name(sra_id)
            output_location = Path(self.geo_downloader.output_dir) / sra_id

            # Unique PBS file name
            pbs_filename = f"fasterqdump_{self.geo_id}_{sra_id}_{i}.pbs"

            # Unique temp path
            temp_location = str(self.temp_path / f"fasterqdump_{self.geo_id}_{sra_id}")

            # Create PBS script
            with open(pbs_filename, "w") as pbs_file:
                pbs_file.write(
                    f"""#!/bin/bash

#PBS -N fasterqdump_{self.geo_id}_{sra_id}_{i}
#PBS -o fasterqdump_{self.geo_id}_{sra_id}_{i}.out
#PBS -e fasterqdump_{self.geo_id}_{sra_id}_{i}.err

#PBS -q workq
# workq - Fila default e sem restrições. Utiliza todos os nós.
# fatq - fila para os fat nodes.
# normq - fila para nodes comuns.
# gpuq - fila para processamento em GPU.
#PBS -V
#PBS -W umask=002

#PBS -l nodes=1:ppn=8
#PBS -l mem=64gb
#PBS -l walltime=48:00:00

# cd $PBS_O_WORKDIR

# Environments
eval \"$(micromamba shell hook --shell bash)\"
micromamba activate bio

# Current Job Parameter
basepath=\"{current_dir}\"
cd $basepath

# Python HEREDOC (EOF Block)

python <<EOF
import sys
from pathlib import Path

# Define the base path of your project
current_path = Path().resolve()
project_root = current_path.parent.parent # Moves 2 levels up from Bloom/bloom/data/ to Bloom

# Add the current PBS job working directory to sys.path
sys.path.insert(0, str(project_root.resolve()))

from bloom.data.geo_sra_downloader import GEODataDownloader

geo_downloader = GEODataDownloader(
    \"{self.geo_id}\",
    \"{self.download_path}\",
    \"{self.email}\",
    \"{self.api_key}\",
    \"{self.ncbi_path}\")

# Create paths if they do not exist already
Path(\"{temp_location}\").mkdir(parents=True, exist_ok=True)
Path(\"{output_location}\").mkdir(parents=True, exist_ok=True)

geo_downloader._process_sra_to_fastq(\"{sra_file_name}\",
                                     threads=\"8\",
                                     bufsize=\"8MB\",
                                     curcache=\"512MB\",
                                     mem=\"64GB\",
                                     include_technical=True,
                                     split_files=True,
                                     split_3=True,
                                     log_level=\"6\",
                                     verbose=True,
                                     temp_directory=\"{temp_location}\",
                                     output_directory=\"{output_location}\")

EOF
"""
                )

    def run_jobs(self):

        # Create PBS scripts
        for i, sra_id in enumerate(self.sra_list):

            # Unique PBS file name
            openpbs_submit = "qsub"
            pbs_filename = f"fasterqdump_{self.geo_id}_{sra_id}_{i}.pbs"

            # Submit job
            result = subprocess.run([openpbs_submit, pbs_filename], capture_output=True, text=True, check=True)
            print("Job submitted:", result.stdout.strip())

    def merge_files(self):

        # Get name of files
        out_file_list = []
        err_file_list = []
        # Create PBS scripts
        for i, sra_id in enumerate(self.sra_list):

            # Unique out and err PBS file name
            out_filename = f"fasterqdump_{self.geo_id}_{sra_id}_{i}.out"
            err_filename = f"fasterqdump_{self.geo_id}_{sra_id}_{i}.err"

            # Append to list only if it exists
            if os.path.exists(out_filename):
                out_file_list.append(out_filename)
            if os.path.exists(err_filename):
                err_file_list.append(err_filename)

        # Output file name
        out_file_merged = self.geo_downloader.output_dir / f"{self.geo_id}_merged.out"
        err_file_merged = self.geo_downloader.output_dir / f"{self.geo_id}_merged.err"

        # Merge only if list is not empty
        if len(out_file_list) >= 1:
            self._merge(out_file_list, out_file_merged)
        if len(err_file_list) >= 1:
            self._merge(err_file_list, err_file_merged)

    def delete_files(self):

        # Submit PBS scripts
        for i, sra_id in enumerate(self.sra_list):

            # Unique PBS, out and err file name
            pbs_filename = f"fasterqdump_{self.geo_id}_{sra_id}_{i}.pbs"
            out_filename = f"fasterqdump_{self.geo_id}_{sra_id}_{i}.out"
            err_filename = f"fasterqdump_{self.geo_id}_{sra_id}_{i}.err"

            # Removing files
            if os.path.exists(pbs_filename):
                os.remove(pbs_filename)
            if os.path.exists(out_filename):
                os.remove(out_filename)
            if os.path.exists(err_filename):
                os.remove(err_filename)

    def delete_SRA_files(self):

        # Create PBS scripts
        for i, sra_id in enumerate(self.sra_list):

            # Parameters
            sra_file_name = self._chech_sra_file_name(sra_id)
            sra_file_name.unlink(missing_ok=True)

    def _chech_sra_file_name(self, sra_id):

        # Check which SRA file exists
        sra_file = self.ncbi_path / "sra" / f"{sra_id}.sra"
        sra_file_lite = self.ncbi_path / "sra" / f"{sra_id}.sralite"
        if sra_file.exists():
            return sra_file
        if sra_file_lite.exists():
            return sra_file_lite
        return None

    def _merge(self, file_list, output_file):
        # Merge files
        with open(output_file, "w") as outfile:
            for file in file_list:
                outfile.write(f"> {file}\n")  # Write the filename as a header
                with open(file, "r") as infile:
                    outfile.write(infile.read())  # Append file content
                outfile.write(
                    ("-" * 50) + "\n\n"
                )  # Add 50 dashes and 2 blank lines at the end

# Entry Point
if __name__ == "__main__":

    # Parameters
    geo_id = "GSE285812"
    api_key = "c5087c87794c22daeb8f52d13fc5a363d108"
    user_email = "eduardogade@gmail.com"
    operation = sys.argv[1]

    # Dataset
    ncbi_path = "/storage2/egusmao/ncbi_sra/"
    download_path = "/storage2/egusmao/projects/Bloom/data/raw/"
    # download_path = "/Users/egg/Projects/Bloom/data/raw/"
    temp_path = "/storage2/egusmao/tmp/"

    # Input handling
    run_pbs = RunPBS(geo_id, ncbi_path, download_path, temp_path, user_email, api_key)

    # Operation
    if operation == "metadata":
        run_pbs.create_metadata()
    elif operation == "download":
        run_pbs.download()
    elif operation == "make":
        run_pbs.create_files()
    elif operation == "run":
        run_pbs.run_jobs()
    elif operation == "merge":
        run_pbs.merge_files()
    elif operation == "delete":
        run_pbs.delete_files()
    elif operation == "delete_all_sra_files":
        run_pbs.delete_SRA_files()

