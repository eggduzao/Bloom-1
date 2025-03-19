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
    
    def __init__(self, geo_id, ncbi_path, download_path, email=None, api_key=None):

        self.geo_id = geo_id
        self.num_workers = num_workers
        self.ncbi_path = Path(ncbi_path).resolve()
        self.download_path = Path(download_path).resolve()
        self.email = email
        self.api_key = api_key
        self.geo_downloader = GEODataDownloader(self.geo_id,
                                                self.download_path,
                                                self.email,
                                                self.api_key)
        self.sra_list = self.geo_downloader.get_sra_id()
        
    def create_metadata(self):

        self.geo_downloader.create_metadata_table()

    def download(self):

        self.geo_downloader.download_raw_data()

    def create_files(self):

        # Create PBS scripts
        for i, sra_id in enumerate(self.sra_list):

            # Parameters
            sra_file_name = self.ncbi_path / f"{sra_id}.sralite"
            output_location = self.geo_downloader.output_dir

            # Unique PBS file name
            pbs_filename = f"{self.geo_id}_{sra_id}_{i}.pbs"

            # Create PBS script
            with open(pbs_filename, "w") as pbs_file:
                pbs_file.write(
                    f"""#!/bin/bash

#PBS -N fastqdump_{self.geo_id}_{sra_id}_{i}
#PBS -o fastqdump_{self.geo_id}_{sra_id}_{i}.out
#PBS -e fastqdump_{self.geo_id}_{sra_id}_{i}.err

#PBS -q workq
# workq - Fila default e sem restrições. Utiliza todos os nós.
# fatq - fila para os fat nodes.
# normq - fila para nodes comuns.
# gpuq - fila para processamento em GPU.
#PBS -V
#PBS -W umask=002

#PBS -l nodes=1:ppn=4
#PBS -l mem=32gb
#PBS -l walltime=12:00:00

# cd $PBS_O_WORKDIR

# Environments
eval "$(micromamba shell hook --shell bash)"
micromamba activate bio

# Current Job Parameter
basepath=\"{output_location}\"
cd $basepath

# Uncompress Control
geo_downloader = GEODataDownloader({self.geo_id}, {self.download_path}, {self.email}, {self.api_key})
geo_downloader._process_sra_to_fastq({sra_file_name},
                                     log_level = \"5\",
                                     gzip_files = True,
                                     split_files = False,
                                     split_3: bool = True,
                                     output_directory={output_location})

"""
                )

    def run_jobs(self):

        # Create PBS scripts
        for i, sra_id in enumerate(self.sra_list):

            # Unique PBS file name
            pbs_filename = f"{self.geo_id}_{sra_id}_{i}.pbs"

            # Submit job
            result = subprocess.run(["qsub", pbs_filename], capture_output=True, text=True, check=True)
            print("Job submitted:", result.stdout.strip())

    def merge_files(self):

        # Get name of files
        out_file_list = []
        err_file_list = []
        # Create PBS scripts
        for i, sra_id in enumerate(self.sra_list):

            # Parameters
            sra_file_name = self.ncbi_path / f"{sra_id}.sralite"

            # Unique out and err PBS file name
            out_filename = f"{self.geo_id}_{sra_id}_{i}.out"
            err_filename = f"{self.geo_id}_{sra_id}_{i}.err"

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
        for i, param in enumerate(["Place", "holder"]):

            # Unique PBS, out and err file name
            pbs_filename = f"{self.geo_id}_{i}.pbs"
            out_filename = f"{self.geo_id}_{i}.out"
            err_filename = f"{self.geo_id}_{i}.err"

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
            sra_file_name = self.ncbi_path / f"{sra_id}.sralite"
            sra_file_name.unlink(missing_ok=True)

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

if __name__ == "__main__":

    # python -c "import os; os.cpu_count()"

    # Parameters
    geo_id = "GSE285812"
    api_key = "c5087c87794c22daeb8f52d13fc5a363d108"
    user_email = "eduardogade@gmail.com"
    operation = sys.argv[1]

    # Dataset
    ncbi_path = "/storage2/egusmao/ncbi_sra/sra"
    download_path = "/storage2/egusmao/projects/Bloom/data/raw/"
    # download_path = "/Users/egg/Projects/Bloom/data/raw/"

    # Input handling
    run_pbs = RunPBS(geo_id, ncbi_path, download_path, user_email, api_key)

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

