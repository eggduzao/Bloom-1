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
    
    def __init__(self, geo_id, num_workers, download_path, email=None, api_key=None):

        self.geo_id = geo_id
        self.num_workers = num_workers
        self.download_path = download_path
        self.email = email
        self.api_key = api_key
        self.geo_downloader = GEODataDownloader(self.geo_id, self.download_path, self.email, self.api_key)

    def create_metadata(self):

        self.geo_downloader.create_metadata_table()

    def download(self):

        self.geo_downloader.download_raw_data(self.num_workers)

if __name__ == "__main__":

    # python -c "import os; os.cpu_count()"

    # Parameters
    geo_id = "GSE285812"
    api_key = "c5087c87794c22daeb8f52d13fc5a363d108"
    user_email = "eduardogade@gmail.com"
    num_workers = 10
    operation = sys.argv[1]

    # Dataset
    download_path = "/storage2/egusmao/projects/Bloom/data/raw/"
    # download_path = "/Users/egg/Projects/Bloom/data/raw/"

    # Input handling
    run_pbs = RunPBS(geo_id, num_workers, download_path, user_email, api_key)

    # Operation
    if operation == "make":
        run_pbs.create_metadata()
    elif operation == "run":
        run_pbs.download()
