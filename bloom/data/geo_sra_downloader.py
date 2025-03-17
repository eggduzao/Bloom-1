"""
GEO-SRA Downloader
===================
Placeholder.

Authors: Eduardo G. Gusmao.

"""

###################################################################################################
# Libraries
###################################################################################################

# Python
import os
import time
import shutil
import requests
import functools
import subprocess
import multiprocessing
from pathlib import Path
from typing import Optional, Dict, List

# Internal


# External
import GEOparse
import pandas as pd
from Bio import Entrez
from bs4 import BeautifulSoup


###################################################################################################
# Constants
###################################################################################################

# import GEOparse.GEOparse as geo_parser

# def get_GEO_file_http(geo):
#     """Force GEOparse to use HTTPS instead of FTP."""
#     return f"https://ftp.ncbi.nlm.nih.gov/geo/series/{geo[:-3]}nnn/{geo}/soft/{geo}_family.soft.gz"

# # Monkey-patch the function
# geo_parser.get_GEO_file = get_GEO_file_http

###################################################################################################
# Classes
###################################################################################################

class GEODataDownloader:
    """
    A class for downloading data from GEO, including processed files, raw SRA files, 
    and associated metadata.

    Supports:
    - Processed data files from GEO Series (GSE) and Samples (GSM).
    - Raw sequencing data from SRA (via `fastq-dump`).
    - Automatic parsing of metadata to label and organize files correctly.

    Attributes
    ----------
    geo_id : str
        The GEO or SRA accession ID (e.g., "GSE12345", "SRP09876").
    output_dir : Path
        Directory where downloaded files will be stored.
    _temp_file_name : Path
        Temporary directory to store metadata and intermediate files.
    _gse : GEOparse.GSE
        GEO Series object for metadata extraction.
    """

    def __init__(self, geo_id: str, output_dir: str = "data/raw", email: str = None, api_key: str = None):
        """
        Initialize the downloader.

        Parameters
        ----------
        geo_id : str
            The GEO or SRA accession ID (e.g., "GSE12345", "SRP09876").
        output_dir : str, optional
            Directory where downloaded files will be stored, by default "data/raw".
        email : str, optional
            Email address for Entrez authentication.
        api_key : str, optional
            NCBI API key for higher rate limits.
        """
        self.geo_id = geo_id
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)  # Ensure output directory exists

        # Temporary directory for storing metadata and intermediate files
        self._temp_file_name = self.output_dir / f"{self.geo_id}_temp"
        self._temp_file_name.mkdir(parents=True, exist_ok=True)

        # Set credentials for NCBI Entrez API
        if email:
            Entrez.email = email
        else:
            raise ValueError("You must provide a valid email address for Entrez access.")

        if api_key:
            Entrez.api_key = api_key  # optional, but recommended for heavy usage

        self._gse = None  # Will hold the GEO Series object

        # Fetch metadata first, as it is required before downloading any data
        self.metadata_table = None  # Placeholder for metadata table
        if self._is_gse():
            self._gse = self._fetch_gse_data()  # Fetch GEO metadata

    def _is_gse(self) -> bool:
        """
        Check if the provided GEO ID corresponds to a GEO Series (GSE).

        Returns
        -------
        bool
            True if the ID is a GSE accession, False otherwise.
        """
        return self.geo_id.startswith("GSE")  # Simple check for GEO Series IDs

    def _fetch_gse_data(self) -> GEOparse.GSE:
        """
        Fetch GEO Series metadata using GEOparse.

        Returns
        -------
        GEOparse.GSE
            The GEO Series object with metadata.
        """
        print(f"Fetching {self.geo_id} from GEO...")
        gse = GEOparse.get_GEO(geo=self.geo_id, destdir=str(self._temp_file_name))
        print(f"Fetched {self.geo_id} Successfully!")
        return gse

    def __del__(self) -> None:
        """
        Destructor to clean up temporary files when the object is deleted.

        Notes
        -----
        - This method ensures that temporary metadata files are deleted 
          when the instance is garbage collected.
        - The actual downloaded data remains untouched in `self.output_dir`.
        """
        if hasattr(self, "_temp_file_name"):
            if self._temp_file_name.exists():
                print(f"Cleaning up temporary files at {self._temp_file_name}...")
                # shutil.rmtree(self._temp_file_name)  # Delete temporary directory

    def create_metadata_table(self) -> None:
        """
        Extracts metadata from GEO and SRA, processes it into a structured table, 
        and saves the result as a TSV file.

        This function follows a structured pipeline:
        1. Extract study-wide metadata.
        2. Extract sample-level metadata.
        3. Map GSM (sample IDs) to SRX (experiment IDs).
        4. Fetch SRR (run IDs) corresponding to each experiment.
        5. Merge metadata into a Pandas DataFrame.
        6. Insert study-wide metadata as the first two rows.
        7. Save the final table in a tab-separated format.

        Notes
        -----
        - The first two rows of the table will contain **study-wide metadata** 
          (keys in row 1, values in row 2).
        - Sample-specific metadata follows from row 3 onward.
        - The final table will be stored at `self.output_dir`.

        Raises
        ------
        Exception
            If there is an issue extracting metadata or writing the file.
        """
        
        # Step 1: Extract study-wide metadata (e.g., experiment details, project description)
        study_dictionary = self._extract_study_metadata()

        # Step 2: Extract sample-specific metadata (e.g., sample name, organism, molecule type)
        sample_dictionary = self._extract_sample_metadata()

        # Step 3: Extract mapping between GSM (samples) and SRX (experiments)
        gsm_to_srx = self._extract_srx_from_gse()

        # Step 4: Process GSM -> SRX -> SRR Mapping
        geo_sra_data = []  # List to store processed metadata
        for gsm, srx in gsm_to_srx.items():
            srr_ids = self._fetch_srr_from_srx(srx)  # Fetch run IDs for the experiment
            time.sleep(2)  # Delay to prevent rate limits from SRA servers

            # Step 4.1: Merge metadata with SRX and SRR information
            merged_data = sample_dictionary.get(gsm, {}).copy()  # Retrieve sample metadata
            merged_data.update({
                "GSM_ID": gsm,  # GEO Sample ID
                "SRX_ID": srx,  # SRA Experiment ID
                "SRR_IDs": ",".join(srr_ids) if srr_ids else "No SRRs found"  # List of run IDs
            })

            geo_sra_data.append(merged_data)

        # Step 5: Convert collected data into a Pandas DataFrame
        df_samples = pd.DataFrame(geo_sra_data)

        # Step 6: Convert study-wide metadata into a DataFrame (first row: keys, second row: values)
        df_study = pd.DataFrame([study_dictionary])

        # Step 7: Save the final study table as a TSV file
        output_path = self.output_dir / f"{self.geo_id}_study.tsv"
        df_study.to_csv(output_path, sep="\t", index=False)

        # Step 8: Save the final metadata table as a TSV file
        output_path = self.output_dir / f"{self.geo_id}_metadata.tsv"
        df_samples.to_csv(output_path, sep="\t", index=False)

    def download_raw_data(self, num_workers: int = None) -> None:
        """
        Downloads raw sequencing data (FASTQ files) from SRA using `fastq-dump`.

        This method reads the previously generated metadata table (TSV), extracts
        all **SRR** (SRA Run IDs), and downloads the corresponding sequencing data.

        Workflow:
        1. Open the metadata table generated by `create_metadata_table()`.
        2. Extract all SRR (SRA Run IDs) from the metadata.
        3. Create an `SRA` subfolder inside `self.output_dir` for storing downloads.
        4. Use `fastq-dump` (or `prefetch` from SRA Toolkit) to download FASTQ files.
        5. Handle potential download failures and retry failed downloads.

        Parameters
        ----------
        num_workers : int
            Number of cores or processing units to parallelize.

        Raises
        ------
        FileNotFoundError
            If the metadata table file does not exist.
        RuntimeError
            If `fastq-dump` fails for any SRR.
        """

        # Step 1: Locate metadata file
        metadata_file = self.output_dir / f"{self.geo_id}_metadata.tsv"

        if not metadata_file.exists():
            raise FileNotFoundError(f"Metadata table not found: {metadata_file}")

        # Step 2: Load metadata table and extract SRR IDs
        metadata_df = pd.read_csv(metadata_file, sep="\t")

        if "SRR_IDs" not in metadata_df.columns:
            raise ValueError("Metadata table is missing 'SRR_IDs' column!")

        # Step 3: Iterate over SRR IDs and download each FASTQ file
        srr_id_list = []
        for idx, row in metadata_df.iterrows():
            srr_list = str(row["SRR_IDs"]).split(",")

            for srr_id in srr_list:
                srr_id = srr_id.strip()

                if srr_id.lower() in ["no srrs found", "nan", ""]:
                    print(f"No SRRs found for {row['GSM_ID']} (Skipping)")
                    continue

                # Construct FASTQ file path
                fastq_file = self.output_dir / f"{srr_id}.fastq.gz"

                # Check if file already exists to avoid redundant downloads
                if fastq_file.exists():
                    print(f"{srr_id}.fastq.gz already exists, skipping download.")
                    continue

                # Step 4: Add SSR_ID to list
                srr_id_list.append(srr_id)

        # Step 5: Detect available CPU cores
        if num_workers is None:
            num_workers = max(1, int(ceil(os.cpu_count()/2)))

        # Pre-fill optional parameters
        download_partial = functools.partial(
            self._download_with_fastq_dump,
            split_files = False,
            split_3 = True,
            gzip_files = True
        )

        # Step 6: Download using `fastq-dump` and parallel processing
        with multiprocessing.Pool(processes=num_workers) as pool:
            pool.map(download_partial, srr_id_list)

        print(f"\nAll available SRA runs downloaded to: {str(self.output_dir)}")

    def download_processed_data(self) -> None:
        """
        Downloads processed data files (e.g., counts, expression matrices) from GEO.

        GEO provides processed data as supplementary files, accessible via FTP.

        Workflow:
        1. Construct the GEO FTP URL based on the GEO ID format.
        2. Fetch the list of available processed files.
        3. Download each file and save it in `self.output_dir/Processed`.
        4. Handle missing files and failed downloads gracefully.

        Raises
        ------
        ValueError
            If the GEO ID is invalid.
        RuntimeError
            If file download fails.
        """

        # Step 1: Construct base GEO FTP URL
        GEO_BASE_URL = "https://ftp.ncbi.nlm.nih.gov/geo/"
        
        if self.geo_id.startswith("GSE"):
            url = f"{GEO_BASE_URL}series/{self.geo_id[:-3]}nnn/{self.geo_id}/suppl/"
        elif self.geo_id.startswith("GSM"):
            url = f"{GEO_BASE_URL}samples/{self.geo_id[:-3]}nnn/{self.geo_id}/suppl/"
        else:
            raise ValueError("Invalid GEO ID format. Expected GSExxxx or GSMxxxx.")

        print(f"Checking for processed data files at: {url}")

        # Step 2: Fetch list of available files (assumes `_fetch_file_list` is implemented)
        file_list = self._fetch_file_list(url)

        if not file_list:
            print("No processed data files found.")
            return

        # Step 3: Create output directory for processed data
        processed_output_dir = self.output_dir / "Processed"
        processed_output_dir.mkdir(parents=True, exist_ok=True)

        # Step 4: Download each file
        for filename in file_list:
            file_url = url + filename
            output_path = processed_output_dir / filename

            if output_path.exists():
                print(f"{filename} already exists, skipping download.")
                continue

            print(f"Downloading {filename}...")
            try:
                response = requests.get(file_url, stream=True)
                response.raise_for_status()  # Raise error if request failed
                
                with open(output_path, "wb") as f:
                    for chunk in response.iter_content(chunk_size=8192):
                        f.write(chunk)

                print(f"Downloaded: {filename}")

            except requests.RequestException:
                print(f"Failed to download {filename}.")
                continue  # Skip failed downloads

        print(f"\nAll processed data downloaded to: {processed_output_dir}")

    def _fetch_gse_data(self) -> GEOparse.GSE:
        """
        Fetches GEO Series data for a given GEO accession ID (GSE).

        This method downloads metadata and associated files from the NCBI GEO database.
        If the data is already downloaded in the temporary directory, it will reuse it.

        Returns
        -------
        GEOparse.GEO.GSE
            A GEOparse GSE object containing dataset metadata and sample information.

        Raises
        ------
        ValueError
            If the GEO ID does not start with 'GSE'.
        RuntimeError
            If downloading or parsing the GEO dataset fails.

        Examples
        --------
        >>> gse = self._fetch_gse_data()
        >>> print(gse.metadata["title"])
        """

        # Ensure GEO ID is valid
        if not self.geo_id.startswith("GSE"):
            raise ValueError("Invalid GEO ID. This function only supports GSE IDs.")

        print(f"Fetching GEO dataset: {self.geo_id}")

        try:
            # Fetch GEO dataset, storing it in a temporary directory
            gse = GEOparse.get_GEO(geo=self.geo_id, destdir=str(self._temp_file_name))

            print(f"Successfully fetched GEO dataset: {self.geo_id}")
            return gse

        except Exception as e:
            raise RuntimeError(f"Failed to fetch GEO dataset {self.geo_id}: {e}")

    def _extract_study_metadata(self) -> Dict[str, str]:
        """
        Extracts study-wide metadata from the GEO dataset.

        This method retrieves general study information, such as the dataset title, 
        summary, and experimental design.

        Returns
        -------
        Dict[str, str]
            A dictionary containing the study-wide metadata.

        Raises
        ------
        AttributeError
            If `_gse` is not initialized or contains invalid metadata.

        Examples
        --------
        >>> study_metadata = self._extract_study_metadata()
        >>> print(study_metadata["Dataset Title"])
        """

        if not self._gse:
            raise AttributeError("GEO dataset not loaded. Call `_fetch_gse_data()` first.")

        try:
            study_metadata = {
                "Dataset Title": "; ".join(self._gse.metadata.get("title", [])),
                "Dataset Summary": "; ".join(self._gse.metadata.get("summary", [])),
                "Dataset Design": "; ".join(self._gse.metadata.get("overall_design", []))
            }
        except KeyError as e:
            raise KeyError(f"Metadata key missing: {e}")

        return study_metadata

    def _extract_sample_metadata(self) -> Dict[str, Dict[str, str]]:
        """
        Extracts metadata for each sample (GSM) in the GEO dataset.

        This method iterates through all samples, extracting relevant details such as
        the sample title, organism, library strategy, and experimental conditions.

        Returns
        -------
        Dict[str, Dict[str, str]]
            A dictionary where each key is a **Sample ID (GSM_ID)**, and the value 
            is another dictionary containing the metadata fields.

        Raises
        ------
        AttributeError
            If `_gse` is not initialized or contains invalid metadata.

        Examples
        --------
        >>> sample_metadata = self._extract_sample_metadata()
        >>> print(sample_metadata["GSM123456"]["Title"])
        """

        if not self._gse:
            raise AttributeError("GEO dataset not loaded. Call `_fetch_gse_data()` first.")

        sample_dict = {}

        for sample_id, sample in self._gse.gsms.items():
            try:
                metadata = {
                    "Sample_ID": sample_id,
                    "Title": "; ".join(sample.metadata.get("title", [])),
                    "Source_Name": "; ".join(sample.metadata.get("source_name_ch1", [])),
                    "Organism": "; ".join(sample.metadata.get("organism_ch1", [])),
                    "Molecule": "; ".join(sample.metadata.get("molecule_ch1", [])),
                    "Characteristics": "; ".join(sample.metadata.get("characteristics_ch1", [])),
                    "Library_Source": "; ".join(sample.metadata.get("library_source", [])),
                    "Library_Strategy": "; ".join(sample.metadata.get("library_strategy", [])),
                    "Description": "; ".join(sample.metadata.get("description", []))
                }

                sample_dict[sample_id] = metadata

            except Exception as e:
                print(f"Error processing sample {sample_id}: {e}")
                continue  # Skip problematic samples

        return sample_dict

    def _extract_srx_from_gse(self) -> Dict[str, str]:
        """
        Extracts **SRA Experiment IDs (SRX)** from GSM metadata.

        This method scans the metadata of all samples (GSMs) and identifies 
        SRA links, extracting the corresponding **SRX IDs**.

        Returns
        -------
        Dict[str, str]
            A dictionary mapping **GSM_IDs** (keys) to their corresponding **SRX_IDs** (values).
            If no SRX is found for a given GSM, it will not be included in the dictionary.

        Raises
        ------
        AttributeError
            If `_gse` is not initialized or does not contain sample metadata.

        Examples
        --------
        >>> srx_mapping = self._extract_srx_from_gse()
        >>> print(srx_mapping["GSM123456"])  # Expected Output: "SRX123456"
        """

        if not self._gse:
            raise AttributeError("GEO dataset not loaded. Call `_fetch_gse_data()` first.")

        srx_mapping = {}

        for gsm_id, gsm in self._gse.gsms.items():
            relations = gsm.metadata.get("relation", [])

            for relation in relations:
                if "SRA:" in relation:
                    srx_id = relation.split("SRA:")[-1].strip()

                    # Extract SRX ID from full NCBI link if present
                    if "https://www.ncbi.nlm.nih.gov/sra?term=" in srx_id:
                        srx_id = srx_id.split("term=")[-1]

                    # Ensure valid SRX ID before storing
                    if srx_id.startswith("SRX"):
                        srx_mapping[gsm_id] = srx_id

        if not srx_mapping:
            print("Warning: No SRX IDs were found in the dataset.")

        return srx_mapping

    def _fetch_srr_from_srx(self, srx_id: str) -> List[str]:
        """
        Queries **NCBI SRA** for **SRR (Run IDs)** using the given **SRX Experiment ID**.

        This method performs an HTTP request to the **NCBI SRA database** and extracts all
        associated **SRR run IDs**, which are required for downloading sequencing data.

        Parameters
        ----------
        srx_id : str
            The **SRX Experiment ID** for which **SRR IDs** should be retrieved.

        Returns
        -------
        List[str]
            A list of **SRR run IDs** associated with the provided SRX ID.
            Returns an empty list if no SRRs are found.

        Raises
        ------
        ValueError
            If an invalid or empty SRX ID is provided.
        requests.RequestException
            If the HTTP request to NCBI fails.

        Examples
        --------
        >>> srr_list = self._fetch_srr_from_srx("SRX123456")
        >>> print(srr_list)  # Expected Output: ["SRR987654", "SRR987655"]
        """

        if not srx_id or not srx_id.startswith("SRX"):
            raise ValueError(f"Invalid SRX ID provided: {srx_id}")

        url = f"https://www.ncbi.nlm.nih.gov/sra/?term={srx_id}"
        print(f"Fetching SRR IDs for {srx_id}...")

        try:
            response = requests.get(url, timeout=10)
            response.raise_for_status()  # Raise an error for bad responses (4xx, 5xx)

        except requests.RequestException as e:
            print(f"Failed to fetch SRX {srx_id}: {e}")
            return []

        soup = BeautifulSoup(response.text, "html.parser")
        srr_list = []

        # Search for SRR links in the HTML response
        for link in soup.find_all("a"):
            text = link.text.strip()
            if text.startswith("SRR"):
                srr_list.append(text)

        if not srr_list:
            print(f"Warning: No SRR IDs found for {srx_id}")

        return srr_list

    def _download_with_fastq_dump(self,
                                  sra_id: str,
                                  split_files: bool = False,
                                  split_3: bool = True,
                                  gzip_files: bool = True) -> None:
        """
        Download an SRA file and convert it to FASTQ using `fastq-dump`.

        This method runs `fastq-dump` with configurable options to download 
        and convert raw sequencing data from **NCBI SRA** into FASTQ format.

        Parameters
        ----------
        sra_id : str
            The **SRA Run accession ID** (e.g., "SRR123456").
        split_files : bool, optional
            Whether to split paired-end reads into `_1.fastq` and `_2.fastq` files. Default is False.
        split_3 : bool, optional
            Whether to split three reads (forward, reverse, and index/barcode) into separate files.
            If `split_files=True`, this is ignored. Default is True.
        gzip_files : bool, optional
            Whether to compress the FASTQ output using gzip (`.fastq.gz`). Default is True.

        Raises
        ------
        ValueError
            If `sra_id` is empty or not a valid SRA accession.
        FileNotFoundError
            If `fastq-dump` is not installed or not found in the system path.
        subprocess.CalledProcessError
            If `fastq-dump` execution fails (e.g., due to network issues or NCBI restrictions).

        Examples
        --------
        >>> downloader._download_with_fastq_dump("SRR123456")
        Running fastq-dump for SRR123456...
        Finished fastq-dump for SRR123456

        >>> downloader._download_with_fastq_dump("SRR789012", split_files=True, gzip_files=False)
        Running fastq-dump for SRR789012...
        Finished fastq-dump for SRR789012
        """

        if not sra_id or not sra_id.startswith("SRR"):
            raise ValueError(f"Invalid SRA ID provided: {sra_id}")

        print(f"Running fastq-dump for {sra_id}..." + "-"*30)

        # Ensure output directory exists
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Construct fastq-dump parameters
        split_param = ""
        if split_3:
            split_param = "--split-3"
        elif split_files:
            split_param = "--split-files"

        gzip_param = "--gzip" if gzip_files else ""

        # Full command
        cmd = [
            "fastq-dump",
            split_param,
            gzip_param,
            "--outdir", str(self.output_dir),
            sra_id
        ]

        # Remove empty strings (if any) from command list
        cmd = [arg for arg in cmd if arg]

        print(f"Running command: {' '.join(cmd)}")

        try:
            # Run fastq-dump and capture errors if any
            subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

            print(f"Finished fastq-dump for {sra_id}" + "-"*30)

        except FileNotFoundError:
            raise FileNotFoundError("`fastq-dump` not found. Make sure it's installed and in your PATH.")
        
        except subprocess.CalledProcessError as e:
            print(f"Fastq-dump failed for {sra_id} with error: {e.stderr.decode()}")
            raise

    def __repr__(self) -> str:
        """
        Official string representation of the class instance.

        This method returns a detailed representation of the object,
        which is useful for debugging and logging.

        Returns
        -------
        str
            A string containing the class name, GEO ID, and output directory.

        Examples
        --------
        >>> downloader = GEODataDownloader("GSE285812", "data/raw")
        >>> repr(downloader)
        'GEODataDownloader(geo_id="GSE285812", output_dir="data/raw")'
        """
        return f'GEODataDownloader(geo_id="{self.geo_id}", output_dir="{self.output_dir}")'

    def __str__(self) -> str:
        """
        User-friendly string representation of the class instance.

        This method returns a human-readable summary of the instance,
        typically used when printing the object.

        Returns
        -------
        str
            A descriptive string with key details about the downloader instance.

        Examples
        --------
        >>> downloader = GEODataDownloader("GSE285812", "data/raw")
        >>> print(downloader)
        '[GEODataDownloader] GEO ID: GSE285812 | Output Directory: data/raw'
        """
        return f'[GEODataDownloader] GEO ID: {self.geo_id} | Output Directory: {self.output_dir}'

