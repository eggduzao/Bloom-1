
"""
core
-------

Parses command-line arguments and performs core processing functions for dimension analysis.
"""

############################################################################################################
### Import
############################################################################################################

import os
import argparse

from . import __version__
from .utils import TableReader
from .metrics import (
    CompletenessMetric,
    ConformanceMetric,
    ConsistencyMetric,
    AgreementMetric,
    RelevanceMetric,
    RepresentativenessMetric,
    ContextualizationMetric,
)

############################################################################################################
### Constants
############################################################################################################

SEED = 1987

############################################################################################################
### Argument Parsing
############################################################################################################

def parse_args():
    """
    Parses command-line arguments.

    Returns
    -------
    argparse.Namespace
        An object containing the parsed command-line arguments.

    Positional Arguments
    --------------------
    input_training_folder : str
        Path to the input training folder.
    output_training_folder : str
        Path to the output training folder where results will be saved.
    root_name : str
        Root directory name used to structure the file paths.

    Optional Arguments
    ------------------
    --severity : float, optional
        Specifies the training severity. Default is 1.0.

    Standard Options
    ----------------
    -h, --help : Show help message and exit.
    -v, --version : Show the version of the Stainalyzer tool and exit.

    Notes
    -----
    - Ensure that input paths are valid directories.
    - The tool assumes images are formatted correctly within the input directory.
    
    """
    parser = argparse.ArgumentParser(
        description="Stainalyzer Tool: A robust tool for DAB-stained image analysis in microscopy.",
        epilog="Example Usage: Stainalyzer /path/to/input/root_folder/ /path/to/output root_folder --severity 0.5"
    )

    # Positional arguments
    parser.add_argument(
        "input_training_folder",
        type=str,
        help="Path to the input training folder."
    )
    parser.add_argument(
        "output_training_folder",
        type=str,
        help="Path to the output training folder where results will be saved."
    )
    parser.add_argument(
        "root_name",
        type=str,
        help="Root directory name used to structure the file paths."
    )

    # Optional arguments
    parser.add_argument(
        "--severity",
        type=float,
        default=0.5,
        help="Specifies the training severity from 0.0 (increase false positives) to 1.0 (decrease true positives) (default: 0.5)."
    )

    # Version and Help
    parser.add_argument(
        "-v", "--version",
        action="version",
        version=f'Stainalyzer {__version__}',
        help="Show program version and exit."
    )

    args = parser.parse_args()

    # Validate arguments
    if args.severity < 0.0 or args.severity > 1.0:
        parser.error("--severity must be between 0.0 and 1.0.")

    return args

############################################################################################################
### Classes
############################################################################################################


def core_function():
    """
    Core processing function for analyzing quality dimensions.

    Placeholder

    Parameters
    ----------
    Placeholder : str
        Placeholder.
    Placeholder : str
        Placeholder.
    Placeholder : str
        Placeholder.
    Placeholder : float, optional
        Placeholder.

    Raises
    ------
    Placeholder
        Placeholder.
    Placeholder
        Placeholder.
    
    Notes
    -----
    Placeholder

    """

    prefix = ("/Users/egg/Projects/Genoma_SUS_Code/genoma_sus_fenotipos/"
              "genoma_sus_fenotipos/projects/estudo_dimensoes_qualidade/")

    input_table_file_name = os.path.join(prefix, "data/raw/sabe_final_2010.xlsx")
    table_yaml_file_name = os.path.join(prefix, "rules/table.yaml")
    global_yaml_file_name = os.path.join(prefix, "rules/global.yaml ")
    output_location = os.path.join(prefix, "data/processed/")

    table_reader = TableReader(path_name=input_table_file_name, first_header=True, row_labels=False)
    table = table_reader.get_pandas_dataframe()

    table_yaml = YamlUtils(input_path=table_yaml_file_name)
    global_yaml = YamlUtils(input_path=global_yaml_file_name)

    missing_values = [None, np.nan, pd.NA, "NA", "N/A"]
    completeness_metric = CompletenessMetric(table, missing_values, global_yaml, table_yaml)

    missing_count, total_missing_count, completeness = calculate_completeness()
    row_result_vector, col_result_vector = calculate_bagging_completeness()

    print(missing_count)

    # calculate_uniqueness()
    # calculate_bagging_uniqueness()















if __name__ == "__main__":
    core_function()



"""
class PhenotypeTable:
    
    A container class for representing a phenotype dataset.

    Attributes:
    -----------
    header : PhenotypeHeader
        An instance of PhenotypeHeader representing the dataset's header.
    entries : list[PhenotypeEntry]
        A list of PhenotypeEntry instances representing individual data entries.
    

    def __init__(self, header, entries):
        "" "
        Initializes a PhenotypeTable with a header and entries.

        Parameters:
        -----------
        header : PhenotypeHeader
            The header for the phenotype dataset.
        entries : list[PhenotypeEntry]
            A list of entries in the phenotype dataset.
        "" "
        self.header = header  # Instance of PhenotypeHeader
        self.entries = entries  # List of PhenotypeEntry instances

    def get_header(self):
        "" "
        Gets the header of the phenotype table.

        Returns:
        --------
        PhenotypeHeader
            The header instance of the phenotype table.
        "" "
        return self.header

    def set_header(self, header):
        "" "
        Sets the header of the phenotype table.

        Parameters:
        -----------
        header : PhenotypeHeader
            The new header to set for the phenotype table.
        "" "
        self.header = header

    def add_entry(self, entry):
        "" "
        Adds an entry to the phenotype table.

        Parameters:
        -----------
        entry : PhenotypeEntry
            The entry to add to the table.
        "" "
        self.entries.append(entry)

    def remove_entry(self, index):
        "" "
        Removes an entry from the phenotype table by index.

        Parameters:
        -----------
        index : int
            The index of the entry to remove.

        Raises:
        -------
        IndexError
            If the index is out of range.
        "" "
        try:
            del self.entries[index]
        except IndexError as e:
            raise IndexError("Invalid index: {0}".format(index)) from e

    def get_entry(self, index):
        "" "
        Gets an entry from the phenotype table by index.

        Parameters:
        -----------
        index : int
            The index of the entry to retrieve.

        Returns:
        --------
        PhenotypeEntry
            The entry at the specified index.

        Raises:
        -------
        IndexError
            If the index is out of range.
        "" "
        try:
            return self.entries[index]
        except IndexError as e:
            raise IndexError("Invalid index: {0}".format(index)) from e

    def calculate_dimension_metric(self, metric_name):
        "" "
        Calculates a dimension metric for the phenotype table.

        Parameters:
        -----------
        metric_name : str
            The name of the dimension metric to calculate. Options include:
            "Completeness", "Conformance", "Consistency", "Agreement",
            "Relevance", "Representativeness", "Contextualization".

        Returns:
        --------
        float
            A placeholder value for the calculated metric.

        Notes:
        ------
        This method will eventually need to interact with other metric classes
        (e.g., CompletenessMetric) to perform the calculations.
        
        TODO:
        -----
        - Integrate this method with metric calculation classes.
        - Handle unsupported metric names.
        "" "
        # TODO - Validate metric_name and handle unsupported metrics
        # TODO - Interact with metric calculation classes once implemented
        return 0.0  # Placeholder return value

    def __len__(self):
        "" "
        Gets the number of entries in the phenotype table.

        Returns:
        --------
        int
            The number of entries in the table.
        "" "
        return len(self.entries)

    def __str__(self):
        "" "
        Returns a string representation of the phenotype table.

        Returns:
        --------
        str
            A summary string of the phenotype table.
        "" "
        return f"PhenotypeTable with {len(self.entries)} entries."

    def search_entries(self, **kwargs):
        "" "
        Searches for entries in the phenotype table matching criteria.

        Parameters:
        -----------
        kwargs : dict
            Key-value pairs of attributes to search for in entries.

        Returns:
        --------
        list[PhenotypeEntry]
            A list of entries matching the search criteria.

        TODO:
        -----
        - Implement filtering logic for search.
        - Define criteria matching rules with PhenotypeEntry attributes.
        "" "
        # TODO - Implement search logic based on PhenotypeEntry attributes
        return []  # Placeholder return value

    def to_dataframe(self):
        "" "
        Converts the phenotype table to a pandas DataFrame.

        Returns:
        --------
        pandas.DataFrame
            The phenotype table as a DataFrame.

        TODO:
        -----
        - Define conversion logic based on PhenotypeEntry and PhenotypeHeader.
        "" "
        # TODO - Implement conversion logic
        raise NotImplementedError("Conversion to DataFrame not implemented yet.")

class PhenotypeHeader:
    "" "
    Represents the header of a phenotype table.

    Attributes:
    -----------
    fields : list[FieldDescriptor]
        An ordered list of FieldDescriptor instances representing the header fields.
    "" "

    def __init__(self):
        "" "
        Initializes an empty PhenotypeHeader with an empty list of fields.
        "" "
        self.fields = []

    def add_field(self, field, field_id=None):
        "" "
        Adds a FieldDescriptor to the header.

        Parameters:
        -----------
        field : FieldDescriptor
            The FieldDescriptor to add.
        field_id : str, optional
            The ID to assign to the field. If not provided, defaults to an
            incremented string ID (e.g., "00000", "00001", etc.).

        Notes:
        ------
        Field IDs must be unique within the PhenotypeHeader.
        
        TODO:
        -----
        - Validate that field_id is unique if provided.
        "" "
        if field_id is None:
            field_id = f"{len(self.fields):05d}"
        # TODO - Ensure field_id uniqueness
        field.id = field_id
        self.fields.append(field)

    def get_field(self, field_id):
        "" "
        Retrieves a FieldDescriptor by its ID.

        Parameters:
        -----------
        field_id : str
            The ID of the field to retrieve.

        Returns:
        --------
        FieldDescriptor
            The field with the specified ID.

        Raises:
        -------
        ValueError
            If no field with the given ID exists.
        "" "
        for field in self.fields:
            if field.id == field_id:
                return field
        raise ValueError(f"Field with ID {field_id} not found.")

    def remove_field(self, field_id):
        "" "
        Removes a FieldDescriptor by its ID.

        Parameters:
        -----------
        field_id : str
            The ID of the field to remove.

        Raises:
        -------
        ValueError
            If no field with the given ID exists.
        "" "
        for index, field in enumerate(self.fields):
            if field.id == field_id:
                del self.fields[index]
                return
        raise ValueError(f"Field with ID {field_id} not found.")

    def list_fields(self):
        "" "
        Lists all FieldDescriptors in the header.

        Returns:
        --------
        list[FieldDescriptor]
            A list of all fields in the header.
        "" "
        return self.fields

    def update_field(self, field_id, new_field):
        "" "
        Updates a FieldDescriptor by its ID.

        Parameters:
        -----------
        field_id : str
            The ID of the field to update.
        new_field : FieldDescriptor
            The new FieldDescriptor to replace the existing one.

        Raises:
        -------
        ValueError
            If no field with the given ID exists.
        "" "
        for index, field in enumerate(self.fields):
            if field.id == field_id:
                self.fields[index] = new_field
                new_field.id = field_id
                return
        raise ValueError(f"Field with ID {field_id} not found.")

    def __len__(self):
        "" "
        Gets the number of fields in the header.

        Returns:
        --------
        int
            The number of fields in the header.
        "" "
        return len(self.fields)

    def __str__(self):
        "" "
        Returns a string representation of the PhenotypeHeader.

        Returns:
        --------
        str
            A summary string of the header.
        "" "
        return f"PhenotypeHeader with {len(self.fields)} fields."

    def to_dict(self):
        "" "
        Converts the header to a dictionary representation.

        Returns:
        --------
        dict
            A dictionary where keys are field IDs and values are field descriptors.

        TODO:
        -----
        - Define the FieldDescriptor serialization logic.
        "" "
        # TODO - Implement FieldDescriptor serialization to dictionary
        return {field.id: None for field in self.fields}  # Placeholder

    def from_dict(self, header_dict):
        "" "
        Populates the PhenotypeHeader from a dictionary representation.

        Parameters:
        -----------
        header_dict : dict
            A dictionary where keys are field IDs and values are field descriptors.

        TODO:
        -----
        - Implement deserialization logic for FieldDescriptor.
        "" "
        # TODO - Implement deserialization logic for FieldDescriptor
        raise NotImplementedError("Deserialization from dictionary not implemented yet.")
"""

