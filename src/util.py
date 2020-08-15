"""
Util
===================
The Util classes contains many utilities needed by other classes such as the paths to input files.
"""

###################################################################################################
# Libraries
###################################################################################################

# Python
import os
import sys
import codecs
import traceback
from optparse import OptionParser, BadOptionError, AmbiguousOptionError
from configparser import ConfigParser

# Internal

# External


###################################################################################################
# Data Path Handling
###################################################################################################

# TODO ALLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL

def get_rgtdata_path():
  return os.path.expanduser(os.getenv("RGTDATA", os.path.join(os.getenv("HOME"), "rgtdata")))

class ConfigurationFile:
  """
  Represent the data path configuration file (data.config). It serves as a superclass to classes that will contain
  default variables (such as paths, parameters to tools, etc.) for a certain purpose (genomic data, motif data, etc.).
  *Variables:*
    - self.config -- Represents the configuration file.
    - self.data_dir -- Represents the root path to data files.
  """

  def __init__(self):
    # Reading config file directory
    data_config_file_name = os.path.join(get_rgtdata_path(), "data.config")

    # Parsing config file
    self.config = ConfigParser()
    # self.config.read(data_config_file_name)
    self.config.read_file(codecs.open(data_config_file_name, "r", "utf8"))

    # Overwriting config using user options
    # self.config.read(data_config_file_name + ".user")
    self.config.read_file(codecs.open(data_config_file_name + ".user", "r", "utf8"))

    # Reading data directory
    self.data_dir = os.path.split(data_config_file_name)[0]

class GenomeData(ConfigurationFile):
  """Represent genomic data. Inherits ConfigurationFile."""

  def __init__(self, organism):
    """Initializes GenomeData.
    *Keyword arguments:*
      - organism -- Organism alias.
    """
    ConfigurationFile.__init__(self)
    self.organism = organism
    self.genome = os.path.join(self.data_dir, self.config.get(organism, 'genome'))
    self.chromosome_sizes = os.path.join(self.data_dir, self.config.get(organism, 'chromosome_sizes'))
    self.genes_gencode = os.path.join(self.data_dir, self.config.get(organism, 'genes_Gencode'))
    self.genes_refseq = os.path.join(self.data_dir, self.config.get(organism, 'genes_RefSeq'))
    self.annotation = os.path.join(self.data_dir, self.config.get(organism, 'annotation'))
    self.annotation_dump_dir = os.path.dirname(os.path.join(self.data_dir, self.annotation))
    self.gene_alias = os.path.join(self.data_dir, self.config.get(organism, 'gene_alias'))
    if organism in ["hg19", "hg38", "mm9"]:
      self.repeat_maskers = os.path.join(self.data_dir, self.config.get(organism, 'repeat_maskers'))
    else:
      self.repeat_maskers = None

  def get_organism(self):
    """Returns the current organism."""
    return self.organism

  def get_genome(self):
    """Returns the current path to the genome fasta file."""
    return self.genome

  def get_chromosome_sizes(self):
    """Returns the current path to the chromosome sizes text file."""
    return self.chromosome_sizes

  def get_gene_regions(self):
    """Returns the current path to the gene_regions BED file."""
    return self.genes_gencode

  def get_genes_gencode(self):
    """Returns the current path to the gene_regions BED file."""
    return self.genes_gencode

  def get_genes_refseq(self):
    """Returns the current path to the gene_regions BED file."""
    return self.genes_refseq

  def get_annotation(self):
    """
    Returns the current path to the gencode annotation gtf file.
    """
    return self.annotation

  def get_annotation_dump_dir(self):
    """Returns the current path to the gencode annotation gtf file."""
    return self.annotation_dump_dir

  def get_gene_alias(self):
    """Returns the current path to the gene alias txt file."""
    return self.gene_alias

  def get_repeat_maskers(self):
    """Returns the current path to directory for repeat maskers."""
    if self.repeat_maskers:
      return self.repeat_maskers
    else:
      print("*** There is no repeat masker data for " + self.organism)

class HmmData(ConfigurationFile):
  """Represent HMM data. Inherits ConfigurationFile."""

  def __init__(self):
    """Error handling.
    
    *Keyword arguments:*
    
      - msg -- String containing the error message.
    
    *Return:*
    
      - return -- An error message to the user.
    """

    ConfigurationFile.__init__(self)
    self.default_hmm_dnase = os.path.join(self.data_dir, self.config.get('HmmData', 'default_hmm_dnase'))
    self.default_hmm_dnase_bc = os.path.join(self.data_dir, self.config.get('HmmData', 'default_hmm_dnase_bc'))
    self.default_hmm_atac_paired = os.path.join(self.data_dir,
                                                    self.config.get('HmmData', 'default_hmm_atac_paired'))
    self.default_hmm_atac_single = os.path.join(self.data_dir,
                                                   self.config.get('HmmData', 'default_hmm_atac_single'))
    self.default_hmm_histone = os.path.join(self.data_dir, self.config.get('HmmData', 'default_hmm_histone'))
    self.default_hmm_dnase_histone = os.path.join(self.data_dir,
                                                      self.config.get('HmmData', 'default_hmm_dnase_histone'))
    self.default_hmm_dnase_histone_bc = os.path.join(self.data_dir,
                                                         self.config.get('HmmData', 'default_hmm_dnase_histone_bc'))
    self.default_hmm_atac_histone = os.path.join(self.data_dir,
                                                     self.config.get('HmmData', 'default_hmm_atac_histone'))
    self.default_hmm_atac_histone_bc = os.path.join(self.data_dir,
                                                        self.config.get('HmmData', 'default_hmm_atac_histone_bc'))
    self.default_bias_table_F_SH = os.path.join(self.data_dir,
                                                    self.config.get('HmmData', 'default_bias_table_F_SH'))
    self.default_bias_table_R_SH = os.path.join(self.data_dir,
                                                    self.config.get('HmmData', 'default_bias_table_R_SH'))
    self.default_bias_table_F_DH = os.path.join(self.data_dir,
                                                    self.config.get('HmmData', 'default_bias_table_F_DH'))
    self.default_bias_table_R_DH = os.path.join(self.data_dir,
                                                    self.config.get('HmmData', 'default_bias_table_R_DH'))
    self.default_bias_table_F_ATAC = os.path.join(self.data_dir,
                                                      self.config.get('HmmData', 'default_bias_table_F_ATAC'))
    self.default_bias_table_R_ATAC = os.path.join(self.data_dir,
                                                      self.config.get('HmmData', 'default_bias_table_R_ATAC'))
    self.dependency_model = os.path.join(self.data_dir, self.config.get('HmmData', 'dependency_model'))
    self.slim_dimont_predictor = os.path.join(self.data_dir, self.config.get('HmmData', 'slim_dimont_predictor'))
    self.default_test_fa = os.path.join(self.data_dir, self.config.get('HmmData', 'default_test_fa'))

  def get_default_hmm_dnase(self):
    """Returns the current default DNase only hmm."""
    return self.default_hmm_dnase

  def get_default_hmm_dnase_bc(self):
    """Returns the current default DNase only hmm."""
    return self.default_hmm_dnase_bc

  def get_default_hmm_atac_paired(self):
    """Returns the current default atac only hmm."""
    return self.default_hmm_atac_paired

  def get_default_hmm_atac_single(self):
    """Returns the current default atac only hmm."""
    return self.default_hmm_atac_single

  def get_default_hmm_histone(self):
    """Returns the current default Histone only hmm."""
    return self.default_hmm_histone

  def get_default_hmm_dnase_histone(self):
    """Returns the current default DNase+histone hmm."""
    return self.default_hmm_dnase_histone

  def get_default_hmm_dnase_histone_bc(self):
    """Returns the current default DNase+histone hmm."""
    return self.default_hmm_dnase_histone_bc

  def get_default_hmm_atac_histone(self):
    """Returns the current default atac+histone hmm."""
    return self.default_hmm_atac_histone

  def get_default_hmm_atac_histone_bc(self):
    """Returns the current default atac+histone hmm."""
    return self.default_hmm_atac_histone_bc

  def get_default_bias_table_F_SH(self):
    """Returns the DNase-seq single hit default bias table for the forward strand."""
    return self.default_bias_table_F_SH

  def get_default_bias_table_R_SH(self):
    """Returns the DNase-seq single hit default bias table for the reverse strand."""
    return self.default_bias_table_R_SH

  def get_default_bias_table_F_DH(self):
    """Returns the DNase-seq double hit default bias table for the forward strand."""
    return self.default_bias_table_F_DH

  def get_default_bias_table_R_DH(self):
    """Returns the DNase-seq double hit default bias table for the reverse strand."""
    return self.default_bias_table_R_DH

  def get_default_bias_table_F_ATAC(self):
    """Returns the ATAC-seq default bias table for the forward strand."""
    return self.default_bias_table_F_ATAC

  def get_default_bias_table_R_ATAC(self):
    """Returns the ATAC-seq default bias table for the reverse strand."""
    return self.default_bias_table_R_ATAC

  def get_dependency_model(self):
    return self.dependency_model

  def get_slim_dimont_predictor(self):
    return self.slim_dimont_predictor

  def get_default_test_fa(self):
    return self.default_test_fa


###################################################################################################
# Argument Parsing
###################################################################################################

# From here.

class HelpfulOptionParser(OptionParser):
  """This class represents an OptionParser that prints full help on errors. Inherits OptionParser.

  *Keyword arguments:*

    - OptionParser -- Inherited OptionParser object.
  """

  def error(self, msg):
    """Error handling.
    
    *Keyword arguments:*
    
      - msg -- String containing the error message.
    
    *Return:*
    
      - return -- An error message to the user.
    """
    self.print_help(sys.stderr)
    self.exit(2, "\n%s: error: %s\n" % (self.get_prog_name(), msg))


class PassThroughOptionParser(HelpfulOptionParser):
  """When unknown arguments are encountered, bundle with largs and try again, until rargs is depleted.
     sys.exit(status) will still be called if a known argument is passed incorrectly (e.g. missing arguments or
     bad argument types, etc.). Inherits HelpfulOptionParser.

  *Keyword arguments:*

    - HelpfulOptionParser -- Inherited HelpfulOptionParser object.
  """

  def _process_args(self, largs, rargs, values):
    """Overwrites _process_args to achieve desired effect.
    
    *Keyword arguments:*
    
      - largs -- The tool's optional arguments.
      - rargs -- The tool's required arguments.
    """
    while rargs:
      try:
        HelpfulOptionParser._process_args(self, largs, rargs, values)
      except (BadOptionError, AmbiguousOptionError):
        pass
        # largs.append(err.opt_str)


###################################################################################################
# Warnings and Errors Handling
###################################################################################################

class ErrorHandler:
  """Handles errors and warnings in a standardized way.

  *Error Dictionary Standard:*

    Each entry consists of a key+list in the form X:[Y,Z,W] where:
      - X -- The key representing the internal error name.
      - Y -- Error number.
      - Z -- Exit status.
      - W -- Error message to be print.

  *Warning Dictionary Standard:*

    Each entry consists of a key+list in the form X:[Y,Z] where:
      - X -- The key representing the internal warning name.
      - Y -- Warning number.
      - Z -- Warning message to be print.
    """

  def __init__(self):
    """Initializes required objects for warning and error handling.
    """

    self.program_name = os.path.basename(sys.argv[0])

    self.error_dictionary = {
      "DEFAULT_ERROR": [0, 0, "Undefined error. Program terminated with exit status 0."],
      "MOTIF_ANALYSIS_OPTION_ERROR": [1, 0,
                                      "You must define one specific analysis. Run '" + self.program_name + " -h' for help."],

      "FP_WRONG_ARGUMENT": [2, 0,
                            "You must provide at least one and no more than one experimental matrix as input argument."],
      "XXXXXXX": [26, 0, "Xxxxxx"]
    }
    self.error_number = 0
    self.exit_status = 1
    self.error_message = 2

    self.warning_dictionary = {
      "DEFAULT_WARNING": [0, "Undefined warning."],
      "FP_ONE_REGION": [1,
                        "There are more than one 'regions' file in the experiment matrix. Only the first will be used."],
      "FP_MANY_DNASE": [2,
                        "There are more than one DNASE or ATAC 'reads' file. Only the first one will be used."],
      "XXXXXXX": [12, "Xxxxxx"]
    }
    self.warning_number = 0
    self.warning_message = 1

  def throw_error(self, error_type, add_msg=""):
    """Throws the specified error type. If the error type does not exist, throws a default error message and exits.
    
    *Keyword arguments:*
    
      - error_type -- Error type.
      - add_msg -- Message to add to the error.
    """

    # Fetching error type
    try:
      error_number = self.error_dictionary[error_type][self.error_number]
      exit_status = self.error_dictionary[error_type][self.exit_status]
      error_message = self.error_dictionary[error_type][self.error_message]
    except (KeyError, IndexError):
      error_number = self.error_dictionary["DEFAULT_ERROR"][self.error_number]
      exit_status = self.error_dictionary["DEFAULT_ERROR"][self.exit_status]
      error_message = self.error_dictionary["DEFAULT_ERROR"][self.error_message]

    # Handling error
    complete_error_message = ("--------------------------------------------------\n"
                              "Error Number: " + str(error_number) + ".\n"
                              "Program: " + self.program_name + ".\n"
                              "Report: " + error_message + " " + add_msg + "\n"
                              "Behaviour: The program will quit with exit status " + str(exit_status) + ".\n"
                              "--------------------------------------------------")
    print(complete_error_message, file=sys.stderr)
    traceback.print_exc()
    sys.exit(exit_status)

  def throw_warning(self, warning_type, add_msg=""):
    """Throws the specified warning type. If the warning type does not exist, throws a default warning message and exits.
    
    *Keyword arguments:*
    
      - warning_type -- Warning type.
      - add_msg -- Message to add to the error.
    """

    # Fetching warning type
    try:
      warning_number = self.warning_dictionary[warning_type][self.warning_number]
      warning_message = self.warning_dictionary[warning_type][self.warning_message]
    except (KeyError, IndexError):
      warning_number = self.warning_dictionary["DEFAULT_WARNING"][self.warning_number]
      warning_message = self.warning_dictionary["DEFAULT_WARNING"][self.warning_message]

    # Handling warning
    complete_warning_message = ("--------------------------------------------------\n"
                                "Warning Number: " + str(warning_number) + ".\n"
                                "Program: " + self.program_name + ".\n"
                                "Report: " + warning_message + " " + add_msg + "\n"
                                "--------------------------------------------------")
    print(complete_warning_message, file=sys.stderr)


###################################################################################################
# Auxiliary Functions as Static Methods
###################################################################################################

class AuxiliaryFunctions:
  """Class for small auxiliary functions.
  """

  @staticmethod
  def string_is_int(s):
    """Verifies whether a string is a representation of a numeric integer.
    
    *Keyword arguments:*
    
      - s -- String to be verified.
    
    *Return:*
    
      - return -- True if it is a numeric integer or false if it is not.
    """
    try:
      int(s)
      return True
    except ValueError:
      return False

  @staticmethod
  def string_is_float(s):
    """Verifies whether a string is a representation of a numeric float.
    
    *Keyword arguments:*
    
      - s -- String to be verified.
    
    *Return:*
    
      - return -- True if it is a numeric float or false if it is not.
    """
    try:
      float(s)
      return True
    except ValueError:
      return False

  @staticmethod
  def which(program):
    """Fetch path to an executable program.
    
    *Keyword arguments:*
    
      - program -- Program to check.
    
    *Return:*
    
      - return -- Returns the path to the executable program or None if it does not exist or it is not executable.
    """

    def is_exe(fpath):
      return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
      if is_exe(program):
        return program
    else:
      for path in os.environ["PATH"].split(os.pathsep):
        path = path.strip('"')
        exe_file = os.path.join(path, program)
        if is_exe(exe_file):
          return exe_file
    return None

  @staticmethod
  def overlap(): # TODOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    pass # TODOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

  @staticmethod
  def revcomp(s):
    """Reverse complement a DNA (A, C, G, T, N) sequence.
    
    *Keyword arguments:*
    
      - s -- DNA sequence string to be reverse complemented.
    
    *Return:*
    
      - return -- The reverse complement of the sequence, if it exists.
    """
    revDict = dict([("A", "T"), ("T", "A"), ("C", "G"), ("G", "C"), ("N", "N")])
    return "".join([revDict[e] for e in s[::-1]])
    # TODO - Add error to unknown alphabet.

