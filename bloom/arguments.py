from __future__ import print_function
"""
Arguments Module
===================
Placeholder.

Authors: Eduardo G. Gusmao.

"""

###################################################################################################
# Libraries
###################################################################################################

# Python
import os
import sys
import optparse
import subprocess
import multiprocessing

# Internal
from bloom.__version__ import __version__
from bloom.io import InputFileType
from bloom.util import PassThroughOptionParser, ErrorHandler, AuxiliaryFunctions

# External

###################################################################################################
# ArgumentParser Class
###################################################################################################

class ArgumentParser():
  """This class represents TODO.

  *Keyword arguments:*

    - argument1 -- Short description. This argument represents a long description. It can be:
      - Possibility 1: A possibility 1.
      - Possibility 2: A possibility 2.

    - argument2 -- Short description. This argument represents a long description. It can be:
      - Possibility 1: A possibility 1.
      - Possibility 2: A possibility 2.
  """

  def __init__(self):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    
    # Initializing class objects
    self.error_handler = ErrorHandler() 
    self.usage_message = None
    self.version_message = None
    self.parser = None
    self.options = None
    self.arguments = None

    # Load options and arguments
    self.load_usage_message()
    self.load_version_message()
    self.load_parser()
    self.load_options()
    self.load_options_and_arguments()
    self.option_argument_validity_check()

  #############################################################################
  # Main Operations
  #############################################################################

  def load_usage_message(self):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Usage message
    self.usage_message = ("%prog [options] <Input Matrix File> <Output Matrix File> <Output Contact File>\n\n"

                          "-----------------------------------------------------------------------------\n"
                          "Bloom is a program designed to reveal occult patterns\n"
                          "in 3C-based data, such as Hi-C.\n\n"

                          "The Contact Matrix can be in one of these formats:\n"
                          "- Format 1.\n" # TODO
                          "  Format 2.\n" # TODO
                          "- Format 3.\n\n" # TODO

                          "Bloom's documentation can be found at:\n"
                          "WEBSITE\n\n" # TODO

                          "For more information, please refer to:\n"
                          "WEBSITE\n\n" # TODO

                          "For further questions or comments please refer to:\n"
                          "ORIGINAL PAPER\n"  # TODO
                          "-----------------------------------------------------------------------------")


  def load_version_message(self):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Version message
    self.version_message = "Bloom - Analysis of 3C-based contact matrices. Version: " + str(__version__)

  def load_parser(self):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Initializing Option Parser
    self.parser = PassThroughOptionParser(usage = self.usage_message, version = self.version_message)

  def add_option(self, option_alias, option_variable, option_type, meta_type, default_option, help_message):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Creating option alias
    option_alias = "-" + option_alias

    # Creating option name based on variable name
    option_name = "--" + option_variable.replace("_", "-")

    # Correct meta_type to uppercase
    meta_type = meta_type.upper()

    # Hidden options from help message
    if(not help_message): help_message = optparse.SUPPRESS_HELP
    
    # Add option
    if(option_type == "bool" and default_option):
      self.parser.add_option(option_alias, option_name, dest = option_variable, action="store_false",
                             default=default_option, help=help_message)
    elif(option_type == "bool"):
      self.parser.add_option(option_alias, option_name, dest = option_variable, action="store_true",
                             default=default_option, help=help_message)
    else:
      self.parser.add_option(option_alias, option_name, dest = option_variable, type = option_type,
                             metavar = meta_type, default=default_option, help=help_message)

  def load_options(self):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Options

    chromosome_help = ("If this option is used, bloom will perform its analysis only on the "
                     "desired chromosome. The chromosome should follow the input file's "
                     "chromosome name conventions.")
    self.add_option("x", "chromosome", "string", "STRING", None, chromosome_help)

    organism_help = ("The organism in which the analysis is being applied to. It can be one of "
                     "the following: hg19, hg38, mm9, mm10. If you are analyzing other organism, "
                     "please check the manual for the necessary files and configuration.")
    self.add_option("o", "organism", "string", "STRING", "hg19", organism_help)

    resolution_help = ("Current resolution of the input file matrix. In case of multiple resolutions "
                       "(e.g. .hic and .mcool files), please select the lowest possible resolution. "
                       "This option can be left blank and Bloom will detect the minimal resolution "
                       "of the input matrix. However, this detection can be time consuming.")
    self.add_option("r", "resolution", "int", "INT", None, resolution_help)

    input_type_help = ("The format of the input file. Bloom currently supports: hic, cool, mcool and "
                     "sparse. Sparse refers to a bedgraph3+3 file. Please refer to the manual for "
                     "details on such format. This option can be left blank and Bloom will detect "
                     "the minimal resolution. However, this detection can be time consuming.")
    self.add_option("i", "input_type", "string", "STRING", None, input_type_help)

    cores_help = ("Number of cores the user would like to dedicate for the parallel steps of "
                  "Bloom pipeline. Default to half the total number of CPU cores.")
    self.add_option("c", "cores", "int", "INT", multiprocessing.cpu_count() / 2, cores_help)

    temporary_location_help = ("A directory (does not need to exist) to write required temporary files. "
                               "Defaults to a new directory named 'btmp' in the current working directory. "
                               "All temporary files will be deleted after Bloom's execution.")
    self.add_option("t", "temporary_location", "string", "PATH", os.path.join(os.getcwd(), "btmp"), temporary_location_help)

    output_type_help = ("The desired format of the output file. Since only one matrix is generated by "
                        "Bloom, possible alternatives are: hic, cool and sparse.")
    self.add_option("y", "output_type", "string", "STRING", "sparse", output_type_help)

    output_resolution_help = ("The desired output resolution (bps). At this moment it is advised to "
                              "choose between: 1000, 2000, 5000, 10000. Currently, it is not advised "
                              "to select a resolution lower than 1000 bps.")
    self.add_option("e", "output_resolution", "int", "INT", 1000, output_resolution_help)

    avoid_distance_help = ("Distance, in bins (given the desired output resolution) in which to "
                           "avoid values of bins from the diagonal. Depending on the experimental "
                           "setup and organism, it is advised to be between 10000 to 150000 bps.")
    self.add_option("a", "avoid_distance", "int", "INT", 50000, avoid_distance_help)

    seed_help = ("Select the same pseudo-random number generator seed to be able to"
                 "replicate an experiment (same output is guaranteed for the same seed).")
    self.add_option("s", "seed", "int", "INT", 123, seed_help)


    # Hidden Options

    pre_min_contig_removed_bins_help = None
    self.add_option("A", "pre_min_contig_removed_bins", "int", "INT", 5, pre_min_contig_removed_bins_help)

    pre_remove_threshold_help = None
    self.add_option("B", "pre_remove_threshold", "float", "FLOAT", 1.0, pre_remove_threshold_help)

    sica_pvalue_threshold_help = None
    self.add_option("C", "sica_pvalue_threshold", "float", "FLOAT", 0.95, sica_pvalue_threshold_help)

    sica_bottom_bin_ext_range_help = None
    self.add_option("D", "sica_bottom_bin_ext_range", "string", "STRING", "3,10", sica_bottom_bin_ext_range_help)

    sica_left_bin_ext_range_help = None
    self.add_option("E", "sica_left_bin_ext_range", "string", "STRING", "3,10", sica_left_bin_ext_range_help)

    sica_right_bin_ext_range_help = None
    self.add_option("F", "sica_right_bin_ext_range", "string", "STRING", "1,4", sica_right_bin_ext_range_help)

    sica_top_bin_ext_range_help = None
    self.add_option("G", "sica_top_bin_ext_range", "string", "STRING", "1,4", sica_top_bin_ext_range_help)

    sica_bonuscrosslb_range_help = None
    self.add_option("H", "sica_bonuscrosslb_range", "string", "STRING", "0.25,0.3", sica_bonuscrosslb_range_help)

    sica_bonuscross_range_help = None
    self.add_option("I", "sica_bonuscross_range", "string", "STRING", "0.1,0.25", sica_bonuscross_range_help)

    sica_bonuslb_range_help = None
    self.add_option("J", "sica_bonuslb_range", "string", "STRING", "0.1,0.25", sica_bonuslb_range_help)

    goba_vertical_multiplier_help = None
    self.add_option("K", "goba_vertical_multiplier", "string", "STRING", "0.5,0.75", goba_vertical_multiplier_help)

    goba_ortogonal_multiplier_help = None
    self.add_option("L", "goba_ortogonal_multiplier", "string", "STRING", "0.1,0.3", goba_ortogonal_multiplier_help)

    goba_filling_frequency_help = None
    self.add_option("M", "goba_filling_frequency", "float", "FLOAT", 0.75, goba_filling_frequency_help)

    goba_banding_value_mult_range_help = None
    self.add_option("N", "goba_banding_value_mult_range", "string", "STRING", "0.4,0.6", goba_banding_value_mult_range_help)

    goba_banding_further_range_help = None
    self.add_option("O", "goba_banding_further_range", "string", "STRING", "0.9,0.99", goba_banding_further_range_help)

    goba_banding_frequency_help = None
    self.add_option("P", "goba_banding_frequency", "float", "FLOAT", 0.5, goba_banding_frequency_help)

    goba_outing_value_mult_range_help = None
    self.add_option("Q", "goba_outing_value_mult_range", "string", "STRING", "0.3,0.5", goba_outing_value_mult_range_help)

    goba_outing_further_range_help = None
    self.add_option("R", "goba_outing_further_range", "string", "STRING", "0.9,0.99", goba_outing_further_range_help)

    goba_outing_frequency_help = None
    self.add_option("S", "goba_outing_frequency", "float", "FLOAT", 0.34, goba_outing_frequency_help)

    dpmm_random_degrade_range_help = None
    self.add_option("T", "dpmm_random_degrade_range", "string", "STRING", "0.01,0.02", dpmm_random_degrade_range_help)

    dpmm_degrade_multiplier_help = None
    self.add_option("U", "dpmm_degrade_multiplier", "float", "FLOAT", 0.05, dpmm_degrade_multiplier_help)

    dpmm_half_length_bin_interval_help = None
    self.add_option("V", "dpmm_half_length_bin_interval", "string", "STRING", "1,5", dpmm_half_length_bin_interval_help)

    dpmm_value_range_help = None
    self.add_option("W", "dpmm_value_range", "string", "STRING", "0.0001,0.001", dpmm_value_range_help)

    dpmm_random_range_help = None
    self.add_option("X", "dpmm_random_range", "string", "STRING", "0.00001,0.0001", dpmm_random_range_help)

    dpmm_iteration_multiplier_help = None
    self.add_option("Y", "dpmm_iteration_multiplier", "int", "INT", 1, dpmm_iteration_multiplier_help)

    ifs_multiplier_help = None
    self.add_option("Z", "ifs_multiplier", "int", "INT", 1000, ifs_multiplier_help)

    ifs_min_matrix_threshold_help = None
    self.add_option("b", "ifs_min_matrix_threshold", "int", "INT", 0, ifs_min_matrix_threshold_help)

    """
    # Examples:
    
    parameter_1_help = ("Parameter 1 does a lot of things like bla bla bla bla bla bla bla bla bl"
                        "bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla ")
    self.add_option("a", "parameter1", "string", "PATH", os.getcwd(), parameter_1_help)

    parameter_2_help = ("Parameter 2 does a lot of things like bla bla bla bla bla bla bla bla bl"
                        "bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla ")
    self.add_option("b", "parameter_2", "string", "string", None, parameter_2_help)

    parameter_3_help = ("Parameter 3 does a lot of things like bla bla bla bla bla bla bla bla bl"
                        "bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla ")
    self.add_option("c", "parameter3", "int", "INT", 1, parameter_3_help)

    parameter_4_help = ("Parameter 4 does a lot of things like bla bla bla bla bla bla bla bla bl"
                        "bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla ")
    self.add_option("d", "parameter4_44_4", "float", "FLOAT", 0.1, parameter_4_help)

    parameter_5_help = ("Parameter 5 does a lot of things like bla bla bla bla bla bla bla bla bl"
                        "bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla ")
    self.add_option("e", "parameter5", "bool", "BOOLFALSE", None, parameter_5_help)

    parameter_6_help = None
    self.add_option("parameter6", "string", "STRING", "STRING", parameter_6_help)
    """

  def load_options_and_arguments(self):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Processing Options
    self.options, self.arguments = self.parser.parse_args()

  def option_argument_validity_check(self):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Verify number of arguments
    if(len(sys.argv) == 1):
      self.parser.print_help(sys.stderr)
      sys.exit(1)
    if(len(self.arguments) != 3):
      pass
      # self.error_handler.throw_error("TODO") # TODO

    # Arguments
    input_matrix = self.arguments[0]
    output_matrix = self.arguments[1]
    output_matrix_directory = os.path.dirname(os.path.abspath(os.path.expanduser(output_matrix)))
    output_contacts = self.arguments[2]
    output_contacts_directory = os.path.dirname(os.path.abspath(os.path.expanduser(output_contacts)))

    # Verify arguments
    if(not os.path.isfile(input_matrix)):
      pass
      # self.error_handler.throw_error("TODO") # TODO
    if(not os.path.isdir(output_matrix_directory)):
      pass
      # self.error_handler.throw_error("TODO") # TODO
    if(not os.path.isdir(output_contacts_directory)):
      pass
      # self.error_handler.throw_error("TODO") # TODO

    # Options
    organism = self.options.organism
    resolution = self.options.resolution
    input_type = self.options.input_type
    cores = self.options.cores
    temporary_location = os.path.abspath(os.path.expanduser(self.options.temporary_location))
    output_type = self.options.output_type
    output_resolution = self.options.output_resolution
    avoid_distance = self.options.avoid_distance
    seed = self.options.seed

    # Verify operational options
    if(not AuxiliaryFunctions.string_is_int(resolution)):
      pass
      # self.error_handler.throw_error("TODO") # TODO
    if(input_type == "hic"):
      self.options.input_type = InputFileType.HIC
    elif(input_type == "cool"):
      self.options.input_type = InputFileType.COOL
    elif(input_type == "mcool"):
      self.options.input_type = InputFileType.MCOOL
    elif(input_type == "sparse"):
      self.options.input_type = InputFileType.SPARSE
    else:
      self.options.input_type = InputFileType.UNKNOWN
    if(not AuxiliaryFunctions.string_is_int(cores)):
      pass
      # self.error_handler.throw_error("TODO") # TODO
    if(not AuxiliaryFunctions.string_is_validpath(temporary_location)):
      pass
      # self.error_handler.throw_error("TODO") # TODO
    else:
      temporary_location = self.create_temporary_directory(input_matrix, temporary_location)
    if(output_type == "hic"):
      self.options.output_type = InputFileType.HIC
    elif(output_type == "cool"):
      self.options.output_type = InputFileType.COOL
    elif(output_type == "mcool"):
      pass
      # self.error_handler.throw_error("TODO") # TODO
    elif(output_type == "sparse"):
      self.options.output_type = InputFileType.SPARSE
    else:
      # self.error_handler.throw_warning("TODO") # TODO - Output file not recognized, printing sparse instead.
      self.options.output_type = InputFileType.SPARSE
    if(not AuxiliaryFunctions.string_is_int(output_resolution)):
      pass
      # self.error_handler.throw_error("TODO") # TODO
    if(not AuxiliaryFunctions.string_is_int(avoid_distance)):
      pass
      # self.error_handler.throw_error("TODO") # TODO
    if(not AuxiliaryFunctions.string_is_int(seed)):
      pass
      # self.error_handler.throw_error("TODO") # TODO

  #############################################################################
  # Auxiliary Operations
  #############################################################################

  def create_temporary_directory(self, input_file_name, temporary_location):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Creating temorary directory as: temporary folder / input file name
    input_contact_matrix_name = os.path.splitext(os.path.basename(input_file_name))[0]
    temporary_directory = None
    try:
      temporary_directory = os.path.join(temporary_location, input_contact_matrix_name)
      temp_creation_command = ["mkdir", "-p", temporary_directory]
      temp_creation_process = subprocess.run(temp_creation_command , stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)
    except Exception:
      raise
      # self.error_handler.throw_error("TODO") # TODO

    # Returning the name of the temporary directory
    return temporary_directory


