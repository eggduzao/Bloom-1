"""
XXX
===================
Placeholder.

# Operational Options
--organism hg19
--resolution 10000
--input-type hic
--ncpu 8
--temporary-location /home/johndoe/

# TODO

# Required Arguments
input-matrix
output-matrix
output-contacts

Authors: Eduardo G. Gusmao.

"""

###################################################################################################
# Libraries
###################################################################################################

# Python
import os
import optparse

# Internal
from bloom.__version__ import __version__
from bloom.contact_map import InputFileType
from bloom.util import PassThroughOptionParser, ErrorHandler, AuxiliaryFunctions

# External

###################################################################################################
# Basic Objects
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
    self.usage_message = ("%prog [options] <Contact Matrix>\n\n"

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

    # Operational Options
    organism_help = ("organism does a lot of things like bla bla bla bla bla bla bla bla bl"
                     "bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla ")
    self.add_option("o", "organism", "string", "string", None, organism_help)

    resolution_help = ("resolution does a lot of things like bla bla bla bla bla bla bla bla bl"
                        "bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla ")
    self.add_option("r", "resolution", "int", "INT", None, resolution_help)

    input_type_help = ("input_type does a lot of things like bla bla bla bla bla bla bla bla bl"
                     "bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla ")
    self.add_option("i", "input_type", "string", "string", None, input_type_help)

    cores_help = ("cores does a lot of things like bla bla bla bla bla bla bla bla bl"
                  "bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla")
    self.add_option("c", "cores", "int", "INT", None, cores_help)

    temporary_location_help = ("temporary_location does a lot of things like bla bla bla bla bla "
                               "bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla b")
    self.add_option("t", "temporary_location", "string", "PATH", os.getcwd(), temporary_location_help)

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
    if(len(self.arguments) != 3):
      self.error_handler.throw_error("TODO") # TODO

    # Arguments
    input_matrix = self.arguments[0]
    output_matrix = self.arguments[1]
    output_matrix_directory = os.path.dirname(os.path.abspath(os.path.expanduser(output_matrix)))
    output_contacts = self.arguments[2]
    output_contacts_directory = os.path.dirname(os.path.abspath(os.path.expanduser(output_contacts)))

    # Verify arguments
    if(not os.path.isfile(input_matrix)):
      self.error_handler.throw_error("TODO") # TODO
    if(not os.path.isdir(output_matrix_directory)):
      self.error_handler.throw_error("TODO") # TODO
    if(not os.path.isdir(output_contacts_directory)):
      self.error_handler.throw_error("TODO") # TODO

    # Operational ptions
    organism = self.options.organism
    resolution = self.options.resolution
    input_type = self.options.input_type
    cores = self.options.cores
    temporary_location = os.path.abspath(os.path.expanduser(self.options.temporary_location))

    # Verify operational options
    if(not AuxiliaryFunctions.string_is_int(resolution)):
      self.error_handler.throw_error("TODO") # TODO
    if(input_type == "hic"): input_type = InputFileType.HIC
    elif(input_type == "cool"): input_type = InputFileType.COOL
    elif(input_type == "mcool"): input_type = InputFileType.MCOOL
    elif(input_type == "sparse"): input_type = InputFileType.SPARSE
    else: input_type = InputFileType.UNKNOWN
    if(not AuxiliaryFunctions.string_is_int(cores)):
      self.error_handler.throw_error("TODO") # TODO
    if(not AuxiliaryFunctions.string_is_validpath(temporary_location)):
      self.error_handler.throw_error("TODO") # TODO
    else:
      temporary_location = self.create_temporary_directory(temporary_location)

  #############################################################################
  # Auxiliary Operations
  #############################################################################

  def create_temporary_directory(self, temporary_location):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Creating temorary directory as: temporary folder / input file name
    input_contact_matrix_name = os.path.splitext(os.path.basename(self.input_file_name))[0]
    temporary_directory = None
    try:
      temporary_directory = os.path.join(temporary_location, input_contact_matrix_name)
      temp_creation_command = ["mkdir", "-p", temporary_directory]
      temp_creation_process = subprocess.run(temp_creation_command , stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)
    except Exception:
      self.error_handler.throw_error("TODO") # TODO

    # Returning the name of the temporary directory
    return temporary_directory


