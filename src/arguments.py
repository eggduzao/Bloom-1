"""
XXX
===================
Placeholder.

Authors: Eduardo G. Gusmao.

"""

###################################################################################################
# Libraries
###################################################################################################

# Python
from os import getcwd
from optparse import SUPPRESS_HELP

# Internal
from src.__version__ import __version__
from src.util import PassThroughOptionParser, ErrorHandler

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

  def add_option(self, option_variable, option_type, meta_type, default_option, help_message):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Creating option name based on variable name
    option_name = "--" + option_variable.replace("_", "-")

    # Correct meta_type to uppercase
    meta_type = meta_type.upper()

    # Hidden options from help message
    if(not help_message): help_message = SUPPRESS_HELP
    
    # Add option
    if(option_type == "bool" and default_option):
      self.parser.add_option(option_name, dest = option_variable, action="store_false",
                             default=default_option, help=help_message)
    elif(option_type == "bool"):
      self.parser.add_option(option_name, dest = option_variable, action="store_true",
                             default=default_option, help=help_message)
    else:
      self.parser.add_option(option_name, dest = option_variable, type = option_type,
                             metavar = meta_type, default=default_option, help=help_message)


  def load_options(self):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    parameter_1_help = ("Parameter 1 does a lot of things like bla bla bla bla bla bla bla bla bl"
                        "bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla ")
    self.add_option("parameter1", "string", "PATH", getcwd(), parameter_1_help)

    parameter_2_help = ("Parameter 2 does a lot of things like bla bla bla bla bla bla bla bla bl"
                        "bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla ")
    self.add_option("parameter_2", "string", "string", None, parameter_2_help)

    parameter_3_help = ("Parameter 3 does a lot of things like bla bla bla bla bla bla bla bla bl"
                        "bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla ")
    self.add_option("parameter3", "int", "INT", 1, parameter_3_help)

    parameter_4_help = ("Parameter 4 does a lot of things like bla bla bla bla bla bla bla bla bl"
                        "bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla ")
    self.add_option("parameter4_44_4", "float", "FLOAT", 0.1, parameter_4_help)

    parameter_5_help = ("Parameter 5 does a lot of things like bla bla bla bla bla bla bla bla bl"
                        "bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla ")
    self.add_option("parameter5", "bool", "BOOLFALSE", None, parameter_5_help)

    parameter_6_help = None
    self.add_option("parameter6", "string", "STRING", "STRING", parameter_6_help)

  def load_options_and_arguments(self):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Processing Options
    self.options, self.arguments = self.parser.parse_args()

  def argument_validity_check(self): # TODO
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    if len(arguments) < 1:
        print(usage_message)
        exit(1)
    # if(not arguments or len(arguments) > 1): error_handler.throw_error("FP_WRONG_ARGUMENT") # TODO

