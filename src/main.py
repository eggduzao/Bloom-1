"""
Main Module
===================
Placeholder.

Authors: Eduardo G. Gusmao.

"""

###################################################################################################
# Libraries
###################################################################################################

# Python


# Internal
from src.__version__ import __version__
from src.util import ErrorHandler
from src.arguments import ArgumentParser



# External







###################################################################################################
# Tool Execution
###################################################################################################

def main():
  """Entry point to the execution of Bloom.
    
  *Keyword arguments:*
    
    - None.
    
  *Return:*
    
    - None.
  """

  ###############################################################################
  # Initialization
  ###############################################################################

  # Initializing ErrorHandler
  error_handler = ErrorHandler()

  ###############################################################################
  # Parameter Handling
  ###############################################################################

  # Initializing argument parser
  argument_parser = ArgumentParser()
  args = argument_parser.arguments
  opts = argument_parser.options

  print(opts.parameter1)

  print("The analysis finished successfully.")


