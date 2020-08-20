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
from bloom.__version__ import __version__
from bloom.util import ErrorHandler
from bloom.arguments import ArgumentParser

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







  # Test
  print(opts.parameter1)
  print("The analysis finished successfully.")


