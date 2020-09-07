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
from time import time

# Internal
from bloom.__version__ import __version__
from bloom.util import ErrorHandler, Juicer
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






  """
  # Test 1
  print(opts.parameter1)
  print("The analysis finished successfully.")
  """

  """
  # Test 2
  t1 = time()
  juicer = Juicer(8)
  input_file_name = "/usr/users/egadegu/Projects/Papantonis_Bloom/Data/5_dn_isHiC_Human_CM/dnHiC/HEART.hic"
  juicer.add_dump("100000", "1:1000000:20000000", "1:1000000:20000000", input_file_name, "/usr/users/egadegu/Projects/Papantonis_Bloom/Bloom/a2.txt")
  juicer.add_dump("100000", "2:1000000:20000000", "2:1000000:20000000", input_file_name, "/usr/users/egadegu/Projects/Papantonis_Bloom/Bloom/b2.txt")
  juicer.add_dump("100000", "3:1000000:20000000", "3:1000000:20000000", input_file_name, "/usr/users/egadegu/Projects/Papantonis_Bloom/Bloom/c2.txt")
  juicer.add_dump("100000", "4:1000000:20000000", "4:1000000:20000000", input_file_name, "/usr/users/egadegu/Projects/Papantonis_Bloom/Bloom/d2.txt")
  juicer.add_dump("100000", "5:1000000:20000000", "5:1000000:20000000", input_file_name, "/usr/users/egadegu/Projects/Papantonis_Bloom/Bloom/e2.txt")
  juicer.add_dump("100000", "6:1000000:20000000", "6:1000000:20000000", input_file_name, "/usr/users/egadegu/Projects/Papantonis_Bloom/Bloom/f2.txt")
  juicer.add_dump("100000", "7:1000000:20000000", "7:1000000:20000000", input_file_name, "/usr/users/egadegu/Projects/Papantonis_Bloom/Bloom/g2.txt")
  juicer.add_dump("100000", "8:1000000:20000000", "8:1000000:20000000", input_file_name, "/usr/users/egadegu/Projects/Papantonis_Bloom/Bloom/h2.txt")
  juicer.run_dump()
  t2 = time()
  print(t2-t1)
  # 1 CPU = 14.55 s
  # 2 CPU = 4.04 s
  # cmp --silent a1.txt a2.txt && echo 'equal' || echo 'different'
  # cmp --silent b1.txt b2.txt && echo 'equal' || echo 'different'
  # cmp --silent c1.txt c2.txt && echo 'equal' || echo 'different'
  # cmp --silent d1.txt d2.txt && echo 'equal' || echo 'different'
  # cmp --silent e1.txt e2.txt && echo 'equal' || echo 'different'
  # cmp --silent f1.txt f2.txt && echo 'equal' || echo 'different'
  # cmp --silent g1.txt g2.txt && echo 'equal' || echo 'different'
  # cmp --silent h1.txt h2.txt && echo 'equal' || echo 'different'
  """

  # Test 3
  t1 = time()
  juicer = Juicer(2)
  input_file_name = "/usr/users/egadegu/Projects/Papantonis_Bloom/Data/5_dn_isHiC_Human_CM/dnHiC/HEART.hic"
  juicer.add_dump("100000", "1:1000000:20000000", "1:1000000:20000000", input_file_name, "/usr/users/egadegu/Projects/Papantonis_Bloom/Bloom/a2.txt")
  juicer.add_dump("100000", "2:1000000:20000000", "2:1000000:20000000", input_file_name, "/usr/users/egadegu/Projects/Papantonis_Bloom/Bloom/b2.txt")
  juicer.run_dump()
  t2 = time()
  print(t2-t1)

  #input_file_name = "/usr/users/egadegu/Projects/Papantonis_Bloom/Data/5_dn_isHiC_Human_CM/dnHiC/HEART.mcool"

