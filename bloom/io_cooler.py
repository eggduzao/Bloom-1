"""
IO Cooler Module
===================
Placeholder.

Authors: Eduardo G. Gusmao.

"""

###################################################################################################
# Libraries
###################################################################################################

# Python
import os
import gc
import sys
import codecs
import traceback
import subprocess
import configparser
import multiprocessing

# Internal
from bloom.util import ConfigurationFile, ChromosomeSizes, ErrorHandler
from bloom.io_bedgraph import Bedgraph

# External


###################################################################################################
# Cooler Auxiliary Class
###################################################################################################

class Cooler(ConfigurationFile):
  """This class represents TODO.

  *Keyword arguments:*

    - argument1 -- Short description. This argument represents a long description. It can be:
      - Possibility 1: A possibility 1.
      - Possibility 2: A possibility 2.

    - argument2 -- Short description. This argument represents a long description. It can be:
      - Possibility 1: A possibility 1.
      - Possibility 2: A possibility 2.
  """

  def __init__(self, organism, ncpu):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Configuration file initialization
    ConfigurationFile.__init__(self)
    self.cooler_command = self.config.get("Cooler", "command")
    self.chromosome_sizes_file_name = os.path.join(self.bloom_data_path, self.config.get("ChromosomeSizes", organism))
    self.bed_graph_handler = Bedgraph(organism, ncpu)

    # Auxiliary Parameters
    self.ncpu = ncpu
    self.organism = organism
    self.process_queue = []
    self.cooler_resolution_list = [str(i*1000) for i in [1, 5, 10, 50, 100, 500, 1000]] # Please check: http://dcic.4dnucleome.org/data%20standards/

    # Error handler
    self.error_handler = ErrorHandler()

  def dump(self, resolution, region1, region2, input_file_name, output_file_name):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
  
    # Execution of Cooler's dump
    if(resolution):
      command = [self.cooler_command, "dump", "-t", "pixels", "--join", "-r", region1, "-r2", region2,
                 "::".join(input_file_name, "resolutions/" + str(resolution)), ">", output_file_name]
    else:
      command = [self.cooler_command, "dump", "-t", "pixels", "--join", "-r", region1, "-r2", region2,
                 input_file_name, ">", output_file_name]
    dump_process = subprocess.run(command , stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)
    return dump_process

  def add_dump(self, resolution, region1, region2, input_file_name, output_file_name):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    self.process_queue.append((resolution, region1, region2, input_file_name, output_file_name))

  def run_dump(self):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    
    pool = multiprocessing.Pool(self.ncpu)
    dump_process_output = pool.starmap(self.dump, [arguments for arguments in self.process_queue])
    pool.close()
    pool.join()
    pool = None
    self.process_queue = None
    gc.collect()
    successful_execution = True
    for cp in dump_process_output:
      try:
        cp.check_returncode()
      except subprocess.CalledProcessError:
        successful_execution = False # TODO - Error: One or more processes didnt execute correctly.

    return successful_execution

  def load(self, genome_id, resolution, sparse_matrix_dictionary, temporary_location, output_file_name):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Sparse matrix to bedgraph
    bedgraph_file_name = temporary_location + "bedgraph_file_name.bg2"
    self.bed_graph_handler.load(resolution, sparse_matrix_dictionary, bedgraph_file_name)
  
    # Execution of Cooler's dump
    command = [self.cooler_command, "load", "-f", "bg2", "--assembly", genome_id, "--count-as-float",
               ":".join(self.chromosome_sizes_file_name, str(resolution)), bedgraph_file_name, output_file_name]
    load_process = subprocess.run(command , stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)

    # Remove the temporary bedgraph file
    remove_command = ["rm", "-rf", bedgraph_file_name]
    remove_process = subprocess.run(remove_command , stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)

    return dump_process

  def add_load(self, genome_id, resolution, sparse_matrix_dictionary, temporary_location, output_file_name):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    self.process_queue.append((genome_id, resolution, sparse_matrix_dictionary, temporary_location, output_file_name))

  def run_load(self):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    
    pool = multiprocessing.Pool(self.ncpu)
    load_process_output = pool.starmap(self.load, [arguments for arguments in self.process_queue])
    pool.close()
    pool.join()
    pool = None
    self.process_queue = None
    gc.collect()
    successful_execution = True
    for cp in load_process_output:
      try:
        cp.check_returncode()
      except subprocess.CalledProcessError:
        successful_execution = False # TODO - Error: One or more processes didnt execute correctly.

  def identify_minimal_resolution(self, input_file_name, temporary_location, region = "chr1:1,000,000-5,000,000"):

    # Current result resolution
    resolution = None

    # Use built-in functions over the possible resolution list
    for res in self.cooler_resolution_list:

      # Add job
      output_file_name = os.path.join(temporary_location, "res_test_" + res + ".txt")
      self.add_dump(res, region, region, input_file_name, output_file_name)

      # Run job
      successful_execution = self.run_dump()

      # Remove temporary files
      remove_command = ["rm", "-rf", output_file_name]
      remove_process = subprocess.run(remove_command , stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)

      # Check resolution
      if(successful_execution):
        resolution = res
        break

    # Return
    return resolution

  def filetype_is_cooler(self, input_file_name, temporary_location, region = "chr1:1,000,000-5,000,000"):

    # Current result resolution
    is_cooler = False

    # Use built-in functions over the possible resolution list
    for res in self.cooler_resolution_list:

      # Add job
      output_file_name = os.path.join(temporary_location, "file_test_" + res + ".txt")
      self.add_dump(res, region, region, input_file_name, output_file_name)

      # Run job
      successful_execution = self.run_dump()

      # Remove temporary files
      remove_command = ["rm", "-rf", output_file_name]
      remove_process = subprocess.run(remove_command , stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)

      # Check resolution
      if(successful_execution):
        is_cooler = True
        break

    # Return
    return is_cooler 

