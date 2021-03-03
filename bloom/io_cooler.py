from __future__ import print_function
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
from bloom.contact_map import ContactMap
from bloom.util import ConfigurationFile, ChromosomeSizes, ErrorHandler
from bloom.io_bedgraph import Bedgraph

# External


###################################################################################################
# Cooler Class
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

    # Auxiliary Parameters
    self.ncpu = ncpu
    self.organism = organism
    self.process_queue = []
    self.cooler_resolution_list = [str(i*1000) for i in [1, 2, 5, 10, 25, 50, 100, 250, 500, 1000, 2500, 5000, 10000]] # Please check: http://dcic.4dnucleome.org/data%20standards/

    # Error handler
    self.error_handler = ErrorHandler()

    # Bedgraph handler
    self.bed_graph_handler = Bedgraph(self.organism, self.ncpu)

  #############################################################################
  # Read: Single Cooler (.cool) -> Bedgraph (.bg)
  #############################################################################

  def dump_single(self, region1, region2, input_file_name, output_file_name):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
  
    # Execution of Cooler's dump
    dump_command = [self.cooler_command, "dump", "-t", "pixels", "--join", "-r", region1, "-r2", region2, "-o", output_file_name, input_file_name]
    dump_process = subprocess.run(dump_command, stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)

    # Return dump process
    return dump_process

  def add_dump_single(self, region1, region2, input_file_name, output_file_name):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Append job to queue
    self.process_queue.append((region1, region2, input_file_name, output_file_name))

  def run_dump_single(self, return_type = "success"):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    
    # Execute job queue
    pool = multiprocessing.Pool(self.ncpu)
    dump_process_output = pool.starmap(self.dump_single, [arguments for arguments in self.process_queue])
    pool.close()
    pool.join()

    # Clean queue
    pool = None
    self.process_queue = []
    gc.collect()

    # Check execution status
    successful_execution = True
    for cp in dump_process_output:
      try:
        cp.check_returncode()
      except subprocess.CalledProcessError:
        successful_execution = False # TODO - Error: One or more processes didnt execute correctly.
        # raise
        # self.error_handler.throw_error("TODO") # TODO

    # Return mode
    if(return_type == "success"):
      return successful_execution
    elif(return_type == "process_out"):
      return dump_process_output
    else:
      return None

  #############################################################################
  # Read: Multiple Cooler (.mcool) -> Bedgraph (.bg)
  #############################################################################

  def dump_multiple(self, resolution, region1, region2, input_file_name, output_file_name):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
  
    # Execution of Cooler's dump
    command = [self.cooler_command, "dump", "-t", "pixels", "--join", "-r", region1, "-r2", region2, "-o", output_file_name, "::".join([input_file_name, "resolutions/" + str(resolution)])]
    dump_process = subprocess.run(command , stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)

    # Return dump process
    return dump_process

  def add_dump_multiple(self, resolution, region1, region2, input_file_name, output_file_name):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Append job to queue
    self.process_queue.append((resolution, region1, region2, input_file_name, output_file_name))

  def run_dump_multiple(self, return_type = "success"):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
 
    # Execute job queue   
    pool = multiprocessing.Pool(self.ncpu)
    dump_process_output = pool.starmap(self.dump_multiple, [arguments for arguments in self.process_queue])
    pool.close()
    pool.join()

    # Clean queue
    pool = None
    self.process_queue = []
    gc.collect()

    # Check execution status
    successful_execution = True
    for cp in dump_process_output:
      try:
        cp.check_returncode()
      except subprocess.CalledProcessError:
        successful_execution = False # TODO - Error: One or more processes didnt execute correctly.
        # raise
        # self.error_handler.throw_error("TODO") # TODO

    # Return mode
    if(return_type == "success"):
      return successful_execution
    elif(return_type == "process_out"):
      return dump_process_output
    else:
      return None

  #############################################################################
  # Write: Contact Map -> Cooler (.(m)cool)
  #############################################################################

  def load(self, contact_map, temporary_location, output_file_name):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Sparse matrix to bedgraph
    bedgraph_file_name = os.path.join(temporary_location, "bedgraph_file_name.bg2")
    self.bed_graph_handler.load(contact_map, bedgraph_file_name)
  
    # Execution of Cooler's dump
    command = [self.cooler_command, "load", "-f", "bg2", "--assembly", contact_map.organism, "--count-as-float",
               ":".join([self.chromosome_sizes_file_name, str(contact_map.resolution)]), bedgraph_file_name, output_file_name]
    load_process = subprocess.run(command, stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)

    # Remove the temporary bedgraph file
    remove_command = ["rm", "-rf", bedgraph_file_name]
    remove_process = subprocess.run(remove_command , stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)

    # Return load process
    return load_process

  def add_load(self, contact_map, temporary_location, output_file_name):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Append job to queue
    self.process_queue.append((contact_map, temporary_location, output_file_name))

  def run_load(self, return_type = "success"):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    
    # Execute job queue
    pool = multiprocessing.Pool(self.ncpu)
    load_process_output = pool.starmap(self.load, [arguments for arguments in self.process_queue])
    pool.close()
    pool.join()

    # Clean queue
    pool = None
    self.process_queue = []
    gc.collect()

    # Check execution status
    successful_execution = True
    for cp in load_process_output:
      try:
        cp.check_returncode()
      except subprocess.CalledProcessError:
        successful_execution = False # TODO - Error: One or more processes didnt execute correctly.
        # raise
        # self.error_handler.throw_error("TODO") # TODO

    # Return mode
    if(return_type == "success"):
      return successful_execution
    elif(return_type == "process_out"):
      return load_process_output
    else:
      return None

  #############################################################################
  # Auxiliary IO identifying methods
  #############################################################################

  def identify_minimal_resolution(self, input_file_name, temporary_location, check_type = "cool", region = "chr1:0-5,000,000"):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Current result resolution
    resolution = None

    # If checking type is single cooler (.cool)
    if(check_type == "cool"):

      # Temporary dump file
      output_file_name = os.path.join(temporary_location, "res_test.txt")

      # Add single dump to queue
      self.add_dump_single(region, region, input_file_name, output_file_name)

      # Run single dump job
      successful_execution = self.run_dump_single()

      # Check resolution
      if(successful_execution):
        resolution = self.bed_graph_handler.identify_minimal_resolution(output_file_name)

      # Remove temporary files
      remove_command = ["rm", "-rf", output_file_name]
      remove_process = subprocess.run(remove_command , stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)

    # If checking type is multiple cooler (.mcool)
    if(check_type == "mcool"):

      # Use built-in functions over the possible resolution list
      for res in self.cooler_resolution_list:

        # Temporary dump file
        output_file_name = os.path.join(temporary_location, "res_test_" + res + ".txt")

        # Add single dump to queue
        self.add_dump_multiple(res, region, region, input_file_name, output_file_name)

        # Run single dump job
        successful_execution = self.run_dump_multiple()

        # Remove temporary files
        remove_command = ["rm", "-rf", output_file_name]
        remove_process = subprocess.run(remove_command , stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)

        # Check resolution
        if(successful_execution):
          resolution = res
          break

    # Return
    return resolution

  def filetype_is_cooler(self, input_file_name, temporary_location, check_type = "cool", region = "chr1:0-5,000,000"):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Current result resolution
    is_cooler = False

    # If checking type is single cooler (.cool)
    if(check_type == "cool"):

      # Temporary dump file
      output_file_name = os.path.join(temporary_location, "res_test.txt")

      # Add single dump to queue
      self.add_dump_single(region, region, input_file_name, output_file_name)

      # Run single dump job
      successful_execution = self.run_dump_single()

      # Remove temporary files
      remove_command = ["rm", "-rf", output_file_name]
      remove_process = subprocess.run(remove_command , stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)

      # Check resolution
      if(successful_execution):
        is_cooler = True

    # If checking type is multiple cooler (.mcool)
    if(check_type == "mcool"):

      # Use built-in functions over the possible resolution list
      for res in self.cooler_resolution_list:

        # Temporary dump file
        output_file_name = os.path.join(temporary_location, "res_test_" + res + ".txt")

        # Add single dump to queue
        self.add_dump_multiple(res, region, region, input_file_name, output_file_name)

        # Run single dump job
        successful_execution = self.run_dump_multiple()

        # Remove temporary files
        remove_command = ["rm", "-rf", output_file_name]
        remove_process = subprocess.run(remove_command , stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)

        # Check resolution
        if(successful_execution):
          is_cooler = True
          break

    # Return
    return is_cooler 

    """
    # Current result
    is_cooler = False

    # Temporary list file
    output_file_name = os.path.join(temporary_location, "iscool_test.txt")

    # Execute list command
    list_command = ["cooler", "ls", input_file_name, ">", output_file_name]
    list_process = subprocess.run(list_command , stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)

    # Number of lines in output file
    try:
      nlines = AuxiliaryFunctions.number_of_lines(output_file_name)
    except Exception:
      pass
      # self.error_handler.throw_error("TODO") # TODO

    # If checking type is single cooler (.cool)
    if(check_type == "cool" and nlines == 1):
      is_cooler = True

    # If checking type is multiple cooler (.mcool)
    if(check_type == "mcool" and nlines > 1):
      is_cooler = True

    # Return
    return is_cooler 

    """

