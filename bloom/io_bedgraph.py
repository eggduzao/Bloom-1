from __future__ import print_function
"""
IO Bedgraph Module
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
from bloom.util import ConfigurationFile, ChromosomeSizes, ErrorHandler, AuxiliaryFunctions

# External
import numpy


###################################################################################################
# Bedgraph Class
###################################################################################################

class Bedgraph(ConfigurationFile):
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

    # Auxiliary Parameters
    self.ncpu = ncpu
    self.organism = organism
    self.process_queue = []

    # Chromosome sizes
    self.chromosome_sizes = ChromosomeSizes(self.organism)

    # Error handler
    self.error_handler = ErrorHandler()

  #############################################################################
  # Read: Bedgraph (.bg) -> Contact Map
  #############################################################################

  def dump(self, chromosome, input_file_name, contact_map, logit = False, pseudocount = 1.0):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
 
    # Open bedgraph input file [chrom1, pos11, pos12, chrom2, po21, pos22, count]
    input_file = codecs.open(input_file_name, "rU", "utf8")
   
    # Iterate through bedgraph file
    for line in input_file:

      # Get file columns
      ll = line.strip().split("\t")
      chrom1 = ll[0]
      chrom2 = ll[3]
      if(chrom1 != chromosome or chrom1 != chrom2): continue
      pos11 = int(ll[1])
      pos21 = int(ll[4])
      count = float(ll[6])

      # Value or log to write
      if(logit): count = numpy.log(count + 1) + pseudocount

      # Updating matrix
      if(count > 0):
        contact_map.add(chromosome, min(pos11, pos21), max(pos11, pos21), count)

    # Close bedgraph
    input_file.close()

    # Successful execution
    return True

  def add_dump(self, chromosome, input_file_name, contact_map, logit = False, pseudocount = 1.0):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Append job to queue
    self.process_queue.append((chromosome, input_file_name, contact_map, logit, pseudocount))

  def run_dump(self, return_type = "success"):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    
    # Execute job queue
    pool = multiprocessing.Pool(self.ncpu)
    dump_process_output = pool.starmap(self.dump, [arguments for arguments in self.process_queue])
    pool.close()
    pool.join()

    # Clean queue
    pool = None
    self.process_queue = []
    gc.collect()

    # Check execution status
    successful_execution = True
    for cp in dump_process_output:
      if(not cp):
        successful_execution = False
        pass
        # self.error_handler.throw_error("TODO") # TODO

    # Return mode
    if(return_type == "success"):
      return successful_execution
    elif(return_type == "process_out"):
      return dump_process_output
    else:
      return None

  #############################################################################
  # Write: Contact Map -> Bedgraph (.bg)
  #############################################################################

  def load(self, contact_map, output_file_name, start_index = 0):
    """Returns TODO.
       Important = UPPER MATRIX TO BEDGRAPH
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Open output bedgraph file
    output_file = codecs.open(output_file_name, "w", "utf8")

    # Iterating on valid chromosomes
    for chromosome in contact_map.valid_chromosome_list:

      # Iterating on matrix
      for key, value in contact_map.matrix[chromosome].items():

        # Write entry to bedgraph file
        entry = [chromosome, str(key[0]), str(key[0] + contact_map.resolution), chromosome, str(key[1]), str(key[1] + contact_map.resolution), str(value)]
        output_file.write("\t".join(entry) + "\n")

    # Close temporary output file
    output_file.close()

    # Successful execution
    return True

  def add_load(self, contact_map, output_file_name, start_index = 0):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Append job to queue
    self.process_queue.append((contact_map, output_file_name, start_index))

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
      if(not cp):
        successful_execution = False
        pass
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

  def identify_minimal_resolution(self, input_file_name, multiple_check = 1000):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Current result resolution
    resolution = None

    # Open bedgraph input file [chrom1, pos11, pos12, chrom2, po21, pos22, count]
    input_file = codecs.open(input_file_name, "rU", "utf8")
   
    # Check resolution
    ll = input_file.readline().strip().split("\t")
    pos11 = int(ll[1])
    pos12 = int(ll[2])
    raw_resolution = max(pos11, pos12) - min(pos11, pos12)
    if(int(raw_resolution % multiple_check) == 0):
      resolution = raw_resolution
    elif(int(raw_resolution % multiple_check) in list(range(1, int(multiple_check / 2)))): 
      resolution = AuxiliaryFunctions.floor_multiple(raw_resolution, multiple_check)
    else:
      resolution = AuxiliaryFunctions.ceil_multiple(raw_resolution, multiple_check)

    # Close bedgraph
    input_file.close()

    # Return
    return resolution

  def filetype_is_bedgraph(self, input_file_name):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Current result resolution
    is_bedgraph = True

    # Open bedgraph input file [chrom1, pos11, pos12, chrom2, po21, pos22, count]
    try:
      input_file = codecs.open(input_file_name, "rU", "utf8")
      line = input_file.readline()
    except Exception:
      is_bedgraph = False

    # Iterate through bedgraph file
    if(is_bedgraph):
      for line in input_file:

        # Get file columns
        try:
          ll = line.strip().split("\t")
          if(len(ll) != 7):
            is_bedgraph = False
            break
        except Exception:
          is_bedgraph = False
          break

    # Close bedgraph
    if(is_bedgraph):
      input_file.close()

    # Return
    return is_bedgraph 


