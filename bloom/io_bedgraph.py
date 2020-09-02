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
from bloom.util import ConfigurationFile, ChromosomeSizes, ErrorHandler, AuxiliaryFunctions

# External


###################################################################################################
# Bedgraph Auxiliary Class
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

  # BEDGRAPH -> UPPER MATRIX DICTIONARY
  def dump(self, chromosome, input_file_name, sparse_matrix_dictionary):
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
      pos11 = ll[1]
      pos21 = ll[4]
      count = ll[6]
      region = ":".join([chrom1, str(min(int(pos11), int(pos21))), str(max(int(pos11), int(pos21)))])
      try:
        sparse_matrix_dictionary[region] += float(count)
      except Exception:
        sparse_matrix_dictionary[region] = float(count)

    # Close bedgraph
    input_file.close()

    # Successful execution
    return True

  def add_dump(self, chromosome, input_file_name, sparse_matrix_dictionary):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    self.process_queue.append((chromosome, input_file_name, sparse_matrix_dictionary))

  def run_dump(self, return_type = "success"):
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
      if(not cp):
        successful_execution = False
        self.error_handler.throw_error("TODO") # TODO - Error: One or more processes didnt execute correctly. Or warning showing the parameters.

    if(return_type == "success"):
      return successful_execution
    elif(return_type == "process_out"):
      return dump_process_output
    else:
      return None

  # UPPER MATRIX DICTIONARY -> BEDGRAPH
  def load(self, resolution, sparse_matrix_dictionary, output_file_name, start_index = 0):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Open output bedgraph file
    output_file = codecs.open(input_file_name, "w", "utf8")

    # Iterate over chromosome list
    for chrom in self.chromosome_sizes.chromosome_sizes_list:

      # Iterate over position 1
      for i in range(start_index, self.chromosome_sizes.chromosome_sizes_dictionary[chrom], resolution):

        # Iterate over position 2
        for j in range(i, self.chromosome_sizes.chromosome_sizes_dictionary[chrom], resolution):

          # General check
          if(((i + resolution - start_index) > self.chromosome_sizes.chromosome_sizes_dictionary[chrom]) or
             ((j + resolution - start_index) > self.chromosome_sizes.chromosome_sizes_dictionary[chrom])): 
            continue

          # Get count
          region = ":".join([chrom, str(i), str(j)])
          try:
            count = sparse_matrix_dictionary[region]
          except Exception: continue

          # Write entry to bedgraph file
          entry = [chrom, str(i), str(i + resolution), chrom, str(j), str(j + resolution), str(count)]
          output_file.write("\t".join(entry) + "\n")

    # Close temporary output file
    output_file.close()

    # Successful execution
    return True

  def add_load(self, resolution, sparse_matrix_dictionary, output_file_name, start_index = 0):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    self.process_queue.append((resolution, sparse_matrix_dictionary, output_file_name, start_index))

  def run_load(self, return_type = "success"):
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
      if(not cp):
        successful_execution = False
        self.error_handler.throw_error("TODO") # TODO - Error: One or more processes didnt execute correctly.

    if(return_type == "success"):
      return successful_execution
    elif(return_type == "process_out"):
      return load_process_output
    else:
      return None

  def identify_minimal_resolution(self, input_file_name):
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
    raw_resolution = 
    if(raw_resolution % 10 == 0):
      pass
    elif(raw_resolution % 10 in [1, 2, 3, 4]): 
      resolution = AuxiliaryFunctions.floor_multiple(raw_resolution, 10)
    else:
      resolution = AuxiliaryFunctions.ceil_multiple(raw_resolution, 10)

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

