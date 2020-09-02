"""
Contact Map Module
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
import optparse
import traceback
import subprocess
import configparser
import multiprocessing

# Internal
from bloom.util import ErrorHandler, ChromosomeSizes, AuxiliaryFunctions
from bloom.io_bedgraph import Bedgraph
from bloom.io_juicer import Juicer
from bloom.io_cooler import Cooler

# External
import numpy


###################################################################################################
# Basic Objects
###################################################################################################

class InputFileType():
  """This class represents TODO.

  *Keyword arguments:*

    - argument1 -- Short description. This argument represents a long description. It can be:
      - Possibility 1: A possibility 1.
      - Possibility 2: A possibility 2.

    - argument2 -- Short description. This argument represents a long description. It can be:
      - Possibility 1: A possibility 1.
      - Possibility 2: A possibility 2.
  """

  UNKNOWN = 0
  SPARSE = 1
  HIC = 2
  COOL = 3
  MCOOL = 4

class ContactMap():
  """This class represents TODO.

  *Keyword arguments:*

    - argument1 -- Short description. This argument represents a long description. It can be:
      - Possibility 1: A possibility 1.
      - Possibility 2: A possibility 2.

    - argument2 -- Short description. This argument represents a long description. It can be:
      - Possibility 1: A possibility 1.
      - Possibility 2: A possibility 2.
  """

  def __init__(self, input_contact_matrix_file_name, temporary_location, input_resolution = None, input_file_type = InputFileType.UNKNOWN):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Main objects
    self.input_file_name = input_contact_matrix_file_name
    self.temporary_location = temporary_location
    self.input_resolution = input_resolution
    self.input_file_type = input_file_type
    self.matrix = dict()
    self.resolution = None
    
    # Auxiliary objects
    self.total_zero = 0
    self.total_nonzero = 0
    self.total_sum_nonzero = 0.0
    self.max = -numpy.inf
    self.min = numpy.inf

    # Utilitary objects
    self.error_handler = ErrorHandler()
    self.chromosome_sizes = ChromosomeSizes()
    self.bedgraph_handler = Bedgraph()
    self.juicer_handler = Juicer()
    self.cooler_handler = Cooler()

    # Load file
    self.load_matrix()

  #############################################################################
  # Input File Loading
  #############################################################################

  def load_matrix(self):

    # Verify if input file exists
    self.verify_input_file()

    # Verify input file type
    self.detect_file_type()

    # Verify input resolution(s)
    self.detect_input_resolutions()
  
    # Load matrix based on file type and resolution
      # Check if file is Juicer (.hic)
      if(self.input_file_type == InputFileType.HIC):
        self.load_matrix_from_hic()

      # Check if file is singular cooler (.cool)
      elif(self.input_file_type == InputFileType.COOL):
        self.load_matrix_from_cool()

      # Check if file is multiple cooler (.mcool)
      elif(self.input_file_type == InputFileType.MCOOL):
        self.load_matrix_from_mcool()

      # Check if file is sparse text bedgraph (.bg2, .bed, .txt)
      elif(self.input_file_type == InputFileType.SPARSE):
        self.load_matrix_from_sparse()

  def verify_input_file(self):

    # Verify if input file exists
    if(not os.path.isfile(self.input_file_name)):
      self.error_handler.throw_error("TODO") # TODO

  def detect_file_type(self):

    # Detect file type if file is unknown
    if(self.input_file_type == InputFileType.UNKNOWN):

      # Check if file is Juicer (.hic)
      if(self.juicer_handler.filetype_is_juicer(self.input_file_name, self.temporary_location)):
        self.input_file_type = InputFileType.HIC

      # Check if file is singular cooler (.cool)
      elif(self.cooler_handler.filetype_is_cooler(self.input_file_name, self.temporary_location, check_type = "cool")):
        self.input_file_type = InputFileType.COOL

      # Check if file is multiple cooler (.mcool)
      elif(self.cooler_handler.filetype_is_cooler(self.input_file_name, self.temporary_location, check_type = "mcool")):
        self.input_file_type = InputFileType.MCOOL

      # Check if file is sparse text bedgraph (.bg2, .bed, .txt)
      elif(self.bedgraph_handler.filetype_is_bedgraph(self.input_file_name)):
        self.input_file_type = InputFileType.SPARSE

      # No recognizable file type
      else:
        self.error_handler.throw_error("TODO") # TODO

  def detect_input_resolutions(self):

    # Detect file type if file is unknown
    if(self.input_resolution == None):

      # Check if file is Juicer (.hic)
      if(self.input_file_type == InputFileType.HIC):
        self.input_resolution = self.juicer_handler.identify_minimal_resolution(self.input_file_name, self.temporary_location)

      # Check if file is singular cooler (.cool)
      elif(self.input_file_type == InputFileType.COOL):
        self.input_resolution = self.cooler_handler.identify_minimal_resolution(self.input_file_name, self.temporary_location, check_type = "cool")

      # Check if file is multiple cooler (.mcool)
      elif(self.input_file_type == InputFileType.MCOOL):
        self.input_resolution = self.cooler_handler.identify_minimal_resolution(self.input_file_name, self.temporary_location, check_type = "mcool")

      # Check if file is sparse text bedgraph (.bg2, .bed, .txt)
      elif(self.input_file_type == InputFileType.SPARSE):
        self.input_resolution = self.bedgraph_handler.identify_minimal_resolution(self.input_file_name)

      # No recognizable file type
      else:
        self.error_handler.throw_error("TODO") # TODO

    # If resolution continues to be unknown it was not detectable
    if(self.input_resolution == None):
      self.error_handler.throw_error("TODO") # TODO

  def load_matrix_from_hic(self):

    # Iterating on chromosomes
    for chrom in self.chromosome_sizes.chromosome_sizes_list: # TODO

self.juicer_handler.add_dump(resolution, region1, region2, input_file_name, output_file_name)

    # Running all jobs
    dump_process_output = self.juicer_handler.run_dump(return_type = "process_out")

    # Verification of dumping processes
    self.load_process_verification(dump_process_output)

    # Bedgraph

  def load_matrix_from_cool(self):

    # Iterating on chromosomes
    for chrom in self.chromosome_sizes.chromosome_sizes_list:

self.cooler_handler.add_dump_single(resolution, region1, region2, input_file_name, output_file_name)

    # Running all jobs
    dump_process_output = self.cooler_handler.run_dump_single(return_type = "process_out")

    # Verification of dumping processes
    self.load_process_verification(dump_process_output)

    # Bedgraph

  def load_matrix_from_mcool(self):

    # Iterating on chromosomes
    for chrom in self.chromosome_sizes.chromosome_sizes_list:

self.cooler_handler.add_dump_multiple(resolution, region1, region2, input_file_name, output_file_name)


    # Running all jobs
    dump_process_output = self.cooler_handler.run_dump_multiple(return_type = "process_out")

    # Verification of dumping processes
    self.load_process_verification(dump_process_output)

    # Bedgraph

  def load_matrix_from_sparse(self):

    # Iterating on chromosomes
    for chrom in self.chromosome_sizes.chromosome_sizes_list:

      # Adding chromosome dump job
      self.bedgraph_handler.add_dump(chrom, self.input_file_name, self.matrix)

    # Running all jobs
    dump_process_output = self.bedgraph_handler.run_dump(return_type = "process_out")

    # Verification of loading processes
    self.load_process_verification(dump_process_output)

  def load_process_verification(self, dump_process_output):

    # Check if all chromosomes executed successfully
    for k in range(0, len(dump_process_output)):

      # Get specific process and chromosome
      cp = dump_process_output[k]
      chrom = self.chromosome_sizes.chromosome_sizes_list[k]

      # Verify if chromosome executed correctly
      try:
        returncode = cp.check_returncode()
        if(returncode > 0): 
          self.error_handler.throw_warning("TODO", chrom) # TODO - the process had an error, and exited with that code
        elif(returncode < 0): 
          self.error_handler.throw_warning("TODO", chrom) # TODO - the process was killed with a signal of -1 * exitcode
      except subprocess.CalledProcessError:
        self.error_handler.throw_error("TODO", chrom) # TODO

  #############################################################################
  # Output File Writing
  #############################################################################

  def write_matrix_as_sparse(self):
    pass

  def write_matrix_as_hic(self):
    pass

  def write_matrix_as_mcool(self):
    pass

  def write_matrix_as_cool(self):
    pass

  #############################################################################
  # Binary Operations
  #############################################################################

  def create_bin_file_from_text_file(self, text_file_name, bin_file_name):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Write bin file line by line
    text_file = open(text_file_name, "r")
    bin_file = open(bin_file_name, "wb")
    for line in text_file: bin_file.write(bytearray(line, "utf-8"))
    text_file.close()
    bin_file.close()

  def create_dictionary_from_bin_file(self, bin_file_name):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Dictionary
    res_dict = dict()

    # Iterate on each byte of the file
    string_array = []
    bin_file = open(bin_file_name, "rb")
    byte = bin_file.read(1)
    while byte != b"":
      string = byte.decode("utf-8")
      if(string == "\n"):
        ss = "".os.path.join(string_array).split("\t")
        res_dict[":".os.path.join(ss[:3])] = float(ss[3])
        string_array = []
      else: string_array.append(string)
      byte = bin_file.read(1)

    # Returning objects
    return res_dict


  #############################################################################
  # Auxiliary Operations
  #############################################################################









