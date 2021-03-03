from __future__ import print_function
"""
IO Module
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
from bloom.util import ChromosomeSizes, ErrorHandler, AuxiliaryFunctions
from bloom.io_bedgraph import Bedgraph
from bloom.io_juicer import Juicer
from bloom.io_cooler import Cooler

# External
import numpy


###################################################################################################
# InputFileType Class
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


###################################################################################################
# IO Class
###################################################################################################

class IO():
  """This class represents TODO.

  *Keyword arguments:*

    - argument1 -- Short description. This argument represents a long description. It can be:
      - Possibility 1: A possibility 1.
      - Possibility 2: A possibility 2.

    - argument2 -- Short description. This argument represents a long description. It can be:
      - Possibility 1: A possibility 1.
      - Possibility 2: A possibility 2.
  """

  def __init__(self, input_contact_map_file_name, temporary_location, organism, ncpu, input_resolution = None, input_file_type = InputFileType.UNKNOWN):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Main objects
    self.input_file_name = input_contact_map_file_name
    self.temporary_location = temporary_location
    self.organism = organism
    self.ncpu = ncpu
    self.input_resolution = input_resolution
    self.input_file_type = input_file_type

    # Utilitary objects
    self.error_handler = ErrorHandler()
    self.chromosome_sizes = ChromosomeSizes(self.organism)
    self.bedgraph_handler = Bedgraph(self.organism, self.ncpu)
    self.juicer_handler = Juicer(self.ncpu)
    self.cooler_handler = Cooler(self.organism, self.ncpu)


  #############################################################################
  # Read: File -> Contact Map (Class)
  #############################################################################

  def read(self):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - contact_map -- A contact map containing all the values from this class' input file.
    """

    # Verify if input file exists and check file type
    self.verify_input_file()

    # Verify input file type
    if(not self.input_file_type or self.input_file_type == InputFileType.UNKNOWN):
      self.detect_file_type()

    # Verify input resolution(s)
    if(not self.input_resolution):
      self.detect_input_resolutions()

    # Create new contact map
    contact_map = ContactMap(self.organism, self.input_resolution)

    # Load contact map based on file type and resolution

    # Check if file is Juicer (.hic)
    if(self.input_file_type == InputFileType.HIC):
      self.load_contact_map_from_hic(contact_map)

    # Check if file is singular cooler (.cool)
    elif(self.input_file_type == InputFileType.COOL):
      self.load_contact_map_from_cool(contact_map)

    # Check if file is multiple cooler (.mcool)
    elif(self.input_file_type == InputFileType.MCOOL):
      self.load_contact_map_from_mcool(contact_map)

    # Check if file is sparse text bedgraph (.bg2, .bed, .txt)
    elif(self.input_file_type == InputFileType.SPARSE):
      self.load_contact_map_from_sparse(contact_map)

    # Input file type not recognized
    else:
      pass
      # self.error_handler.throw_error("TODO") # TODO

    # Return new contact map
    return contact_map

  def verify_input_file(self):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Verify if input file exists
    if(not os.path.isfile(self.input_file_name)):
      pass
      # self.error_handler.throw_error("TODO") # TODO

    # Possible file types
    possible_file_type_list = [InputFileType.UNKNOWN, InputFileType.HIC, InputFileType.COOL, InputFileType.MCOOL, InputFileType.SPARSE]

    # Check if file is Juicer (.hic)
    if(self.input_file_type not in possible_file_type_list):
      self.input_file_type = InputFileType.UNKNOWN

  def detect_file_type(self):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

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
        pass
        # self.error_handler.throw_error("TODO") # TODO

    # If file type continues to be unknown it was not detectable
    if(self.input_file_type == InputFileType.UNKNOWN):
      pass
      # self.error_handler.throw_error("TODO") # TODO

  def detect_input_resolutions(self):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

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
        pass
        # self.error_handler.throw_error("TODO") # TODO

    # If resolution continues to be unknown it was not detectable
    if(self.input_resolution == None):
      pass
      # self.error_handler.throw_error("TODO") # TODO

  def load_contact_map_from_hic(self, contact_map):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # List of temporary files to remove
    list_files_to_remove = []

    # Iterating on chromosomes
    for chrom in self.chromosome_sizes.chromosome_sizes_list:

      # Regions
      chrom_wo_chr = chrom.split("chr")[-1]
      start = "1"
      end = str(contact_map.total_1d_bp[chrom])
      region = ":".join([chrom_wo_chr, start, end])

      # Temporary bedgraph file
      bedgraph_file_name = os.path.join(self.temporary_location, "bedgraph_file_name_" + chrom + ".bg2")
      list_files_to_remove.append(bedgraph_file_name)

      # Adding juicer dump job
      self.juicer_handler.dump(self.input_resolution, region, region, self.input_file_name, self.temporary_location, bedgraph_file_name, output_type = "bedgraph")
      # self.juicer_handler.add_dump(self.input_resolution, region, region, self.input_file_name, self.temporary_location, bedgraph_file_name, output_type = "bedgraph")

      # Adding bedgraph dump job
      self.bedgraph_handler.dump(chrom, bedgraph_file_name, contact_map)
      # self.bedgraph_handler.add_dump(chrom, bedgraph_file_name, contact_map)

    # Running juicer jobs
    # dump_process_output = self.juicer_handler.run_dump(return_type = "process_out")

    # Verification of juicer dumping processes
    #self.load_process_verification(dump_process_output)

    # Running bedgraph jobs
    # dump_process_output = self.bedgraph_handler.run_dump(return_type = "process_out")

    # Verification of bedgraph loading processes
    #self.load_process_verification(dump_process_output)

    # Removing temporary files
    remove_command = ["rm", "-rf"] + list_files_to_remove
    remove_process = subprocess.run(remove_command, stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)

  def load_contact_map_from_cool(self, contact_map):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # List of temporary files to remove
    list_files_to_remove = []

    # Iterating on chromosomes
    for chrom in self.chromosome_sizes.chromosome_sizes_list:

      # Regions
      start = "1"
      end = '{:,}'.format(contact_map.total_1d_bp[chrom])
      region = chrom + ":" + start + "-" + end

      # Temporary bedgraph file
      bedgraph_file_name = os.path.join(self.temporary_location, "bedgraph_file_name_" + chrom + ".bg2")
      list_files_to_remove.append(bedgraph_file_name)

      # Adding chromosome dump job
      self.cooler_handler.dump_single(region, region, self.input_file_name, bedgraph_file_name)
      # self.cooler_handler.add_dump_single(region, region, self.input_file_name, bedgraph_file_name)

      # Adding bedgraph dump job
      self.bedgraph_handler.dump(chrom, bedgraph_file_name, contact_map)
      # self.bedgraph_handler.add_dump(chrom, bedgraph_file_name, contact_map)

    # Running cooler jobs
    # dump_process_output = self.cooler_handler.run_dump_single(return_type = "process_out")

    # Verification of cooler dumping processes
    #self.load_process_verification(dump_process_output)

    # Running bedgraph jobs
    # dump_process_output = self.bedgraph_handler.run_dump(return_type = "process_out")

    # Verification of bedgraph loading processes
    #self.load_process_verification(dump_process_output)

    # Removing temporary files
    remove_command = ["rm", "-rf"] + list_files_to_remove
    remove_process = subprocess.run(remove_command , stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)

  def load_contact_map_from_mcool(self, contact_map):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # List of temporary files to remove
    list_files_to_remove = []

    # Iterating on chromosomes
    for chrom in self.chromosome_sizes.chromosome_sizes_list:

      # Regions
      start = "1"
      end = '{:,}'.format(contact_map.total_1d_bp[chrom])
      region = chrom + ":" + start + "-" + end

      # Temporary bedgraph file
      bedgraph_file_name = os.path.join(self.temporary_location, "bedgraph_file_name_" + chrom + ".bg2")
      list_files_to_remove.append(bedgraph_file_name)

      # Adding chromosome dump job
      self.cooler_handler.dump_multiple(self.input_resolution, region, region, self.input_file_name, bedgraph_file_name)
      # self.cooler_handler.add_dump_multiple(self.input_resolution, region, region, self.input_file_name, bedgraph_file_name)

      # Adding bedgraph dump job
      self.bedgraph_handler.dump(chrom, bedgraph_file_name, contact_map)
      # self.bedgraph_handler.add_dump(chrom, bedgraph_file_name, contact_map)

    # Running cooler jobs
    # dump_process_output = self.cooler_handler.run_dump_multiple(return_type = "process_out")

    # Verification of cooler dumping processes
    #self.load_process_verification(dump_process_output)

    # Running bedgraph jobs
    # dump_process_output = self.bedgraph_handler.run_dump(return_type = "process_out")

    # Verification of bedgraph loading processes
    #self.load_process_verification(dump_process_output)

    # Removing temporary files
    remove_command = ["rm", "-rf"] + list_files_to_remove
    remove_process = subprocess.run(remove_command , stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)

  def load_contact_map_from_sparse(self, contact_map):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Iterating on chromosomes
    for chrom in self.chromosome_sizes.chromosome_sizes_list:

      # Adding chromosome dump job
      self.bedgraph_handler.dump(chrom, self.input_file_name, contact_map)
      # self.bedgraph_handler.add_dump(chrom, self.input_file_name, contact_map)

    # Running all jobs
    # dump_process_output = self.bedgraph_handler.run_dump(return_type = "process_out")

    # Verification of loading processes
    #self.load_process_verification(dump_process_output)

  def load_process_verification(self, dump_process_output):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Check if all chromosomes executed successfully
    for k in range(0, len(dump_process_output)):

      # Get specific process and chromosome
      cp = dump_process_output[k]
      chrom = self.chromosome_sizes.chromosome_sizes_list[k]

      # Verify if chromosome executed correctly
      try:
        returncode = cp.check_returncode() # TODO - AttributeError: 'bool' object has no attribute 'check_returncode'
        if(returncode > 0):
          pass
          # self.error_handler.throw_warning("TODO", chrom) # TODO - the process had an error, and exited with that code
        elif(returncode < 0): 
          pass
          # self.error_handler.throw_warning("TODO", chrom) # TODO - the process was killed with a signal of -1 * exitcode
      except subprocess.CalledProcessError:
        raise
        # self.error_handler.throw_error("TODO") # TODO


  #############################################################################
  # Write: Contact Matrix -> File
  #############################################################################

  def write(self, contact_map, output_file_name, output_format = InputFileType.SPARSE):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Create output location if it does not exist
    output_location = os.path.dirname(output_file_name)
    if(not os.path.isdir(output_location)):
      try:
        os.mkdir(output_location)
      except OSError:
        raise
        # self.error_handler.throw_error("TODO") # TODO

    # Write file as Juicer (.hic)
    if(output_format == InputFileType.HIC):
      self.write_contact_map_as_hic(contact_map, output_file_name)

    # Write file as singular cooler (.cool)
    elif(output_format == InputFileType.COOL):
      self.write_contact_map_as_cool(contact_map, output_file_name)

    # Write file as sparse text bedgraph (.bg2, .bed, .txt)
    elif(output_format == InputFileType.SPARSE):
      self.write_contact_map_as_sparse(contact_map, output_file_name)

    # No recognizable output format
    else:
      pass
      # self.error_handler.throw_error("TODO") # TODO

  def write_contact_map_as_hic(self, contact_map, output_file_name):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Adding load job
    self.juicer_handler.load(contact_map, self.temporary_location, output_file_name)
    # self.juicer_handler.add_load(contact_map, self.temporary_location, output_file_name)

    # Running load job
    # self.juicer_handler.run_load(return_type = "success")

  def write_contact_map_as_cool(self, contact_map, output_file_name):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Adding load job
    self.cooler_handler.load(contact_map, self.temporary_location, output_file_name)
    # self.cooler_handler.add_load(contact_map, self.temporary_location, output_file_name)

    # Running load job
    # self.cooler_handler.run_load(return_type = "success")

  def write_contact_map_as_sparse(self, contact_map, output_file_name):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Adding load job
    self.bedgraph_handler.load(contact_map, output_file_name, start_index = 0)
    # self.bedgraph_handler.add_load(contact_map, output_file_name, start_index = 0)

    # Running load job
    # self.bedgraph_handler.run_load(return_type = "success")


  def write_loop_list(self, loop_list, output_file_name):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Write loop list
    output_file = open(output_file_name, "w")
    for loop in loop_list:
      output_file.write(loop)
    output_file.close()

