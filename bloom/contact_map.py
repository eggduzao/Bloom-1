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

  def __init__(self, input_contact_matrix_file_name, temporary_location, organism, ncpu, input_resolution = None, input_file_type = InputFileType.UNKNOWN):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Main objects
    self.input_file_name = input_contact_matrix_file_name
    self.temporary_location = temporary_location
    self.organism = organism
    self.ncpu = ncpu
    self.input_resolution = input_resolution
    self.input_file_type = input_file_type
    self.matrix = dict()
    self.resolution = None
    
    # Auxiliary objects
    self.total_bins = dict() # per chromosome
    self.total_1d_bins = dict() # per chromosome
    self.total_zero_bins = dict() # per chromosome
    self.total_nonzero_bins = dict() # per chromosome
    self.total_nonzero_value = dict() # per chromosome
    self.max_value = dict() # per chromosome
    self.min_value = dict() # per chromosome

    # Utilitary objects
    self.error_handler = ErrorHandler()
    self.chromosome_sizes = ChromosomeSizes(self.organism)
    self.bedgraph_handler = Bedgraph(self.organism, self.ncpu)
    self.juicer_handler = Juicer(self.ncpu)
    self.cooler_handler = Cooler(self.organism, self.ncpu)

    # Load file
    self.load_matrix()

  #############################################################################
  # Input File Loading
  #############################################################################

  def load_matrix(self):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

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
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Verify if input file exists
    if(not os.path.isfile(self.input_file_name)):
      self.error_handler.throw_error("TODO") # TODO

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
        self.error_handler.throw_error("TODO") # TODO

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
        self.error_handler.throw_error("TODO") # TODO

    # If resolution continues to be unknown it was not detectable
    if(self.input_resolution == None):
      self.error_handler.throw_error("TODO") # TODO

  def load_matrix_from_hic(self):
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
      chrom_wo_chr = chrom.split("chr")[0]
      start = "1"
      end = str(self.chromosome_sizes.chromosome_sizes_dictionary[chrom])
      region = ":".join([chrom_wo_chr, start, end])

      # Temporary bedgraph file
      bedgraph_file_name = os.path.join(self.temporary_location, "bedgraph_file_name_" + chrom + ".bg2")
      list_files_to_remove.append(bedgraph_file_name)

      # Adding juicer dump job
      self.juicer_handler.add_dump(self.resolution, region, region, self.input_file_name, self.temporary_location, bedgraph_file_name, output_type = "bedgraph")

      # Adding bedgraph dump job
      self.bedgraph_handler.add_dump(chrom, bedgraph_file_name, self.matrix)

    # Running juicer jobs
    dump_process_output = self.juicer_handler.run_dump(return_type = "process_out")

    # Verification of juicer dumping processes
    self.load_process_verification(dump_process_output)

    # Running bedgraph jobs
    dump_process_output = self.bedgraph_handler.run_dump(return_type = "process_out")

    # Verification of bedgraph loading processes
    self.load_process_verification(dump_process_output)

    # Removing temporary files
    remove_command = ["rm", "-rf"] + list_files_to_remove
    remove_process = subprocess.run(remove_command , stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)

  def load_matrix_from_cool(self):
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
      end = '{:,}'.format(self.chromosome_sizes.chromosome_sizes_dictionary[chrom])
      region = chrom + ":" + start + "-" + end

      # Temporary bedgraph file
      bedgraph_file_name = os.path.join(self.temporary_location, "bedgraph_file_name_" + chrom + ".bg2")
      list_files_to_remove.append(bedgraph_file_name)

      # Adding chromosome dump job
      self.cooler_handler.add_dump_single(region, region, self.input_file_name, bedgraph_file_name)

      # Adding bedgraph dump job
      self.bedgraph_handler.add_dump(chrom, bedgraph_file_name, self.matrix)

    # Running cooler jobs
    dump_process_output = self.cooler_handler.run_dump_single(return_type = "process_out")

    # Verification of cooler dumping processes
    self.load_process_verification(dump_process_output)

    # Running bedgraph jobs
    dump_process_output = self.bedgraph_handler.run_dump(return_type = "process_out")

    # Verification of bedgraph loading processes
    self.load_process_verification(dump_process_output)

    # Removing temporary files
    remove_command = ["rm", "-rf"] + list_files_to_remove
    remove_process = subprocess.run(remove_command , stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)

  def load_matrix_from_mcool(self):
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
      end = '{:,}'.format(self.chromosome_sizes.chromosome_sizes_dictionary[chrom])
      region = chrom + ":" + start + "-" + end

      # Temporary bedgraph file
      bedgraph_file_name = os.path.join(self.temporary_location, "bedgraph_file_name_" + chrom + ".bg2")
      list_files_to_remove.append(bedgraph_file_name)

      # Adding chromosome dump job
      self.cooler_handler.add_dump_multiple(self.resolution, region, region, self.input_file_name, bedgraph_file_name)

      # Adding bedgraph dump job
      self.bedgraph_handler.add_dump(chrom, bedgraph_file_name, self.matrix)

    # Running cooler jobs
    dump_process_output = self.cooler_handler.run_dump_multiple(return_type = "process_out")

    # Verification of cooler dumping processes
    self.load_process_verification(dump_process_output)

    # Running bedgraph jobs
    dump_process_output = self.bedgraph_handler.run_dump(return_type = "process_out")

    # Verification of bedgraph loading processes
    self.load_process_verification(dump_process_output)

    # Removing temporary files
    remove_command = ["rm", "-rf"] + list_files_to_remove
    remove_process = subprocess.run(remove_command , stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)

  def load_matrix_from_sparse(self, input_file_name):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Iterating on chromosomes
    for chrom in self.chromosome_sizes.chromosome_sizes_list:

      # Adding chromosome dump job
      self.bedgraph_handler.add_dump(chrom, self.input_file_name, self.matrix)

    # Running all jobs
    dump_process_output = self.bedgraph_handler.run_dump(return_type = "process_out")

    # Verification of loading processes
    self.load_process_verification(dump_process_output)

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

  def write_matrix_as_hic(self, output_file_name):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Adding load job
    self.juicer_handler.add_load(self.organism, self.resolution, self.matrix, self.temporary_location, output_file_name)

    # Running load job
    self.juicer_handler.run_load(return_type = "success")

  def write_matrix_as_cool(self, output_file_name):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Adding load job
    self.cooler_handler.add_load(self.organism, self.resolution, self.matrix, self.temporary_location, output_file_name)

    # Running load job
    self.cooler_handler.run_load(return_type = "success")

  def write_matrix_as_sparse(self, output_file_name):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Adding load job
    self.bedgraph_handler.add_load(self.resolution, self.matrix, output_file_name, start_index = 0)

    # Running load job
    self.bedgraph_handler.run_load(return_type = "success")


  #############################################################################
  # Operational Operations
  #############################################################################

  def set(self, chrom i, j, value):

    # Set matrix value
    key = ":".join([chrom, str(i), str(j)])
    self.matrix[key] = value

  def add(self, chrom, i, j, value):

    # Set matrix value
    key = ":".join([chrom, str(i), str(j)])
    try:
      self.matrix[key] += value
    except Exception:
      self.matrix[key] = value

  def get(self, chrom, i, j):

    # Return value
    key = ":".join([chrom, str(i), str(j)])
    try:
      return self.matrix[key]
    except Exception:
      return 0


  #############################################################################
  # Sparsity Operations
  #############################################################################

  def total_upper_matrix(self):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Total bins = Summation of the total number of bins of each chromosome map upper triangular matrix without the diagonal bins
    self.total_bins = 0

    # Get valid chromosome list
    valid_chromosome_dict = dict()
    for key, value in self.matrix.iteritems():
      chrom = key.split(":")[0]
      valid_chromosome_dict[chrom] = True
    valid_chromosome_list = sorted(valid_chromosome_dict.keys())

    # Iterate on valid chromosome list
    for chrom in valid_chromosome_list:

      if(not valid_chromosome_list[chrom]): continue

      chrom_size = self.chromosome_sizes.chromosome_sizes_dictionary[chrom]
      total_1d_bins = AuxiliaryFunctions.floor_multiple(chrom_size, self.resolution) / self.resolution
      total_bins = ((total_1d_bins * (total_1d_bins-1))/2) + total_1d_bins

      try:
        self.total_1d_bins[chrom] = total_1d_bins
      except Exception:
        self.error_handler.throw_error("TODO") # TODO - Error: One or more processes didnt execute correctly.

      try:
        self.total_bins[chrom] = total_bins
      except Exception:
        self.error_handler.throw_error("TODO") # TODO - Error: One or more processes didnt execute correctly.

  def update_sparsity(self):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Iterate over matrix
    for key, value in self.matrix.iteritems():
  
      # Get chromosome
      chrom = key.split(":")[0]
      
      # Update values
      try:
        self.total_nonzero_bins[chrom] += 1
      except Exception:
        self.total_nonzero_bins[chrom] = 1

      try:
        self.total_nonzero_value[chrom] += value
      except Exception:
        self.total_nonzero_value[chrom] = value

      try:
        if(value > self.max_value[chrom]):
          self.max_value[chrom] = value
      except Exception:
        self.max_value[chrom] = value

      try:
        if(value < self.min_value[chrom]):
          self.min_value[chrom] = value
      except Exception:
        self.min_value[chrom] = value
   
    # Update total bins = 0
    self.total_zero_bins = self.total_bins - self.total_nonzero_bins

  def get_sparsity(self, chromosome):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Return sparsity level
    return self.total_nonzero_bins[chromosome] / self.total_bins[chromosome]

  def get_sparseity_weighted_sum(self, chromosome):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Return sparsity weighted sum
    return self.total_nonzero_value[chromosome] * (self.total_nonzero_bins[chromosome] / self.total_bins[chromosome])


  #############################################################################
  # Binary Operations
  #############################################################################




  #############################################################################
  # Auxiliary Operations
  #############################################################################

  def compare_matrices(self, matrix):

    # Get input keys and compare
    mykeys = sorted(self.matrix.keys())
    cpkeys = sorted(matrix)
    if(mykeys != cpkeys):
      return False

    # Compare each element
    for k in mykeys:
      myvalue = self.matrix[k]
      cpvalue = matrix[k]
      if((myvalue > cpvalue + (0.1 * cpvalue)) or (myvalue < cpvalue - (0.1 * cpvalue))):
        return False

    # Return
    return True

  def distance_from_diagonal(self, key):

    # Fetching rows and columns
    kk = key.split(":")
    i_bin = float(int(kk[1]) / self.resolution)
    j_bin = float(int(kk[2]) / self.resolution)
    
    # Calculating standardized distance from diagonal
    std_distance = ((j_bin - i_bin) / 2. ) / (self.total_1d_bins[chrom] / 2)

    # Returning distance
    return std_distance












