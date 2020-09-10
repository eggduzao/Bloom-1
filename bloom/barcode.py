"""
Barcode Module
===================
Placeholder.

Authors: Eduardo G. Gusmao.

"""

# DUMP FILES (OUTPUT OF THIS CLASS) ARE ALREADY ON BG2 FORMAT


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
from bloom.contact_map import ContactMap, InputFileType
from bloom.util import ChromosomeSizes, BarcodeFiles, ErrorHandler
from bloom.io_bedgraph import Bedgraph
from bloom.io_juicer import Juicer

# External

###################################################################################################
# Basic Objects
###################################################################################################

class Barcode():
  """This class represents TODO.

  *Keyword arguments:*

    - argument1 -- Short description. This argument represents a long description. It can be:
      - Possibility 1: A possibility 1.
      - Possibility 2: A possibility 2.

    - argument2 -- Short description. This argument represents a long description. It can be:
      - Possibility 1: A possibility 1.
      - Possibility 2: A possibility 2.
    """

  def __init__(self, ncpu, organism):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Initialization
    self.ncpu = ncpu
    self.organism = organism

    # Class objects
    self.cfile = None

    # Utilitary objects
    self.error_handler = ErrorHandler()
    self.barcode_handler = BarcodeFiles(self.organism)
    self.chromosome_sizes = ChromosomeSizes(self.organism)
    self.bedgraph_handler = Bedgraph(self.organism, self.ncpu)
    self.juicer_handler = Juicer(self.ncpu)

  def check_best_barcode_to_input(self, input_file_name, temporary_location, input_file_type = InputFileType.UNKNOWN): # This will assign a cfile matrix (output) to the input contact matrix, based on all bfiles matrices barcodes

    # Return status
    return_status = False
    bfile_name = None

    # Get matrix from input file
    try:
      input_matrix = ContactMap(input_file_name, temporary_location, self.organism, self.ncpu, input_resolution = 500000, input_file_type = input_file_type)
    except Exception:
      return return_status

    # Iterating on barcode dictionary
    for barcode_file_name in self.barcode_handler.barcode_file_list:

      # Dump barcode file
      bfile_name = os.path.join(temporary_location, "bfile_name.dump")
      self.binfile_to_dumpfile(barcode_file_name, bfile_name)

      # Convert dump file to matrix
      current_barcode = ContactMap(bfile_name, temporary_location, self.organism, self.ncpu, input_resolution = 500000, input_file_type = InputFileType.SPARSE)

      # Compare input_matrix with current_barcode with 1% error margin
      comp = input_matrix.compare_matrices(current_barcode)
      if(comp):

        # Dump cfile
        return_status = True
        current_cfile = self.barcode_handler.barcode_file_dictionary[barcode_file_name]
        self.cfile = os.path.join(temporary_location, "cfile_name.dump")
        self.binfile_to_dumpfile(current_cfile, self.cfile)
        
        # Stop search
        break

      elif(not comp):
        continue

    # Remove temporary file
    if(bfile_name):
      remove_command = ["rm", "-rf", bfile_name]
      remove_process = subprocess.run(remove_command , stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)

    # Return
    return return_status

  def binfile_to_dumpfile(self, bin_file_name, dump_file_name):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Iterate on each byte of the file
    string_array = []
    bin_file = codecs.open(bin_file_name, "rb", "utf8")
    dump_file = codecs.open(dump_file_name, "w", "utf8")
    byte = bin_file.read(1)
    while byte != b"":
      string = byte.decode("utf-8")
      string_array.append(string)
      if(string == "\n"):
        ss = "".join(string_array)
        dump_file.write(ss)
        string_array = []
      byte = bin_file.read(1)

    # Closing files
    bin_file.close()
    dump_file.close()

