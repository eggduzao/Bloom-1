from __future__ import print_function
"""
Barcode Test
===================
Placeholder.

Authors: Eduardo G. Gusmao.

"""

###################################################################################################
# Libraries
###################################################################################################

# Python
import os
import sys
import subprocess

# Internal
from bloom.contact_map import ContactMap
from bloom.barcode import Barcode
from bloom.io import InputFileType, IO

# External
import numpy as np


###################################################################################################
# BarcodeTest Class
###################################################################################################

class BarcodeTest():
  """Placeholder.

  *Keyword arguments:*

    - argument1 -- Short description. This argument represents a long description. It can be:
      - Possibility 1: A possibility 1.
      - Possibility 2: A possibility 2.

    - argument2 -- Short description. This argument represents a long description. It can be:
      - Possibility 1: A possibility 1.
      - Possibility 2: A possibility 2.
  """

  def __init__(self, input_location, temporary_location, output_location):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Class objects
    self.input_location = input_location
    self.temporary_location = temporary_location
    self.output_location = output_location

    # Creating folders
    for folder_name in [self.input_location, self.temporary_location, self.output_location]:
      if(not os.path.exists(folder_name)):
        folder_creation_command = ["mkdir", "-p", folder_name]
        folder_creation_process = subprocess.run(folder_creation_command , stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)

  def write_list(self, name, list_to_write, output_file):
    """Write list_to_write to output_file
    """
    output_file.write("### " + name + "\n")
    if(list_to_write):
      for i in range(0, len(list_to_write)):
        output_file.write(str(i) + ": " + str(list_to_write[i]) + "\n")
    else:
      output_file.write("None\n")
    output_file.write("\n")

  def write_dict(self, name, dict_to_write, output_file):
    """Write dict_to_write to output_file
    """
    output_file.write("### " + name + "\n")
    if(dict_to_write):
      keys = sorted(dict_to_write.keys())
      if(keys):
        for k in keys:
          output_file.write(str(k) + ": " + str(dict_to_write[k]) + "\n")
      else:
        output_file.write("None\n")
    else:
      output_file.write("None\n")
    output_file.write("\n")

  def write_contact_map(self, name, contact_map, output_file):
    """Write contact_map matrix to output_file
    """
    output_file.write("### " + name + "\n")
    if(contact_map and contact_map.matrix):
      keys = sorted(contact_map.matrix.keys())
      if(keys):
        for k in keys:
          regions = sorted(contact_map.matrix[k].keys(), key = lambda x: x[0])
          if(regions):
            for r in regions:
              output_file.write("\t".join([str(k), str(r[0]), str(r[1]), str(contact_map.matrix[k][r])]) + "\n")
          else:
            output_file.write(str(k) + ": " + "None\n")
      else:
        output_file.write("None\n")
    else:
      output_file.write("None\n")
    output_file.write("\n")

  def write_matrix(self, name, matrix, output_file):
    """Write dict_to_write to output_file
    """
    output_file.write("### " + name + "\n")
    if(matrix.size != 0):
      counter = 0
      for vec in matrix:
        output_file.write("[" + str(counter) + "]\t".join([str(e) for e in vec]) + "\n")
        counter += 1
    else:
      output_file.write("None\n")
    output_file.write("\n")

  def main_guide(self):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Input file
    input_file_name_bg = os.path.join(self.input_location, "bedgraph_matrix", "MESC_500000.bg") # OK
    input_file_name_hic = os.path.join(self.input_location, "juicer_matrix", "SCONT_500000.hic") # OK
    input_file_name_co = os.path.join(self.input_location, "cooler_matrix", "MESC_500000.cool") # Error because cooler trim values with chrom size > genome
    input_file_name_mco = os.path.join(self.input_location, "cooler_matrix", "E14R1.mcool") # Error because it is another dataset

    # Output file
    output_file_name_bg_map1 = os.path.join(self.output_location, "main_guide_bg_map1.txt")
    output_file_name_bg_map2 = os.path.join(self.output_location, "main_guide_bg_map2.txt")
    output_file_contact_bg = os.path.join(self.output_location, "main_contact_bg.txt")
    output_file_name_hic_map1 = os.path.join(self.output_location, "main_guide_hic_map1.txt")
    output_file_name_hic_map2 = os.path.join(self.output_location, "main_guide_hic_map2.txt")
    output_file_contact_hic = os.path.join(self.output_location, "main_contact_hic.txt")
    output_file_name_co_map1 = os.path.join(self.output_location, "main_guide_co_map1.txt")
    output_file_name_co_map2 = os.path.join(self.output_location, "main_guide_co_map2.txt")
    output_file_contact_co = os.path.join(self.output_location, "main_contact_co.txt")
    output_file_name_mco_map1 = os.path.join(self.output_location, "main_guide_mco_map1.txt")
    output_file_name_mco_map2 = os.path.join(self.output_location, "main_guide_mco_map2.txt")
    output_file_contact_mco = os.path.join(self.output_location, "main_contact_mco.txt")


    # Read bedgraph
    io = IO(input_file_name_bg, self.temporary_location, "mm9", 4, input_resolution = 500000, input_file_type = InputFileType.SPARSE)
    contact_map = io.read()
    contact_map.update_valid_chromosome_list()
    io.write(contact_map, output_file_name_bg_map1, InputFileType.SPARSE)

    # Barcode
    barcode_instance = Barcode(4, "mm9", contact_map, self.temporary_location)
    final_contact_map, final_loop_list = barcode_instance.main_guide()
    if(final_contact_map):
      final_contact_map.update_valid_chromosome_list()
      io.write(final_contact_map, output_file_name_bg_map2, InputFileType.SPARSE)
      io.write_loop_list(final_loop_list, output_file_contact_bg)

    # Read bedgraph
    io = IO(input_file_name_hic, self.temporary_location, "hg19", 4, input_resolution = 500000, input_file_type = InputFileType.HIC)
    contact_map = io.read()
    contact_map.update_valid_chromosome_list()
    io.write(contact_map, output_file_name_hic_map1, InputFileType.SPARSE)

    # Barcode
    barcode_instance = Barcode(4, "hg19", contact_map, self.temporary_location)
    final_contact_map, final_loop_list = barcode_instance.main_guide()
    if(final_contact_map):
      final_contact_map.update_valid_chromosome_list()
      io.write(final_contact_map, output_file_name_hic_map2, InputFileType.SPARSE)
      io.write_loop_list(final_loop_list, output_file_contact_hic)

    # Read bedgraph
    io = IO(input_file_name_co, self.temporary_location, "mm9", 4, input_resolution = 500000, input_file_type = InputFileType.COOL)
    contact_map = io.read()
    contact_map.update_valid_chromosome_list()
    io.write(contact_map, output_file_name_co_map1, InputFileType.SPARSE)

    # Barcode
    barcode_instance = Barcode(4, "mm9", contact_map, self.temporary_location)
    final_contact_map, final_loop_list = barcode_instance.main_guide()
    if(final_contact_map):
      final_contact_map.update_valid_chromosome_list()
      io.write(final_contact_map, output_file_name_co_map2, InputFileType.SPARSE)
      io.write_loop_list(final_loop_list, output_file_contact_co)

    # Read bedgraph
    io = IO(input_file_name_mco, self.temporary_location, "mm9", 4, input_resolution = 500000, input_file_type = InputFileType.MCOOL)
    contact_map = io.read()
    contact_map.update_valid_chromosome_list()
    io.write(contact_map, output_file_name_mco_map1, InputFileType.SPARSE)

    # Barcode
    barcode_instance = Barcode(4, "mm9", contact_map, self.temporary_location)
    final_contact_map, final_loop_list = barcode_instance.main_guide()
    if(final_contact_map):
      final_contact_map.update_valid_chromosome_list()
      io.write(final_contact_map, output_file_name_mco_map2, InputFileType.SPARSE)
      io.write_loop_list(final_loop_list, output_file_contact_mco)

###################################################################################################
# Execution
###################################################################################################

if __name__ == "__main__":

  # Fixed Parameters
  current_name = os.path.splitext(os.path.basename(__file__))[0]
  current_path = os.path.dirname(os.path.abspath(__file__))
  input_location = os.path.join(current_path, "input")
  temporary_location = os.path.join(current_path, "temp", current_name)
  output_location = os.path.join(current_path, "output", current_name)
  
  # Creating test
  test = BarcodeTest(input_location, temporary_location, output_location)

  # Tests
  test.main_guide()

