from __future__ import print_function
"""
Cooler Test
===================
Placeholder.

Authors: Eduardo G. Gusmao.

"""

# cp /usr/users/egadegu/Projects/Papantonis_Bloom/Data/2_isHiC_Mouse_ES14/isHiC/E14R1.mcool /usr/users/egadegu/Projects/Papantonis_Bloom/Bloom/test/input/cooler_matrix

"""

cooler info /home/uni08/egadegu/Projects/Papantonis_Bloom/Bloom/test/output/test_5_io_cooler/load.cool
cooler ls /home/uni08/egadegu/Projects/Papantonis_Bloom/Bloom/test/input/cooler_matrix/MESC_500000.cool
cooler ls /home/uni08/egadegu/Projects/Papantonis_Bloom/Bloom/test/input/cooler_matrix/E14R1.mcool

cooler dump -t pixels --join -r chr1:0-100,000,000 -r2 chr1:0-100,000,000 /home/uni08/egadegu/Projects/Papantonis_Bloom/Bloom/test/output/test_5_io_cooler/load.cool > /home/uni08/egadegu/Projects/Papantonis_Bloom/Bloom/test/output/test_5_io_cooler/dump_single.txt

cooler dump -t pixels --join -r chr1:0,000,000-5,000,000 -r2 chr1:0,000,000-5,000,000 /home/uni08/egadegu/Projects/Papantonis_Bloom/Bloom/test/input/cooler_matrix/E14R1.mcool::resolutions/500000 > /home/uni08/egadegu/Projects/Papantonis_Bloom/Bloom/test/output/test_5_io_cooler/test.cool

["cooler", "dump", "-t", "pixels", "--join", "-r", "chr1:0-100,000,000", "-r2", "chr1:0-100,000,000", "/home/uni08/egadegu/Projects/Papantonis_Bloom/Bloom/test/input/cooler_matrix/MESC_500000.cool", ">", "/home/uni08/egadegu/Projects/Papantonis_Bloom/Bloom/test/output/test_5_io_cooler/dump_single.txt"

subprocess.run(command)

 stdout = "/home/uni08/egadegu/Projects/Papantonis_Bloom/Bloom/test/output/test_5_io_cooler/out.txt", stderr = "/home/uni08/egadegu/Projects/Papantonis_Bloom/Bloom/test/output/test_5_io_cooler/err.txt")

cooler dump -t pixels --join -r chr1:0-100,000,000 -r2 chr1:0-100,000,000 -o /home/uni08/egadegu/Projects/Papantonis_Bloom/Bloom/test/output/test_5_io_cooler/dump_single.txt /home/uni08/egadegu/Projects/Papantonis_Bloom/Bloom/test/input/cooler_matrix/MESC_500000.cool

cooler load -f bg2 --assembly mm9 --count-as-float /usr/users/egadegu/bloom_data/chrom_sizes/chrom.sizes.mm9:500000 /home/uni08/egadegu/Projects/Papantonis_Bloom/Bloom/test/temp/test_5_io_cooler/bedgraph_file_name.bg2 /home/uni08/egadegu/Projects/Papantonis_Bloom/Bloom/test/output/test_5_io_cooler/load.cool

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
from bloom.io_cooler import Cooler
from bloom.io_bedgraph import Bedgraph

# External
import numpy as np


###################################################################################################
# CoolerTest Class
###################################################################################################

class CoolerTest():
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


  def dump_single(self):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Input file
    input_file_name = os.path.join(self.input_location, "cooler_matrix", "MESC_500000.cool")

    # Output file
    output_file_name = os.path.join(self.output_location, "dump_single.txt")

    # Dump juicer / bedgraph
    cooler = Cooler("mm9", 4)
    region = "chr1:0-100,000,000"
    cooler.dump_single(region, region, input_file_name, output_file_name)


  def dump_multiple(self):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Input file
    input_file_name = os.path.join(self.input_location, "cooler_matrix", "E14R1.mcool")

    # Output file
    output_file_name = os.path.join(self.output_location, "dump_multiple.txt")

    # Dump juicer / bedgraph
    cooler = Cooler("mm9", 4)
    region = "chr1:0-100,000,000"
    cooler.dump_multiple("500000", region, region, input_file_name, output_file_name)


  def load(self):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Input file
    input_file_name = os.path.join(self.input_location, "bedgraph_matrix", "MESC_500000.bg")

    # Output file
    output_file_name = os.path.join(self.output_location, "load.cool")
    output_file_name_cm = os.path.join(self.output_location, "load.txt")
    output_file = open(output_file_name_cm, "w")

    # Contact map
    contact_map = ContactMap("mm9", 500000)
    contact_map.update_valid_chromosome_list()
    
    # Load contact map from bedgraph
    bedgraph = Bedgraph("mm9", 4)
    bedgraph.dump("chr1", input_file_name, contact_map, logit = False, pseudocount = 0.0)
    contact_map.update_valid_chromosome_list()
    self.write_contact_map("contact map", contact_map, output_file)

    # Load
    cooler = Cooler("mm9", 4)
    cooler.load(contact_map, self.temporary_location, output_file_name)

    # Termination
    output_file.close()


  def identify_minimal_resolution(self):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Input file
    input_file_name_1 = os.path.join(self.input_location, "cooler_matrix", "E14R1.mcool")
    input_file_name_2 = os.path.join(self.input_location, "cooler_matrix", "MESC_500000.cool")

    # Output file
    output_file_name = os.path.join(self.output_location, "identify_minimal_resolution.txt")
    output_file = open(output_file_name, "w")

    # Identify resolution
    cooler = Cooler("mm9", 4)
    resolution_1 = cooler.identify_minimal_resolution(input_file_name_1, self.temporary_location, check_type = "mcool", region = "chr1:0-5,000,000")
    resolution_2 = cooler.identify_minimal_resolution(input_file_name_2, self.temporary_location, check_type = "cool", region = "chr1:0-5,000,000")
    self.write_list("resolution", [resolution_1, resolution_2], output_file)

    # Termination
    output_file.close()


  def filetype_is_cooler(self):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Input file
    input_file_name_1 = os.path.join(self.input_location, "cooler_matrix", "E14R1.mcool")
    input_file_name_2 = os.path.join(self.input_location, "cooler_matrix", "MESC_500000.cool")
    input_file_name_bedgraph = os.path.join(self.input_location, "bedgraph_matrix", "MESC_500000.bg")

    # Output file
    output_file_name = os.path.join(self.output_location, "filetype_is_cooler.txt")
    output_file = open(output_file_name, "w")

    # Identify filetype
    cooler = Cooler("mm9", 4)
    is_cool_1 = cooler.filetype_is_cooler(input_file_name_1, self.temporary_location, check_type = "cool", region = "chr1:0-5,000,000")
    is_cool_2 = cooler.filetype_is_cooler(input_file_name_2, self.temporary_location, check_type = "cool", region = "chr1:0-5,000,000")
    is_cool_3 = cooler.filetype_is_cooler(input_file_name_bedgraph, self.temporary_location, check_type = "cool", region = "chr1:0-5,000,000")
    is_mcool_1 = cooler.filetype_is_cooler(input_file_name_1, self.temporary_location, check_type = "mcool", region = "chr1:0-5,000,000")
    is_mcool_2 = cooler.filetype_is_cooler(input_file_name_2, self.temporary_location, check_type = "mcool", region = "chr1:0-5,000,000")
    is_mcool_3 = cooler.filetype_is_cooler(input_file_name_bedgraph, self.temporary_location, check_type = "mcool", region = "chr1:0-5,000,000")
    self.write_list("file type is cool", [is_cool_1, is_cool_2, is_cool_3, is_mcool_1, is_mcool_2, is_mcool_3], output_file)

    # Termination
    output_file.close()


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
  test = CoolerTest(input_location, temporary_location, output_location)

  # Tests
  test.dump_single()
  test.dump_multiple()
  test.load()
  test.identify_minimal_resolution()
  test.filetype_is_cooler()

