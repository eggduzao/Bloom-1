from __future__ import print_function
"""
IO Test
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
from bloom.io import InputFileType, IO

# External
import numpy as np

"""

cooler dump -t pixels --join -r chr1:1-197,000,000 -r2 chr1:1-197,000,000 -o /home/uni08/egadegu/Projects/Papantonis_Bloom/Bloom/test/temp/test_6_io/bedgraph_file_name_chr1.bg2 /home/uni08/egadegu/Projects/Papantonis_Bloom/Bloom/test/input/cooler_matrix/MESC_500000.cool

pip3.9 install --user cooler

java -Djava.awt.headless=true -Xmx32000m -jar /usr/users/egadegu/bloom_data/bin/juicer_tools_1.22.01.jar dump observed NONE /usr/users/egadegu/Projects/Papantonis_Bloom/Bloom/test/output/test_6_io/write_hic.hic chr1:0:249000000 chr1:0:249000000 BP 500000 /usr/users/egadegu/Projects/Papantonis_Bloom/Bloom/test/output/test_6_io/write_hic.txt

cooler info /usr/users/egadegu/Projects/Papantonis_Bloom/Bloom/test/output/test_6_io/write_co.cool

"""

###################################################################################################
# IOTest Class
###################################################################################################

class IOTest():
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


  def read(self):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Input file
    input_file_name_bg = os.path.join(self.input_location, "bedgraph_matrix", "MESC_500000.bg")
    input_file_name_hic = os.path.join(self.input_location, "juicer_matrix", "SCONT_500000.hic")
    input_file_name_co = os.path.join(self.input_location, "cooler_matrix", "MESC_500000.cool")
    input_file_name_mco = os.path.join(self.input_location, "cooler_matrix", "E14R1.mcool")

    # Output file
    output_file_name_bg = os.path.join(self.output_location, "read_bg.txt")
    output_file_bg = open(output_file_name_bg, "w")
    output_file_name_hic = os.path.join(self.output_location, "read_hic.txt")
    output_file_hic = open(output_file_name_hic, "w")
    output_file_name_co = os.path.join(self.output_location, "read_co.txt")
    output_file_co = open(output_file_name_co, "w")
    output_file_name_mco = os.path.join(self.output_location, "read_mco.txt")
    output_file_mco = open(output_file_name_mco, "w")
    output_file_name_un = os.path.join(self.output_location, "read_un.txt")
    output_file_un = open(output_file_name_un, "w")

    # Read bedgraph
    io = IO(input_file_name_bg, self.temporary_location, "mm9", 4, input_resolution = 500000, input_file_type = InputFileType.SPARSE)
    contact_map = io.read()
    contact_map.update_valid_chromosome_list()
    self.write_contact_map("contact map", contact_map, output_file_bg)

    io = IO(input_file_name_hic, self.temporary_location, "hg19", 4, input_resolution = 500000, input_file_type = InputFileType.HIC)
    contact_map = io.read()
    contact_map.update_valid_chromosome_list()
    self.write_contact_map("contact map", contact_map, output_file_hic)

    io = IO(input_file_name_co, self.temporary_location, "mm9", 4, input_resolution = 500000, input_file_type = InputFileType.COOL)
    contact_map = io.read()
    contact_map.update_valid_chromosome_list()
    self.write_contact_map("contact map", contact_map, output_file_co)

    io = IO(input_file_name_mco, self.temporary_location, "mm9", 4, input_resolution = 1000000, input_file_type = InputFileType.MCOOL)
    contact_map = io.read()
    contact_map.update_valid_chromosome_list()
    self.write_contact_map("contact map", contact_map, output_file_mco)

    io = IO(input_file_name_bg, self.temporary_location, "mm9", 4)
    contact_map = io.read()
    contact_map.update_valid_chromosome_list()
    self.write_contact_map("contact map", contact_map, output_file_un)

    # Termination
    output_file_bg.close()
    output_file_hic.close()
    output_file_co.close()
    output_file_mco.close()
    output_file_un.close()

    return contact_map

  def write(self, contact_map):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Input file
    input_file_name_bg = os.path.join(self.input_location, "bedgraph_matrix", "MESC_500000.bg")
    input_file_name_hic = os.path.join(self.input_location, "juicer_matrix", "SCONT_500000.hic")
    input_file_name_co = os.path.join(self.input_location, "cooler_matrix", "MESC_500000.cool")

    # Output file
    output_file_name_bg = os.path.join(self.output_location, "write_bg.bg")
    output_file_name_hic = os.path.join(self.output_location, "write_hic.hic")
    output_file_name_co = os.path.join(self.output_location, "write_co.cool")

    # Write bedgraph
    io = IO(input_file_name_bg, self.temporary_location, "mm9", 4, input_resolution = 500000, input_file_type = InputFileType.SPARSE)
    io.write(contact_map, output_file_name_bg, InputFileType.SPARSE)

    # Write hic
    io = IO(input_file_name_bg, self.temporary_location, "mm9", 4, input_resolution = 500000, input_file_type = InputFileType.SPARSE)
    io.write(contact_map, output_file_name_hic, InputFileType.HIC)

    # Write cool
    io = IO(input_file_name_bg, self.temporary_location, "mm9", 4, input_resolution = 500000, input_file_type = InputFileType.SPARSE)
    io.write(contact_map, output_file_name_co, InputFileType.COOL)


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
  test = IOTest(input_location, temporary_location, output_location)

  # Tests
  contact_map = test.read()
  test.write(contact_map)

