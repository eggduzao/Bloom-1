from __future__ import print_function
"""
Juicer Test
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
from bloom.io_juicer import Juicer
from bloom.io_bedgraph import Bedgraph

# External
import numpy as np

"""

java -Djava.awt.headless=true -Xmx32000m -jar /usr/users/egadegu/bloom_data/bin/juicer_tools_1.22.01.jar dump observed NONE /home/uni08/egadegu/Projects/Papantonis_Bloom/Bloom/test/output/test_4_io_juicer/load.hic 1:0:249000000 1:0:249000000 BP 500000 /home/uni08/egadegu/Projects/Papantonis_Bloom/Bloom/test/output/test_4_io_juicer/load_dump.txt

"""


###################################################################################################
# JuicerTest Class
###################################################################################################

class JuicerTest():
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

  def dump(self):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Input file
    input_file_name = os.path.join(self.input_location, "juicer_matrix", "SCONT_500000.hic")

    # Output file
    output_file_name_juicer = os.path.join(self.output_location, "dump_juicer.txt")
    output_file_name_bedgraph = os.path.join(self.output_location, "dump_bedgraph.txt")

    # Contact map
    contact_map = ContactMap("hg19", 500000)
    contact_map.update_valid_chromosome_list()

    # Dump juicer / bedgraph
    juicer = Juicer(4)
    region = ":".join(["1", "0", str(contact_map.total_1d_bp["chr1"])])
    juicer.dump(contact_map.resolution, region, region, input_file_name, self.temporary_location, output_file_name_juicer, output_type = "juicer")
    juicer.dump(contact_map.resolution, region, region, input_file_name, self.temporary_location, output_file_name_bedgraph, output_type = "bedgraph")

    return output_file_name_bedgraph


  def load(self, bedgraph_file_name):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Output file
    output_file_name = os.path.join(self.output_location, "load.hic")

    # Contact map
    contact_map = ContactMap("hg19", 500000)
    contact_map.update_valid_chromosome_list()
    
    # Load contact map from bedgraph
    bedgraph = Bedgraph("hg19", 4)
    bedgraph.dump("chr1", bedgraph_file_name, contact_map, logit = False, pseudocount = 0.0)
    contact_map.update_valid_chromosome_list()

    # Load
    juicer = Juicer(4)
    juicer.load(contact_map, self.temporary_location, output_file_name)


  def identify_minimal_resolution(self):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Input file
    input_file_name_1 = os.path.join(self.input_location, "juicer_matrix", "SCONT_25000.hic")
    input_file_name_2 = os.path.join(self.input_location, "juicer_matrix", "CM_unbloom.hic")

    # Output file
    output_file_name = os.path.join(self.output_location, "identify_minimal_resolution.txt")
    output_file = open(output_file_name, "w")

    # Identify resolution - single resolution
    juicer = Juicer(4)
    resolution_1 = juicer.identify_minimal_resolution(input_file_name_1, self.temporary_location, region = "1:1000000:5000000")

    # Identify resolution - multiple resolutions
    resolution_2 = juicer.identify_minimal_resolution(input_file_name_2, self.temporary_location, region = "1:1000000:5000000")
    self.write_list("resolution", [resolution_1, resolution_2], output_file)

    # Termination
    output_file.close()


  def filetype_is_juicer(self):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Input file
    input_file_name_1 = os.path.join(self.input_location, "juicer_matrix", "SCONT_25000.hic")
    input_file_name_2 = os.path.join(self.input_location, "juicer_matrix", "CM_unbloom.hic")
    input_file_name_bedgraph = os.path.join(self.input_location, "bedgraph_matrix", "MESC_500000.bg")

    # Output file
    output_file_name = os.path.join(self.output_location, "filetype_is_juicer.txt")
    output_file = open(output_file_name, "w")

    # Identify filetype
    juicer = Juicer(4)
    is_juicer_1 = juicer.filetype_is_juicer(input_file_name_1, self.temporary_location, region = "1:1000000:5000000")
    is_juicer_2 = juicer.filetype_is_juicer(input_file_name_2, self.temporary_location, region = "1:1000000:5000000")
    is_juicer_3 = juicer.filetype_is_juicer(input_file_name_bedgraph, self.temporary_location, region = "1:1000000:5000000")
    self.write_list("file type is juicer", [is_juicer_1, is_juicer_2, is_juicer_3], output_file)

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
  test = JuicerTest(input_location, temporary_location, output_location)

  # Tests
  bedgraph_file_name = test.dump()
  test.load(bedgraph_file_name)
  test.identify_minimal_resolution()
  test.filetype_is_juicer()

