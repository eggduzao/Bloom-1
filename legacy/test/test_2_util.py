from __future__ import print_function
"""
Util Test
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
from bloom.util import ConfigurationFile, ChromosomeSizes, BarcodeFiles, ExcList, AuxiliaryFunctions

# External
import numpy as np


###################################################################################################
# UtilTest Class
###################################################################################################

class UtilTest():
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


  def configuration(self):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Output file
    output_file_name = os.path.join(self.output_location, "configuration.txt")
    output_file = open(output_file_name, "w")

    # Configuration file
    configuration_file = ConfigurationFile()
    self.write_list("ConfigurationFile", [configuration_file.bloom_data_path, configuration_file.bloom_config_file_name], output_file)

    # ChromosomeSizes
    chromosome_sizes = ChromosomeSizes("hg19")
    self.write_list("ChromosomeSizes", [chromosome_sizes.organism, chromosome_sizes.chromosome_sizes_file_name], output_file)
    self.write_list("chromosome_sizes_list", chromosome_sizes.chromosome_sizes_list, output_file)
    self.write_dict("chromosome_sizes_dictionary", chromosome_sizes.chromosome_sizes_dictionary, output_file)

    # BarcodeFiles
    barcode_files = BarcodeFiles()
    self.write_list("BarcodeFiles", [barcode_files.bcf], output_file)
    self.write_list("barcode_number_list", barcode_files.barcode_number_list, output_file)
    self.write_dict("barcode_file_dictionary", barcode_files.barcode_file_dictionary, output_file)

    # ExcList
    exclist = ExcList("hg19", 1000)
    self.write_list("ExcList", [exclist.exclist_file_name], output_file)
    self.write_dict("exclude_dictionary", exclist.exclude_dictionary, output_file)

    # Termination
    output_file.close()


  def auxiliary(self):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Output file
    output_file_name = os.path.join(self.output_location, "auxiliary.txt")
    output_file = open(output_file_name, "w")

    # string_is_int
    a = AuxiliaryFunctions.string_is_int("1")
    b = AuxiliaryFunctions.string_is_int("1.0")
    c = AuxiliaryFunctions.string_is_int("1e10")
    d = AuxiliaryFunctions.string_is_int(np.inf)
    self.write_list("string_is_int", [a, b, c, d], output_file)

    # string_is_float
    a = AuxiliaryFunctions.string_is_float("1")
    b = AuxiliaryFunctions.string_is_float("1.0")
    c = AuxiliaryFunctions.string_is_float("1e10")
    d = AuxiliaryFunctions.string_is_float(np.inf)
    self.write_list("string_is_float", [a, b, c, d], output_file)

    # string_is_validfile
    a = AuxiliaryFunctions.string_is_validfile(self.output_location)
    b = AuxiliaryFunctions.string_is_validfile(output_file_name)
    c = AuxiliaryFunctions.string_is_validfile("1")
    self.write_list("string_is_validfile", [a, b, c], output_file)

    # string_is_validpath
    a = AuxiliaryFunctions.string_is_validpath(self.output_location)
    b = AuxiliaryFunctions.string_is_validpath(output_file_name)
    c = AuxiliaryFunctions.string_is_validpath("1")
    self.write_list("string_is_validpath", [a, b, c], output_file)

    # overlap_count
    a = AuxiliaryFunctions.overlap_count([1, 2], [5, 8])
    b = AuxiliaryFunctions.overlap_count([-1, -2], [-5, -8])
    c = AuxiliaryFunctions.overlap_count([1, 2], [1, 2])
    d = AuxiliaryFunctions.overlap_count((0, 20), (15, 30))
    e = AuxiliaryFunctions.overlap_count([5, 8], (1, 2))
    self.write_list("overlap_count", [a, b, c, d, e], output_file)

    # ceil_multiple
    a = AuxiliaryFunctions.ceil_multiple(31, 15)
    b = AuxiliaryFunctions.ceil_multiple(-31, 15)
    self.write_list("ceil_multiple", [a, b], output_file)

    # floor_multiple
    a = AuxiliaryFunctions.floor_multiple(31, 15)
    b = AuxiliaryFunctions.floor_multiple(-31, 15)
    self.write_list("floor_multiple", [a, b], output_file)

    # shorten_integer
    a = AuxiliaryFunctions.shorten_integer("1000000000")
    b = AuxiliaryFunctions.shorten_integer("100000000")
    c = AuxiliaryFunctions.shorten_integer("100000")
    d = AuxiliaryFunctions.shorten_integer("123456789")
    self.write_list("shorten_integer", [a, b, c, d], output_file)

    # expand_integer
    a = AuxiliaryFunctions.expand_integer(a)
    b = AuxiliaryFunctions.expand_integer(b)
    c = AuxiliaryFunctions.expand_integer(c)
    d = AuxiliaryFunctions.expand_integer(d)
    self.write_list("expand_integer", [a, b, c, d], output_file)

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
  test = UtilTest(input_location, temporary_location, output_location)

  # Tests
  test.configuration()
  test.auxiliary()

