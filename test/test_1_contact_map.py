from __future__ import print_function
"""
Contact Map Test
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
#import unittest

# Internal
from bloom.contact_map import ContactMap

# External
import numpy as np


###################################################################################################
# ContactMapTest Class
###################################################################################################

class ContactMapTest():
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

  def load_contact_map(self):
    """Returns TODO.

    *Keyword arguments:*

      - argument -- An argument.

    *Return:*

      - return -- A return.
    """

    # Output file
    output_file_name = os.path.join(self.output_location, "1_load_contact_map.txt")
    output_file = open(output_file_name, "w")

    # Create blank contact map
    contact_map = ContactMap("hg19", 500000)
    self.write_list("organism", [contact_map.organism], output_file)
    self.write_list("resolution", [contact_map.resolution], output_file)
    self.write_contact_map("matrix", contact_map, output_file)
    ###
    self.write_dict("min_value_diagonal", contact_map.min_value_diagonal, output_file)
    self.write_dict("max_value_diagonal", contact_map.max_value_diagonal, output_file)
    self.write_dict("total_value_diagonal", contact_map.total_value_diagonal, output_file)
    self.write_dict("min_value_no_diagonal", contact_map.min_value_no_diagonal, output_file)
    self.write_dict("max_value_no_diagonal", contact_map.max_value_no_diagonal, output_file)
    self.write_dict("total_value_no_diagonal", contact_map.total_value_no_diagonal, output_file)
    ###
    self.write_dict("total_bins", contact_map.total_bins, output_file)
    self.write_dict("total_nonzero_bins", contact_map.total_nonzero_bins, output_file)
    self.write_dict("total_zero_bins", contact_map.total_zero_bins, output_file)
    ###
    self.write_dict("total_1d_bp", contact_map.total_1d_bp, output_file)
    self.write_dict("total_1d_bins", contact_map.total_1d_bins, output_file)
    ###
    self.write_dict("total_bins_triangle", contact_map.total_bins_triangle, output_file)
    self.write_dict("total_nonzero_bins_triangle", contact_map.total_nonzero_bins_triangle, output_file)
    self.write_dict("total_zero_bins_triangle", contact_map.total_zero_bins_triangle, output_file)
    ###
    self.write_list("valid_chromosome_list", [contact_map.valid_chromosome_list], output_file)

    # Termination
    output_file.close()

    return contact_map

  def matrix_operations(self):
    """Returns TODO.

    *Keyword arguments:*

      - argument -- An argument.

    *Return:*

      - return -- A return.
    """

    # Output file
    output_file_name = os.path.join(self.output_location, "2_matrix_operations.txt")
    output_file = open(output_file_name, "w")
    output_file_name_matrix1 = os.path.join(self.output_location, "2_matrix_operations_matrix1.txt")
    output_file_matrix1 = open(output_file_name_matrix1, "w")
    output_file_name_matrix2 = os.path.join(self.output_location, "2_matrix_operations_matrix2.txt")
    output_file_matrix2 = open(output_file_name_matrix2, "w")

    # Contact map
    contact_map = ContactMap("hg19", 500000)

    # Set
    contact_map.set("chr1", 10 * contact_map.resolution, 10 * contact_map.resolution, 1.0)
    contact_map.set("chr1", 5 * contact_map.resolution, 10 * contact_map.resolution, 2.0)
    contact_map.set("chr1", 6 * contact_map.resolution, 11 * contact_map.resolution, 3.0)
    contact_map.set("chr1", 9 * contact_map.resolution, 10 * contact_map.resolution, 4.2)
    contact_map.set("chr10", 10 * contact_map.resolution, 20 * contact_map.resolution, 10.0)
    self.write_contact_map("set", contact_map, output_file)

    # Add
    contact_map.add("chr1", 10 * contact_map.resolution, 10 * contact_map.resolution, 9.0)
    contact_map.add("chr1", 5 * contact_map.resolution, 10 * contact_map.resolution, 10.0)
    contact_map.add("chr3", 5 * contact_map.resolution, 10 * contact_map.resolution, 10.0)
    self.write_contact_map("add", contact_map, output_file)

    # Get
    value1 = contact_map.get("chr1", 10 * contact_map.resolution, 10 * contact_map.resolution)
    value2 = contact_map.get("chr3", 5 * contact_map.resolution, 10 * contact_map.resolution)
    value3 = contact_map.get("chr3", 6 * contact_map.resolution, 10 * contact_map.resolution)
    self.write_list("get", [value1, value2, value3], output_file)

    # Delete strips
    contact_map.delete_strips("chr1", 3 * contact_map.resolution, 8 * contact_map.resolution)
    value1 = contact_map.get("chr1", 5 * contact_map.resolution, 10 * contact_map.resolution)
    value2 = contact_map.get("chr1", 6 * contact_map.resolution, 11 * contact_map.resolution)
    value3 = contact_map.get("chr1", 10 * contact_map.resolution, 10 * contact_map.resolution)
    assert value1 == 0, "should be 0"
    assert value2 == 0, "should be 0"
    assert value3 == 10.0, "should be 10.0"

    # Set from matrix
    numpy_matrix = np.random.rand(contact_map.total_1d_bins["chr2"],contact_map.total_1d_bins["chr2"])
    contact_map.set_from_matrix("chr2", numpy_matrix, matrix_type = "numpy_array", storage_type = "upper_triangle")
    self.write_contact_map("set from matrix", contact_map, output_file_matrix1)
    value1 = contact_map.get("chr2", 1 * contact_map.resolution, 1 * contact_map.resolution)
    value2 = contact_map.get("chr2", 4 * contact_map.resolution, 8 * contact_map.resolution)
    value3 = contact_map.get("chr2", 5 * contact_map.resolution, 20 * contact_map.resolution)
    self.write_list("set from matrix", [value1, value2, value3], output_file)

    # Get full matrix
    new_numpy_matrix = contact_map.get_full_matrix("chr2", symmetric = True, return_type = "numpy_array")
    self.write_matrix("get full matrix", new_numpy_matrix, output_file_matrix2)
    value1 = contact_map.get("chr2", 1 * contact_map.resolution, 1 * contact_map.resolution)
    value2 = contact_map.get("chr2", 4 * contact_map.resolution, 8 * contact_map.resolution)
    value3 = contact_map.get("chr2", 5 * contact_map.resolution, 20 * contact_map.resolution)
    self.write_list("get from matrix", [value1, value2, value3], output_file)

    # Termination
    output_file.close()
    output_file_matrix1.close()
    output_file_matrix2.close()

  def resolution_bin_bp_operations(self):
    """Returns TODO.

    *Keyword arguments:*

      - argument -- An argument.

    *Return:*

      - return -- A return.
    """

    # Output file
    output_file_name = os.path.join(self.output_location, "3_resolution_bin_bp_operations.txt")
    output_file = open(output_file_name, "w")

    # Contact map
    contact_map = ContactMap("hg19", 500000)

    # BP to bin
    bin11 = contact_map.bp_to_bin(10 * contact_map.resolution)
    bin21, bin22 = contact_map.bp_to_bin(5 * contact_map.resolution, 7 * contact_map.resolution)
    self.write_list("bp to bin", [bin11, bin21, bin22], output_file)

    # bin to BP
    bp11 = contact_map.bin_to_bp(10)
    bp21, bp22 = contact_map.bin_to_bp(5, 7)
    self.write_list("bin to bp", [bp11, bp21, bp22], output_file)

    # Ceil BP
    c1 = contact_map.ceil_bp(1000)
    (c2, c3) = contact_map.ceil_bp((123456, 654321))
    c4 = contact_map.ceil_bp(600000.0)
    self.write_list("ceil bp", [c1, c2, c3, c4], output_file)

    # Floor BP
    c1 = contact_map.floor_bp(1000)
    (c2, c3) = contact_map.floor_bp((123456, 654321))
    c4 = contact_map.floor_bp(600000.0)
    self.write_list("ceil bp", [c1, c2, c3, c4], output_file)

    # Termination
    output_file.close()

  def auxiliary_vectors_operations(self):
    """Returns TODO.

    *Keyword arguments:*

      - argument -- An argument.

    *Return:*

      - return -- A return.
    """

    # Output file
    output_file_name = os.path.join(self.output_location, "4_auxiliary_vectors_operations.txt")
    output_file = open(output_file_name, "w")

    # Empty chromosome list
    contact_map = ContactMap("hg19", 500000)
    contact_map.update_valid_chromosome_list()
    self.write_list("empty - valid chromosome list", [contact_map.valid_chromosome_list], output_file)

    # All statistics
    contact_map = ContactMap("hg19", 500000)
    contact_map.set("chr1", 1 * contact_map.resolution, 1 * contact_map.resolution, 1.0)
    contact_map.set("chr1", 2 * contact_map.resolution, 2 * contact_map.resolution, 3.0)
    contact_map.set("chr1", 1 * contact_map.resolution, 2 * contact_map.resolution, 5.0)
    contact_map.set("chr1", 3 * contact_map.resolution, 10 * contact_map.resolution, 2.0)
    contact_map.set("chr1", 10 * contact_map.resolution, 20 * contact_map.resolution, 4.0)
    contact_map.set("chr2", 1 * contact_map.resolution, 1 * contact_map.resolution, 1.0)
    contact_map.set("chr3", 1 * contact_map.resolution, 1 * contact_map.resolution, 1.0)
    contact_map.set("chr4", 1 * contact_map.resolution, 1 * contact_map.resolution, 1.0)
    contact_map.set("chr5", 1 * contact_map.resolution, 1 * contact_map.resolution, 1.0)
    contact_map.set("chr10", 1 * contact_map.resolution, 1 * contact_map.resolution, 1.0)
    contact_map.update_valid_chromosome_list()
    self.write_list("all - valid_chromosome_list", [contact_map.valid_chromosome_list], output_file)
    contact_map.calculate_all_statistics()
    ###
    self.write_list("all - organism", [contact_map.organism], output_file)
    self.write_list("all - resolution", [contact_map.resolution], output_file)
    self.write_contact_map("all - matrix", contact_map, output_file)
    self.write_dict("all - min_value_diagonal", contact_map.min_value_diagonal, output_file)
    self.write_dict("all - max_value_diagonal", contact_map.max_value_diagonal, output_file)
    self.write_dict("all - total_value_diagonal", contact_map.total_value_diagonal, output_file)
    self.write_dict("all - min_value_no_diagonal", contact_map.min_value_no_diagonal, output_file)
    self.write_dict("all - max_value_no_diagonal", contact_map.max_value_no_diagonal, output_file)
    self.write_dict("all - total_value_no_diagonal", contact_map.total_value_no_diagonal, output_file)
    self.write_dict("all - total_bins", contact_map.total_bins, output_file)
    self.write_dict("all - total_nonzero_bins", contact_map.total_nonzero_bins, output_file)
    self.write_dict("all - total_zero_bins", contact_map.total_zero_bins, output_file)
    self.write_dict("all - total_1d_bp", contact_map.total_1d_bp, output_file)
    self.write_dict("all - total_1d_bins", contact_map.total_1d_bins, output_file)
    self.write_dict("all - total_bins_triangle", contact_map.total_bins_triangle, output_file)
    self.write_dict("all - total_nonzero_bins_triangle", contact_map.total_nonzero_bins_triangle, output_file)
    self.write_dict("all - total_zero_bins_triangle", contact_map.total_zero_bins_triangle, output_file)

    # Statistic value
    contact_map = ContactMap("hg19", 500000)
    contact_map.set("chr1", 1 * contact_map.resolution, 1 * contact_map.resolution, 1.0)
    contact_map.set("chr1", 2 * contact_map.resolution, 2 * contact_map.resolution, 3.0)
    contact_map.set("chr1", 1 * contact_map.resolution, 2 * contact_map.resolution, 5.0)
    contact_map.set("chr1", 3 * contact_map.resolution, 10 * contact_map.resolution, 2.0)
    contact_map.set("chr1", 10 * contact_map.resolution, 20 * contact_map.resolution, 4.0)
    contact_map.set("chr2", 1 * contact_map.resolution, 1 * contact_map.resolution, 1.0)
    contact_map.set("chr3", 1 * contact_map.resolution, 1 * contact_map.resolution, 1.0)
    contact_map.set("chr4", 1 * contact_map.resolution, 1 * contact_map.resolution, 1.0)
    contact_map.set("chr5", 1 * contact_map.resolution, 1 * contact_map.resolution, 1.0)
    contact_map.set("chr10", 1 * contact_map.resolution, 1 * contact_map.resolution, 1.0)
    contact_map.update_valid_chromosome_list()
    contact_map.calculate_statistics_value()
    self.write_list("value - organism", [contact_map.organism], output_file)
    self.write_list("value - resolution", [contact_map.resolution], output_file)
    self.write_contact_map("value - matrix", contact_map, output_file)
    self.write_dict("value - min_value_diagonal", contact_map.min_value_diagonal, output_file)
    self.write_dict("value - max_value_diagonal", contact_map.max_value_diagonal, output_file)
    self.write_dict("value - total_value_diagonal", contact_map.total_value_diagonal, output_file)
    self.write_dict("value - min_value_no_diagonal", contact_map.min_value_no_diagonal, output_file)
    self.write_dict("value - max_value_no_diagonal", contact_map.max_value_no_diagonal, output_file)
    self.write_dict("value - total_value_no_diagonal", contact_map.total_value_no_diagonal, output_file)
    self.write_dict("value - total_bins", contact_map.total_bins, output_file)
    self.write_dict("value - total_nonzero_bins", contact_map.total_nonzero_bins, output_file)
    self.write_dict("value - total_zero_bins", contact_map.total_zero_bins, output_file)
    self.write_dict("value - total_1d_bp", contact_map.total_1d_bp, output_file)
    self.write_dict("value - total_1d_bins", contact_map.total_1d_bins, output_file)
    self.write_dict("value - total_bins_triangle", contact_map.total_bins_triangle, output_file)
    self.write_dict("value - total_nonzero_bins_triangle", contact_map.total_nonzero_bins_triangle, output_file)
    self.write_dict("value - total_zero_bins_triangle", contact_map.total_zero_bins_triangle, output_file)
    self.write_list("value - valid_chromosome_list", [contact_map.valid_chromosome_list], output_file)

    # Statistic non-value
    contact_map = ContactMap("hg19", 500000)
    contact_map.set("chr1", 1 * contact_map.resolution, 1 * contact_map.resolution, 1.0)
    contact_map.set("chr1", 2 * contact_map.resolution, 2 * contact_map.resolution, 3.0)
    contact_map.set("chr1", 1 * contact_map.resolution, 2 * contact_map.resolution, 5.0)
    contact_map.set("chr1", 3 * contact_map.resolution, 10 * contact_map.resolution, 2.0)
    contact_map.set("chr1", 10 * contact_map.resolution, 20 * contact_map.resolution, 4.0)
    contact_map.set("chr2", 1 * contact_map.resolution, 1 * contact_map.resolution, 1.0)
    contact_map.set("chr3", 1 * contact_map.resolution, 1 * contact_map.resolution, 1.0)
    contact_map.set("chr4", 1 * contact_map.resolution, 1 * contact_map.resolution, 1.0)
    contact_map.set("chr5", 1 * contact_map.resolution, 1 * contact_map.resolution, 1.0)
    contact_map.set("chr10", 1 * contact_map.resolution, 1 * contact_map.resolution, 1.0)
    contact_map.update_valid_chromosome_list()
    contact_map.calculate_all_non_value_statistics()
    self.write_list("nonvalue - organism", [contact_map.organism], output_file)
    self.write_list("nonvalue - resolution", [contact_map.resolution], output_file)
    self.write_contact_map("nonvalue - matrix", contact_map, output_file)
    self.write_dict("nonvalue - min_value_diagonal", contact_map.min_value_diagonal, output_file)
    self.write_dict("nonvalue - max_value_diagonal", contact_map.max_value_diagonal, output_file)
    self.write_dict("nonvalue - total_value_diagonal", contact_map.total_value_diagonal, output_file)
    self.write_dict("nonvalue - min_value_no_diagonal", contact_map.min_value_no_diagonal, output_file)
    self.write_dict("nonvalue - max_value_no_diagonal", contact_map.max_value_no_diagonal, output_file)
    self.write_dict("nonvalue - total_value_no_diagonal", contact_map.total_value_no_diagonal, output_file)
    self.write_dict("nonvalue - total_bins", contact_map.total_bins, output_file)
    self.write_dict("nonvalue - total_nonzero_bins", contact_map.total_nonzero_bins, output_file)
    self.write_dict("nonvalue - total_zero_bins", contact_map.total_zero_bins, output_file)
    self.write_dict("nonvalue - total_1d_bp", contact_map.total_1d_bp, output_file)
    self.write_dict("nonvalue - total_1d_bins", contact_map.total_1d_bins, output_file)
    self.write_dict("nonvalue - total_bins_triangle", contact_map.total_bins_triangle, output_file)
    self.write_dict("nonvalue - total_nonzero_bins_triangle", contact_map.total_nonzero_bins_triangle, output_file)
    self.write_dict("nonvalue - total_zero_bins_triangle", contact_map.total_zero_bins_triangle, output_file)
    self.write_list("nonvalue - valid_chromosome_list", [contact_map.valid_chromosome_list], output_file)

    # Termination
    output_file.close()

  def distance_operations(self):
    """Returns TODO.

    *Keyword arguments:*

      - argument -- An argument.

    *Return:*

      - return -- A return.
    """

    # Output file
    output_file_name = os.path.join(self.output_location, "5_distance_operations.txt")
    output_file = open(output_file_name, "w")

    # Contact map
    contact_map = ContactMap("hg19", 500000)

    # Manhattan
    dist_1 = contact_map.bin_distance_from_diagonal_manhattan((1, 1))
    dist_2 = contact_map.bin_distance_from_diagonal_manhattan((1, 2))
    dist_3 = contact_map.bin_distance_from_diagonal_manhattan((1, 4))
    self.write_list("manhattan", [dist_1, dist_2, dist_3], output_file)

    # Euclidean
    dist_1 = contact_map.bin_distance_from_diagonal_euclidean((1, 1))
    dist_2 = contact_map.bin_distance_from_diagonal_euclidean((1, 2))
    dist_3 = contact_map.bin_distance_from_diagonal_euclidean((1, 4))
    self.write_list("euclidean", [dist_1, dist_2, dist_3], output_file)

    # Termination
    output_file.close()

  def sparsity_statistics(self):
    """Returns TODO.

    *Keyword arguments:*

      - argument -- An argument.

    *Return:*

      - return -- A return.
    """

    # Output file
    output_file_name = os.path.join(self.output_location, "6_sparsity_statistics.txt")
    output_file = open(output_file_name, "w")

    # TODO

    # Termination
    output_file.close()

  def multimatrix_operations(self):
    """Returns TODO.

    *Keyword arguments:*

      - argument -- An argument.

    *Return:*

      - return -- A return.
    """

    # Output file
    output_file_name = os.path.join(self.output_location, "7_multimatrix_operations.txt")
    output_file = open(output_file_name, "w")

    # Original
    contact_map_1 = ContactMap("hg19", 500000)
    contact_map_1.set("chr1", 1 * contact_map_1.resolution, 10 * contact_map_1.resolution, 1.0)
    contact_map_1.set("chr1", 2 * contact_map_1.resolution, 10 * contact_map_1.resolution, 2.0)
    contact_map_1.set("chr1", 3 * contact_map_1.resolution, 10 * contact_map_1.resolution, 3.0)

    # 2 - True - same entries
    contact_map_2 = ContactMap("hg19", 500000)
    contact_map_2.set("chr1", 1 * contact_map_2.resolution, 10 * contact_map_2.resolution, 1.0)
    contact_map_2.set("chr1", 2 * contact_map_2.resolution, 10 * contact_map_2.resolution, 2.0)
    contact_map_2.set("chr1", 3 * contact_map_2.resolution, 10 * contact_map_2.resolution, 3.0)

    # 3 - True - subset of entries
    contact_map_3 = ContactMap("hg19", 500000)
    contact_map_3.set("chr1", 1 * contact_map_3.resolution, 10 * contact_map_3.resolution, 1.0)
    contact_map_3.set("chr1", 2 * contact_map_3.resolution, 10 * contact_map_3.resolution, 2.0)

    # 4 - False - completely different
    contact_map_4 = ContactMap("hg19", 500000)
    contact_map_4.set("chr1", 4 * contact_map_4.resolution, 10 * contact_map_4.resolution, 1.0)
    contact_map_4.set("chr1", 5 * contact_map_4.resolution, 10 * contact_map_4.resolution, 2.0)
    contact_map_4.set("chr1", 6 * contact_map_4.resolution, 10 * contact_map_4.resolution, 3.0)

    # 5 - True - slightly different values
    contact_map_5 = ContactMap("hg19", 500000)
    contact_map_5.set("chr1", 1 * contact_map_5.resolution, 10 * contact_map_5.resolution, 1.1)
    contact_map_5.set("chr1", 2 * contact_map_5.resolution, 10 * contact_map_5.resolution, 1.9)
    contact_map_5.set("chr1", 3 * contact_map_5.resolution, 10 * contact_map_5.resolution, 3.1)

    # 6 - False - very different values
    contact_map_6 = ContactMap("hg19", 500000)
    contact_map_6.set("chr1", 1 * contact_map_6.resolution, 10 * contact_map_6.resolution, 1.2)
    contact_map_6.set("chr1", 2 * contact_map_6.resolution, 10 * contact_map_6.resolution, 2.0)
    contact_map_6.set("chr1", 3 * contact_map_6.resolution, 10 * contact_map_6.resolution, 3.0)

    # 7 - False - equal entries, but bigger
    contact_map_7 = ContactMap("hg19", 500000)
    contact_map_7.set("chr1", 1 * contact_map_7.resolution, 10 * contact_map_7.resolution, 1.0)
    contact_map_7.set("chr1", 2 * contact_map_7.resolution, 10 * contact_map_7.resolution, 2.0)
    contact_map_7.set("chr1", 3 * contact_map_7.resolution, 10 * contact_map_7.resolution, 3.0)
    contact_map_7.set("chr1", 3 * contact_map_7.resolution, 11 * contact_map_7.resolution, 3.0)

    # 8 - False - Empty
    contact_map_8 = ContactMap("hg19", 500000)

    # Comparison
    result_12 = contact_map_1.match_subset("chr1", contact_map_2, similarity_degree = 0.1)
    result_13 = contact_map_1.match_subset("chr1", contact_map_3, similarity_degree = 0.1)
    result_14 = contact_map_1.match_subset("chr1", contact_map_4, similarity_degree = 0.1)
    result_15 = contact_map_1.match_subset("chr1", contact_map_5, similarity_degree = 0.1)
    result_16 = contact_map_1.match_subset("chr1", contact_map_6, similarity_degree = 0.1)
    result_17 = contact_map_1.match_subset("chr1", contact_map_7, similarity_degree = 0.1)
    result_18 = contact_map_1.match_subset("chr1", contact_map_8, similarity_degree = 0.1)
    self.write_list("match subset", [result_12, result_13, result_14, result_15, result_16, result_17, result_18], output_file)

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
  test = ContactMapTest(input_location, temporary_location, output_location)

  # Tests
  test.load_contact_map()
  test.matrix_operations()
  test.resolution_bin_bp_operations()
  test.auxiliary_vectors_operations()
  test.distance_operations()
  test.sparsity_statistics()
  test.multimatrix_operations()

# cat test/output/test_1_contact_map/1_load_contact_map.txt
# cat test/output/test_1_contact_map/2_matrix_operations.txt
# cat test/output/test_1_contact_map/3_resolution_bin_bp_operations.txt
# cat test/output/test_1_contact_map/4_auxiliary_vectors_operations.txt
# cat test/output/test_1_contact_map/5_distance_operations.txt
# cat test/output/test_1_contact_map/6_sparsity_statistics.txt
# cat test/output/test_1_contact_map/7_multimatrix_operations.txt
