from __future__ import print_function
"""
Pipeline Test
===================
Placeholder.

Authors: Eduardo G. Gusmao.

# Issues:

1. In the preprocess module, when we augment the resolution, the diagonal gets broader, therefore the number of bins to
   void (remove) has to be bigger somehow based on the new resolution.

2. Create diagonal_borderline in SICA.

3. Verify the formulas for star, goba fill, dpmm degrade.

4. Dpmm fill shape:
   4.1. Check perfect min and max edges
   4.2. Check why the value tends to go to the bottom right of the contact map!

5. Why IFS is reporting values slightly larger than 1.

"""

###################################################################################################
# Libraries
###################################################################################################

# Python
import os
import sys
import subprocess
from time import time

# Internal
#from bloom.contact_map import ContactMap
from bloom.io import InputFileType, IO
from bloom.preprocess import Preprocess
from bloom.sica import Sica
from bloom.goba import Goba
from bloom.dpmm import Dpmm
from bloom.ifs import Ifs

# External
import numpy as np


###################################################################################################
# PipelineTest Class
###################################################################################################

class PipelineTest():
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

  def write_dict_of_dict(self, name, dict_to_write, output_file):
    """Write dict_to_write to output_file
    """
    output_file.write("### " + name + "\n")
    if(dict_to_write):
      keys = sorted(dict_to_write.keys())
      if(keys):
        for k in keys:
          output_file.write("#" + str(k) + "\n")
          for kk in dict_to_write[k]:
            output_file.write(str(kk) + "\t" + str(dict_to_write[k][kk]) + "\n")
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



  ###################################################################################################
  # 8. Preprocess
  ###################################################################################################

  def convert_to_minimal_resolution(self):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Input file
    #input_file_name_bg = os.path.join(self.input_location, "bedgraph_matrix", "MESC_500000.bg")
    input_file_name_bg = os.path.join(self.input_location, "test_matrix", "MESC_100000_chrX_first5M.bg")

    # Output file
    output_initial_matrix_bg = os.path.join(self.output_location, "1_initial_matrix.bg")
    output_initial_matrix_hc = os.path.join(self.output_location, "1_initial_matrix.hic")

    output_minres_matrix_bg = os.path.join(self.output_location, "2_PRE_minres_matrix.bg")
    output_minres_matrix_hc = os.path.join(self.output_location, "2_PRE_minres_matrix.hic")

    # Read bedgraph
    io = IO(input_file_name_bg, self.temporary_location, "mm9", 4, input_resolution = 100000, input_file_type = InputFileType.SPARSE)
    contact_map = io.read()
    #contact_map.update_valid_chromosome_list()
    contact_map.valid_chromosome_list = ["chrX"] # Force to play only with chrX
    
    # Write bedgraph and hic
    io.write(contact_map, output_initial_matrix_bg, InputFileType.SPARSE)
    io.write(contact_map, output_initial_matrix_hc, InputFileType.HIC)

    # Convert to minimal resolution
    preprocess = Preprocess(4, contact_map, minimal_resolution = 50000, min_contig_removed_bins = 2, remove_threshold = 1)
    minimal_res_contact_map = preprocess.convert_to_minimal_resolution(recalculate_statistics = True)

    # Write bedgraph and hic
    io = IO(input_file_name_bg, self.temporary_location, "mm9", 4, input_resolution = 50000, input_file_type = InputFileType.SPARSE)
    io.write(minimal_res_contact_map, output_minres_matrix_bg, InputFileType.SPARSE)
    io.write(minimal_res_contact_map, output_minres_matrix_hc, InputFileType.HIC)

    # Return objects
    return io, preprocess, minimal_res_contact_map

  def main_remove_blacklist(self, io, preprocess, contact_map):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Output file
    output_blackrem_matrix_bg = os.path.join(self.output_location, "3_PRE_blackrem_matrix.bg")
    output_blackrem_matrix_hc = os.path.join(self.output_location, "3_PRE_blackrem_matrix.hic")
    output_removed_blacklist_file_name = os.path.join(self.output_location, "3_PRE_removed_blacklist_file_name.txt")
    output_removed_blacklist_file = open(output_removed_blacklist_file_name, "w")

    # Remove blacklisted regions
    preprocess.main_remove_blacklist(contact_map)

    # Write bedgraph and hic
    io.write(contact_map, output_blackrem_matrix_bg, InputFileType.SPARSE)
    io.write(contact_map, output_blackrem_matrix_hc, InputFileType.HIC)

    # Write removed dictionary
    self.write_dict_of_dict("output_removed_blacklist", preprocess.removed_dict, output_removed_blacklist_file)

    # Return objects
    output_removed_blacklist_file.close()
    return io, preprocess, contact_map


  def main_void_statistics(self, io, preprocess, contact_map):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Output file
    output_voidrem_matrix_bg = os.path.join(self.output_location, "4_PRE_voidrem_matrix.bg")
    output_voidrem_matrix_hc = os.path.join(self.output_location, "4_PRE_voidrem_matrix.hic")
    output_voidrem_file_name = os.path.join(self.output_location, "4_PRE_voidrem_file_name.txt")
    output_voidrem_file = open(output_voidrem_file_name, "w")

    # Removing void spaces
    preprocess.main_void_statistics(contact_map)

    # Write bedgraph and hic
    io.write(contact_map, output_voidrem_matrix_bg, InputFileType.SPARSE)
    io.write(contact_map, output_voidrem_matrix_hc, InputFileType.HIC)

    # Write removed dictionary
    self.write_dict_of_dict("output_voidrem_file", preprocess.removed_dict, output_voidrem_file)

    # Termination
    output_voidrem_file.close()
    return io, preprocess, contact_map



  ###################################################################################################
  # 9. SICA
  ###################################################################################################

  def main_calculate_distributions(self, io, preprocess, contact_map):

    # Output file
    annotdict_output_file_name = os.path.join(self.output_location, "5_1_SICA_annotdict.txt")
    distdict_output_file_name = os.path.join(self.output_location, "5_2_SICA_distdict.txt")
    distdiag_output_file_name = os.path.join(self.output_location, "5_3_SICA_distdiag.txt")
    pvalue_output_file_name = os.path.join(self.output_location, "5_4_SICA_pvalue.txt")
    sigvalue_output_file_name = os.path.join(self.output_location, "5_5_SICA_sigvalue.txt")
    annotdict_output_file = open(annotdict_output_file_name, "w")
    distdict_output_file = open(distdict_output_file_name, "w")
    distdiag_output_file = open(distdiag_output_file_name, "w")
    pvalue_output_file = open(pvalue_output_file_name, "w")
    sigvalue_output_file = open(sigvalue_output_file_name, "w")

    # Create SICA object
    sica = Sica(4, contact_map, avoid_distance = 100000, removed_dict = preprocess.removed_dict, pvalue_threshold = 0.95)

    # Calculating distributions and pvalues
    sica.main_calculate_distributions()

    # Writing
    self.write_dict_of_dict("annotation_dictionary", sica.annotation_dictionary, annotdict_output_file)
    self.write_dict_of_dict("distribution_dictionary", sica.distribution_dictionary, distdict_output_file)
    self.write_dict_of_dict("dist_to_diag_dictionary", sica.dist_to_diag_dictionary, distdiag_output_file)
    self.write_dict_of_dict("pvalue_dictionary", sica.pvalue_dictionary, pvalue_output_file)
    self.write_dict("significant_values_dictionary", sica.significant_values_dictionary, sigvalue_output_file)

    # Termination
    annotdict_output_file.close()
    distdict_output_file.close()
    distdiag_output_file.close()
    pvalue_output_file.close()
    sigvalue_output_file.close()
    return io, preprocess, sica, contact_map

    # Annotation dictionary
    #self.annotation_dictionary = dict() # Same as contact_map matrix, but with an extra annotation flag:
    # R = Points removed because they fall into a 0 contig or blacklist.
    # A = Points to avoid because they fall in 0-avoid_distance
    # T1, T2, T3, T4, T5 = Tad hierarchy levels depending on point's distance to diagonal.
    # C1, C2, C3 = Compartment hierarchy levels depending on point's distance to diagonal.
    # O = Points farther from diagonal than the biggest compartment length (C1)
    # Upper case letters are important points. Lower case letters are less important points.

    # Auxiliary dictionaries
    #self.distribution_dictionary = dict() # per chromosome / per SicaDist -> Contains all the matrix's signal within that specific SicaDist
    #self.dist_to_diag_dictionary = dict() # per chromosome / per regular matrix key -> Manhattan distance to the diagonal
    #self.pvalue_dictionary = dict() # per chromosome / per SicaDist -> Contains [name of the fitted distribution (FD), parameters of FD, value at pvalue_threshold (given FD)]
    #self.significant_values_dictionary = dict() # per chromosome -> All significant (peaks) values of the matrix, i.e. >= value at pvalue_threshold (given FD)


  def main_star_contacts(self, io, preprocess, sica, contact_map):

    # Output matrix
    output_star_matrix_bg = os.path.join(self.output_location, "6_SICA_star_matrix.bg")
    output_star_matrix_hc = os.path.join(self.output_location, "6_SICA_star_matrix.hic")

    # Output files
    annotdict_output_file_name = os.path.join(self.output_location, "6_1_SICA_annotdict.txt")
    distdict_output_file_name = os.path.join(self.output_location, "6_2_SICA_distdict.txt")
    distdiag_output_file_name = os.path.join(self.output_location, "6_3_SICA_distdiag.txt")
    pvalue_output_file_name = os.path.join(self.output_location, "6_4_SICA_pvalue.txt")
    sigvalue_output_file_name = os.path.join(self.output_location, "6_5_SICA_sigvalue.txt")
    annotdict_output_file = open(annotdict_output_file_name, "w")
    distdict_output_file = open(distdict_output_file_name, "w")
    distdiag_output_file = open(distdiag_output_file_name, "w")
    pvalue_output_file = open(pvalue_output_file_name, "w")
    sigvalue_output_file = open(sigvalue_output_file_name, "w")

    # Star significant interactions
    sica.main_star_contacts(bottom_bin_ext_range = [2,4], left_bin_ext_range = [2,4], right_bin_ext_range = [1,2], top_bin_ext_range = [1,2])

    # Writing
    io.write(contact_map, output_star_matrix_bg, InputFileType.SPARSE)
    io.write(contact_map, output_star_matrix_hc, InputFileType.HIC)
    self.write_dict_of_dict("annotdict_output_file", sica.annotation_dictionary, annotdict_output_file)
    self.write_dict_of_dict("distdict_output_file", sica.distribution_dictionary, distdict_output_file)
    self.write_dict_of_dict("distdiag_output_file", sica.dist_to_diag_dictionary, distdiag_output_file)
    self.write_dict_of_dict("pvalue_output_file", sica.pvalue_dictionary, pvalue_output_file)
    self.write_dict("sigvalue_output_file", sica.significant_values_dictionary, sigvalue_output_file)

    # Termination
    annotdict_output_file.close()
    distdict_output_file.close()
    distdiag_output_file.close()
    pvalue_output_file.close()
    sigvalue_output_file.close()
    return io, preprocess, sica, contact_map



  ###################################################################################################
  # 10. GOBA
  ###################################################################################################

  def main_fill(self, io, preprocess, sica, contact_map):

    # Output file
    output_goba_matrix_bg = os.path.join(self.output_location, "7_GOBA_fill_matrix.bg")
    output_goba_matrix_hc = os.path.join(self.output_location, "7_GOBA_fill_matrix.hic")

    # Goba
    goba = Goba(contact_map, sica)

    # Performing fill
    goba.main_fill()

    # Writing
    io.write(contact_map, output_goba_matrix_bg, InputFileType.SPARSE)
    io.write(contact_map, output_goba_matrix_hc, InputFileType.HIC)

    # Termination
    return io, preprocess, sica, goba, contact_map



  ###################################################################################################
  # 11. DPMM
  ###################################################################################################

  def diagonal_degrade(self, io, preprocess, sica, goba, contact_map):

    # Output file
    output_degrade_matrix_bg = os.path.join(self.output_location, "8_DPMM_degrade_matrix.bg")
    output_degrade_matrix_hc = os.path.join(self.output_location, "8_DPMM_degrade_matrix.hic")

    # Creating dpmm
    dpmm = Dpmm(4, contact_map, sica, random_degrade_range = [0.01, 0.02], degrade_multiplier = 0.05,
               half_length_bin_interval = [1, 3], value_range = [10e-4, 10e-3], random_range = [10e-8, 10e-7], iteration_multiplier = 10, seed = 123)

    # Performing diagonal degrade
    dpmm.diagonal_degrade()

    # Writing
    io.write(contact_map, output_degrade_matrix_bg, InputFileType.SPARSE)
    io.write(contact_map, output_degrade_matrix_hc, InputFileType.HIC)

    # Termination
    return io, preprocess, sica, goba, dpmm, contact_map


  def introduce_shapes(self, io, preprocess, sica, goba, dpmm, contact_map):

    # Output file
    output_shape_matrix_bg = os.path.join(self.output_location, "9_DPMM_shape_matrix.bg")
    output_shape_matrix_hc = os.path.join(self.output_location, "9_DPMM_shape_matrix.hic")

    # Performing shapes
    dpmm.introduce_shapes()

    # Writing
    io.write(contact_map, output_shape_matrix_bg, InputFileType.SPARSE)
    io.write(contact_map, output_shape_matrix_hc, InputFileType.HIC)

    # Termination
    return io, preprocess, sica, goba, dpmm, contact_map



  ###################################################################################################
  # 12. IFS
  ###################################################################################################

  def main_calculate_ifs(self, io, preprocess, sica, goba, dpmm, contact_map):

    # Output file
    output_ifs_matrix_bg = os.path.join(self.output_location, "10_IFS_calcifs_matrix.bg")
    output_ifs_matrix_hc = os.path.join(self.output_location, "10_IFS_calcifs_matrix.hic")
    output_loop_file_name = os.path.join(self.output_location, "12_IFS_loop_file_name.bg")
    output_matrix_file_name = os.path.join(self.output_location, "12_IFS_matrix_file_name.hic")
    matrix_output_format = InputFileType.HIC

    # Create Ifs
    ifs = Ifs(contact_map, sica, goba, dpmm, io, output_loop_file_name, output_matrix_file_name, matrix_output_format)

    # Calculate Ifs
    ifs.main_calculate_ifs()

    # Writing
    io.write(contact_map, output_ifs_matrix_bg, InputFileType.SPARSE)
    io.write(contact_map, output_ifs_matrix_hc, InputFileType.HIC)

    # Termination
    return io, preprocess, sica, goba, dpmm, ifs, contact_map

  def main_fix_matrix(self, io, preprocess, sica, goba, dpmm, ifs, contact_map):

    # Output file
    output_finalfix_matrix_bg = os.path.join(self.output_location, "11_IFS_finalfix_matrix.bg")
    output_finalfix_matrix_hc = os.path.join(self.output_location, "11_IFS_finalfix_matrix.hic")

    # Calculate Ifs
    ifs.main_fix_matrix(multiplier = 10, min_matrix_threshold = 0.01)

    # Writing
    io.write(contact_map, output_finalfix_matrix_bg, InputFileType.SPARSE)
    io.write(contact_map, output_finalfix_matrix_hc, InputFileType.HIC)

    # Termination
    return io, preprocess, sica, goba, dpmm, ifs, contact_map




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
  test = PipelineTest(input_location, temporary_location, output_location)

  # Preprocessing tests
  print("Preprocessing tests")
  pt1 = time()
  print("Runing convert_to_minimal_resolution")
  io, preprocess, contact_map = test.convert_to_minimal_resolution()
  pt2 = time()
  print("Finished - Time = " + str(round((pt2 - pt1)/60, 2)) + "m")
  print("Runing minimal_res_contact_map")
  io, preprocess, contact_map = test.main_remove_blacklist(io, preprocess, contact_map)
  pt3 = time()
  print("Finished - Time = " + str(round((pt3 - pt2)/60, 2)) + "m")
  print("Runing blackrem_contact_map")
  io, preprocess, contact_map = test.main_void_statistics(io, preprocess, contact_map)
  pt4 = time()
  print("Finished - Time = " + str(round((pt4 - pt3)/60, 2)) + "m\n")

  # Sica tests
  print("Sica tests")
  pt1 = time()
  print("Runing main_calculate_distributions")
  io, preprocess, sica, contact_map = test.main_calculate_distributions(io, preprocess, contact_map)
  pt2 = time()
  print("Finished - Time = " + str(round((pt2 - pt1)/60, 2)) + "m")
  print("Runing main_star_contacts")
  io, preprocess, sica, contact_map = test.main_star_contacts(io, preprocess, sica, contact_map)
  pt3 = time()
  print("Finished - Time = " + str(round((pt3 - pt2)/60, 2)) + "m\n")

  # Goba tests
  print("Goba tests")
  pt1 = time()
  print("Runing main_fill")
  io, preprocess, sica, goba, contact_map = test.main_fill(io, preprocess, sica, contact_map)
  pt2 = time()
  print("Finished - Time = " + str(round((pt2 - pt1)/60, 2)) + "m\n")

  # Dpmm tests
  print("Dpmm tests")
  pt1 = time()
  print("Runing diagonal_degrade")
  io, preprocess, sica, goba, dpmm, contact_map = test.diagonal_degrade(io, preprocess, sica, goba, contact_map)
  pt2 = time()
  print("Finished - Time = " + str(round((pt2 - pt1)/60, 2)) + "m")
  print("Runing introduce_shapes")
  io, preprocess, sica, goba, dpmm, contact_map = test.introduce_shapes(io, preprocess, sica, goba, dpmm, contact_map)
  pt3 = time()
  print("Finished - Time = " + str(round((pt3 - pt2)/60, 2)) + "m\n")

  # Ifs tests
  print("Ifs tests")
  pt1 = time()
  print("Runing main_calculate_ifs")
  io, preprocess, sica, goba, dpmm, ifs, contact_map = test.main_calculate_ifs(io, preprocess, sica, goba, dpmm, contact_map)
  pt2 = time()
  print("Finished - Time = " + str(round((pt2 - pt1)/60, 2)) + "m")
  print("Runing main_fix_matrix")
  io, preprocess, sica, goba, dpmm, ifs, contact_map = test.main_fix_matrix(io, preprocess, sica, goba, dpmm, ifs, contact_map)
  pt3 = time()
  print("Finished - Time = " + str(round((pt3 - pt2)/60, 2)) + "m\n")


