from __future__ import print_function
"""
Main Module
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
import time

# Internal
from bloom.util import ErrorHandler, AuxiliaryFunctions
from bloom.arguments import ArgumentParser
from bloom.io import InputFileType, IO
from bloom.contact_map import ContactMap
from bloom.preprocess import Preprocess
from bloom.sica import Sica
from bloom.goba import Goba
from bloom.dpmm import Dpmm
from bloom.ifs import Ifs
from bloom.barcode import Barcode

# External


###################################################################################################
# Tool Execution
###################################################################################################

def main():
  """Entry point to the execution of Bloom.
    
  *Keyword arguments:*
    
    - None.
    
  *Return:*
    
    - None.
  """

  ###############################################################################
  # Initialization
  ###############################################################################

  # Initializing ErrorHandler
  error_handler = ErrorHandler()


  ###############################################################################
  # Parameter Handling
  ###############################################################################

  # Initializing argument parser
  argument_parser = ArgumentParser()
  args = argument_parser.arguments
  opts = argument_parser.options

  # Arguments
  input_matrix = os.path.abspath(os.path.expanduser(args[0]))
  output_matrix = os.path.abspath(os.path.expanduser(args[1]))
  output_contacts = os.path.abspath(os.path.expanduser(args[2]))

  # Options
  organism = opts.organism
  resolution = opts.resolution
  input_type = opts.input_type
  cores = opts.cores
  temporary_location = os.path.abspath(os.path.expanduser(opts.temporary_location))
  output_type = opts.output_type
  output_resolution = opts.output_resolution
  avoid_distance = opts.avoid_distance
  seed = opts.seed

  # Hidden Options
  pre_min_contig_removed_bins = opts.pre_min_contig_removed_bins
  pre_remove_threshold = opts.pre_remove_threshold
  sica_pvalue_threshold = opts.sica_pvalue_threshold
  sica_bottom_bin_ext_range = opts.sica_bottom_bin_ext_range.split(",")
  sica_left_bin_ext_range = opts.sica_left_bin_ext_range.split(",")
  sica_right_bin_ext_range = opts.sica_right_bin_ext_range.split(",")
  sica_top_bin_ext_range = opts.sica_top_bin_ext_range.split(",")
  dpmm_random_degrade_range = opts.dpmm_random_degrade_range.split(",")
  dpmm_degrade_multiplier = opts.dpmm_degrade_multiplier
  dpmm_half_length_bin_interval = opts.dpmm_half_length_bin_interval.split(",")
  dpmm_value_range = opts.dpmm_value_range.split(",")
  dpmm_random_range = opts.dpmm_random_range.split(",")
  dpmm_iteration_multiplier = opts.dpmm_iteration_multiplier


  ###############################################################################
  # Reading user's contact map
  ###############################################################################

  # Create contact map with input file
  io_instance = IO(input_matrix, temporary_location, organism, cores, input_resolution = resolution, input_file_type = input_type)
  contact_map = io_instance.read()
  contact_map.calculate_all_statistics()


  ###############################################################################
  # Barcode
  ###############################################################################

  # Verify barcode
  barcode_instance = Barcode(cores, organism, contact_map, temporary_location)
  final_contact_map, final_loop_list = barcode_instance.main_guide()

  # Print matrix and loop list
  if(final_matrix and final_loop):
    io_instance.write(final_contact_map, output_matrix, output_format = output_type)
    io_instance.write_loop_list(final_loop_list, output_contacts)
    sys.exit(0)


  ###############################################################################
  # Preprocess
  ###############################################################################

  # Preprocess
  preprocess_instance = Preprocess(cores, contact_map, minimal_resolution = output_resolution,
                                   min_contig_removed_bins = pre_min_contig_removed_bins, remove_threshold = pre_remove_threshold)
  contact_map = preprocess_instance.convert_to_minimal_resolution(recalculate_statistics = True)
  preprocess_instance.check_sparsity()
  preprocess_instance.remove_blacklist()
  preprocess_instance.main_void_statistics()
  contact_map.calculate_all_statistics()


  ###############################################################################
  # ICA / SICA
  ###############################################################################

  # SICA
  sica_instance = Sica(cores, contact_map, avoid_distance, removed_dict = preprocess_instance.removed_dict, pvalue_threshold = sica_pvalue_threshold)
  sica_instance.main_calculate_distributions()
  sica_instance.main_star_contacts(bottom_bin_ext_range = bottom_bin_ext_range, left_bin_ext_range = left_bin_ext_range,
                                   right_bin_ext_range = right_bin_ext_range, top_bin_ext_range = top_bin_ext_range)
  contact_map.calculate_all_statistics()


  ###############################################################################
  # GOBA
  ###############################################################################

  # GOBA
  goba_instance = Goba(contact_map, sica_instance)
  goba_instance.main_fill()
  contact_map.calculate_all_statistics()


  ###############################################################################
  # DPMM
  ###############################################################################

  # DPMM
  dpmm_instance = Dpmm(cores, contact_map, avoid_distance, preprocess_instance.removed_dict, random_degrade_range = dpmm_random_degrade_range,
                       degrade_multiplier = dpmm_degrade_multiplier, half_length_bin_interval = dpmm_half_length_bin_interval, 
                       value_range = dpmm_value_range, random_range = dpmm_random_range, iteration_multiplier = dpmm_iteration_multiplier, seed = seed)
  dpmm_instance.diagonal_degrade()
  dpmm_instance.introduce_shapes()
  contact_map.calculate_all_statistics()


  ###############################################################################
  # IFS
  ###############################################################################

  # IFS
  ifs_instance = Ifs(contact_map, sica_instance, goba_instance, dpmm_instance, io_instance, output_contacts, output_matrix, output_type)
  loop_list = ifs_instance.calculate_ifs()
  ifs_instance.fix_matrix()


  ###############################################################################
  # Termination
  ###############################################################################

  # Removing temporary file
  AuxiliaryFunctions.remove_folder_files(temporary_location)

  # Write final matrix and loop list
  # io_instance.write(contact_map, output_matrix, output_format = output_type)
  # io_instance.write_loop_list(loop_list, output_contacts)



