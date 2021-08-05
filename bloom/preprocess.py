"""
Preprocess Module
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
import random
import codecs
import optparse
import traceback
import subprocess
import configparser
import multiprocessing

# Internal
from bloom.contact_map import ContactMap
from bloom.util import ErrorHandler, ExcList, AuxiliaryFunctions

# External
import numpy as np
import scipy.interpolate
import scipy.ndimage

###################################################################################################
# Basic Objects
###################################################################################################

class Preprocess():
  """This class represents TODO.

  *Keyword arguments:*

    - argument1 -- Short description. This argument represents a long description. It can be:
      - Possibility 1: A possibility 1.
      - Possibility 2: A possibility 2.

    - argument2 -- Short description. This argument represents a long description. It can be:
      - Possibility 1: A possibility 1.
      - Possibility 2: A possibility 2.
  """

  def __init__(self, ncpu, input_contact_map, minimal_resolution = 1000, min_contig_removed_bins = 5, remove_threshold = 1, seed = None):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Seed
    if(seed):
      random.seed(seed)

    # Main objects
    self.input_contact_map = input_contact_map
    self.minimal_resolution = minimal_resolution
    self.min_contig_removed_bins = min_contig_removed_bins
    self.remove_threshold = remove_threshold
    self.seed = seed

    # Statistics dictionaries
    self.removed_dict = dict() # Points removed because falls into a 0 contig or blacklist. For every chromosome -> for every row/col bin -> True or error.
    self.colsum_dict = dict() # For every chromosome -> for every col bin -> total sum of that col signal

    # Auxiliary objects
    self.ncpu = ncpu
    self.process_queue = []

    # Utilitary objects
    self.error_handler = ErrorHandler()
    self.exclist_handler = ExcList(self.input_contact_map.organism, self.minimal_resolution)


  #############################################################################
  # Reshape
  #############################################################################

  def convert_to_minimal_resolution(self, recalculate_statistics = True):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Create new contact map
    new_contact_map = ContactMap(organism = self.input_contact_map.organism, resolution = self.minimal_resolution, matrix = None, seed = self.seed)
    new_contact_map.valid_chromosome_list = [chromosome for chromosome in self.input_contact_map.valid_chromosome_list]

    # Iterating on valid chromosomes
    for chromosome in new_contact_map.valid_chromosome_list:

      # Add reshape process to list
      self.reshape(chromosome, new_contact_map)
      # self.add_reshape(chromosome, new_contact_map)

    # Execute reshape
    # self.run_reshape()

    # Recalculate statistics
    if(recalculate_statistics):
      new_contact_map.calculate_all_statistics()

    # Return new contact map
    return new_contact_map

  def congrid(self, a, newdims, method = 'linear', centre = False, minusone = False): # TODO - Re-style
    """Arbitrary resampling of source array to new dimension sizes.
    Currently only supports maintaining the same number of dimensions.
    To use 1-D arrays, first promote them to shape (x,1).
    
    Uses the same parameters and creates the same co-ordinate lookup points
    as IDL''s congrid routine, which apparently originally came from a VAX/VMS
    routine of the same name.

    method:
    neighbour - closest value from original data
    nearest and linear - uses n x 1-D interpolations using
                         scipy.interpolate.interp1d
    (see Numerical Recipes for validity of use of n 1-D interpolations)
    spline - uses ndimage.map_coordinates

    centre:
    True - interpolation points are at the centres of the bins
    False - points are at the front edge of the bin

    minusone:
    For example- inarray.shape = (i,j) & new dimensions = (x,y)
    False - inarray is resampled by factors of (i/x) * (j/y)
    True - inarray is resampled by(i-1)/(x-1) * (j-1)/(y-1)
    This prevents extrapolation one element beyond bounds of input array.
    """

    if not a.dtype in [np.float64, np.float32]:
      a = np.cast[float](a)

    m1 = np.cast[int](minusone)
    ofs = np.cast[int](centre) * 0.5
    old = np.array( a.shape )
    ndims = len( a.shape )
    if len( newdims ) != ndims:
      print("[congrid] dimensions error. This routine currently only support rebinning to the same number of dimensions.")
      return None
    newdims = np.asarray( newdims, dtype=int )
    dimlist = []

    if method == 'neighbour':
      for i in range( ndims ):
        base = np.indices(newdims)[i]
        dimlist.append( (old[i] - m1) / (newdims[i] - m1) * (base + ofs) - ofs )
      cd = np.array( dimlist ).round().astype(int)
      newa = a[tuple( cd )]
      return newa

    elif method in ['nearest','linear']:
      # calculate new dims
      for i in range( ndims ):
        base = np.arange( newdims[i] )
        dimlist.append( (old[i] - m1) / (newdims[i] - m1) * (base + ofs) - ofs )

      # specify old dims
      olddims = [np.arange(i, dtype = np.float) for i in list( a.shape )]

      # first interpolation - for ndims = any
      mint = scipy.interpolate.interp1d( olddims[-1], a, kind=method, fill_value="extrapolate" )
      newa = mint( dimlist[-1] )

      trorder = [ndims - 1] + list(range( ndims - 1 ))
      for i in range( ndims - 2, -1, -1 ):
        newa = newa.transpose( trorder )
        mint = scipy.interpolate.interp1d( olddims[i], newa, kind=method, fill_value="extrapolate" )
        newa = mint( dimlist[i] )

      if ndims > 1:
        # need one more transpose to return to original dimensions
        newa = newa.transpose( trorder )

      return newa
    elif method in ['spline']:
      oslices = [ slice(0,j) for j in old ]
      oldcoords = np.ogrid[oslices]
      nslices = [ slice(0,j) for j in list(newdims) ]
      newcoords = np.mgrid[nslices]

      newcoords_dims = list(range(newcoords.ndim))
      #make first index last
      newcoords_dims.append(newcoords_dims.pop(0))
      newcoords_tr = newcoords.transpose(newcoords_dims)
      # makes a view that affects newcoords

      newcoords_tr = newcoords_tr + ofs

      deltas = (np.asarray(old) - m1) / (newdims - m1)
      newcoords_tr = newcoords_tr * deltas

      newcoords_tr = newcoords_tr - ofs

      newa = scipy.ndimage.map_coordinates(a, newcoords)
      return newa
    else:
      print("Congrid error: Unrecognized interpolation type. Currently only \'neighbour\', \'nearest\',\'linear\', and \'spline\' are supported.")
      return None

  def reshape(self, chromosome, new_contact_map):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
 
    # Fetching chromosome's old matrix
    old_matrix = self.input_contact_map.get_full_matrix(chromosome, symmetric = True, return_type = "numpy_array")

    # Reshaping using congrid
    total_chromosome_matrix_bins = new_contact_map.total_1d_bins[chromosome]
    new_matrix = self.congrid(old_matrix, [total_chromosome_matrix_bins, total_chromosome_matrix_bins], method = 'linear', centre = True, minusone = False)

    # Loading new matrix into the new contact map
    new_contact_map.set_from_matrix(chromosome, new_matrix, matrix_type = "numpy_array", storage_type = "upper_triangle")

    # Successful execution
    return True

  def add_reshape(self, chromosome, new_contact_map):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Append job to queue
    self.process_queue.append((chromosome, new_contact_map))

  def run_reshape(self):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Execute job queue
    pool = multiprocessing.Pool(self.ncpu)
    reshape_process_output = pool.starmap(self.reshape, [arguments for arguments in self.process_queue])
    pool.close()
    pool.join()

    # Clean queue
    pool = None
    self.process_queue = []
    gc.collect()

    # Check execution status
    successful_execution = True
    for cp in dump_process_output:
      if(not cp):
        successful_execution = False
        self.error_handler.throw_error("TODO") # TODO - Error: One or more processes didnt execute correctly.

    # Return success status
    return successful_execution

  #############################################################################
  # Sparsity
  #############################################################################

  def check_sparsity(self):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    pass
    # Placeholder
    # Future - TODO

  #############################################################################
  # Blacklist
  #############################################################################

  def main_remove_blacklist(self, contact_map):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Get valid chromosome list
    valid_chromosome_list = contact_map.valid_chromosome_list

    # Iterating on valid chromosomes
    for chromosome in valid_chromosome_list:

      # Add chromosome to dictionaries
      try:
        self.removed_dict[chromosome]
      except Exception:
        self.removed_dict[chromosome] = dict()

      # Add remove_from_map job to the queue
      self.remove_blacklist(chromosome, contact_map)
      # self.add_remove_blacklist(chromosome, contact_map)

    # Run remove_from_map
    # self.run_remove_blacklist()

  def remove_blacklist(self, chromosome, contact_map):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Verify if chromosome exists
    try:
      self.exclist_handler.exclude_dictionary[chromosome]
    except Exception:
      return

    # Iterating on matrix
    for key in self.exclist_handler.exclude_dictionary[chromosome].keys():

      # Binary version of key
      kbin = contact_map.bp_to_bin(key)

      # Putting key on removed_dict
      self.removed_dict[chromosome][kbin] = True

  def add_remove_blacklist(self, chromosome, contact_map):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Append job to queue
    self.process_queue.append((chromosome, contact_map))

  def run_remove_blacklist(self):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    
    # Execute job queue
    pool = multiprocessing.Pool(self.ncpu)
    pool.starmap(self.remove_blacklist, [arguments for arguments in self.process_queue])
    pool.close()
    pool.join()

    # Clean queue
    pool = None
    self.process_queue = []
    gc.collect()


  #############################################################################
  # Row/Col Iteration to remove void intervals
  #############################################################################

  def main_void_statistics(self, contact_map):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Get valid chromosome list
    valid_chromosome_list = contact_map.valid_chromosome_list
    
    # Iterating on valid chromosomes
    for chromosome in valid_chromosome_list:

      # Add chromosome to dictionaries
      self.colsum_dict[chromosome] = [0.0] * contact_map.total_1d_bins[chromosome]

      # Add colsum job to the queue
      self.colsum(chromosome, contact_map)
      # self.add_colsum(chromosome, contact_map)

    # Run colsum
    #self.run_colsum()

    # Iterating on valid chromosomes
    for chromosome in valid_chromosome_list:

      # Add chromosome to dictionaries
      try:
        self.removed_dict[chromosome]
      except Exception:
        self.removed_dict[chromosome] = dict()

      # Add remove rowcol job to the queue
      self.remove_rowcol(chromosome)
      # self.add_remove_rowcol(chromosome)

    # Run remove_rowcol
    # self.run_remove_rowcol()

    # Iterating on valid chromosomes
    for chromosome in valid_chromosome_list:

      # Add remove_from_map job to the queue
      self.remove_from_map(chromosome, contact_map)
      # self.add_remove_from_map(chromosome, contact_map)

    # Run remove_from_map
    # self.run_remove_from_map()

  def colsum(self, chromosome, contact_map):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Iterating on matrix
    for key, value in contact_map.matrix[chromosome].items():

      # Binary version of keys
      brow = contact_map.bp_to_bin(key[0])
      bcol = contact_map.bp_to_bin(key[1])

      # Heuristic rule of #bins away from diagonal == # contig bins lower threshold to remove
      #if(abs(bcol - brow) <= self.min_contig_removed_bins): continue

      # Check if value is above threshold
      if(value > self.remove_threshold):

        # Add amount to dictionary
        self.colsum_dict[chromosome][brow] += value
        self.colsum_dict[chromosome][bcol] += value

  def add_colsum(self, chromosome, contact_map):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Append job to queue
    self.process_queue.append((chromosome, contact_map))

  def run_colsum(self):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    
    # Execute job queue
    pool = multiprocessing.Pool(self.ncpu)
    pool.starmap(self.colsum, [arguments for arguments in self.process_queue])
    pool.close()
    pool.join()

    # Clean queue
    pool = None
    self.process_queue = []
    gc.collect()

  def remove_rowcol(self, chromosome):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Help flag
    zero_track = []

    # Iterating on colsum dictionary
    for i in range(0, len(self.colsum_dict[chromosome])):

      # Fetching column count
      colcount = self.colsum_dict[chromosome][i]

      # Column with signal found
      if(colcount > 0):

        # Enough continuous cols
        if(len(zero_track) >= self.min_contig_removed_bins):

          # Put all contiguous columns in removed dictionary
          for z in zero_track:
            self.removed_dict[chromosome][z] = True

        # Re-set zero_track
        zero_track = []

      # Column with signal not found
      else:
        zero_track.append(i)

    # LAST VERIFICATION -> Enough continuous cols
    if(len(zero_track) >= self.min_contig_removed_bins):

      # LAST VERIFICATION -> Put all contiguous columns in removed dictionary
      for z in zero_track:
        self.removed_dict[chromosome][z] = True

  def add_remove_rowcol(self, chromosome):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Append job to queue
    self.process_queue.append((chromosome))

  def run_remove_rowcol(self):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    
    # Execute job queue
    pool = multiprocessing.Pool(self.ncpu)
    pool.starmap(self.remove_rowcol, [arguments for arguments in self.process_queue])
    pool.close()
    pool.join()

    # Clean queue
    pool = None
    self.process_queue = []
    gc.collect()


  #############################################################################
  # Remove from contact map
  #############################################################################

  def remove_from_map(self, chromosome, contact_map):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # List of keys to remove from dictionary
    removed_keys = []

    # Iterating on matrix
    for key, value in contact_map.matrix[chromosome].items():

      # Binary version of keys
      brow = contact_map.bp_to_bin(key[0])
      bcol = contact_map.bp_to_bin(key[1])

      # Check if entry needs to be removed
      try:
        self.removed_dict[chromosome][brow]
        removed_keys.append(key)
        continue
      except Exception:
        pass
      try:
        self.removed_dict[chromosome][bcol]
        removed_keys.append(key)
        continue
      except Exception:
        pass

    # Remove values from dictionary
    for key in removed_keys:
      contact_map.matrix[chromosome].pop(key)

  def add_remove_from_map(self, chromosome, contact_map):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Append job to queue
    self.process_queue.append((chromosome, contact_map))

  def run_remove_from_map(self):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    
    # Execute job queue
    pool = multiprocessing.Pool(self.ncpu)
    pool.starmap(self.remove_from_map, [arguments for arguments in self.process_queue])
    pool.close()
    pool.join()

    # Clean queue
    pool = None
    self.process_queue = []
    gc.collect()


  #############################################################################
  # OLD
  #############################################################################

  """
  def remove_rowcol(self, chromosome, contact_map):

    # Main loop corresponding to the current row/col being analyzed
    for k in range(self.avoid_distance + 1, contact_map.max_bin[chromosome] + 1):

      # Summation of the row/col
      rowcol_k_summation = 0.0

      # Row loop until avoid region is reached
      row = 0
      for r in range(0, contact_map.max_bin[chromosome] - self.avoid_distance):

        # Row walk stop criteria - when it reaches the avoided region
        if((k - r) <= self.avoid_distance):
          break
        else:
          row = r

        # TODO - Update summation
        rowcol_k_summation += #  if contact_map[row, k] has value, put +1.0

        # TODO - Statistics

      # Col loop until the end of the matrix
      for col in range(k + 1, contact_map.max_bin[chromosome] + 1):

        # TODO - Update summation
        rowcol_k_summation += # if contact_map[row, k] has value, put +1.0

        # TODO - Statistics
      
      # Calculation percentage and verifying whether to remove the row/col
      rowcol_k_summation = rowcol_k_summation / (contact_map.max_bin[chromosome])
      if(rowcol_k_summation <= remove_threshold):
        self.removed_dict[chromosome][k] = True

      # Placeholder

    # Placeholder

  """


