from __future__ import print_function
"""
Contact Map Module
===================
This module contains the contact map class and has all operations needed for a cohesive matrix
functioning. From small operations such as binning coersions to large statistical flags.

Authors: Eduardo G. Gusmao.

"""

###################################################################################################
# Libraries
###################################################################################################

# Python
import random

# Internal
from bloom.util import AuxiliaryFunctions

# External
import numpy


###################################################################################################
# ContactMap Class
###################################################################################################

class ContactMap():
  """This class represents a contact map.

  EXAMPLE OF MATRIX:

  chr1_matrix = matrix[chr1]
  chr1_matrix[(pos1, pos2)] = 2.0

  *Keyword arguments:*

    - argument1 -- Short description. This argument represents a long description. It can be:
      - Possibility 1: A possibility 1.
      - Possibility 2: A possibility 2.

    - argument2 -- Short description. This argument represents a long description. It can be:
      - Possibility 1: A possibility 1.
      - Possibility 2: A possibility 2.
  """

  def __init__(self, organism, resolution, error_handler, chromosome_sizes, matrix = None, seed = None):
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
    self.organism = organism
    self.resolution = resolution
    self.matrix = matrix # per chromosome.
    self.seed = seed

    # Auxiliary statistics value dictionaries
    self.min_value_diagonal = dict() # per chromosome = Minimum value (> 0) of diagonal.
    self.max_value_diagonal = dict() # per chromosome = Maximum value of diagonal.
    self.total_value_diagonal = dict() # per chromosome = Total number of counts of diagonal.
    self.min_value_no_diagonal = dict() # per chromosome = Minimum value (> 0) excluding diagonal (upper triangle).
    self.max_value_no_diagonal = dict() # per chromosome = Maximum value excluding diagonal (upper triangle).
    self.total_value_no_diagonal = dict() # per chromosome = Total number of counts excluding diagonal (upper triangle).

    # Auxiliary statistics full matrix bin dictionaries
    self.total_bins = dict() # per chromosome = Total number of bins.
    self.total_nonzero_bins = dict() # per chromosome = Total number of bins with counts > 0.
    self.total_zero_bins = dict() # per chromosome = Total number of bins with counts = 0.

    # Auxiliary statistics 1D (row=col) bp/bin dictionaries
    self.total_1d_bp = dict() # per chromosome = Total number of bp of each chromosome = maximum permitted value.
    self.total_1d_bins = dict() # per chromosome = Total number of bins of each chromosome = maximum permitted value.

    # Auxiliary upper triangular matrix dictionaries
    self.total_bins_triangle = dict() # per chromosome = Total number of bins of each map upper triangular matrix without the diagonal bins.
    self.total_nonzero_bins_triangle = dict() # per chromosome = Total number of bins of each map upper triangular matrix without the diagonal bins with counts > 0.
    self.total_zero_bins_triangle = dict() # per chromosome = Total number of bins of each map upper triangular matrix without the diagonal bins with counts = 0.

    # Auxiliary vectors
    self.valid_chromosome_list = []

    # Utilitary objects
    self.error_handler = error_handler
    self.chromosome_sizes = chromosome_sizes

    # Loading blank matrix
    if(self.matrix == None):
      self.load_blank_matrix()


  #############################################################################
  # Load Blank Matrix
  #############################################################################

  def load_blank_matrix(self): # OK
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Initial matrix is a dictionary
    self.matrix = dict()

    # Create a value dictionary per chromosome
    for chromosome in self.chromosome_sizes.chromosome_sizes_list:
      self.matrix[chromosome] = dict()

      # Include valid chromosome
      self.valid_chromosome_list.append(chromosome)

    self.calculate_all_non_value_statistics(empty_matrix = True)


  #############################################################################
  # Resolution & bin/bp Operations
  #############################################################################

  def bp_to_bin(self, i, j = None): # OK
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    new_i = int(numpy.floor(AuxiliaryFunctions.floor_multiple(i, self.resolution) / self.resolution))
    if(j):
      new_j = int(numpy.floor(AuxiliaryFunctions.floor_multiple(j, self.resolution) / self.resolution))
      return new_i, new_j
    else:
      return new_i

  def bin_to_bp(self, i, j = None):  # OK
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    new_i = int(numpy.floor(i * self.resolution))
    if(j):
      new_j = int(numpy.floor(j * self.resolution))
      return new_i, new_j
    else:
      return new_i

  def ceil_bp(self, bp_obj):  # OK
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    if(isinstance(bp_obj, int)):
      return AuxiliaryFunctions.ceil_multiple(bp_obj, self.resolution)
    elif(isinstance(bp_obj, tuple)):
      return ( AuxiliaryFunctions.ceil_multiple(bp_obj[0], self.resolution), AuxiliaryFunctions.ceil_multiple(bp_obj[1], self.resolution) )
    elif(isinstance(bp_obj, float)):
      return int(AuxiliaryFunctions.ceil_multiple(numpy.ceil(bp_obj), self.resolution))
    else:
      pass
      # self.error_handler.throw_error("TODO") # TODO

  def floor_bp(self, bp_obj):  # OK
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    if(isinstance(bp_obj, int)):
      return AuxiliaryFunctions.floor_multiple(bp_obj, self.resolution)
    elif(isinstance(bp_obj, tuple)):
      return ( AuxiliaryFunctions.floor_multiple(bp_obj[0], self.resolution), AuxiliaryFunctions.floor_multiple(bp_obj[1], self.resolution) )
    elif(isinstance(bp_obj, float)):
      return int(AuxiliaryFunctions.floor_multiple(numpy.floor(bp_obj), self.resolution))
    else:
      pass
      # self.error_handler.throw_error("TODO") # TODO


  #############################################################################
  # Matrix Operations
  #############################################################################

  def set(self, chrom, i, j, value): # OK
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    self.matrix[chrom][(i, j)] = value

  def set_from_matrix(self, chromosome, matrix, matrix_type = "numpy_array", storage_type = "upper_triangle"): # OK
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Iterate through matrix's rows
    for row in range(0, self.total_1d_bins[chromosome]):

      # Iterate through matrix's columns
      for col in range(row, self.total_1d_bins[chromosome]):

        # Set value
        if(matrix[row, col] > 0):
          row_bp = self.bin_to_bp(row)
          col_bp = self.bin_to_bp(col)
          self.set(chromosome, row_bp, col_bp, matrix[row, col])

  def add(self, chrom, i, j, value): # OK
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    try:
      self.matrix[chrom][(i, j)] += value
    except Exception:
      self.matrix[chrom][(i, j)] = value

  def get(self, chrom, i, j): # OK
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    try:
      return self.matrix[chrom][(i, j)]
    except Exception:
      return None
      # self.error_handler.throw_error("TODO") # TODO

  def get_full_matrix(self, chromosome, symmetric = True, return_type = "numpy_array"): # OK
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Creating empty matrix
    max_bin = self.total_1d_bins[chromosome]
    full_matrix = numpy.zeros((max_bin, max_bin))

    # Iterating on internal matrix
    for key, value in self.matrix[chromosome].items():
    
      # Binned locations
      row_bin = self.bp_to_bin(key[0])
      col_bin = self.bp_to_bin(key[1])

      # Writing value on cell
      #try:
      full_matrix[row_bin, col_bin] = value
      if(symmetric):
        full_matrix[col_bin, row_bin] = value
      #except IndexError: # Error - TODO
      #  continue

    return full_matrix

  
  #############################################################################
  # Auxiliary Vectors and Dictionaries Operations
  #############################################################################

  def update_valid_chromosome_list(self): # OK
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # New list of valid chromosomes
    new_valid_chromosome_list = []

    # Iterating on each chromosome to check whether there are values for that chromosome
    for chromosome in self.chromosome_sizes.chromosome_sizes_list:

       # Verify if chromosome's matrix exists
       try:
         chrom_matrix = self.matrix[chromosome]

         # Verify if chromosome's matrix is not empty
         if(chrom_matrix):
           new_valid_chromosome_list.append(chromosome)
       except Exception:
         continue

    # Update the list of valid chromosomes
    self.valid_chromosome_list = new_valid_chromosome_list

  def calculate_all_statistics(self): # OK
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # New value dictionaries
    new_min_value_diagonal = dict()
    new_max_value_diagonal = dict()
    new_total_value_diagonal = dict()
    new_min_value_no_diagonal = dict()
    new_max_value_no_diagonal = dict()
    new_total_value_no_diagonal = dict()

    # New bin dictionaries
    new_total_bins = dict()
    new_total_nonzero_bins = dict()
    new_total_zero_bins = dict()
    new_total_1d_bp = dict()
    new_total_1d_bins = dict()
    new_total_bins_triangle = dict()
    new_total_nonzero_bins_triangle = dict()
    new_total_zero_bins_triangle = dict()

    # Iterating on valid chromosomes
    for chromosome in self.valid_chromosome_list:

      # Initialize all dictionaries
      new_min_value_diagonal[chromosome] = numpy.inf
      new_max_value_diagonal[chromosome] = -numpy.inf
      new_total_value_diagonal[chromosome] = 0
      new_min_value_no_diagonal[chromosome] = numpy.inf
      new_max_value_no_diagonal[chromosome] = -numpy.inf
      new_total_value_no_diagonal[chromosome] = 0
      new_total_bins[chromosome] = 0
      new_total_nonzero_bins[chromosome] = 0
      new_total_zero_bins[chromosome] = 0
      new_total_1d_bp[chromosome] = 0
      new_total_1d_bins[chromosome] = 0
      new_total_bins_triangle[chromosome] = 0
      new_total_nonzero_bins_triangle[chromosome] = 0
      new_total_zero_bins_triangle[chromosome] = 0

      # Updating new_total_1d_bp
      new_total_1d_bp[chromosome] = self.floor_bp(self.chromosome_sizes.chromosome_sizes_dictionary[chromosome])

      # Updating new_total_1d_bins
      new_total_1d_bins[chromosome] = self.bp_to_bin(self.chromosome_sizes.chromosome_sizes_dictionary[chromosome])

      # Updating new_total_bins
      new_total_bins[chromosome] = int(new_total_1d_bins[chromosome] ** 2)

      # Updating new_total_bins_triangle
      new_total_bins_triangle[chromosome] = int((new_total_1d_bins[chromosome] * (new_total_1d_bins[chromosome]-1)) / 2)

      # Iterating on matrix
      for key, value in self.matrix[chromosome].items():

        # Only upper triangle and no 0 values
        if(key[0] > key[1] or value <= 0): continue

        # Updating new_total_nonzero_bins
        try:
          new_total_nonzero_bins[chromosome] += 1
        except Exception:
          new_total_nonzero_bins[chromosome] = 1

        # Diagonal vs upper triangle w/o diagonal
        if(key[0] == key[1]):

          # Minimum and maximum diagonal value
          try:
            if(value < new_min_value_diagonal[chromosome]):
              new_min_value_diagonal[chromosome] = value
          except Exception:
            new_min_value_diagonal[chromosome] = value
          try:
            if(value > new_max_value_diagonal[chromosome]):
              new_max_value_diagonal[chromosome] = value
          except Exception:
            new_max_value_diagonal[chromosome] = value

          # Total value of diagonal
          try:
            new_total_value_diagonal[chromosome] += value
          except Exception:
            new_total_value_diagonal[chromosome] = value

        else:

          # Minimum and maximum non-diagonal value
          try:
            if(value < new_min_value_no_diagonal[chromosome]):
              new_min_value_no_diagonal[chromosome] = value
          except Exception:
            new_min_value_no_diagonal[chromosome] = value
          try:
            if(value > new_max_value_no_diagonal[chromosome]):
              new_max_value_no_diagonal[chromosome] = value
          except Exception:
            new_max_value_no_diagonal[chromosome] = value

          # Total non-diagonal value
          try:
            new_total_value_no_diagonal[chromosome] += value
          except Exception:
            new_total_value_no_diagonal[chromosome] = value

          # Updating new_total_nonzero_bins_triangle
          try:
            new_total_nonzero_bins_triangle[chromosome] += 1
          except Exception:
            new_total_nonzero_bins_triangle[chromosome] = 1

      # Correcting new_total_nonzero_bins for full matrix
      new_total_nonzero_bins[chromosome] = int(new_total_nonzero_bins[chromosome] * 2)

      # Updating new_total_zero_bins
      new_total_zero_bins[chromosome] = new_total_bins[chromosome] - new_total_nonzero_bins[chromosome]

      # Updating new_total_zero_bins_triangle
      new_total_zero_bins_triangle[chromosome] = new_total_bins_triangle[chromosome] - new_total_nonzero_bins_triangle[chromosome]
        
    # Update value dictionaries
    self.min_value_diagonal = new_min_value_diagonal
    self.max_value_diagonal = new_max_value_diagonal
    self.total_value_diagonal = new_total_value_diagonal
    self.min_value_no_diagonal = new_min_value_no_diagonal
    self.max_value_no_diagonal = new_max_value_no_diagonal
    self.total_value_no_diagonal = new_total_value_no_diagonal

    # Update bin dictionaries
    self.total_bins = new_total_bins
    self.total_nonzero_bins = new_total_nonzero_bins
    self.total_zero_bins = new_total_zero_bins
    self.total_1d_bp = new_total_1d_bp
    self.total_1d_bins = new_total_1d_bins
    self.total_bins_triangle = new_total_bins_triangle
    self.total_nonzero_bins_triangle = new_total_nonzero_bins_triangle
    self.total_zero_bins_triangle = new_total_zero_bins_triangle

  def calculate_statistics_value(self): # OK
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # New value dictionaries
    new_min_value_diagonal = dict()
    new_max_value_diagonal = dict()
    new_total_value_diagonal = dict()
    new_min_value_no_diagonal = dict()
    new_max_value_no_diagonal = dict()
    new_total_value_no_diagonal = dict()

    # Iterating on valid chromosomes
    for chromosome in self.valid_chromosome_list:

      # Initialize all dictionaries
      new_min_value_diagonal[chromosome] = numpy.inf
      new_max_value_diagonal[chromosome] = -numpy.inf
      new_total_value_diagonal[chromosome] = 0
      new_min_value_no_diagonal[chromosome] = numpy.inf
      new_max_value_no_diagonal[chromosome] = -numpy.inf
      new_total_value_no_diagonal[chromosome] = 0

      # Iterating on matrix
      for key, value in self.matrix[chromosome].items():

        # Only upper triangle and no 0 values
        if(key[0] > key[1] or value <= 0): continue

        # Diagonal vs upper triangle w/o diagonal
        if(key[0] == key[1]):

          # Minimum and maximum diagonal value
          try:
            if(value < new_min_value_diagonal[chromosome]):
              new_min_value_diagonal[chromosome] = value
          except Exception:
            new_min_value_diagonal[chromosome] = value
          try:
            if(value > new_max_value_diagonal[chromosome]):
              new_max_value_diagonal[chromosome] = value
          except Exception:
            new_max_value_diagonal[chromosome] = value

          # Total value of diagonal
          try:
            new_total_value_diagonal[chromosome] += value
          except Exception:
            new_total_value_diagonal[chromosome] = value

        else:

          # Minimum and maximum non-diagonal value
          try:
            if(value < new_min_value_no_diagonal[chromosome]):
              new_min_value_no_diagonal[chromosome] = value
          except Exception:
            new_min_value_no_diagonal[chromosome] = value
          try:
            if(value > new_max_value_no_diagonal[chromosome]):
              new_max_value_no_diagonal[chromosome] = value
          except Exception:
            new_max_value_no_diagonal[chromosome] = value

          # Total non-diagonal value
          try:
            new_total_value_no_diagonal[chromosome] += value
          except Exception:
            new_total_value_no_diagonal[chromosome] = value
        
    # Update value dictionaries
    self.min_value_diagonal = new_min_value_diagonal
    self.max_value_diagonal = new_max_value_diagonal
    self.total_value_diagonal = new_total_value_diagonal
    self.min_value_no_diagonal = new_min_value_no_diagonal
    self.max_value_no_diagonal = new_max_value_no_diagonal
    self.total_value_no_diagonal = new_total_value_no_diagonal

  def calculate_all_non_value_statistics(self, empty_matrix = False): # OK
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # New bin dictionaries
    new_total_bins = dict()
    new_total_nonzero_bins = dict()
    new_total_zero_bins = dict()
    new_total_1d_bp = dict()
    new_total_1d_bins = dict()
    new_total_bins_triangle = dict()
    new_total_nonzero_bins_triangle = dict()
    new_total_zero_bins_triangle = dict()

    # Iterating on valid chromosomes
    for chromosome in self.valid_chromosome_list:

      # Initialize all dictionaries
      new_total_bins[chromosome] = 0
      new_total_nonzero_bins[chromosome] = 0
      new_total_zero_bins[chromosome] = 0
      new_total_1d_bp[chromosome] = 0
      new_total_1d_bins[chromosome] = 0
      new_total_bins_triangle[chromosome] = 0
      new_total_nonzero_bins_triangle[chromosome] = 0
      new_total_zero_bins_triangle[chromosome] = 0

      # Updating new_total_1d_bp
      new_total_1d_bp[chromosome] = self.floor_bp(self.chromosome_sizes.chromosome_sizes_dictionary[chromosome])

      # Updating new_total_1d_bins
      new_total_1d_bins[chromosome] = self.bp_to_bin(self.chromosome_sizes.chromosome_sizes_dictionary[chromosome])

      # Updating new_total_bins
      new_total_bins[chromosome] = int(new_total_1d_bins[chromosome] ** 2)

      # Updating new_total_bins_triangle
      new_total_bins_triangle[chromosome] = int((new_total_1d_bins[chromosome] * (new_total_1d_bins[chromosome]-1)) / 2)

      # Calculate value-based bins only if chromosome is not empty
      if(not empty_matrix):

        # Iterating on matrix
        for key, value in self.matrix[chromosome].items():

          # Only upper triangle and no 0 values
          if(key[0] >= key[1] or value <= 0): continue

          # Updating new_total_nonzero_bins
          try:
            new_total_nonzero_bins[chromosome] += 1
          except Exception:
            new_total_nonzero_bins[chromosome] = 1

          # Updating new_total_nonzero_bins_triangle
          try:
            new_total_nonzero_bins_triangle[chromosome] += 1
          except Exception:
            new_total_nonzero_bins_triangle[chromosome] = 1

        # Correcting new_total_nonzero_bins for full matrix
        new_total_nonzero_bins[chromosome] = int(new_total_nonzero_bins[chromosome] * 2)

        # Updating new_total_zero_bins
        new_total_zero_bins[chromosome] = new_total_bins[chromosome] - new_total_nonzero_bins[chromosome]

        # Updating new_total_zero_bins_triangle
        new_total_zero_bins_triangle[chromosome] = new_total_bins_triangle[chromosome] - new_total_nonzero_bins_triangle[chromosome]

    # Update bin dictionaries
    self.total_bins = new_total_bins
    self.total_nonzero_bins = new_total_nonzero_bins
    self.total_zero_bins = new_total_zero_bins
    self.total_1d_bp = new_total_1d_bp
    self.total_1d_bins = new_total_1d_bins
    self.total_bins_triangle = new_total_bins_triangle
    self.total_nonzero_bins_triangle = new_total_nonzero_bins_triangle
    self.total_zero_bins_triangle = new_total_zero_bins_triangle

  def calculate_statistics_full(self, empty_matrix = False): # OK
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # New full matrix bin dictionaries
    new_total_bins = dict()
    new_total_nonzero_bins = dict()
    new_total_zero_bins = dict()

    # Iterating on valid chromosomes
    for chromosome in self.valid_chromosome_list:

      # Initialize all dictionaries
      new_total_bins[chromosome] = 0
      new_total_nonzero_bins[chromosome] = 0
      new_total_zero_bins[chromosome] = 0

      # Total bins
      new_total_bins[chromosome] = int(self.bp_to_bin(self.chromosome_sizes.chromosome_sizes_dictionary[chromosome]) ** 2)

      # Calculate value-based bins only if chromosome is not empty
      if(not empty_matrix):

        # Iterating on matrix
        for key, value in self.matrix[chromosome].items():

          # Only upper triangle and no 0 values
          if(key[0] >= key[1] or value <= 0): continue

          # Updating new_total_nonzero_bins
          try:
            new_total_nonzero_bins[chromosome] += 1
          except Exception:
            new_total_nonzero_bins[chromosome] = 1

        # Correcting new_total_nonzero_bins for full matrix
        new_total_nonzero_bins[chromosome] = int(new_total_nonzero_bins[chromosome] * 2)

        # Updating new_total_zero_bins
        new_total_zero_bins[chromosome] = new_total_bins[chromosome] - new_total_nonzero_bins[chromosome]

    # Update total bin dictionaries
    self.total_bins = new_total_bins
    self.total_nonzero_bins = new_total_nonzero_bins
    self.total_zero_bins = new_total_zero_bins

  def calculate_statistics_1D(self): # OK
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # New 1D bin dictionaries
    new_total_1d_bp = dict()
    new_total_1d_bins = dict()

    # Iterating on valid chromosomes
    for chromosome in self.valid_chromosome_list:

      # Initialize all dictionaries
      new_total_1d_bp[chromosome] = 0
      new_total_1d_bins[chromosome] = 0

      # Updating total 1D bp
      new_total_1d_bp[chromosome] = self.floor_bp(self.chromosome_sizes.chromosome_sizes_dictionary[chromosome])

      # Updating total 1D bins
      new_total_1d_bins[chromosome] = self.bp_to_bin(self.chromosome_sizes.chromosome_sizes_dictionary[chromosome])

    # Update 1D bin dictionaries
    self.total_1d_bp = new_total_1d_bp
    self.total_1d_bins = new_total_1d_bins

  def calculate_statistics_upper_triangle(self, empty_matrix = False): # OK
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # New full matrix bin dictionaries
    new_total_bins_triangle = dict()
    new_total_nonzero_bins_triangle = dict()
    new_total_zero_bins_triangle = dict()

    # Iterating on valid chromosomes
    for chromosome in self.valid_chromosome_list:

      # Initialize all dictionaries
      new_total_bins_triangle[chromosome] = 0
      new_total_nonzero_bins_triangle[chromosome] = 0
      new_total_zero_bins_triangle[chromosome] = 0

      # Total bins in upper triangle
      total_1d_bp = self.bp_to_bin(self.chromosome_sizes.chromosome_sizes_dictionary[chromosome])
      new_total_bins_triangle[chromosome] = int((total_1d_bp * (total_1d_bp-1)) / 2)

      # Calculate value-based bins only if chromosome is not empty
      if(not empty_matrix):

        # Iterating on matrix
        for key, value in self.matrix[chromosome].items():

          # Only upper triangle and no 0 values
          if(key[0] >= key[1] or value <= 0): continue

          # Updating new_total_nonzero_bins_triangle
          try:
            new_total_nonzero_bins_triangle[chromosome] += 1
          except Exception:
            new_total_nonzero_bins_triangle[chromosome] = 1

        # Updating new_total_zero_bins_triangle
        new_total_zero_bins_triangle[chromosome] = new_total_bins_triangle[chromosome] - new_total_nonzero_bins_triangle[chromosome]

    # Update triangular bin dictionaries
    self.total_bins_triangle = new_total_bins_triangle
    self.total_nonzero_bins_triangle = new_total_nonzero_bins_triangle
    self.total_zero_bins_triangle = new_total_zero_bins_triangle

  """
  def calculate_non_basic_statistics(self):

    # Update the non-basic vetor's values

    # Total bins = Summation of the total number of bins of each chromosome map upper triangular matrix without the diagonal bins
    self.total_bins_triangle = 0

    # Get valid chromosome list
    valid_chromosome_dict = dict()
    for key, value in self.matrix.iteritems():
      chrom = key.split(":")[0]
      valid_chromosome_dict[chrom] = True
    valid_chromosome_list = sorted(valid_chromosome_dict.keys())

    # Iterate on valid chromosome list
    for chrom in valid_chromosome_list:

      if(not valid_chromosome_list[chrom]): continue

      chrom_size = self.chromosome_sizes.chromosome_sizes_dictionary[chrom]
      total_1d_bins = AuxiliaryFunctions.floor_multiple(chrom_size, self.resolution) / self.resolution
      total_bins_triangle = ((total_1d_bins * (total_1d_bins-1))/2) + total_1d_bins

      # Update total number of 1D bins (numer of rows = number of columns)
      try:
        self.total_1d_bins[chrom] = total_1d_bins
      except Exception:
        self.error_handler.throw_error("TODO") # TODO - Error: One or more processes didnt execute correctly.

      # Update total number of bins
      try:
        self.total_bins_triangle[chrom] = total_bins_triangle
      except Exception:
        self.error_handler.throw_error("TODO") # TODO - Error: One or more processes didnt execute correctly.

      # Update total bins = 0
      try:
        self.total_zero_bins[chrom] = self.total_bins_triangle[chrom] - self.total_nonzero_bins[chrom]
      except Exception:
        self.error_handler.throw_error("TODO") # TODO - Error: One or more processes didnt execute correctly.
  """


  #############################################################################
  # Distance Operations
  #############################################################################

  def bin_distance_from_diagonal_manhattan(self, bin_point): # OK
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    return max(bin_point) - min(bin_point)

  def bin_distance_from_diagonal_euclidean(self, bin_point): # OK
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    return int(numpy.ceil((max(bin_point) - min(bin_point)) / 2))

  #############################################################################
  # Sparsity & Statistical Operations
  #############################################################################

  def get_sparsity(self, chromosome): # TODO
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    pass
    # Return sparsity level
    #return self.total_nonzero_bins[chromosome] / self.total_bins[chromosome]
    # Placeholder
    # Future - TODO

  def get_sparsity_weighted_sum(self, chromosome): # TODO
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    pass
    # Return sparsity weighted sum
    #return self.total_value[chromosome] * (self.total_nonzero_bins[chromosome] / self.total_bins[chromosome])
    # Placeholder
    # Future - TODO

  def standardize(self): # TODO
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    pass 

    # Iterate over matrix
    #for key, value in self.matrix.iteritems(): # TODO
  
    #  # Get chromosome
    #  chrom = key.split(":")[0]
    #  newvalue = (float(value) - float(self.min_value_diagonal[chrom])) / (float(self.max_value_diagonal[chrom]) - float(self.min_value_diagonal[chrom]))
      
    #  # Update values
    #  try:
    #    self.matrix[key] = newvalue
    #  except Exception:
    #    raise
    #    # self.error_handler.throw_error("TODO") # TODO


  #############################################################################
  # Multi-matrix Operations
  #############################################################################

  def compare_matrices(self, chromosome, contact_map, similarity_degree = 0.1):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Get input keys and compare
    mykeys = sorted(self.matrix[chromosome].keys())
    cpkeys = sorted(contact_map.matrix[chromosome].keys())
    if(mykeys != cpkeys):
      return False

    # Compare each element
    for k in mykeys:
      myvalue = self.matrix[chromosome][k]
      cpvalue = contact_map.matrix[chromosome][k]
      if((myvalue > (cpvalue + (similarity_degree * cpvalue))) or (myvalue < (cpvalue - (similarity_degree * cpvalue)))):
        return False

    # Return
    return True

  def match_subset(self, chromosome, contact_map, similarity_degree = 0.1):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Match object
    match = True

    # Iterating on matrix
    for key, value in contact_map.matrix[chromosome].items():
      try:
        myvalue = self.matrix[chromosome][key]
      except Exception:
        #print("No value")
        #print("\t".join([str(e) for e in [chromosome, key[0], key[0] + self.resolution, key[1], key[1] + self.resolution, value]]))
        match = False
        break
      if((myvalue > (value + (similarity_degree * value))) or (myvalue < (value - (similarity_degree * value)))):
        #print("No match")
        #print("\t".join([str(e) for e in [chromosome, key[0], key[0] + self.resolution, key[1], key[1] + self.resolution, value]]))
        match = False
        break

    # Return
    return match

  #############################################################################
  # Drop-in / Drop-out Operations
  #############################################################################

  def delete_strips(self, chromosome, start, end):
    keys = []
    for key in self.matrix[chromosome]:
      if start <= key[0] < end or start <= key[1] < end:
        keys.append(key)
    for key in keys:
      self.matrix[chromosome].pop(key)


  def delete_signal(self,chromosome,i,j):
    try:
      self.matrix[chromosome].pop((i,j))
    except Exception:
      print ("Please specify a valid coordinate divisible by the resolution!")

  def delete_signal_randomly_select_bins(self,chromosome,percentage):
    all_bin_pairs = list(self.matrix[chromosome].keys())
    ## get a list of N bin pairs to delete
    N = int(self.total_bins_triangle[chromosome] * percentage)
    bin_pairs = random.sample(all_bin_pairs, N)
    for bin_pair in bin_pairs:
      self.delete_signal(chromosome,*bin_pair) 

  def add_signal_randomly_select_bins(self,chromosome,percentage,value):
    all_bin_pairs = list(self.matrix[chromosome].keys())
    ## get a list of N bin pairs to delete
    N = int(self.total_bins_triangle[chromosome] * percentage)
    bin_pairs = random.sample(all_bin_pairs, N)
    for bin_pair in bin_pairs:
      self.add(chromosome,*bin_pair, value)

  def main_delete_signal_randomly_select_bins(self,chromosomes,percentage):
    for chromosome in chromosomes:
      self.delete_signal_randomly_select_bins(chromosome,percentage)

  def main_add_signal_randomly_select_bins(self,chromosomes,percentage,value):
    for chromosome in chromosomes:
      self.add_signal_randomly_select_bins(chromosome,percentage,value)

