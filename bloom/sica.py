"""
SICA Module
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
from bloom.util import ChromosomeSizes, AuxiliaryFunctions

# External

###################################################################################################
# Basic Objects
###################################################################################################

class Sica():
  """This class represents TODO.

  *Keyword arguments:*

    - argument1 -- Short description. This argument represents a long description. It can be:
      - Possibility 1: A possibility 1.
      - Possibility 2: A possibility 2.

    - argument2 -- Short description. This argument represents a long description. It can be:
      - Possibility 1: A possibility 1.
      - Possibility 2: A possibility 2.
    """

  def __init__(self, matrix, diagonal_avoid_bins):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    
    # Class objects
    self.matrix = matrix
    self.diagonal_avoid_bins = diagonal_avoid_bins

    # Stars
    self.star_level = 1
    self.star_coord_list = []
    self.star_coord_dict = dict()

    # Tads
    self.tad_boundaries_list = []
    self.tad_boundaries_dict = dict()










  def star(self, chromosome, diag_dist_scale, neighbor_dist_scale, max_neighbors, diag_dist_rand_range = [0.0, 0.1], neighbor_dist_rand_range = [0.0, 0.1]):

    # Initialization
    star_score_threshold = self.matrix.get_sparsity_weighted_sum(chrom) # TODO - Check for a better threshold

    # Iterate through matrix
    for key in self.matrix.keys():

      # Check chromosome
      kk = key.split(":")
      chrom = kk[0]
      if(chrom != chromosome): continue

      # Initialization
      p1 = int(kk[1])
      p2 = int(kk[2])
      star_score_to_add = self.matrix.get_sparsity_weighted_sum(chrom) # TODO - Check for a better value to add

      




    # What to do with the diagonal?

    # diag_dist_scale = a number that is going to be divided by the distance to the diagonal. The closest to the diagonal, the higher the value to add in the star

    # neighbor_dist_scale = a number to be divided by the distance to the neighbor. The farther the neighbor, the less signal to add.

    # max_neighbors = maximum 1D neighbor distance to sum to have the value to add

    self.total_bins = dict() # per chromosome
    self.total_1d_bins = dict() # per chromosome
    self.total_zero_bins = dict() # per chromosome
    self.total_nonzero_bins = dict() # per chromosome
    self.total_nonzero_value = dict() # per chromosome
    self.max_value = dict() # per chromosome
    self.min_value = dict() # per chromosome



distance_from_diagonal(self, key)
random.uniform(diag_dist_rand_range[0], diag_dist_rand_range[1])
random.uniform(neighbor_dist_rand_range[0], neighbor_dist_rand_range[1])

  def pretad(self):

    # Create the pre-TADs with different intensities given a distance to the diagonal.
    # The more distant, the less pre-TADs

    # Also make the pre-TADs relying on the existing starred signal density.

















