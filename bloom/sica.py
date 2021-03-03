from __future__ import print_function
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
import warnings
import traceback
import subprocess
import configparser
import multiprocessing

# Internal
from bloom.contact_map import ContactMap
from bloom.util import ErrorHandler, AuxiliaryFunctions

# External
import numpy as np
import scipy
import scipy.stats as st

###################################################################################################
# Sica Class
###################################################################################################

class SicaDist():
  """This class represents TODO.

  *Keyword arguments:*

    - argument1 -- Short description. This argument represents a long description. It can be:
      - Possibility 1: A possibility 1.
      - Possibility 2: A possibility 2.

    - argument2 -- Short description. This argument represents a long description. It can be:
      - Possibility 1: A possibility 1.
      - Possibility 2: A possibility 2.
  """

  def __init__(self, contact_map, avoid_distance):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    
    # Fetching distance in bins given resolution
    avoid_distance_bp_min, avoid_distance_bp_max = contact_map.bp_to_bin(0, avoid_distance)
    t5_distance_bp_min, t5_distance_bp_max = contact_map.bp_to_bin(avoid_distance, 250000)
    t4_distance_bp_min, t4_distance_bp_max = contact_map.bp_to_bin(250000, 500000)
    t3_distance_bp_min, t3_distance_bp_max = contact_map.bp_to_bin(500000, 1000000)
    t2_distance_bp_min, t2_distance_bp_max = contact_map.bp_to_bin(1000000, 2000000)
    t1_distance_bp_min, t1_distance_bp_max = contact_map.bp_to_bin(2000000, 3000000)
    c3_distance_bp_min, c3_distance_bp_max = contact_map.bp_to_bin(3000000, 5000000)
    c2_distance_bp_min, c2_distance_bp_max = contact_map.bp_to_bin(5000000, 10000000)
    c1_distance_bp_min, c1_distance_bp_max = contact_map.bp_to_bin(10000000, 20000000)

    # Fixed distances
    self.A = (avoid_distance_bp_min, avoid_distance_bp_max)
    self.T5 = (t5_distance_bp_min, t5_distance_bp_max) # T5 TAD - Minimum
    self.T4 = (t4_distance_bp_min, t4_distance_bp_max) # T4 TAD
    self.T3 = (t3_distance_bp_min, t3_distance_bp_max) # T3 TAD - Average
    self.T2 = (t2_distance_bp_min, t2_distance_bp_max) # T2 TAD -> Small compartment
    self.T1 = (t1_distance_bp_min, t1_distance_bp_max) # T1 TAD - Maximum -> Small compartment
    self.C3 = (c3_distance_bp_min, c3_distance_bp_max) # C1 - Medium compartment
    self.C2 = (c2_distance_bp_min, c2_distance_bp_max) # C2 - Large compartment
    self.C1 = (c1_distance_bp_min, c1_distance_bp_max) # C3 - Very large compartment

    # Distance list
    self.sica_dist_dict = dict([("a", self.A), ("t5", self.T5), ("t4", self.T4), ("t3", self.T3), ("t2", self.T2), 
                                ("t1", self.T1), ("c3", self.C3), ("c2", self.C2), ("c1", self.C1)])


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

  def __init__(self, ncpu, contact_map, avoid_distance, removed_dict = None, pvalue_threshold = 0.95):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    
    # Class objects
    self.contact_map = contact_map
    self.avoid_distance = avoid_distance
    self.pvalue_threshold = pvalue_threshold
    self.removed_dict = removed_dict

    # Annotation dictionary
    self.annotation_dictionary = dict() # Same as contact_map matrix, but with an extra annotation flag:
    # R = Points removed because they fall into a 0 contig or blacklist.
    # A = Points to avoid because they fall in 0-avoid_distance
    # T1, T2, T3, T4, T5 = Tad hierarchy levels depending on point's distance to diagonal.
    # C1, C2, C3 = Compartment hierarchy levels depending on point's distance to diagonal.
    # O = Points farther from diagonal than the biggest compartment length (C1)
    # Upper case letters are important points. Lower case letters are less important points.

    # Auxiliary dictionaries
    self.distribution_dictionary = dict() # per chromosome / per SicaDist -> Contains all the matrix's signal within that specific SicaDist
    self.dist_to_diag_dictionary = dict() # per chromosome / per regular matrix key -> Manhattan distance to the diagonal
    self.pvalue_dictionary = dict() # per chromosome / per SicaDist -> Contains [name of the fitted distribution (FD), parameters of FD, value at pvalue_threshold (given FD)]
    self.significant_values_dictionary = dict() # per chromosome -> All significant (peaks) values of the matrix, i.e. >= value at pvalue_threshold (given FD)

    # Auxiliary distributions
    self.distribution_list = [st.alpha, st.beta, st.burr, st.dgamma, st.dweibull, st.erlang, st.expon, st.exponpow, st.genexpon,
                              st.gilbrat, st.gumbel_r, st.invweibull, st.kstwobign, st.levy, st.ncx2, st.wald, st.weibull_min]

    # Auxiliary parameters
    self.ncpu = ncpu
    self.process_queue = []

    # Utilitary objects
    self.error_handler = ErrorHandler()
    self.dist_handler = SicaDist(contact_map, avoid_distance)


  #############################################################################
  # Annotation
  #############################################################################

  # Each region distance is fit into a distribution, then a p-value is calculated
  # and a threshold is set to mark a point as a peak if p-value >= threshold.
  # Then, depending on the distance from diagonal, annotate each point with
  # a annotation letter: R, A, O, Tx, Cx.

  def main_calculate_distributions(self):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Iterating on valid chromosomes - Calculate histograms
    for chromosome in self.contact_map.valid_chromosome_list:

      # Put chromosome in distribution_dictionary, annotation dictionary and dist_to_diag_dictionary
      self.distribution_dictionary[chromosome] = dict()
      self.annotation_dictionary[chromosome] = dict()
      self.dist_to_diag_dictionary[chromosome] = dict()

      # Add histogram calculation job to queue
      self.add_calculate_histogram(chromosome)

    # Run histogram calculation jobs
    self.run_calculate_histogram()

    # Iterating on valid chromosomes - Calculate pvalues
    for chromosome in self.contact_map.valid_chromosome_list:

      # Put chromosome in pvalue_dictionary
      self.pvalue_dictionary[chromosome] = dict()

      # Iterating on SicaDist
      for skey, svalue in self.dist_handler.sica_dist_dict.items():

        # Add p-value calculation job to queue
        self.add_calculate_pvalues(chromosome, svalue)

    # Run p-value calculation jobs
    self.run_calculate_pvalues()

  def calculate_histogram(self, chromosome):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Iterating on matrix
    for key, value in self.contact_map.matrix[chromosome].items():

      # Check distance to diagonal
      min_bin_key, max_bin_key = self.contact_map.bp_to_bin(key[0], key[1])
      bin_key = (min_bin_key, max_bin_key)
      distance_to_diag = self.contact_map.bin_distance_from_diagonal_manhattan(bin_key)

      # Add distance to diagonal to dictionary
      self.dist_to_diag_dictionary[chromosome][key] = distance_to_diag

      # Default annotation - Outside any distance
      annotation = "o"

      # Attribute distance - Removed dict
      try:
        self.removed_dict[chromosome][bin_key[0]]
        annotation = "r"
      except Exception:
        pass
      try:
        self.removed_dict[chromosome][bin_key[1]]
        annotation = "r"
      except Exception:
        pass

      # Attribute distance - Sica dist dict
      else:

        # Iterating on SicaDist distances
        for skey, svalue in self.dist_handler.sica_dist_dict.items():

          if(svalue[0] <= distance_to_diag < svalue[1]):
            annotation = skey
            try:
              self.distribution_dictionary[chromosome][svalue].append(value)
            except Exception:
              self.distribution_dictionary[chromosome][svalue] = [value]
          break

      # Put value in the annotation dictionary
      self.annotation_dictionary[chromosome][key] = annotation


  def add_calculate_histogram(self, chromosome):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Append job to queue
    self.process_queue.append((chromosome))

  def run_calculate_histogram(self):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Execute job queue
    pool = multiprocessing.Pool(self.ncpu)
    pool.starmap(self.calculate_histogram, [arguments for arguments in self.process_queue])
    pool.close()
    pool.join()

    # Clean queue
    pool = None
    self.process_queue = []
    gc.collect()

  def best_fit_distribution(data, bins = 100):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Histogram from original data
    y, x = np.histogram(data, bins=bins, density=True)
    x = (x + np.roll(x, -1))[:-1] / 2.0

    # Best parameters initialization
    best_distribution = st.norm
    best_params = (0.0, 1.0)
    best_sse = np.inf
    best_pvalue = 0.9

    # Estimate distribution parameters from data
    for distribution in self.distribution_list:

      # Try to fit the distribution
      try:

        # Ignore warnings from data that can't be fit
        with warnings.catch_warnings():
          warnings.filterwarnings('ignore')

          # Fit distribution to data
          params = distribution.fit(data)

          # Separate parts of parameters
          arg = params[:-2]
          loc = params[-2]
          scale = params[-1]

          # Calculate fitted PDF, error with fit in distribution and value at pvalue_threshold
          pdf = distribution.pdf(x, loc=loc, scale=scale, *arg)
          sse = np.sum(np.power(y - pdf, 2.0))
          value_at_pvalue = distribution.ppf(self.pvalue_threshold, *arg, loc=loc, scale=scale) if arg else distribution.ppf(self.pvalue_threshold, loc=loc, scale=scale)

          # Update best distribution
          if(best_sse > sse > 0):
            best_distribution = distribution
            best_params = params
            best_sse = sse
            best_pvalue = value_at_pvalue

      except Exception as e:
        # PASS - TODO WARNING - IF NO DIST CAN BE FIT, SEND ERROR
        pass

    return best_distribution.name, best_params, best_pvalue

  def calculate_pvalues(self, chromosome, sica_dist):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Calculating best p-value given the current distribution
    data = self.distribution_dictionary[chromosome][sica_dist]
    best_distribution, best_params, best_pvalue = self.best_fit_distribution(data)
    self.pvalue_dictionary[chromosome][sica_dist] = [best_distribution, best_params, best_pvalue]

    # Iterating on matrix to update annotation dictionary
    for key, value in self.contact_map.matrix[chromosome].items():

      # Add distance to diagonal to dictionary
      distance_to_diag = self.dist_to_diag_dictionary[chromosome][key]

      # Check if current bin is inside the sica dist
      if(sica_dist[0] <= distance_to_diag < sica_dist[1]):

        # Update annotation disctionary if value is bigger than or equal threshold value at p-value
        if(value >= best_pvalue):
          self.annotation_dictionary[chromosome][key] = self.annotation_dictionary[chromosome][key].upper()
          try:
            self.significant_values_dictionary[chromosome].append(value)
          except Exception:
            self.significant_values_dictionary[chromosome] = [value]

  def add_calculate_pvalues(self, chromosome, sica_dist):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Append job to queue
    self.process_queue.append((chromosome, sica_dist))

  def run_calculate_pvalues(self):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Execute job queue
    pool = multiprocessing.Pool(self.ncpu)
    pool.starmap(self.calculate_pvalues, [arguments for arguments in self.process_queue])
    pool.close()
    pool.join()

    # Clean queue
    pool = None
    self.process_queue = []
    gc.collect()


  #############################################################################
  # Star Dictionaries
  #############################################################################

  def main_star_contacts(self):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Iterating on valid chromosomes - Starring contacts
    for chromosome in self.contact_map.valid_chromosome_list:

      # Add histogram calculation job to queue
      self.add_star_contacts(chromosome)

    # Run histogram calculation jobs
    self.run_star_contacts()

  def star_contacts(self, chromosome, bottom_bin_ext_range = [3,10], left_bin_ext_range = [3,10], right_bin_ext_range = [1,4], top_bin_ext_range = [1,4]):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Iterating on matrix
    for key, value in self.contact_map.matrix[chromosome].items():

      # Bin is a valid peak
      if(self.annotation_dictionary[chromosome][key].isupper()):

        # Key point
        key_row = key[0]
        key_col = key[1]

        # Selecting lengths
        bottom = random.randint(bottom_bin_ext_range[0], bottom_bin_ext_range[1])
        left = random.randint(left_bin_ext_range[0], left_bin_ext_range[1])
        right = random.randint(right_bin_ext_range[0], right_bin_ext_range[1])
        top = random.randint(top_bin_ext_range[0], top_bin_ext_range[1])

        # Iterating on point
        for i in range(-top, bottom):
          for j in range(-left, right):

            # Absolute i and j
            absi = abs(i)
            absj = abs(j)

            # If point is valid to be modified
            if(((i <= 0) and (j <= 0) and ((absi + absj) <= min(left, top))) or
               ((i != j) and (i >= 0) and (j >= 0) and ((absi + absj) <= min(right, bottom))) or
               ((i < 0) and (j > 0) and ((absi + absj) <= min(right, top))) or
               ((i > 0) and (j < 0) and ((absi + absj) <= min(left, bottom)))):

              # Decrease
              decrease = (absi + absj) / value

              # Bonuscrosslb
              bonuscrosslb = 0
              if(((i == 0) and (j <= 0)) or ((j == 0) and (i >= 0))):
                bonuscrosslb = random.uniform(0.25, 0.3) * value

              # Bonuscross
              bonuscross = 0
              if((i == 0) or (j == 0)):
                bonuscross = random.uniform(0.1, 0.25) * value

              # Bonuslb
              bonuslb = 0
              if((i >= 0) and (j >= 0)):
                bonuslb = random.uniform(0.1, 0.25) * value

              # Final score
              final_score = value - decrease + bonuscrosslb + bonuscross + bonuslb
              bpi, bpj = self.contact_map.matrix.bin_to_bp(i, j)
              self.contact_map.matrix.add(chromosome, key_row + bpi, key_col + bpj, final_score)

  def add_star_contacts(self, chromosome):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Append job to queue
    self.process_queue.append((chromosome))

  def run_star_contacts(self):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Execute job queue
    pool = multiprocessing.Pool(self.ncpu)
    pool.starmap(self.star_contacts, [arguments for arguments in self.process_queue])
    pool.close()
    pool.join()

    # Clean queue
    pool = None
    self.process_queue = []
    gc.collect()

