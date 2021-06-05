from __future__ import print_function
"""
DPMM Module
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
import traceback
import subprocess
import configparser
import multiprocessing

# Internal
from bloom.sica import SicaDist
from bloom.contact_map import ContactMap
from bloom.util import ErrorHandler

# External
import numpy as np

###################################################################################################
# Dpmm Class
###################################################################################################

class Dpmm():
  """This class represents TODO.

  *Keyword arguments:*

    - argument1 -- Short description. This argument represents a long description. It can be:
      - Possibility 1: A possibility 1.
      - Possibility 2: A possibility 2.

    - argument2 -- Short description. This argument represents a long description. It can be:
      - Possibility 1: A possibility 1.
      - Possibility 2: A possibility 2.
  """

  def __init__(self, ncpu, contact_map, sica_instance, random_degrade_range = [0.01, 0.02], degrade_multiplier = 0.05,
               half_length_bin_interval = [1, 5], value_range = [10e-4, 10e-3], random_range = [10e-8, 10e-7], iteration_multiplier = 1000, seed = 123):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Main objects
    self.contact_map = contact_map

    # Degrade objects
    self.sica_instance = sica_instance
    self.random_degrade_range = random_degrade_range
    self.degrade_multiplier = degrade_multiplier

    # Shape objects
    self.half_length_bin_interval = half_length_bin_interval
    self.value_range = value_range
    self.random_range = random_range
    self.iteration_multiplier = iteration_multiplier

    # Auxiliary parameters
    self.ncpu = ncpu
    self.process_queue = []
    random.seed(seed)

    # Utilitary objects
    self.error_handler = ErrorHandler()
    self.sica_dist_handler = SicaDist(self.contact_map, self.sica_instance.avoid_distance)

  #############################################################################
  # Diagonal Degrade
  #############################################################################

  def diagonal_degrade(self):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Get valid chromosome list
    valid_chromosome_list = self.contact_map.valid_chromosome_list
    
    # Iterating on valid chromosomes
    for chromosome in valid_chromosome_list:

      # Add introduce_squares to the queue
      #self.add_degrade(self.contact_map, chromosome)
      self.degrade(self.contact_map, chromosome)

    # Run introduce squares
    #self.run_degrade()

  def degrade(self, contact_map, chromosome):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    # Vector of elements to add
    elements_to_add = []

    # Get highest matrix value
    #highest_value = contact_map.max_value[chromosome]
    avoid_numpy_array = np.array(self.sica_instance.distribution_dictionary[chromosome][self.sica_dist_handler.A])
    average_avoided_value = np.mean(avoid_numpy_array)

    # Iterate on the rows of this chromosome's contact map
    for row in range(0, contact_map.total_1d_bins[chromosome]):

      # Calculate row's random multiplier and exponential multiplier
      row_random_r = (random.uniform(self.random_degrade_range[0] * average_avoided_value, self.random_degrade_range[1] * average_avoided_value) + average_avoided_value)
      row_expone_l = self.degrade_multiplier * average_avoided_value # TODO - Decide this formula as a random range based also on the sparsity level and the avoid_distance

      # Iterate on columns backwards
      counter = 1
      avoid_distance_in_bins = self.contact_map.bp_to_bin(self.sica_instance.avoid_distance)
      for col in range(row + avoid_distance_in_bins, row - 1, -1):

        # Calculate value to add in the matrix
        value_to_add = row_random_r + ((1 + row_expone_l) * counter)
        bprow = self.contact_map.bin_to_bp(row)
        bpcol = self.contact_map.bin_to_bp(col)
        elements_to_add.append((chromosome, bprow, bpcol, value_to_add))

        # Update counter
        counter += 1

    # Adding elements
    for element in elements_to_add:
      self.contact_map.add(element[0], element[1], element[2], element[3])

  def add_degrade(self, contact_map, chromosome):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Append job to queue
    self.process_queue.append((genome_id, contact_map, chromosome))

  def run_degrade(self):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    
    # Execute job queue
    pool = multiprocessing.Pool(self.ncpu)
    pool.starmap(self.degrade, [arguments for arguments in self.process_queue])
    pool.close()
    pool.join()

    # Clean queue
    pool = None
    self.process_queue = []
    gc.collect()


  #############################################################################
  # Shapes
  #############################################################################

  def introduce_shapes(self):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Get valid chromosome list
    valid_chromosome_list = self.contact_map.valid_chromosome_list
    
    # Iterating on valid chromosomes
    for chromosome in valid_chromosome_list:

      # Get number of iterations based on the size of the matrix
      iterations = self.contact_map.total_bins[chromosome] * self.iteration_multiplier

      # Add introduce_squares to the queue
      #self.add_introduce_squares(self.contact_map, chromosome, iterations)
      self.introduce_squares(self.contact_map, chromosome, iterations)

    # Run introduce squares
    #self.run_introduce_squares()

    # Iterating on valid chromosomes
    for chromosome in valid_chromosome_list:

      # Get number of iterations based on the size of the matrix
      iterations = self.contact_map.total_bins[chromosome] * self.iteration_multiplier

      # Add introduce_circles to the queue
      #self.add_introduce_circles(self.contact_map, chromosome, iterations)
      self.introduce_circles(self.contact_map, chromosome, iterations)

    # Run introduce circles
    #self.run_introduce_circles()

  def introduce_squares(self, contact_map, chromosome, iterations):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Maximum bp
    maximum_bp = contact_map.total_1d_bp[chromosome]
  
    # Performing the total number of iterations
    for main_iteration in range(0, iterations):

      # Square middle point
      chosen_length = random.randint(self.half_length_bin_interval[0], self.half_length_bin_interval[1]) * contact_map.resolution
      middle_bp_i = random.randrange(chosen_length, maximum_bp - (2 * chosen_length), contact_map.resolution)
      middle_bp_j = random.randrange(middle_bp_i + chosen_length, maximum_bp - chosen_length, contact_map.resolution)

      # Base value to add
      base_value = random.uniform(self.value_range[0], self.value_range[1])

      # Add a square
      for i in range(middle_bp_i - chosen_length, middle_bp_i + chosen_length + 1, contact_map.resolution):
        for j in range(middle_bp_j - chosen_length, middle_bp_j + chosen_length + 1, contact_map.resolution):

          # Value to add
          value_to_add = base_value + random.uniform(self.random_range[0], self.random_range[1])
          contact_map.add(chromosome, i, j, value_to_add)

  def add_introduce_squares(self, contact_map, chromosome, iterations):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Append job to queue
    self.process_queue.append((genome_id, contact_map, chromosome, iterations))

  def run_introduce_squares(self):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    
    # Execute job queue
    pool = multiprocessing.Pool(self.ncpu)
    pool.starmap(self.introduce_squares, [arguments for arguments in self.process_queue])
    pool.close()
    pool.join()

    # Clean queue
    pool = None
    self.process_queue = []
    gc.collect()

  def introduce_circles(self, contact_map, chromosome, iterations):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Maximum bp
    maximum_bp = contact_map.total_1d_bp[chromosome]
  
    # Performing the total number of iterations
    for main_iteration in range(0, iterations):

      # Square middle point
      base_chosen_length = random.randint(self.half_length_bin_interval[0], self.half_length_bin_interval[1])
      chosen_length = base_chosen_length * contact_map.resolution
      middle_bp_i = random.randrange(chosen_length, maximum_bp - (2 * chosen_length), contact_map.resolution)
      middle_bp_j = random.randrange(middle_bp_i + chosen_length, maximum_bp - chosen_length, contact_map.resolution)

      # Base value to add
      base_value = random.uniform(self.value_range[0], self.value_range[1])

      # Add a circle
      counterI = - base_chosen_length
      counterJ = - base_chosen_length
      for i in range(middle_bp_i - chosen_length, middle_bp_i + chosen_length + 1, contact_map.resolution):
        for j in range(middle_bp_j - chosen_length, middle_bp_j + chosen_length + 1, contact_map.resolution):

          # Circle constraint
          if((abs(counterI) + abs(counterJ)) > base_chosen_length):
            continue

          # Value to add
          value_to_add = base_value + random.uniform(self.random_range[0], self.random_range[1])
          contact_map.add(chromosome, i, j, value_to_add)

          # Updating counter j
          counterJ += 1

          # Updating counter i
        counterI += 1

  def add_introduce_circles(self, contact_map, chromosome, iterations):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Append job to queue
    self.process_queue.append((genome_id, contact_map, chromosome, iterations))

  def run_introduce_circles(self):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    
    # Execute job queue
    pool = multiprocessing.Pool(self.ncpu)
    pool.starmap(self.introduce_circles, [arguments for arguments in self.process_queue])
    pool.close()
    pool.join()

    # Clean queue
    pool = None
    self.process_queue = []
    gc.collect()

"""
introduce_shapes(self, chromosome, diag_dist_scale, neighbor_dist_scale, max_neighbors, diag_dist_rand_range = [0.0, 0.1], neighbor_dist_rand_range = [0.0, 0.1])
    
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

#distance_from_diagonal(self, key)
#random.uniform(diag_dist_rand_range[0], diag_dist_rand_range[1])
#random.uniform(neighbor_dist_rand_range[0], neighbor_dist_rand_range[1])
"""


###################################################################################################
# DPMM Main Class
###################################################################################################

class DPMM():
  """This class represents TODO.

  *Keyword arguments:*

    - argument1 -- Short description. This argument represents a long description. It can be:
      - Possibility 1: A possibility 1.
      - Possibility 2: A possibility 2.

    - argument2 -- Short description. This argument represents a long description. It can be:
      - Possibility 1: A possibility 1.
      - Possibility 2: A possibility 2.
  """

  def __init__(self, prior, alpha, D, manip=None, phi=None, label=None):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Main objects
    self.prior = prior
    self.alpha = alpha
    self._D = D  # data
    if manip is None:
      manip = NullManip()
    self.manip = manip

    # Auxiliary latency
    self._initD()
    self.manip.init(self.D)
    self.n = len(self.D)

    # Initialize r_i array
    self.p = self.alpha * self.prior.pred(self.mD)[:, np.newaxis]

    # Dirichlet parameterization
    if phi is None:
      self.init_phi()
    else:
      self.phi = phi
      self.label = label
      self.nphi = [np.sum(label == i) for i in xrange(label.max())]

  def init_phi(self):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Parameter initialization
    self.label = np.zeros((self.n), dtype=int)
    self.phi = []
    self.nphi = []
    for i in xrange(self.n):
      self.update_c_i(i)
    self.update_phi()

  @property
  def mD(self):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Latent update
    if self.manip_needs_update:
      self._mD = self.manip(self.D)
      self.manip_needs_update = False
    return self._mD

  def _initD(self):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Latent vector
    if isinstance(self._D, PseudoMarginalData):
      self.D = np.mean(self._D.data, axis=1)
    else:
      self.D = self._D
    self.manip_needs_update = True

  def draw_new_label(self, i):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Do not do yield explicitely
    picked = pick_discrete(self.p[i]*np.append([1], self.nphi)) - 1
    return picked

  def del_c_i(self, i):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Decrement mixture instance
    label = self.label[i]
    self.nphi[label] -= 1

    # If mix is deleted, delete params
    if self.nphi[label] == 0:
      del self.phi[label]
      del self.nphi[label]
      self.label[self.label >= label] -= 1
      self.p = np.delete(self.p, label+1, axis=1)

  def update_c_i(self, i):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    label = self.draw_new_label(i)

    # Draw
    if label == -1:
      new_phi = self.prior.post(self.mD[i]).sample()
      self.phi.append(new_phi)
      self.nphi.append(1)
      self.label[i] = len(self.phi)-1
      self.p = np.append(self.p, np.zeros((self.n, 1), dtype=float), axis=1)
      self.p[i+1:, -1] = self.prior.like1(self.mD[i+1:], new_phi)
    else:
      self.label[i] = label
      self.nphi[label] += 1

  def update_c(self):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Update drawing
    for i in xrange(self.n):
      self.del_c_i(i)
      self.update_c_i(i)

  def update_phi(self):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Update parameter
    tot = 0
    for i in xrange(len(self.phi)):
      index = self.label == i
      tot += sum(index)
      data = self.mD[index]
      new_phi = self.prior.post(data).sample()
      self.phi[i] = new_phi
    self.p[:, 1:] = self.prior.like1(self.mD[:, np.newaxis], np.array(self.phi))

  def update_latent_data(self):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Iteration to calculate weights -> selecting a representative sample
    if isinstance(self._D, PseudoMarginalData):
      for i, ph in enumerate(self.phi):

        index = np.nonzero(self.label == i)[0]
        data = self._D[index]
        
        ps = self.prior.like1(self.manip(data.data), ph) / data.interim_prior
        ps /= np.sum(ps, axis=1)[:, np.newaxis]

        for j, p in enumerate(ps):
          self.D[index[j]] = data.data[j, pick_discrete(p)]

      self.p[:, 0] = self.alpha * self.prior.pred(self.mD)
      self.manip_needs_update = True
    else:
      pass # TODO - Report

  def update(self, n=1):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # General update
    for j in xrange(n):
      self.update_c()
      self.update_latent_data()
      self.update_phi()
      self.manip.update(self.D, self.phi, self.label, self.prior)
      self.manip_needs_update = True


###################################################################################################
# GaussND Class
###################################################################################################

class GaussND():
  """This class represents TODO.

  *Keyword arguments:*

    - argument1 -- Short description. This argument represents a long description. It can be:
      - Possibility 1: A possibility 1.
      - Possibility 2: A possibility 2.

    - argument2 -- Short description. This argument represents a long description. It can be:
      - Possibility 1: A possibility 1.
      - Possibility 2: A possibility 2.
  """

  def __init__(self, mu, Sig):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # GND Parameters
    self.mu = np.atleast_1d(mu)
    self.Sig = np.atleast_2d(Sig)
    self.d = len(self.mu)

  def cond(self, x):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Params of cond
    fixed = np.nonzero([x_ is not None for x_ in x])
    nonfixed = np.nonzero([x_ is None for x_ in x])
    mu1 = self.mu[nonfixed]
    mu2 = self.mu[fixed]
    Sig11 = self.Sig[nonfixed, nonfixed]
    Sig12 = self.Sig[fixed, nonfixed]
    Sig22 = self.Sig[fixed, fixed]

    # Calculation of cond
    new_mu = mu1 + np.dot(Sig12, np.dot(np.linalg.inv(Sig22), x[fixed[0]] - mu2))
    new_Sig = Sig11 - np.dot(Sig12, np.dot(np.linalg.inv(Sig22), Sig12.T))

    # Return object
    return GaussND(new_mu, new_Sig)

  def sample(self, size=None):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Sample drawn as a MN
    if self.d == 1:
      return np.random.normal(self.mu, scale=np.sqrt(self.Sig), size=size)
    else:
      return np.random.multivariate_normal(self.mu, self.Sig, size=size)


###################################################################################################
# GMM Class
###################################################################################################

class GMM():
  """This class represents TODO.

  *Keyword arguments:*

    - argument1 -- Short description. This argument represents a long description. It can be:
      - Possibility 1: A possibility 1.
      - Possibility 2: A possibility 2.

    - argument2 -- Short description. This argument represents a long description. It can be:
      - Possibility 1: A possibility 1.
      - Possibility 2: A possibility 2.
  """

  def __init__(self, components, proportions):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # GMM
    self.components = components
    self.proportions = proportions
    self.d = self.components[0].d

  def cond(self, x):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Conditional GMM
    components = [c.cond(x) for c in self.components]
    return GMM(components, self.proportions)

  def sample(self, size=None):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Draw from a multinomial
    if size is None:
      nums = np.random.multinomial(1, self.proportions)
      c = nums.index(1) # which class got picked
      return self.components[c].sample()
    else:
      n = np.prod(size)
      if self.d == 1:
        out = np.empty((n,), dtype=float)
        nums = np.random.multinomial(n, self.proportions)
        i = 0
        for component, num in zip(self.components, nums):
          out[i:i+num] = component.sample(size=num)
          i += num
        out = out.reshape(size)
      else:
        out = np.empty((n, self.d), dtype=float)
        nums = np.random.multinomial(n, self.proportions)
        i = 0
        for component, num in zip(self.components, nums):
          out[i:i+num] = component.sample(size=num)
          i += num
        if isinstance(size, int):
          out = out.reshape((size, self.d))
        else:
          out = out.reshape(size+(self.d,))
      return out


###################################################################################################
# Prior Class
###################################################################################################

class Prior():
  """This class represents TODO.

  *Keyword arguments:*

    - argument1 -- Short description. This argument represents a long description. It can be:
      - Possibility 1: A possibility 1.
      - Possibility 2: A possibility 2.

    - argument2 -- Short description. This argument represents a long description. It can be:
      - Possibility 1: A possibility 1.
      - Possibility 2: A possibility 2.
  """

  def __init__(self, post=None, *args, **kwargs):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    if post is None:
      post = type(self)
    self._post = post

  def sample(self, size=None):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    raise NotImplementedError

  def like1(self, x, *args, **kwargs):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    raise NotImplementedError

  def likelihood(self, D, *args, **kwargs):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    return np.prod(self.like1(D, *args, **kwargs))

  def lnlikelihood(self, D, *args, **kwargs):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    return np.log(self.likelihood(D, *args, **kwargs))

  def __call__(self, *args):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    raise NotImplementedError

  def _post_params(self, D):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    raise NotImplementedError

  def post(self, D):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    return self._post(*self._post_params(D))

  def pred(self, x):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    raise NotImplementedError


###################################################################################################
# NormInvWish Class
###################################################################################################

class NormInvWish(Prior):
  """This class represents TODO.

  *Keyword arguments:*

    - argument1 -- Short description. This argument represents a long description. It can be:
      - Possibility 1: A possibility 1.
      - Possibility 2: A possibility 2.

    - argument2 -- Short description. This argument represents a long description. It can be:
      - Possibility 1: A possibility 1.
      - Possibility 2: A possibility 2.
  """

  def __init__(self, mu_0, kappa_0, Lam_0, nu_0):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    self.mu_0 = np.array(mu_0, dtype=float)
    self.kappa_0 = float(kappa_0)
    self.Lam_0 = np.array(Lam_0, dtype=float)
    self.nu_0 = int(nu_0)
    self.d = len(mu_0)
    self.model_dtype = np.dtype([('mu', float, self.d), ('Sig', float, (self.d, self.d))])
    super(NormInvWish, self).__init__()

  def _S(self, D):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    Dbar = np.mean(D, axis=0)
    return np.dot((D-Dbar).T, (D-Dbar))

  def sample(self, size=None):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    Sig = random_invwish(dof=self.nu_0, invS=self.Lam_0, size=size)
    if size is None:
      ret = np.zeros(1, dtype=self.model_dtype)
      ret['Sig'] = Sig
      ret['mu'] = np.random.multivariate_normal(self.mu_0, Sig/self.kappa_0)
      return ret[0]
    else:
      ret = np.zeros(size, dtype=self.model_dtype)
      ret['Sig'] = Sig
      for r in ret.ravel():
        r['mu'] = np.random.multivariate_normal(self.mu_0, r['Sig']/self.kappa_0)
      return ret

  def like1(self, *args):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    if len(args) == 2:
      x, theta = args
      mu = theta['mu']
      Sig = theta['Sig']
    elif len(args) == 3:
      x, mu, Sig = args
    assert x.shape[-1] == self.d
    assert mu.shape[-1] == self.d
    assert Sig.shape[-1] == Sig.shape[-2] == self.d
    norm = np.sqrt((2*np.pi)**self.d * np.linalg.det(Sig))
    einsum = np.einsum("...i,...ij,...j", x-mu, np.linalg.inv(Sig), x-mu)
    return np.exp(-0.5*einsum)/norm

  def __call__(self, *args):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    if len(args) == 1:
      mu = args[0]['mu']
      Sig = args[0]['Sig']
    elif len(args) == 2:
      mu, Sig = args
    nu_0, d = self.nu_0, self.d
    Z = (2.0**(nu_0*d/2.0) * gammad(d, nu_0/2.0) * \
        (2.0*np.pi/self.kappa_0)**(d/2.0) / np.linalg.det(self.Lam_0)**(nu_0/2.0))
    detSig = np.linalg.det(Sig)
    invSig = np.linalg.inv(Sig)
    einsum = np.einsum("...i,...ij,...j", mu-self.mu_0, invSig, mu-self.mu_0)
    return 1./Z * detSig**(-((nu_0+d)/2.0+1.0)) * \
           np.exp(-0.5*np.trace(np.einsum("...ij,...jk->...ik", self.Lam_0, invSig), axis1=-2, axis2=-1) - self.kappa_0/2.0*einsum)

  def _post_params(self, D):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    shape = D.shape
    if len(shape) == 2:
      n = shape[0]
      Dbar = np.mean(D, axis=0)
    elif len(shape) == 1:
      n = 1
      Dbar = np.mean(D)
    kappa_n = self.kappa_0 + n
    nu_n = self.nu_0 + n
    mu_n = (self.kappa_0 * self.mu_0 + n * Dbar) / kappa_n
    x = (Dbar-self.mu_0)[:, np.newaxis]
    Lam_n = (self.Lam_0 + self._S(D) + self.kappa_0*n/kappa_n*np.dot(x, x.T))
    return mu_n, kappa_n, Lam_n, nu_n

  def pred(self, x):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    return multivariate_t_density(self.nu_0-self.d+1, self.mu_0, self.Lam_0*(self.kappa_0+1)/(self.kappa_0 - self.d + 1), x)

  def evidence(self, D):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    shape = D.shape
    if len(shape) == 2:
      n, d = shape
    elif len(shape) == 1:
      n, d = 1, shape[0]
    assert d == self.d
    mu_n, kappa_n, Lam_n, nu_n = self._post_params(D)
    detLam0 = np.linalg.det(self.Lam_0)
    detLamn = np.linalg.det(Lam_n)
    num = gammad(d, nu_n/2.0) * detLam0**(self.nu_0/2.0)
    den = np.pi**(n*d/2.0) * gammad(d, self.nu_0/2.0) * detLamn**(nu_n/2.0)
    return num/den * (self.kappa_0/kappa_n)**(d/2.0)

