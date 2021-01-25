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
from bloom.contact_map import ContactMap
from bloom.util import ErrorHandler

# External

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

  def __init__(self, ncpu, contact_map, avoid_distance, removed_dict, random_degrade_range = [0.01, 0.02], degrade_multiplier = 0.05,
               half_length_bin_interval = [1, 5], value_range = [10e-3, 10e-4], random_range = [10e-9, 10e-10], iteration_multiplier = 1000, seed = 123):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Main objects
    self.contact_map = contact_map

    # Degrade objects
    self.avoid_distance = avoid_distance
    self.removed_dict = removed_dict
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
      self.add_degrade(self.contact_map, chromosome)

    # Run introduce squares
    self.run_degrade()

  def degrade(self, contact_map, chromosome):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Get highest matrix value
    highest_value = contact_map.max_value[chromosome]

    # Iterate on the rows of this chromosome's contact map
    for row in range(0, contact_map.max_bin[chromosome] + 1):

      # Check whether row is to be avoided (removed due to blacklist or zeros)
      try:
        avoid = self.removed_dict[row]
        continue
      except Exception: pass

      # Calculate row's random multiplier and exponential multiplier
      row_random_r = (random.uniform(self.random_degrade_range[0] * highest_value, self.random_degrade_range[1] * highest_value) + highest_value)
      row_expone_l = self.degrade_multiplier * highest_value # TODO - Decide this formula as a random range based also on the sparsity level and the avoid_distance

      # Iterate on "diagonal-avoided" columns backwards
      counter = 1
      for col in range(row + self.avoid_distance, row - 1, -1):

        # Calculate value to add in the matrix
        value_to_add = row_random_r + ((1 + row_expone_l) ** counter)

        # TODO - ADD VALUE TO ROW / COL

        # Update counter
        counter += 1


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
      iterations = self.contact_map.total_bins * self.iteration_multiplier

      # Add introduce_squares to the queue
      self.add_introduce_squares(contact_map, chromosome, iterations)

    # Run introduce squares
    self.run_introduce_squares()

    # Iterating on valid chromosomes
    for chromosome in valid_chromosome_list:

      # Get number of iterations based on the size of the matrix
      iterations = self.contact_map.total_bins * self.iteration_multiplier

      # Add introduce_circles to the queue
      self.add_introduce_circles(contact_map, chromosome, iterations)

    # Run introduce circles
    self.run_introduce_circles()

  def introduce_squares(self, contact_map, chromosome, iterations):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Maximum bp
    maximum_bp = contact_map.max_bp[chromosome] + 1
  
    # Performing the total number of iterations
    for main_iteration in range(0, iterations):

      # Square middle point
      middle_bp = random.randrange(0, maximum_bp, contact_map.resolution)
      chosen_length = random.randint(self.half_length_bin_interval[0], self.half_length_bin_interval[1]) * contact_map.resolution

      # Base value to add
      base_value = random.uniform(self.value_range[0], self.value_range[1])

      # Add a square
      for i in range(middle_bp - chosen_length, middle_bp + chosen_length + 1, contact_map.resolution):
        for j in range(middle_bp - chosen_length, middle_bp + chosen_length + 1, contact_map.resolution):

          # Minimum and maximum coordinates
          min_pos = min(i, j)
          max_pos = max(i, j)

          # Value to add
          value_to_add = base_value + random.uniform(self.random_range[0], self.random_range[1])

          # Add value to matrix
          # TODO

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
    maximum_bp = contact_map.max_bp[chromosome] + 1
  
    # Performing the total number of iterations
    for main_iteration in range(0, iterations):

      # Square middle point
      middle_bp = random.randrange(0, maximum_bp, contact_map.resolution)
      base_chosen_length = random.randint(self.half_length_bin_interval[0], self.half_length_bin_interval[1])
      chosen_length = base_chosen_length * contact_map.resolution

      # Base value to add
      base_value = random.uniform(self.value_range[0], self.value_range[1])

      # Add a circle
      counterI = - base_chosen_length
      counterJ = - base_chosen_length
      for i in range(middle_bp - chosen_length, middle_bp + chosen_length + 1, contact_map.resolution):
        for j in range(middle_bp - chosen_length, middle_bp + chosen_length + 1, contact_map.resolution):

          # Circle constraint
          if((abs(counterI) + abs(counterJ)) > base_chosen_length):
            continue

          # Minimum and maximum coordinates
          min_pos = min(i, j)
          max_pos = max(i, j)

          # Value to add
          value_to_add = base_value + random.uniform(self.random_range[0], self.random_range[1])

          # Add value to matrix
          # TODO

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







































