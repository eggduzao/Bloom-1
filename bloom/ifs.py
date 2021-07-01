from __future__ import print_function
"""
IFS Module
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
from bloom.sica import Sica, SicaDist
from bloom.util import ErrorHandler, AuxiliaryFunctions

# External
import numpy as np
import scipy
import scipy.stats as st

###################################################################################################
# Ifs Class
###################################################################################################

class Ifs():
  """This class represents TODO.

  *Keyword arguments:*

    - argument1 -- Short description. This argument represents a long description. It can be:
      - Possibility 1: A possibility 1.
      - Possibility 2: A possibility 2.

    - argument2 -- Short description. This argument represents a long description. It can be:
      - Possibility 1: A possibility 1.
      - Possibility 2: A possibility 2.
  """

  def __init__(self, contact_map, sica_instance, goba_instance, dpmm_instance, io_instance,
               output_loop_file_name, output_matrix_file_name, matrix_output_format, seed = None):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Seed
    if(seed):
      random.seed(seed)

    # Class objects
    self.contact_map = contact_map
    self.sica_instance = sica_instance
    self.goba_instance = goba_instance
    self.dpmm_instance = dpmm_instance
    self.io_instance = io_instance
    self.output_loop_file_name = output_loop_file_name
    self.output_matrix_file_name = output_matrix_file_name
    self.matrix_output_format = matrix_output_format
    self.seed = seed

    # Auxiliary objects
    self.ifs_list = []

    # Auxiliary parameters
    self.ncpu = self.sica_instance.ncpu
    self.process_queue = []

    # Utilitary objects
    self.error_handler = ErrorHandler()


  #############################################################################
  # IFS
  #############################################################################

  def main_calculate_ifs(self, min_to_zero = True):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Iterating on valid chromosomes - Calculate histograms
    for chromosome in self.contact_map.valid_chromosome_list:

      # Add histogram calculation job to queue
      #self.add_calculate_ifs(chromosome)
      self.calculate_ifs(chromosome)

    # Run histogram calculation jobs
    #self.run_calculate_ifs()

    # Sort and write IFS list
    self.sort_ifs_list()
    self.standardize_ifs_list(min_to_zero = min_to_zero)
    self.write_ifs()

  def calculate_ifs(self, chromosome):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Allowed distances
    convoluted_dict = self.sica_instance.dist_handler.get_key_to_bin_dict(self.sica_instance.dist_handler.tbin_dist_dict, upper = True)

    # Iterating on matrix
    for key, value in self.sica_instance.annotation_dictionary[chromosome].items():

      # Check if convoluted to t
      try:
        if(not convoluted_dict[value.upper()]):
          continue
      except Exception:
        continue

      # Check if value is real after convolution
      try:
        newvalue = self.contact_map.matrix[chromosome][key]
      except Exception:
        continue

      # Check distance to diagonal
      newvalue = newvalue + (random.uniform(0, 0.1) * newvalue)
      newvalue = np.log2(newvalue)
      self.ifs_list.append([chromosome] + [key[0], key[1]] + [newvalue])

  def add_calculate_ifs(self, chromosome):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Append job to queue
    self.process_queue.append((chromosome))

  def run_calculate_ifs(self):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Execute job queue
    pool = multiprocessing.Pool(self.ncpu)
    pool.starmap(self.calculate_ifs, [arguments for arguments in self.process_queue])
    pool.close()
    pool.join()

    # Clean queue
    pool = None
    self.process_queue = []
    gc.collect()

  def sort_ifs_list(self):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Sort IFS list by: chromosome, pos1, pos2
    self.ifs_list = sorted(self.ifs_list, key = lambda x: (x[0], x[1], x[2]))

  def standardize_ifs_list(self, min_to_zero = True):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Calculate minimum and maximum
    minV = np.inf
    maxV = -np.inf
    for v in self.ifs_list:
      value = float(v[3])
      if(value < minV): minV = value
      if(value > maxV): maxV = value

    # Minimum is 0
    if(min_to_zero):
      minV = 0.0
 
    # Standardize list
    for v in self.ifs_list:
      v[3] = (float(v[3]) - minV) / (maxV - minV)
      v[3] = min(1.0, v[3])

  def write_ifs(self):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Open output matrix file
    output_loop_file = open(self.output_loop_file_name, "w")

    # Write IFS list into bedgraph
    for v in self.ifs_list:
      chromosome = str(v[0])
      p11 = str(v[1])
      p12 = str(v[1] + self.contact_map.resolution)
      p21 = str(v[2])
      p22 = str(v[2] + self.contact_map.resolution)
      value = str(round(v[3], 6))
      output_loop_file.write("\t".join([chromosome, p11, p12, chromosome, p21, p22, value]) + "\n")

    # Closing file
    output_loop_file.close()


  #############################################################################
  # Fix Matrix
  #############################################################################

  def main_fix_matrix(self, multiplier = 1000, min_matrix_threshold = 1):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Iterating on valid chromosomes - Calculate histograms
    for chromosome in self.contact_map.valid_chromosome_list:

      # Add histogram calculation job to queue
      #self.add_fix_matrix(chromosome, multiplier, min_matrix_threshold)
      self.fix_matrix(chromosome, multiplier, min_matrix_threshold)

    # Run histogram calculation jobs
    #self.run_fix_matrix()

    # Write matrix
    self.write_matrix()

  def fix_matrix(self, chromosome, multiplier, min_matrix_threshold):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # List of keys to remove from dictionary
    maxV = -np.inf
    removed_keys = []

    # Iterating on matrix
    for key, value in self.contact_map.matrix[chromosome].items():

      # Binary version of keys
      brow = self.contact_map.bp_to_bin(key[0])
      bcol = self.contact_map.bp_to_bin(key[1])

      # 1. Remove intervals outside the chrom dict boundaries (min / max // bin / bp)
      if( (key[0] < 0) or (key[0] >= self.contact_map.total_1d_bp[chromosome]) or (key[0] % self.contact_map.resolution) or 
          (key[1] < 0) or (key[1] >= self.contact_map.total_1d_bp[chromosome]) or (key[1] % self.contact_map.resolution)):
        removed_keys.append(key)
        continue

      # 2. Make all sica.removed_dict rows/columns = 0 (remove the entry)
      try:
        self.sica_instance.removed_dict[chromosome][brow]
        removed_keys.append(key)
        continue
      except Exception:
        pass
      try:
        self.sica_instance.removed_dict[chromosome][bcol]
        removed_keys.append(key)
        continue
      except Exception:
        pass

      # 3. Remove all values <= min_matrix_threshold
      if(value <= min_matrix_threshold):
        removed_keys.append(key)
        continue

      # 4. Standardization: calculate the maximum value
      if(value > maxV):
        maxV = value

    # Remove values from dictionary
    for key in removed_keys:
      self.contact_map.matrix[chromosome].pop(key)

    # Iterating on matrix
    for key, value in self.contact_map.matrix[chromosome].items():

      # 4. Standardization: perform stanrdization and adding multiplier
      self.contact_map.matrix[chromosome][key] = min(((value + random.random()) - 0) / (maxV - 0), 1) * multiplier

  def add_fix_matrix(self, chromosome, multiplier, min_matrix_threshold):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Append job to queue
    self.process_queue.append((chromosome, multiplier, min_matrix_threshold))

  def run_fix_matrix(self):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Execute job queue
    pool = multiprocessing.Pool(self.ncpu)
    pool.starmap(self.fix_matrix, [arguments for arguments in self.process_queue])
    pool.close()
    pool.join()

    # Clean queue
    pool = None
    self.process_queue = []
    gc.collect()


  def write_matrix(self):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Writing matrix using IO instance
    self.io_instance.write(self.contact_map, self.output_matrix_file_name, output_format = self.matrix_output_format)


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
# GaussianMeanKnownVariance Class
###################################################################################################

class GaussianMeanKnownVariance(Prior):
  """This class represents TODO.

  *Keyword arguments:*

    - argument1 -- Short description. This argument represents a long description. It can be:
      - Possibility 1: A possibility 1.
      - Possibility 2: A possibility 2.

    - argument2 -- Short description. This argument represents a long description. It can be:
      - Possibility 1: A possibility 1.
      - Possibility 2: A possibility 2.
  """

  def __init__(self, mu_0, sigsqr_0, sigsqr):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    self.mu_0 = mu_0
    self.sigsqr_0 = sigsqr_0
    self.sigsqr = sigsqr
    self._norm1 = np.sqrt(2*np.pi*self.sigsqr)
    self._norm2 = np.sqrt(2*np.pi*self.sigsqr_0)
    super(GaussianMeanKnownVariance, self).__init__()

  def sample(self, size=None):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    if size is None:
      return np.random.normal(self.mu_0, np.sqrt(self.sigsqr_0))
    else:
      return np.random.normal(self.mu_0, np.sqrt(self.sigsqr_0), size=size)

  def like1(self, x, mu):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    return np.exp(-0.5*(x-mu)**2/self.sigsqr) / self._norm1

  def __call__(self, mu):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    return np.exp(-0.5*(mu-self.mu_0)**2/self.sigsqr_0) / self._norm2

  def _post_params(self, D):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    try:
      n = len(D)
    except TypeError:
      n = 1
    Dbar = np.mean(D)
    sigsqr_n = 1./(n/self.sigsqr + 1./self.sigsqr_0)
    mu_n = sigsqr_n * (self.mu_0/self.sigsqr_0 + n*Dbar/self.sigsqr)
    return mu_n, sigsqr_n, self.sigsqr

  def pred(self, x):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    sigsqr = self.sigsqr + self.sigsqr_0
    return np.exp(-0.5*(x-self.mu_0)**2/sigsqr) / np.sqrt(2*np.pi*sigsqr)


###################################################################################################
# InvGamma Class
###################################################################################################

class InvGamma(Prior):
  """This class represents TODO.

  *Keyword arguments:*

    - argument1 -- Short description. This argument represents a long description. It can be:
      - Possibility 1: A possibility 1.
      - Possibility 2: A possibility 2.

    - argument2 -- Short description. This argument represents a long description. It can be:
      - Possibility 1: A possibility 1.
      - Possibility 2: A possibility 2.
  """

  def __init__(self, alpha, beta, mu):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    self.alpha = alpha
    self.beta = beta
    self.mu = mu
    super(InvGamma, self).__init__()

  def sample(self, size=None):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    return 1./np.random.gamma(self.alpha, scale=self.beta, size=size)

  def like1(self, x, var):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    return np.exp(-0.5*(x-self.mu)**2/var) / np.sqrt(2*np.pi*var)

  def __call__(self, var):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    al, be = self.alpha, self.beta
    return be**(-al)/gamma(al) * var**(-1.-al) * np.exp(-1./(be*var))

  def _post_params(self, D):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    try:
      n = len(D)
    except TypeError:
      n = 1
    al_n = self.alpha + n/2.0
    be_n = 1./(1./self.beta + 0.5*np.sum((np.array(D)-self.mu)**2))
    return al_n, be_n, self.mu

  def pred(self, x):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    return t_density(2*self.alpha, self.mu, 1./self.beta/self.alpha, x)

  def evidence(self, D):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    raise NotImplementedError


###################################################################################################
# InvGamma2D Class
###################################################################################################

class InvGamma2D(Prior):
  """This class represents TODO.

  *Keyword arguments:*

    - argument1 -- Short description. This argument represents a long description. It can be:
      - Possibility 1: A possibility 1.
      - Possibility 2: A possibility 2.

    - argument2 -- Short description. This argument represents a long description. It can be:
      - Possibility 1: A possibility 1.
      - Possibility 2: A possibility 2.
  """

  def __init__(self, alpha, beta, mu):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    self.alpha = alpha
    self.beta = beta
    self.mu = np.array(mu)
    assert len(mu) == 2
    super(InvGamma2D, self).__init__()

  def sample(self, size=None):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    return 1./np.random.gamma(self.alpha, scale=self.beta, size=size)

  def like1(self, x, var):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    assert isinstance(x, np.ndarray)
    assert x.shape[-1] == 2
    return np.exp(-0.5*np.sum((x-self.mu)**2, axis=-1)/var) / (2*np.pi*var)

  def lnlikelihood(self, D, var):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    return -0.5*np.sum((D-self.mu)**2)/var - D.shape[0]*np.log(2*np.pi*var)

  def __call__(self, var):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    al, be = self.alpha, self.beta
    return be**(-al)/gamma(al) * var**(-1.-al) * np.exp(-1./(be*var))

  def _post_params(self, D):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    try:
      n = len(D)
    except TypeError:
      n = 1
    al_n = self.alpha + n
    be_n = 1./(1./self.beta + 0.5*np.sum((np.array(D)-self.mu)**2))
    return al_n, be_n, self.mu

  def pred(self, x):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    assert isinstance(x, np.ndarray)
    assert x.shape[-1] == 2
    return multivariate_t_density(2*self.alpha, self.mu, 1./self.beta/self.alpha*np.eye(2), x)

  def evidence(self, D):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    raise NotImplementedError


###################################################################################################
# NormInvChi2 Class
###################################################################################################

class NormInvChi2(Prior):
  """This class represents TODO.

  *Keyword arguments:*

    - argument1 -- Short description. This argument represents a long description. It can be:
      - Possibility 1: A possibility 1.
      - Possibility 2: A possibility 2.

    - argument2 -- Short description. This argument represents a long description. It can be:
      - Possibility 1: A possibility 1.
      - Possibility 2: A possibility 2.
  """

  def __init__(self, mu_0, kappa_0, sigsqr_0, nu_0):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    self.mu_0 = float(mu_0)
    self.kappa_0 = float(kappa_0)
    self.sigsqr_0 = float(sigsqr_0)
    self.nu_0 = float(nu_0)
    super(NormInvChi2, self).__init__()
    model_dtype = np.dtype([('mu', float), ('var', float)])

  def sample(self, size=None):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    if size is None:
      var = 1./np.random.chisquare(df=self.nu_0)*self.nu_0*self.sigsqr_0
      ret = np.zeros(1, dtype=self.model_dtype)
      ret['mu'] = np.random.normal(self.mu_0, np.sqrt(var/self.kappa_0))
      ret['var'] = var
      return ret[0]
    else:
      var = 1./np.random.chisquare(df=self.nu_0, size=size)*self.nu_0*self.sigsqr_0
      ret = np.zeros(size, dtype=self.model_dtype)
      ret['mu'] = (np.random.normal(self.mu_0, np.sqrt(1./self.kappa_0), size=size) * np.sqrt(var))
      ret['var'] = var
      return ret

  def like1(self, *args):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    if len(args) == 3:
      x, mu, var = args
    elif len(args) == 2:
      x, theta = args
      mu = theta['mu']
      var = theta['var']
    return np.exp(-0.5*(x-mu)**2/var) / np.sqrt(2*np.pi*var)

  def __call__(self, *args):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    if len(args) == 2:
      mu, var = args
    elif len(args) == 1:
      mu = args[0]['mu']
      var = args[0]['var']
    return (normal_density(self.mu_0, var/self.kappa_0, mu) * scaled_IX_density(self.nu_0, self.sigsqr_0, var))

  def _post_params(self, D):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    try:
      n = len(D)
    except TypeError:
      n = 1
    Dbar = np.mean(D)
    kappa_n = self.kappa_0 + n
    mu_n = (self.kappa_0*self.mu_0 + n*Dbar)/kappa_n
    nu_n = self.nu_0 + n
    sigsqr_n = ((self.nu_0*self.sigsqr_0 + np.sum((D-Dbar)**2) + n*self.kappa_0/(self.kappa_0+n)*(self.mu_0-Dbar)**2)/nu_n)
    return mu_n, kappa_n, sigsqr_n, nu_n

  def pred(self, x):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    return t_density(self.nu_0, self.mu_0, (1.+self.kappa_0)*self.sigsqr_0/self.kappa_0, x)

  def evidence(self, D):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    mu_n, kappa_n, sigsqr_n, nu_n = self._post_params(D)
    try:
      n = len(D)
    except:
      n = 1
    return (gamma(nu_n/2.0)/gamma(self.nu_0/2.0) * np.sqrt(self.kappa_0/kappa_n) * \
           (self.nu_0*self.sigsqr_0)**(self.nu_0/2.0) / (nu_n*sigsqr_n)**(nu_n/2.0) / np.pi**(n/2.0))

  def marginal_var(self, var):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    return scaled_IX_density(self.nu_0, self.sigsqr_0, var)

  def marginal_mu(self, mu):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    return t_density(self.nu_0, self.mu_0, self.sigsqr_0/self.kappa_0, mu)


###################################################################################################
# NormInvGamma Class
###################################################################################################

class NormInvGamma(Prior):
  """This class represents TODO.

  *Keyword arguments:*

    - argument1 -- Short description. This argument represents a long description. It can be:
      - Possibility 1: A possibility 1.
      - Possibility 2: A possibility 2.

    - argument2 -- Short description. This argument represents a long description. It can be:
      - Possibility 1: A possibility 1.
      - Possibility 2: A possibility 2.
  """

  def __init__(self, m_0, V_0, a_0, b_0):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    self.m_0 = float(m_0)
    self.V_0 = float(V_0)
    self.a_0 = float(a_0)
    self.b_0 = float(b_0)
    super(NormInvGamma, self).__init__()
    model_dtype = np.dtype([('mu', float), ('var', float)])

  def sample(self, size=None):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    if size is None:
      var = 1./np.random.gamma(self.a_0, scale=1./self.b_0)
      ret = np.zeros(1, dtype=self.model_dtype)
      ret['mu'] = np.random.normal(self.m_0, np.sqrt(self.V_0*var))
      ret['var'] = var
      return ret[0]
    else:
      var = 1./np.random.gamma(self.a_0, scale=1./self.b_0, size=size)
      ret = np.zeros(size, dtype=self.model_dtype)
      ret['mu'] = np.random.normal(self.m_0, np.sqrt(self.V_0), size=size)*np.sqrt(var)
      ret['var'] = var
      return ret

  def like1(self, *args):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    if len(args) == 3:
      x, mu, var = args
    elif len(args) == 2:
      x, theta = args
      mu = theta['mu']
      var = theta['var']
    return np.exp(-0.5*(x-mu)**2/var) / np.sqrt(2*np.pi*var)

  def __call__(self, *args):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    if len(args) == 1:
      mu = args[0]['mu']
      var = args[0]['var']
    elif len(args) == 2:
      mu, var = args
    normal = np.exp(-0.5*(self.m_0-mu)**2/(var*self.V_0))/np.sqrt(2*np.pi*var*self.V_0)
    ig = self.b_0**self.a_0/gamma(self.a_0)*var**(-(self.a_0+1))*np.exp(-self.b_0/var)
    return normal*ig

  def _post_params(self, D):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    try:
      n = len(D)
    except TypeError:
      n = 1
    Dbar = np.mean(D)
    invV_0 = 1./self.V_0
    V_n = 1./(invV_0 + n)
    m_n = V_n*(invV_0*self.m_0 + n*Dbar)
    a_n = self.a_0 + n/2.0
    b_n = self.b_0 + 0.5*(np.sum((D-Dbar)**2)+n / (1.0+n*self.V_0)*(self.m_0-Dbar)**2)
    return m_n, V_n, a_n, b_n

  def pred(self, x):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    return t_density(2.0*self.a_0, self.m_0, self.b_0*(1.0+self.V_0)/self.a_0, x)

  def evidence(self, D):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    m_n, V_n, a_n, b_n = self._post_params(D)
    try:
      n = len(D)
    except:
      n = 1
    return (np.sqrt(np.abs(V_n/self.V_0)) * (self.b_0**self.a_0)/(b_n**a_n) * \
            gamma(a_n)/gamma(self.a_0) / (np.pi**(n/2.0)*2.0**(n/2.0)))

  def marginal_var(self, var):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    nu_0 = 2*self.a_0
    sigsqr_0 = 2*self.b_0/nu_0
    return scaled_IX_density(nu_0, sigsqr_0, var)

  def marginal_mu(self, mu):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    mu_0 = self.m_0
    kappa_0 = 1./self.V_0
    nu_0 = 2*self.a_0
    sigsqr_0 = 2*self.b_0/nu_0
    return t_density(nu_0, mu_0, sigsqr_0/kappa_0, mu)

