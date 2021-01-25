"""
GOBA Module
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
# Goba Class
###################################################################################################

class Goba():
  """This class represents TODO.

  *Keyword arguments:*

    - argument1 -- Short description. This argument represents a long description. It can be:
      - Possibility 1: A possibility 1.
      - Possibility 2: A possibility 2.

    - argument2 -- Short description. This argument represents a long description. It can be:
      - Possibility 1: A possibility 1.
      - Possibility 2: A possibility 2.
  """

  def __init__(self):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    pass

  #############################################################################
  # Delineating
  #############################################################################

  def delineate_tads(self):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Delineate TADs according to random, distance to diagonal and distance to important points

  def delineate_compartments(self):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Delineate compartments according to random, distance to diagonal and distance to important points




  #############################################################################
  # Filling
  #############################################################################

  def fill_tads(self):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Fill TADs according to random, distance to diagonal and distance to important points

  def fill_compartments(self):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Fill compartments according to random, distance to diagonal and distance to important points


  #############################################################################
  # Writing
  #############################################################################

  def write_tads(self):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Write the TADs based on the apexDict

  def write_compartments(self):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Write compartments based on the compDict
















class GOBA():
  """This class represents TODO.

  *Keyword arguments:*

    - argument1 -- Short description. This argument represents a long description. It can be:
      - Possibility 1: A possibility 1.
      - Possibility 2: A possibility 2.

    - argument2 -- Short description. This argument represents a long description. It can be:
      - Possibility 1: A possibility 1.
      - Possibility 2: A possibility 2.
  """

  def __init__(self, ncpu, contact_map):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Class objects
    self.contact_map = contact_map

    # Auxiliary parameters
    self.ncpu = ncpu
    self.process_queue = []

    # Utilitary objects
    self.error_handler = ErrorHandler()

  def unshear(D, g):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
   
    # Real and imaginary parts
    D1, D2 = D[..., 0], D[..., 1]
    g1, g2 = g[0], g[1]
    a, b = D1 - g1, D2 - g2
    c, d = (1.0 - g1*D1 - g2*D2), g2*D1 - g1*D2

    # Shear out division
    out = np.empty_like(D)
    den = c**2 + d**2
    out[..., 0] = (a*c + b*d)/den
    out[..., 1] = (b*c - a*d)/den

    # Return objects
    return out

  def draw_g_1d_weak_shear(D, phi, label):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # 1D with Gaussian (canonical)
    Lam = 0.0  
    eta = 0.0
    for i, ph in enumerate(phi):
    index = np.nonzero(label == i)
    Lam += len(index[0])/ph
    eta += np.sum(D[index]/ph)
    var = 1./Lam
    mu = eta*var

    # Return objects
    return np.random.normal(loc=mu, scale=np.sqrt(var))


  def draw_g_2d_weak_shear(D, phi, label):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # 2D with Gaussian (canonical)
    Lam = 0.0  
    eta = 0.0
    for i, ph in enumerate(phi):
    index = np.nonzero(label == i)
    Lam += len(index[0])/ph
    eta += np.sum(D[index]/ph, axis=0)
    var = 1./Lam
    mu = eta*var

    # Return objects
    return np.random.multivariate_normal(mean=mu, cov=var*np.eye(2))



class Linear1DShear():
  """This class represents TODO.

  *Keyword arguments:*

    - argument1 -- Short description. This argument represents a long description. It can be:
      - Possibility 1: A possibility 1.
      - Possibility 2: A possibility 2.

    - argument2 -- Short description. This argument represents a long description. It can be:
      - Possibility 1: A possibility 1.
      - Possibility 2: A possibility 2.
  """

  def __init__(self, g):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    self.g = g

  def init(self, D):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    self.g = np.mean(D)

  def __call__(self, D):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    return D - self.g

  def unmanip(self, D):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    return D + self.g

  def update(self, D, phi, label, prior):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    self.g = draw_g_1d_weak_shear(D, phi, label)



class WeakShear():
  """This class represents TODO.

  *Keyword arguments:*

    - argument1 -- Short description. This argument represents a long description. It can be:
      - Possibility 1: A possibility 1.
      - Possibility 2: A possibility 2.

    - argument2 -- Short description. This argument represents a long description. It can be:
      - Possibility 1: A possibility 1.
      - Possibility 2: A possibility 2.
  """

  def __init__(self, g):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    self.g = g

  def init(self, D):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    self.g = np.mean(D, axis=0)

  def __call__(self, D):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    return D - self.g

  def unmanip(self, D):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    return D + self.g

  def update(self, D, phi, label, prior):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    self.g = draw_g_2d_weak_shear(D, phi, label)



class Shear():
  """This class represents TODO.

  *Keyword arguments:*

    - argument1 -- Short description. This argument represents a long description. It can be:
      - Possibility 1: A possibility 1.
      - Possibility 2: A possibility 2.

    - argument2 -- Short description. This argument represents a long description. It can be:
      - Possibility 1: A possibility 1.
      - Possibility 2: A possibility 2.
  """

  def __init__(self, g):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    self.g = g
    self.Nproposals = 0
    self.Nacceptances = 0

  def init(self, D):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    self.g = np.mean(D, axis=0)

  def __call__(self, D):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    return unshear(D, self.g)

  def unmanip(self, D):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    return unshear(D, -self.g)

  def update(self, D, phi, label, prior):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    prop_g = np.random.multivariate_normal(mean=self.g, cov=np.eye(2)*0.003**2)
    current_e_int = unshear(D, self.g)
    prop_e_int = unshear(D, prop_g)
    current_lnlike = 0.0
    prop_lnlike = 0.0
    for i, ph in enumerate(phi):
      index = label == i
      current_lnlike += prior.lnlikelihood(current_e_int[index], ph)
      prop_lnlike += prior.lnlikelihood(prop_e_int[index], ph)

    if prop_lnlike > current_lnlike:
      self.g = prop_g
      self.Nacceptances += 1
    else:
      u = np.random.uniform()
      if u < np.exp(prop_lnlike - current_lnlike):
        self.g = prop_g
        self.Nacceptances += 1
        self.Nproposals += 1







