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

# Internal

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

  def __init__(self, matrix):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    
    # Class objects
    self.matrix = matrix

  def star(self, diag_dist_scale, diag_dist_rand_range, neighbor_dist_scale, neighbor_dist_rand_range, max_neighbors):

    # diag_dist_scale = a number that is going to be divided by the distance to the diagonal. The closest to the diagonal, the higher the value to add in the star

    # neighbor_dist_scale = a number to be divided by the distance to the neighbor. The farther the neighbor, the less signal to add.

    # max_neighbors = maximum 1D neighbor distance to sum to have the value to add

  def pretad(self):

    # Create the pre-TADs with different intensities given a distance to the diagonal.
    # The more distant, the less pre-TADs

    # Also make the pre-TADs relying on the existing starred signal density.

    # Pre-tad shape = circle?
















