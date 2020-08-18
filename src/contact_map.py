"""
Contact Map Module
===================
Placeholder.

Authors: Eduardo G. Gusmao.

"""

###################################################################################################
# Libraries
###################################################################################################

# Python
from os import getcwd
from os.path import isfile, isdir, abspath, expanduser, join, splitext, basename
from subprocess import run, DEVNULL

# Internal
from src.util import ErrorHandler, ChromosomeSizes

# External
from numpy import inf


###################################################################################################
# Basic Objects
###################################################################################################

class InputFileType():
  """This class represents TODO.

  *Keyword arguments:*

    - argument1 -- Short description. This argument represents a long description. It can be:
      - Possibility 1: A possibility 1.
      - Possibility 2: A possibility 2.

    - argument2 -- Short description. This argument represents a long description. It can be:
      - Possibility 1: A possibility 1.
      - Possibility 2: A possibility 2.
  """

  SPARSE = 0
  HIC = 1
  COOL = 2
  MCOOL = 3

class ContactMap():
  """This class represents TODO.

  *Keyword arguments:*

    - argument1 -- Short description. This argument represents a long description. It can be:
      - Possibility 1: A possibility 1.
      - Possibility 2: A possibility 2.

    - argument2 -- Short description. This argument represents a long description. It can be:
      - Possibility 1: A possibility 1.
      - Possibility 2: A possibility 2.
  """

  def __init__(self, input_contact_matrix_file_name, temporary_directory = getcwd(), input_resolution = None, input_file_type = None):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Main objects
    self.input_file_name = input_contact_matrix_file_name
    self.temporary_directory = temporary_directory
    self.input_resolution = input_resolution
    self.input_file_type = input_file_type
    self.matrix = dict()
    self.resolution = None
    
    # Auxiliary objects
    self.total_zero = 0
    self.total_nonzero = 0
    self.total_sum_nonzero = 0.0
    self.max = -inf
    self.min = inf

    # Utilitary objects
    self.error_handler = ErrorHandler()
    self.chromosome_sizes = ChromosomeSizes()

    # Load file

  #############################################################################
  # Input File Loading
  #############################################################################

  def verify_input_file(self):

    # Verify if input file exists
    if(not isfile(self.input_file_name)): pass # TODO error - input file does not exist

  def create_temporary_directory(self):

    # Creating temorary directory as: temporary folder / input file name
    if(isdir(self.temporary_directory)):
      input_contact_matrix_name = splitext(basename(self.input_file_name))[0]
      try:
        self.temporary_directory = join(self.temporary_directory, input_contact_matrix_name)
        temporary_creation_output = run(["mkdir", "-p", self.temporary_directory])
      except Exception: pass # TODO error - temporary path could not be created
    else: pass # TODO error - temporary_directory must be a path


  def load_matrix(self):

    # Verify input file type


    # Verify input resolution(s)

  
    # Load matrix based on file type and resolution
    pass


  def detect_file_type(self):
    pass

  def detect_input_resolutions(self):
    pass

  def load_matrix_from_sparse(self):
    pass

  def load_matrix_from_hic(self):
    pass

  def load_matrix_from_mcool(self):

    """
    subprocess.run(["ls", "-l"], stdout = XXX, stderr = XXX)

    time cooler dump -t pixels --header --join -r chr1 -r2 chr1 /usr/users/egadegu/Projects/Papantonis_Bloom/Data/5_dn_isHiC_Human_CM/isHiC/H9.mcool::resolutions/1000 &> 5M.txt
    """
    # Create file for each chromosome
    pass



  def load_matrix_from_cool(self):
    pass


  #############################################################################
  # Output File Writing
  #############################################################################

  def write_matrix_as_sparse(self):
    pass

  def write_matrix_as_hic(self):
    pass

  def write_matrix_as_mcool(self):
    pass

  def write_matrix_as_cool(self):
    pass

  #############################################################################
  # Binary Operations
  #############################################################################

  def create_bin_file_from_text_file(self, text_file_name, bin_file_name):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Write bin file line by line
    text_file = open(text_file_name, "r")
    bin_file = open(bin_file_name, "wb")
    for line in text_file: bin_file.write(bytearray(line, "utf-8"))
    text_file.close()
    bin_file.close()

  def create_dictionary_from_bin_file(self, bin_file_name):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Dictionary
    res_dict = dict()

    # Iterate on each byte of the file
    string_array = []
    bin_file = open(bin_file_name, "rb")
    byte = bin_file.read(1)
    while byte != b"":
      string = byte.decode("utf-8")
      if(string == "\n"):
        ss = "".join(string_array).split("\t")
        res_dict[":".join(ss[:3])] = float(ss[3])
        string_array = []
      else: string_array.append(string)
      byte = bin_file.read(1)

    # Returning objects
    return res_dict


  #############################################################################
  # Auxiliary Operations
  #############################################################################









