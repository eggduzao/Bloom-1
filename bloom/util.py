"""
Util Module
===================
The Util classes contains many utilities needed by other classes such as the paths to input files.
"""

###################################################################################################
# Libraries
###################################################################################################

# Python
import os
import gc
import sys
import codecs
import optparse
import traceback
import subprocess
import configparser
import multiprocessing

# Internal

# External




"""

=================================

os.path.isfile - file exists
os.path.isdir - directory exists
os.path.abspath - absolute path
os.path.expanduser - expand user (~)
os.path.join - join paths
os.path.splitext - removes extension of file
os.path.basename - return the last thing after a "/"

=================================

def mysleep(number):
    command = ["sleep", str(number) + "s"]
    out_process = subprocess.run(command, stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)
    return number

--

pool = mp.Pool(mp.cpu_count())

results = pool.starmap(mysleep, [(seconds,) for i in list(range(0, times))])

pool.close()
pool.join()

=================================

def worker(f1, f2):
    os.system("run program x on f1 and f2")

def test_run(pool):
     filelist = os.listdir(files_dir)
     for f1 in filelist:
          for f2 in filelist:
               pool.apply_async(worker, args=(f1, f2))

--

if __name__ == "__main__":
     import multiprocessing as mp
     pool = mp.Pool(NUM_CPUS)
     test_run(pool)
     pool.close()
     pool.join()

=================================

  def wrapper_dump(self, args):

    # Wrapper
    self.juicer_dump(*args)
    
=================================

"""

###################################################################################################
# Configuration File Handling
###################################################################################################

class ConfigurationFile:
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

    # Fetching bloom_data_path
    self.bloom_data_path = os.path.expanduser(os.getenv("RGTDATA", default = os.path.join(os.getenv("HOME"), "bloom_data")))

    # Reading config file directory
    self.bloom_config_file_name = os.path.join(self.bloom_data_path, "data.config")

    # Parsing config file
    self.config = configparser.ConfigParser()
    self.config.read_file(codecs.open(self.bloom_config_file_name, "rU", "utf8"))


class ChromosomeSizes(ConfigurationFile):
  """This class represents TODO.

  *Keyword arguments:*

    - argument1 -- Short description. This argument represents a long description. It can be:
      - Possibility 1: A possibility 1.
      - Possibility 2: A possibility 2.

    - argument2 -- Short description. This argument represents a long description. It can be:
      - Possibility 1: A possibility 1.
      - Possibility 2: A possibility 2.
  """

  def __init__(self, organism):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """

    # Configuration file initialization
    ConfigurationFile.__init__(self)
    self.organism = organism
    self.chromosome_sizes_file_name = os.path.join(self.bloom_data_path, self.config.get("ChromosomeSizes", organism))

    # Creating chromosome sizes dictionary and chromosome list
    self.chromosome_sizes_dictionary = dict()
    chrom_sizes_file = codecs.open(self.chromosome_sizes_file_name, "rU", "utf8")
    for line in chrom_sizes_file:
      ll = line.strip().split("\t")
      self.chromosome_sizes_dictionary[ll[0]] = int(ll[1])
    chrom_sizes_file.close()
    self.chromosome_sizes_list = sorted(self.chromosome_sizes_dictionary.keys())


###################################################################################################
# Argument Parsing
###################################################################################################

class HelpfulOptionParser(optparse.OptionParser):
  """This class represents an OptionParser that prints full help on errors. Inherits OptionParser.

  *Keyword arguments:*

    - OptionParser -- Inherited OptionParser object.
  """

  def error(self, msg):
    """Error handling.
    
    *Keyword arguments:*
    
      - msg -- String containing the error message.
    
    *Return:*
    
      - return -- An error message to the user.
    """
    self.print_help(sys.stderr)
    self.exit(2, "\n%s: error: %s\n" % (self.get_prog_name(), msg))


class PassThroughOptionParser(HelpfulOptionParser):
  """When unknown arguments are encountered, bundle with largs and try again, until rargs is depleted.
     sys.exit(status) will still be called if a known argument is passed incorrectly (e.g. missing arguments or
     bad argument types, etc.). Inherits HelpfulOptionParser.

  *Keyword arguments:*

    - HelpfulOptionParser -- Inherited HelpfulOptionParser object.
  """

  def _process_args(self, largs, rargs, values):
    """Overwrites _process_args to achieve desired effect.
    
    *Keyword arguments:*
    
      - largs -- The tool's optional arguments.
      - rargs -- The tool's required arguments.
    """
    while rargs:
      try:
        HelpfulOptionParser._process_args(self, largs, rargs, values)
      except (optparse.BadOptionError, optparse.AmbiguousOptionError):
        pass
        # largs.append(err.opt_str)


###################################################################################################
# Warnings and Errors Handling
###################################################################################################

class ErrorHandler:
  """Handles errors and warnings in a standardized way.

  *Error Dictionary Standard:*

    Each entry consists of a key+list in the form X:[Y,Z,W] where:
      - X -- The key representing the internal error name.
      - Y -- Error number.
      - Z -- Exit status.
      - W -- Error message to be print.

  *Warning Dictionary Standard:*

    Each entry consists of a key+list in the form X:[Y,Z] where:
      - X -- The key representing the internal warning name.
      - Y -- Warning number.
      - Z -- Warning message to be print.
  """

  def __init__(self):
    """Initializes required objects for warning and error handling.
    """

    self.program_name = os.path.basename(sys.argv[0])

    self.error_dictionary = {
      "DEFAULT_ERROR": [0, 0, "Undefined error. Program terminated with exit status 0."],
      "PLACEHOLDER": [1, 0, "You must define one specific analysis. Run '" + self.program_name + " -h' for help."]
    }
    self.error_number = 0
    self.exit_status = 1
    self.error_message = 2

    self.warning_dictionary = {
      "DEFAULT_WARNING": [0, "Undefined warning."],
      "PLACEHOLDER": [1, "Placeholder."]
    }
    self.warning_number = 0
    self.warning_message = 1

  def throw_error(self, error_type, add_msg=""):
    """Throws the specified error type. If the error type does not exist, throws a default error message and exits.
    
    *Keyword arguments:*
    
      - error_type -- Error type.
      - add_msg -- Message to add to the error.
    """

    # Fetching error type
    try:
      error_number = self.error_dictionary[error_type][self.error_number]
      exit_status = self.error_dictionary[error_type][self.exit_status]
      error_message = self.error_dictionary[error_type][self.error_message]
    except (KeyError, IndexError):
      error_number = self.error_dictionary["DEFAULT_ERROR"][self.error_number]
      exit_status = self.error_dictionary["DEFAULT_ERROR"][self.exit_status]
      error_message = self.error_dictionary["DEFAULT_ERROR"][self.error_message]

    # Handling error
    complete_error_message = ("--------------------------------------------------\n"
                              "Error Number: " + str(error_number) + ".\n"
                              "Program: " + self.program_name + ".\n"
                              "Report: " + error_message + " " + add_msg + "\n"
                              "Behaviour: The program will quit with exit status " + str(exit_status) + ".\n"
                              "--------------------------------------------------")
    print(complete_error_message, file=sys.stderr)
    traceback.print_exc()
    sys.exit(exit_status)

  def throw_warning(self, warning_type, add_msg=""):
    """Throws the specified warning type. If the warning type does not exist, throws a default warning message and exits.
    
    *Keyword arguments:*
    
      - warning_type -- Warning type.
      - add_msg -- Message to add to the error.
    """

    # Fetching warning type
    try:
      warning_number = self.warning_dictionary[warning_type][self.warning_number]
      warning_message = self.warning_dictionary[warning_type][self.warning_message]
    except (KeyError, IndexError):
      warning_number = self.warning_dictionary["DEFAULT_WARNING"][self.warning_number]
      warning_message = self.warning_dictionary["DEFAULT_WARNING"][self.warning_message]

    # Handling warning
    complete_warning_message = ("--------------------------------------------------\n"
                                "Warning Number: " + str(warning_number) + ".\n"
                                "Program: " + self.program_name + ".\n"
                                "Report: " + warning_message + " " + add_msg + "\n"
                                "--------------------------------------------------")
    print(complete_warning_message, file=sys.stderr)


###################################################################################################
# Auxiliary Functions as Static Methods
###################################################################################################

class AuxiliaryFunctions:
  """Class for small auxiliary functions.
  """


# Add check int, float, string, path, file to Util

  @staticmethod
  def string_is_int(s):
    """Verifies whether a string is a representation of a numeric integer.
    
    *Keyword arguments:*
    
      - s -- String to be verified.
    
    *Return:*
    
      - return -- True if it is a numeric integer or false if it is not.
    """
    try:
      int(s)
      return True
    except ValueError:
      return False

  @staticmethod
  def string_is_float(s):
    """Verifies whether a string is a representation of a numeric float.
    
    *Keyword arguments:*
    
      - s -- String to be verified.
    
    *Return:*
    
      - return -- True if it is a numeric float or false if it is not.
    """
    try:
      float(s)
      return True
    except ValueError:
      return False

  @staticmethod
  def string_is_validfile(s):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    # Initialize validity flag
    is_valid = False

    # Check whether file exists and is accessible
    if(os.path.isfile(s)):
      try:
        input_file = codecs.open(s, "rU", "utf8")
        input_file.close()
        is_valid = True
      except Exception:
        is_valid = False

    # Return validity flag
    return is_valid

  @staticmethod
  def string_is_validpath(s):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    # Initialize validity flag
    is_valid = False

    # Check whether path exists and is accessible
    if(os.path.isdir(s)):
      try:
        temp_file_name = os.path.join(s, "temp_file_name.txt")
        input_file = codecs.open(temp_file_name, "w", "utf8")
        input_file.close()
        remove_command = ["rm", "-rf", temp_file_name]
        remove_process = subprocess.run(remove_command , stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)
        is_valid = True
      except Exception:
        is_valid = False

    # Return validity flag
    return is_valid

  @staticmethod
  def overlap_count(interval1, interval2):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    return max(0, min(interval1[1], interval2[1]) - max(interval1[0], interval2[0]))

  @staticmethod
  def revcomp(s):
    """Reverse complement a DNA (A, C, G, T, N) sequence.
    
    *Keyword arguments:*
    
      - s -- DNA sequence string to be reverse complemented.
    
    *Return:*
    
      - return -- The reverse complement of the sequence, if it exists.
    """
    revDict = dict([("A", "T"), ("T", "A"), ("C", "G"), ("G", "C"), ("N", "N")])
    return "".join([revDict[e] for e in s[::-1]])
    # TODO - Add error to unknown alphabet.

  @staticmethod
  def ceil_multiple(num, divisor):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    return num + (divisor - (num%divisor))   

  @staticmethod
  def floor_multiple(num, divisor):
    """Returns TODO.
    
    *Keyword arguments:*
    
      - argument -- An argument.
    
    *Return:*
    
      - return -- A return.
    """
    return num - (num%divisor)
    
