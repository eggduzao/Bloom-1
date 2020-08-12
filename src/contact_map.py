"""
XXX
===================
Placeholder.

"""

###################################################################################################
# Libraries
###################################################################################################

# Python
import os
import sys

# Internal

# External

###################################################################################################
# Basic Objects
###################################################################################################

class ContactMap():
  """Class description.
  """
  def __init__(self): pass

  def create_bin_file_from_text_file(text_file_name, bin_file_name):

    # Write bin file line by line
    text_file = open(text_file_name, "r")
    bin_file = open(bin_file_name, "wb")
    for line in text_file: bin_file.write(bytearray(line, "utf-8"))
    text_file.close()
    bin_file.close()

  def create_dictionary_from_bin_file(bin_file_name):

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



