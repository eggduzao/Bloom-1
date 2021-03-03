"""
Setup Environment
===================
Placeholder.

Authors: Eduardo G. Gusmao.

"""

###################################################################################################
# Libraries
###################################################################################################

# Python
import os
import sys
import subprocess

# Internal

# External

###################################################################################################
# Functions
###################################################################################################

def download(url, output_location):

  # Linux -> wget
  if(sys.platform in ["linux", "linux2"]):
    download_command = ["wget", "-c", "--tries=0", "--read-timeout=30", url, "-P", output_location]
    download_process = subprocess.run(download_command , stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)

  # MAC OS -> curl
  elif(sys.platform in ["darwin"]):
    download_command = ["cd", output_location, "&&", "{", "curl", "-O", url, ";", "cd", "-;", "}"]
    download_process = subprocess.run(download_command , stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)

  # Return download process
  return download_process


###################################################################################################
# Execution
###################################################################################################

# Platform
supported_platforms = ["linux", "linux2", "darwin"]
if(sys.platform not in supported_platforms):
  print("ERROR: This package currently supports only unix-based systems (Linux and MAC OS X).")
  sys.exit(0)

# Environment output location
output_directory = os.path.join(os.path.dirname(os.path.realpath(__file__)), "barcode")
if(not os.path.isdir(output_directory)):
  mkdir_command = ["mkdir", "-p", output_directory]
  mkdir_process = subprocess.run(mkdir_command , stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)

# Environment url list
base_url = "http://wwwuser.gwdg.de/~egadegu/Argyris_Papantonis/Bloom/Auxiliary_Data"
file_name_list = [
  "barcode1bua", "barcode1bub", "barcode1buc", "barcode1bud", "barcode1bue", "barcode1buf", "barcode1bug", "barcode1buh", "barcode1mb", "barcode1lb", 
  "barcode2bua", "barcode2bub", "barcode2buc", "barcode2bud", "barcode2bue", "barcode2buf", "barcode2bug", "barcode2buh", "barcode2mb", "barcode2lb"
]

# Performing all downloads
for file_name in file_name_list:

  full_url = os.path.join(base_url, file_name)
  download(full_url, output_directory)

