from __future__ import print_function
"""
Setup Module
===================
Placeholder.

Authors: Eduardo G. Gusmao.

Install Command: pip3 install --user .
                 pip3 install --user . --upgrade
"""

# Python
import os
import sys
import shutil
import setuptools

###################################################################################################
# TODO Big Unsolved TODOs TODO
###################################################################################################

# Make chromosomes general based on the name provided by the user.

###################################################################################################
# Unsupported Platforms
###################################################################################################

supported_platforms = ["linux", "linux2", "darwin"]
if sys.platform not in supported_platforms:
    print("ERROR: This package currently supports only unix-based systems (Linux and MAC OS X).")
    sys.exit(0)


###################################################################################################
# Parameters
###################################################################################################

"""
Tools Dictionary:
  * Insert in the dictionary bellow a key+tuple X: (Y,Z,W,K) representing:
    - X: A string representing the name the user should provide for this script in order to install the tool.
    - Y: A string representing the name of the program which the user should type in the terminal in order to execute the tool.
    - Z: A string representing the path to the main function that executes the program.
    - W: A list with the package requirements for that program.
    - K: A list with external binary files that will be copied by the setup installation function to the user's bin folder.
Tools Dictionary Standard:
  * All programs should start with "bk-" followed by the program name.
  * The main function called within the script must be termed "main".
"""

# Common dependencies.
common_deps = ["cython>=0.29.0",
               "numpy>=1.4.0",
               "scipy>=1.0.0",
               "pysam>=0.15.0"]

tools_dictionary = {
  "bloom": (
    "bloom",
    "bloom.main:main",
    [],
    []
  )
}


###################################################################################################
# Auxiliary Functions/Classes
###################################################################################################

# TODO

###################################################################################################
# Processing Input Arguments
###################################################################################################

# Defining entry points
current_entry_points = {"console_scripts": []}
for tool_option in tools_dictionary.keys():
  current_entry_points["console_scripts"].append(" = ".join(tools_dictionary[tool_option][:2]))

# Defining install requirements
current_install_requires = common_deps
for tool_option in tools_dictionary.keys():
  current_install_requires += tools_dictionary[tool_option][2]


###################################################################################################
# Creating Data Path
###################################################################################################

# Default bloom_data_path
current_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "bloom_data")
bloom_data_path = os.path.join(os.path.expanduser("~"), "bloom_data")
if os.path.exists(bloom_data_path):
  shutil.rmtree(bloom_data_path)
shutil.copytree(current_path, bloom_data_path)


###################################################################################################
# Setup Function
###################################################################################################

# Parameters
package_name = "Bloom"
package_version = "0.0.1"
short_description = "Computational Framework to reveal occult patterns in 3C contact matrices"
classifiers_list = ["Topic :: Scientific/Engineering :: Bio-Informatics",
                    "Topic :: Scientific/Engineering :: Artificial Intelligence"]
keywords_list = ["NGS"]
author_list = ["Eduardo G. Gusmao"]
corresponding_mail = "eduardo.gusmao@rwth-aachen.de"
license_type = "GPL"

# Fetching additional structural files
readme_file_name = os.path.join(os.path.dirname(os.path.abspath(__file__)), "README.md")

# Fetching long description
readme_file = open(readme_file_name, "rU", encoding="utf-8")
long_description = readme_file.read()
readme_file.close()

# Setup Function
setuptools.setup(name = package_name,
                 version = package_version,
                 description = short_description,
                 long_description = long_description,
                 classifiers = classifiers_list,
                 keywords = ", ".join(keywords_list),
                 author = ", ".join(author_list),
                 author_email = corresponding_mail,
                 license = license_type,
                 packages = setuptools.find_packages(),
                 entry_points = current_entry_points,
                 install_requires = current_install_requires
)


###################################################################################################
# Termination
###################################################################################################

# TODO


