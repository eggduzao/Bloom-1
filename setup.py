from __future__ import print_function
"""
Installs Biokit.

Authors: Eduardo G. Gusmao.

Install Command: pip install --user .
                 pip install --user . --upgrade
"""

# Python
from os.path import join, dirname, abspath
from setuptools import setup, find_packages

###################################################################################################
# Unsupported Platforms
###################################################################################################

# TODO

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
    "src.main:main",
    [],
    []
  )
}

###################################################################################################
# Auxiliary Functions/Classes
###################################################################################################

# TODO mackay

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

# TODO

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
readme_file_name = join(dirname(abspath(__file__)), "README.md")

# Fetching long description
readme_file = open(readme_file_name, "rU")
long_description = readme_file.read()
readme_file.close()

# Setup Function
setup(name=package_name,
      version=package_version,
      description=short_description,
      long_description=long_description,
      classifiers=classifiers_list,
      keywords=", ".join(keywords_list),
      author=", ".join(author_list),
      author_email=corresponding_mail,
      license=license_type,
      packages=find_packages(),
      entry_points=current_entry_points,
      install_requires=current_install_requires
)

###################################################################################################
# Termination
###################################################################################################

# TODO


