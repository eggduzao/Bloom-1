from __future__ import print_function
"""
Setup Module
===================
Placeholder.

Authors: Eduardo G. Gusmao. 

Install Command: pip3.9 install --user .
                 pip3.9 install --user . --upgrade
"""

# Python
import os
import io
import re
import sys
import pwd
import shutil
import optparse
import setuptools


###################################################################################################
# TODO Big Unsolved TODOs TODO
###################################################################################################

# Make chromosomes general based on the name provided by the user.
# Optimization at cpu level.
# Optimization at gpu level.

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
common_deps = ["numpy>=1.12.0",
               "scipy>=1.4.0"]

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

# PassThroughOptionParser Class
class PassThroughOptionParser(optparse.OptionParser):
  """
  An 'unknown option' pass-through implementation of OptionParser.
  When unknown arguments are encountered, bundle with largs and try again,
  until rargs is depleted.
  sys.exit(status) will still be called if a known argument is passed
  incorrectly (e.g. missing arguments or bad argument types, etc.)
  """

  def _process_args(self, largs, rargs, values):
    while rargs:
      try:
        optparse.OptionParser._process_args(self, largs, rargs, values)
      except (optparse.BadOptionError, optparse.AmbiguousOptionError) as err:
        largs.append(err.opt_str)

# recursive_chown_chmod Function
def recursive_chown_chmod(path_to_walk, uid, gid, file_permission, path_permission):
  """
  Recursively applies chown from path.
  """
  for root_dir, directory_list, file_list in os.walk(path_to_walk):
    os.chown(root_dir, uid, gid)
    os.chmod(root_dir, path_permission)
    for f in file_list:
      current_complete_file = os.path.join(root_dir, f)
      os.chown(current_complete_file, uid, gid)
      os.chmod(current_complete_file, file_permission)

def read(*names, **kwargs):
  with io.open(
    os.path.join(os.path.dirname(__file__), *names),
    encoding=kwargs.get("encoding", "utf8")
  ) as fp:
    return fp.read()

def find_version(*file_paths):
  version_file = read(*file_paths)
  version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]", version_file, re.M)
  if version_match:
    return version_match.group(1)
  raise RuntimeError("Unable to find version string.")


###################################################################################################
# Processing Input Arguments
###################################################################################################

# Parameters
current_version = find_version("bloom", "__version__.py")
usage_message = "python setup.py install [python options] [Bloom options]"
version_message = "Bloom. Version: " + str(current_version)

# Initializing Option Parser
parser = PassThroughOptionParser(usage = usage_message, version = version_message)

# Parameter: Copy bloom data folder
param_copy_bloom_data_name = "--copy-bloom-data"
param_copy_bloom_data_help = "Explain here."
parser.add_option(param_copy_bloom_data_name, dest="param_copy_bloom_data", action="store_true", default=False, help=param_copy_bloom_data_help)

# Processing Options
options, arguments = parser.parse_args()
param_copy_bloom_data = options.param_copy_bloom_data

# Manually Removing Additional Options from sys.argv
new_sys_argv = []
for e in sys.argv:
  if param_copy_bloom_data_name == e[:len(param_copy_bloom_data_name)]:
    continue
  new_sys_argv.append(e)
sys.argv = new_sys_argv

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

# Creating data path # TODO
if(param_copy_bloom_data):
  if(os.path.exists(bloom_data_path)):
    shutil.rmtree(bloom_data_path)
  shutil.copytree(current_path, bloom_data_path)


###################################################################################################
# Setup Function
###################################################################################################

# Parameters
package_name = "Bloom"
package_version = str(current_version)
short_description = "Computational Framework to reveal occult patterns in 3C contact matrices"
classifiers_list = ["Topic :: Scientific/Engineering :: Bio-Informatics",
                    "Topic :: Scientific/Engineering :: Artificial Intelligence"]
keywords_list = ["NGS", "Hi-C", "DPMM"]
author_list = ["Eduardo G. Gusmao"]
corresponding_mail = "eduardo.gusmao@med.uni-goettingen.de"
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

# Modifying Permissions when Running Superuser/Admin
if(param_copy_bloom_data):

  # Get current user and set default permissions for files to be visible and binaries executable
  current_user = os.getenv("SUDO_USER")
  default_file_permission = 0o644
  default_path_permission = 0o755

  # Setting the permissions
  if current_user:
    current_user_uid = pwd.getpwnam(current_user).pw_uid
    current_user_gid = pwd.getpwnam(current_user).pw_gid
    recursive_chown_chmod(bloom_data_path, current_user_uid, current_user_gid, default_file_permission,
                          default_path_permission)


