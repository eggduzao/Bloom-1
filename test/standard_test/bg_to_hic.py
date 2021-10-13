from __future__ import print_function

# Python
import os
import sys

# Internal
from bloom.io import InputFileType, IO
from bloom.util import ErrorHandler, JuicerCommand, CoolerCommand, ChromosomeSizes

# Input
chromosome = sys.argv[1]
organism = sys.argv[2]
resolution = int(sys.argv[3])
bg_file_name = sys.argv[4]
temporary_location = sys.argv[5]
hic_file_name = sys.argv[6]

# Initialization of handlers
error_handler = ErrorHandler()
juicer_command = JuicerCommand()
cooler_command = CoolerCommand()
chromosome_sizes = ChromosomeSizes(organism)

# Reading
io = IO(bg_file_name, temporary_location, organism, 4, error_handler, chromosome_sizes, juicer_command, cooler_command, resolution, input_file_type = InputFileType.SPARSE)
contact_map = io.read()
contact_map.valid_chromosome_list = [chromosome]

#Writing
io.write(contact_map, hic_file_name, InputFileType.HIC)

