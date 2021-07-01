from __future__ import print_function

# Python
import os
import sys

# Internal
from bloom.io import InputFileType, IO

# Input
bg_file_name = "./output_matrix.bg"
temporary_location = "./temp/"
hic_file_name = "./output_matrix.hic"

# Parameters
chromosome = "chr14"
organism = "mm9"
resolution = 10000

# Reading
io = IO(bg_file_name, temporary_location, organism, 4, resolution, input_file_type = InputFileType.SPARSE)
contact_map = io.read()
contact_map.valid_chromosome_list = ["chr14"]

#Writing
io.write(resolution, hic_file_name, InputFileType.HIC)

