from __future__ import print_function

# Python
import os
import sys

# Internal
from bloom.io import InputFileType, IO

# Input
bg_file_name = "/usr/users/egadegu/Projects/Papantonis_Bloom/Bloom/test/input/test_matrix/hic_bg_conv/Mariano_Sim_3M.bg"
temporary_location = "/usr/users/egadegu/Projects/Papantonis_Bloom/Bloom/test/temp/hic_to_bg/"
hic_file_name = "/usr/users/egadegu/Projects/Papantonis_Bloom/Bloom/test/input/test_matrix/hic_bg_conv/Mariano_Sim_3M.hic"

# Parameters
chromosome = "chr14"
organism = "hg19"
resolution = 3000

# Reading
io = IO(hic_file_name, temporary_location, organism, 4, resolution, input_file_type = InputFileType.HIC)
contact_map = io.read()
contact_map.valid_chromosome_list = [chromosome]

#Writing
io.write(contact_map, bg_file_name, InputFileType.SPARSE)

