from __future__ import print_function

# Python
import os
import sys

# Internal
from bloom.io import InputFileType, IO

# Input
bg_file_name = "/usr/users/egadegu/Projects/Papantonis_Bloom/Bloom/test/input/test_matrix/test3/wrong_matriano_matrix/Mariano_Sim_3M.bg"
bg_file_name = "/usr/users/egadegu/Projects/Papantonis_Bloom/Bloom/test/input/test_matrix/test3/Mariano_long100-500mod_56_3M.bg"
temporary_location = "/usr/users/egadegu/Projects/Papantonis_Bloom/Bloom/test/temp/hic_to_bg/"

# Parameters
chromosome = "chr14"
organism = "hg19"
resolution = 3000

# Reading
io = IO(bg_file_name, temporary_location, organism, 4, resolution, input_file_type = InputFileType.SPARSE)
contact_map = io.read()
contact_map.valid_chromosome_list = [chromosome]
contact_map.calculate_statistics_full()

# Simple sparsity
sparsity = contact_map.total_zero_bins[chromosome] / contact_map.total_bins[chromosome]
print(sparsity)

