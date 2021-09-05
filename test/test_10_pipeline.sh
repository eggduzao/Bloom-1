#!/bin/bash

# Arguments
input_matrix="/usr/users/egadegu/Projects/Papantonis_Bloom/Bloom/test/input/test_matrix/test3/wrong_matriano_matrix/Mariano_Sim_3M.bg"
output_matrix="/usr/users/egadegu/Projects/Papantonis_Bloom/Bloom/test/output/test_10_pipeline/output_matrix.bg"
output_contacts="/usr/users/egadegu/Projects/Papantonis_Bloom/Bloom/test/output/test_10_pipeline/output_contacts.bg"

# Options
chromosome="--chromosome=chr14"
organism="--organism=hg19"
resolution="--resolution=3000"
input_type="--input-type=sparse"
cores="--cores=4"
temporary_location="--temporary-location=/usr/users/egadegu/Projects/Papantonis_Bloom/Bloom/test/temp/test_10_pipeline/" # Change to your path (created)
output_type="--output-type=hic"
output_resolution="--output-resolution=3000"
avoid_distance="--avoid-distance=6000"
fitting_method="--fitting-method=pvalue" # percentile
seed="--seed=123"

# Hidden Options
pre_min_contig_removed_bins="--pre-min-contig-removed-bins=5"
pre_remove_threshold="--pre-remove-threshold=0.0"
sica_pvalue_threshold="--sica-pvalue-threshold=0.99"
sica_bottom_bin_ext_range="--sica-bottom-bin-ext-range=3,10"
sica_left_bin_ext_range="--sica-left-bin-ext-range=3,10"
sica_right_bin_ext_range="--sica-right-bin-ext-range=1,4"
sica_top_bin_ext_range="--sica-top-bin-ext-range=1,4"
sica_bonuscrosslb_range="--sica-bonuscrosslb-range=0.0,0.05"
sica_bonuscross_range="--sica-bonuscross-range=0.0,0.05"
sica_bonuslb_range="--sica-bonuslb-range=0.0,0.05"
goba_vertical_multiplier="--goba-vertical-multiplier=0.5,1.0"
goba_ortogonal_multiplier="--goba-ortogonal-multiplier=0.25,0.75"
goba_filling_frequency="--goba-filling-frequency=1.0"
goba_banding_value_mult_range="--goba-banding-value-mult-range=0.4,0.6"
goba_banding_further_range="--goba-banding-further-range=0.9,0.99"
goba_banding_frequency="--goba-banding-frequency=0.5"
goba_outing_value_mult_range="--goba-outing-value-mult-range=0.3,0.5"
goba_outing_further_range="--goba-outing-further-range=0.9,0.99"
goba_outing_frequency="--goba-outing-frequency=0.34"
goba_eppoch_resolution="--goba-eppoch-resolution=100000"
goba_eppoch_threshold="--goba-eppoch-threshold=5"
dpmm_random_degrade_range="--dpmm-random-degrade-range=0.01,0.02"
dpmm_degrade_multiplier="--dpmm-degrade-multiplier=0.05"
dpmm_half_length_bin_interval="--dpmm-half-length-bin-interval=1,4"
dpmm_value_range="--dpmm-value-range=0.00001,0.0001"
dpmm_random_range="--dpmm-random-range=0.000001,0.00001"
dpmm_iteration_multiplier="--dpmm-iteration-multiplier=1.0"
dpmm_em_significant_threshold="--dpmm-em-significant-threshold=10"
dpmm_em_signal_threshold="--dpmm-em-signal-threshold=1.0"
dpmm_em_avoid_distance="--dpmm-em-avoid-distance=5"
dpmm_ur_square_size="--dpmm-ur-square-size=1000000"
dpmm_ur_delete_size="--dpmm-ur-delete-size=10000000"
ifs_multiplier="--ifs-multiplier=1000"
ifs_min_matrix_threshold="--ifs-min-matrix-threshold=0"

bloom $chromosome $organism $resolution $input_type $cores $temporary_location $output_type $output_resolution $avoid_distance $fitting_method $seed $pre_min_contig_removed_bins $pre_remove_threshold $sica_pvalue_threshold $sica_bottom_bin_ext_range $sica_left_bin_ext_range $sica_right_bin_ext_range $sica_top_bin_ext_range $sica_bonuscrosslb_range $sica_bonuscross_range $sica_bonuslb_range $goba_vertical_multiplier $goba_ortogonal_multiplier $goba_filling_frequency $goba_banding_value_mult_range $goba_banding_further_range $goba_banding_frequency $goba_outing_value_mult_range $goba_outing_further_range $goba_outing_frequency $goba_eppoch_resolution $goba_eppoch_threshold $dpmm_random_degrade_range $dpmm_degrade_multiplier $dpmm_half_length_bin_interval $dpmm_value_range $dpmm_random_range $dpmm_iteration_multiplier $dpmm_em_significant_threshold $dpmm_em_signal_threshold $dpmm_em_avoid_distance $dpmm_ur_square_size $dpmm_ur_delete_size $ifs_multiplier $ifs_min_matrix_threshold $input_matrix $output_matrix $output_contacts

