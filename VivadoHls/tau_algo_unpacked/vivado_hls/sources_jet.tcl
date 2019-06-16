## Set the top level module
set_top algo_layer2
##
#### Add source code
add_files src/algo_jet_layer2.cpp -cflags "-DTESTMP7"
add_files src/bitonic_sort.cpp
add_files src/sorting_network.cpp
add_files src/GlobalCorrelator_HLS/firmware/simple_fullpfalgo.cpp -cflags "-DTESTMP7 -DHLS_pipeline_II=3"
add_files src/GlobalCorrelator_HLS/puppi/firmware/simple_puppi.cpp -cflags "-DTESTMP7 -DHLS_pipeline_II=3"
#
### Add testbed files
add_files -tb src/algo_jet_layer2_tb.cpp

### Add test input files
add_files -tb data/test1_inp.txt
add_files -tb data/test1_out_ref.txt
add_files -tb data/ttbar_pu140_inp.txt
add_files -tb data/ttbar_pu140_out_ref.txt
add_files -tb data/ttbar_pu140_comb_inp.txt
add_files -tb data/ttbar_pu140_comb_out_ref.txt

