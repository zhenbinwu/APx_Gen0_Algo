## Set the top level module
set_top algo_unpacked
##
#### Add source code
add_files src/algo_unpacked.cpp
add_files src/GlobalCorrelator_HLS/firmware/simple_fullpfalgo.cpp -cflags "-DTESTMP7 -DHLS_pipeline_II=3"
add_files src/GlobalCorrelator_HLS/puppi/firmware/simple_puppi.cpp -cflags "-DTESTMP7 -DHLS_pipeline_II=3"
#
### Add testbed files
add_files -tb src/algo_unpacked_tb.cpp

### Add test input files
add_files -tb data/ttbar_pu140_inp.txt
add_files -tb data/ttbar_pu140_out_ref.txt

