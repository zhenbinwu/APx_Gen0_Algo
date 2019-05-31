############################################################
## This file is generated automatically by Vivado HLS.
## Please DO NOT edit it.
## Copyright (C) 1986-2018 Xilinx, Inc. All Rights Reserved.
############################################################
open_project bitonic_sort
set_top bitonic_sort
add_files ../hls/swap.hpp
add_files ../hls/swap.cpp
add_files ../hls/sorting_network.hpp
add_files ../hls/sorting_network.cpp
add_files ../hls/bitonic_sort.hpp
add_files ../hls/bitonic_sort.cpp
add_files -tb ../hls/tb_bitonic_sort.cpp -cflags "-Wno-unknown-pragmas"
add_files ../../src/GlobalCorrelator_HLS/firmware/simple_fullpfalgo.cpp -cflags "-DTESTMP7 -DHLS_pipeline_II=3"
add_files ../../src/GlobalCorrelator_HLS/puppi/firmware/simple_puppi.cpp -cflags "-DTESTMP7 -DHLS_pipeline_II=3"

open_solution "solution1"
set_part {xcvu9p-flgb2104-2-i}
create_clock -period 300MHz -name default
#source "directives.tcl"

csim_design
csynth_design
cosim_design
export_design -format ip_catalog

exit
