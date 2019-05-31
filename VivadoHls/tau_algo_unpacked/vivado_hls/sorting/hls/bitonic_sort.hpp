#ifndef __BITONIC_SORT_HPP__
#define __BITONIC_SORT_HPP__

// avoid Cosim Error
//#include "/opt/Xilinx/Vivado/2018.2/include/gmp.h"
//#include "/opt/Xilinx/Vivado/2018.2/include/mpfr.h"

#include "../../src/GlobalCorrelator_HLS/firmware/simple_fullpfalgo.h"
#include "ap_int.h"
#include "hls_stream.h"


#define DATA_SIZE 128
#define DATA_W 14

void bitonic_sort(hls::stream<PFChargedObj > data_in[], hls::stream<PFChargedObj > data_out[]);
//void bitonic_sort(PFChargedObj datas1[2][DATA_SIZE]);//,PFChargedObj datas2[DATA_SIZE]);

#endif
