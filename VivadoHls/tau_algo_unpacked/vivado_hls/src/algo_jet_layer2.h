#ifndef ALGO_LAYER2_H
#define ALGO_LAYER2_H

#define AP_INT_MAX_W 1600
#include <ap_int.h>
#include <hls_stream.h>
#include <ap_axi_sdata.h>
#include "GlobalCorrelator_HLS/firmware/simple_fullpfalgo.h"
#include "firmware/parameters.h"
#include "GlobalCorrelator_HLS/puppi/firmware/simple_puppi.h"

#define DATA_SIZE   128
#define NJET        12
#define NPUPPI      256
#define NPUPPI_REG  50
#define NPHIREGIONS 9



// PUPPI Value, from Javier3
// https://www.dropbox.com/sh/sbyea89yl0nawnm/AACA3cq5FpGhAVOYTKk5qABba?dl=0
#define NTRACK_PUPPI 20
#define NCALO_PUPPI  10
#define NMU_PUPPI 2
#define NEMCALO_PUPPI 12 

#define NREGIONS 36
#define NPART 25
#define DEPTH NREGIONS*2
#define NTAUPARTS  10
#define DRCONE 8410
#define DR2MAX 10000

// 0.4/ETAPHI_SCALE 
#define DETAPHI 46
#define DETAPHI_HEX 62

static float PT_SCALE = 4.0;     // quantize in units of 0.25 GeV (can be changed)
static float ETAPHI_FACTOR = 4;  // size of an ecal crystal in phi in integer units (our choice)
static float ETAPHI_SCALE = ETAPHI_FACTOR*(180./M_PI);  // M_PI/180 is the size of an ECal crystal; we make a grid that is 4 times that size
static int16_t PHI_WRAP = 360*ETAPHI_FACTOR;            // what is 3.14 in integer
static int16_t dR_WRAP = 360*ETAPHI_FACTOR;            // what is 3.14 in integer

typedef ap_axis <64*NPART,1,1,1> axi_t;
typedef hls::stream<axi_t> stream_t;
typedef ap_uint<8> count_t;

void algo_layer2(hls::stream<axi_t> &ch_link_in,hls::stream<axi_t> &ne_link_in, hls::stream<axi_t> &em_link_in, hls::stream<axi_t> &mu_link_in, hls::stream<axi_t> &link_out);

#endif
