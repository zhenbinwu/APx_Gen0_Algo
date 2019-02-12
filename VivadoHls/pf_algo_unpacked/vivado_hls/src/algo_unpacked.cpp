#include "algo_unpacked.h"
#include "GlobalCorrelator_HLS/firmware/simple_fullpfalgo.h"

/*
 * algo_unpacked interface exposes fully unpacked input and output link data.
 * This version assumes use of 10G 8b10b links, and thus providing 192bits/BX/link,
 * arranged as an array of 3x 64 bits.
 *
 * !!! N.B. Do NOT use the first bytes of link_in/out[][0] (i.e. link_in/out[][0].range(7,0)
 * as this portion is reserved for transmission of 8b10b input/output link alignment markers.
 *
 * The remaining 184 bits are available for algorithm use.
 *
 */

#define NCHANN (N_CH_OUT+N_CH_OUT+N_CH_OUT+N_CH_OUT+N_CH_OUT)

void algo_unpacked(ap_uint<192> link_in[N_CH_IN], ap_uint<192> link_out[N_CH_OUT])
{

// !!! Retain these 4 #pragma directives below in your algo_unpacked implementation !!!
#pragma HLS ARRAY_PARTITION variable=link_in complete dim=0
#pragma HLS ARRAY_PARTITION variable=link_out complete dim=0
#pragma HLS PIPELINE II=3
#pragma HLS INTERFACE ap_ctrl_hs port=return

        MP7DataWord data_in[NCHANN];
        MP7DataWord data_out[NCHANN];
#pragma HLS ARRAY_PARTITION variable=data_in complete dim=0
#pragma HLS ARRAY_PARTITION variable=data_out complete dim=0
        z0_t vertex = link_in[0].range(31, 22);

	for (int lnk = 0; lnk < N_CH_IN; lnk++) {
#pragma HLS UNROLL
                
		//data_in[lnk+N_CH_IN*0] = link_in[lnk].range( 63,  32);
		//data_in[lnk+N_CH_IN*1] = link_in[lnk].range( 95,  64);
		//data_in[lnk+N_CH_IN*2] = link_in[lnk].range(127,  96);
		//data_in[lnk+N_CH_IN*3] = link_in[lnk].range(159, 128);
		//data_in[lnk+N_CH_IN*4] = link_in[lnk].range(191, 160);
		data_in[lnk*5+0] = link_in[lnk].range( 63,  32);
		data_in[lnk*5+1] = link_in[lnk].range( 95,  64);
		data_in[lnk*5+2] = link_in[lnk].range(127,  96);
		data_in[lnk*5+3] = link_in[lnk].range(159, 128);
		data_in[lnk*5+4] = link_in[lnk].range(191, 160);
                data_out[lnk+N_CH_IN*0] = 0;
                data_out[lnk+N_CH_IN*1] = 0;
                data_out[lnk+N_CH_IN*2] = 0;
                data_out[lnk+N_CH_IN*3] = 0;
                data_out[lnk+N_CH_IN*4] = 0;
	}

        mp7wrapped_pfalgo3_full(data_in, data_out, vertex);

	for (int lnk = 0; lnk < N_CH_OUT; lnk++) {
#pragma HLS UNROLL
                
		link_out[lnk].range(31 ,   0) = 0;
		link_out[lnk].range(63 ,  32) = data_out[lnk*5+0];
		link_out[lnk].range(95 ,  64) = data_out[lnk*5+1];
		link_out[lnk].range(127,  96) = data_out[lnk*5+2];
		link_out[lnk].range(159, 128) = data_out[lnk*5+3];
		link_out[lnk].range(191, 160) = data_out[lnk*5+4];
	}
}
