#include "algo_unpacked.h"
#include "firmware/tau_nn.h"
#include "GlobalCorrelator_HLS/firmware/simple_fullpfalgo.h"
#include "GlobalCorrelator_HLS/puppi/firmware/simple_puppi.h"

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
#define NREGIONS 2
#define NTAU 2
#define N_CH_L1 24
#define MINPT 5

template<int LENGTH>
unsigned int ones_lut_(ap_uint<LENGTH> ornum) {
#pragma HLS function_instantiate variable=ornum
#pragma HLS inline
    ap_uint<LENGTH> count = 0;
    for (int i = 0; i<LENGTH; i++) {
    #pragma HLS UNROLL
        if (ornum > (1<<i)-1) count++;
    }
    return (ornum == 0 ? (unsigned int)(LENGTH-1) : (unsigned int)(LENGTH-count));
}

template<int LENGTH>
unsigned int get_first1(ap_uint<LENGTH> bits) {
    ap_uint<LENGTH> ones = 0;
    for (unsigned int i = 0; i < LENGTH; i++) {
    #pragma HLS UNROLL
        if (bits[i] == 1) ones |= (ap_uint<LENGTH>)(1<<(LENGTH-i)-1);
    }
    return ones_lut_(ones);
}

template<int NPFCH>
void make_inputs(input_t nn_data[NPFCH*8], PFChargedObj pfch[NREGIONS][NTRACK]) {
    unsigned int curind[NREGIONS];
    ap_uint<NTRACK> indr[NREGIONS];
    #pragma HLS ARRAY_PARTITION variable=indr complete
    for (int i = 0; i < NREGIONS; i++) {
        #pragma HLS UNROLL
        curind[i] = 0;
        for (int j = 0; j < NTRACK; j++) {
        #pragma HLS UNROLL
            if (pfch[i][j].hwPt>MINPT) {
                indr[i][j] = 1;
            } else {
                indr[i][j] = 0;
            }
        }
    }
    int theind[NPFCH];
    for (int i = 0; i < NPFCH; i++) {
#pragma HLS PIPELINE II=1
        for (int r = 0; r < NREGIONS; r++) {
            curind[r] = get_first1<NTRACK>(indr[r]);
        }
        theind[i] = (pfch[0][curind[0]].hwPt > pfch[1][curind[1]].hwPt ? 0 : 1); //this only works for NREGIONS = 2...
        nn_data[i*8+0] = input_t(pfch[theind[i]][curind[theind[i]]].hwPt);
        nn_data[i*8+1] = input_t(pfch[theind[i]][curind[theind[i]]].hwEta);
        nn_data[i*8+2] = input_t(pfch[theind[i]][curind[theind[i]]].hwPhi);
        nn_data[i*8+3] = input_t(0);
        nn_data[i*8+4] = input_t(pfch[theind[i]][curind[theind[i]]].hwId == 3 ? 1 : 0);
        nn_data[i*8+5] = input_t(0);
        nn_data[i*8+6] = input_t(0);
        nn_data[i*8+7] = input_t(pfch[theind[i]][curind[theind[i]]].hwId == 0 ? 1 : 0);
        indr[theind[i]][curind[theind[i]]]=0;
    }
}

template<int NPFPHO>
void make_inputs(input_t nn_data[NPFPHO*8], PFNeutralObj pfpho[NREGIONS][NPHOTON]) {
    unsigned int curind[NREGIONS];
    ap_uint<NPHOTON> indr[NREGIONS];
    #pragma HLS ARRAY_PARTITION variable=indr complete
    for (int i = 0; i < NREGIONS; i++) {
        #pragma HLS UNROLL
        curind[i] = 0;
        for (int j = 0; j < NPHOTON; j++) {
        #pragma HLS UNROLL
            if (pfpho[i][j].hwPt>MINPT) {
                indr[i][j] = 1;
            } else {
                indr[i][j] = 0;
            }
        }
    }
    int theind[NPFPHO];
    for (int i = 0; i < NPFPHO; i++) {
#pragma HLS PIPELINE II=1
        for (int r = 0; r < NREGIONS; r++) {
            curind[r] = get_first1<NPHOTON>(indr[r]);
        }
        theind[i] = (pfpho[0][curind[0]].hwPt > pfpho[1][curind[1]].hwPt ? 0 : 1); //this only works for NREGIONS = 2...
        nn_data[i*8+0] = input_t(pfpho[theind[i]][curind[theind[i]]].hwPt);
        nn_data[i*8+1] = input_t(pfpho[theind[i]][curind[theind[i]]].hwEta);
        nn_data[i*8+2] = input_t(pfpho[theind[i]][curind[theind[i]]].hwPhi);
        nn_data[i*8+3] = input_t(1);
        nn_data[i*8+4] = input_t(0);
        nn_data[i*8+5] = input_t(0);
        nn_data[i*8+6] = input_t(0);
        nn_data[i*8+7] = input_t(0);
        indr[theind[i]][curind[theind[i]]] = 0;
    }
}

template<int NPFNE>
void make_inputs(input_t nn_data[NPFNE*8], PFNeutralObj pfne[NREGIONS][NSELCALO]) {
    unsigned int curind[NREGIONS];
    ap_uint<NSELCALO> indr[NREGIONS];
    #pragma HLS ARRAY_PARTITION variable=indr complete
    for (int i = 0; i < NREGIONS; i++) {
        #pragma HLS UNROLL
        curind[i] = 0;
        for (int j = 0; j < NSELCALO; j++) {
        #pragma HLS UNROLL
            if (pfne[i][j].hwPt>MINPT) {
                indr[i][j] = 1;
            } else {
                indr[i][j] = 0;
            }
        }
    }
    int theind[NPFNE];
    for (int i = 0; i < NPFNE; i++) {
#pragma HLS PIPELINE II=1
        for (int r = 0; r < NREGIONS; r++) {
            curind[r] = get_first1<NSELCALO>(indr[r]);
        }
        theind[i] = (pfne[0][curind[0]].hwPt > pfne[1][curind[1]].hwPt ? 0 : 1); //this only works for NREGIONS = 2...
        nn_data[i*8+0] = input_t(pfne[theind[i]][curind[theind[i]]].hwPt);
        nn_data[i*8+1] = input_t(pfne[theind[i]][curind[theind[i]]].hwEta);
        nn_data[i*8+2] = input_t(pfne[theind[i]][curind[theind[i]]].hwPhi);
        nn_data[i*8+3] = input_t(0);
        nn_data[i*8+4] = input_t(0);
        nn_data[i*8+5] = input_t(0);
        nn_data[i*8+6] = input_t(1);
        nn_data[i*8+7] = input_t(0);
        indr[theind[i]][curind[theind[i]]] = 0;
    }
}

template<int NPFMU>
void make_inputs(input_t nn_data[NPFMU*8], PFChargedObj pfmu[NREGIONS][NMU]) {
    unsigned int curind[NREGIONS];
    ap_uint<NMU> indr[NREGIONS];
    #pragma HLS ARRAY_PARTITION variable=indr complete
    for (int i = 0; i < NREGIONS; i++) {
        #pragma HLS UNROLL
        curind[i] = 0;
        for (int j = 0; j < NMU; j++) {
        #pragma HLS UNROLL
            if (pfmu[i][j].hwPt>MINPT) {
                indr[i][j] = 1;
            } else {
                indr[i][j] = 0;
            }
        }
    }
    int theind[NPFMU];
    for (int i = 0; i < NPFMU; i++) {
#pragma HLS PIPELINE II=1
        for (int r = 0; r < NREGIONS; r++) {
            curind[r] = get_first1<NMU>(indr[r]);
        }
        theind[i] = (pfmu[0][curind[0]].hwPt > pfmu[1][curind[1]].hwPt ? 0 : 1); //this only works for NREGIONS = 2...
        nn_data[i*8+0] = input_t(pfmu[theind[i]][curind[theind[i]]].hwPt);
        nn_data[i*8+1] = input_t(pfmu[theind[i]][curind[theind[i]]].hwEta);
        nn_data[i*8+2] = input_t(pfmu[theind[i]][curind[theind[i]]].hwPhi);
        nn_data[i*8+3] = input_t(0);
        nn_data[i*8+4] = input_t(0);
        nn_data[i*8+5] = input_t(1);
        nn_data[i*8+6] = input_t(0);
        nn_data[i*8+7] = input_t(0);
        indr[theind[i]][curind[theind[i]]] = 0;
    }
}

void prep_data(input_t nn_data[N_INPUTS], PFChargedObj pfch[NREGIONS][NTRACK], PFNeutralObj pfpho[NREGIONS][NPHOTON], PFNeutralObj pfne[NREGIONS][NSELCALO], PFChargedObj pfmu[NREGIONS][NMU]) {
#pragma HLS PIPELINE II=3

    input_t tmpch[3*8];
    input_t tmppho[3*8];
    input_t tmpne[3*8];
    input_t tmpmu[1*8];
    #pragma HLS ARRAY_PARTITION variable=tmpch complete
    #pragma HLS ARRAY_PARTITION variable=tmppho complete
    #pragma HLS ARRAY_PARTITION variable=tmpne complete
    #pragma HLS ARRAY_PARTITION variable=tmpmu complete
    make_inputs<3>(tmpch, pfch);
    make_inputs<3>(tmppho, pfpho);
    make_inputs<3>(tmpne, pfne);
    make_inputs<1>(tmpmu, pfmu);

    for (unsigned int itmp = 0; itmp < 3*8; itmp++) {
    #pragma HLS UNROLL
        nn_data[itmp] = tmpch[itmp];
    }
    for (unsigned int itmp = 0; itmp < 3*8; itmp++) {
    #pragma HLS UNROLL
        nn_data[itmp+3*8] = tmppho[itmp];
    }
    for (unsigned int itmp = 0; itmp < 3*8; itmp++) {
    #pragma HLS UNROLL
        nn_data[itmp+6*8] = tmpne[itmp];
    }
    for (unsigned int itmp = 0; itmp < 1*8; itmp++) {
    #pragma HLS UNROLL
        nn_data[itmp+9*8] = tmpmu[itmp];
    }
}

void sel_data(PFChargedObj pfch_sel[NTAU][NREGIONS][NTRACK], PFNeutralObj pfpho_sel[NTAU][NREGIONS][NPHOTON], PFNeutralObj pfne_sel[NTAU][NREGIONS][NSELCALO], PFChargedObj pfmu_sel[NTAU][NREGIONS][NMU], PFChargedObj pfch[NREGIONS][NTRACK], PFNeutralObj pfpho[NREGIONS][NPHOTON], PFNeutralObj pfne[NREGIONS][NSELCALO], PFChargedObj pfmu[NREGIONS][NMU], etaphi_t taueta[NTAU], etaphi_t tauphi[NTAU]) {

    PFChargedObj dummyc; dummyc.hwPt = 0; dummyc.hwEta = 0; dummyc.hwPhi = 0; dummyc.hwId = 0; dummyc.hwZ0 = 0;
    PFNeutralObj dummyn; dummyn.hwPt = 0; dummyn.hwEta = 0; dummyn.hwPhi = 0; dummyn.hwId = 0; dummyn.hwPtPuppi = 0;

    for (int t = 0; t < NTAU; t++) {
        etaphi_t theeta = taueta[t];
        etaphi_t thephi = tauphi[t];
        for (int i = 0; i < NREGIONS; i++) {
            #pragma HLS UNROLL
            for (int j = 0; j < NTRACK; j++) {
            #pragma HLS UNROLL
                etaphi_t deta = theeta - pfch[i][j].hwEta;
                etaphi_t dphi = thephi - pfch[i][j].hwPhi;
                if (deta<92 && deta>-92 && dphi<92 && dphi>-92) {
                    pfch_sel[t][i][j] = pfch[i][j];
                } else {
                    pfch_sel[t][i][j] = dummyc;
                }
            }
            for (int j = 0; j < NPHOTON; j++) {
            #pragma HLS UNROLL
                etaphi_t deta = theeta - pfpho[i][j].hwEta;
                etaphi_t dphi = thephi - pfpho[i][j].hwPhi;
                if (deta<92 && deta>-92 && dphi<92 && dphi>-92) {
                    pfpho_sel[t][i][j] = pfpho[i][j];
                } else {
                    pfpho_sel[t][i][j] = dummyn;
                }
            }
            for (int j = 0; j < NSELCALO; j++) {
            #pragma HLS UNROLL
                etaphi_t deta = theeta - pfne[i][j].hwEta;
                etaphi_t dphi = thephi - pfne[i][j].hwPhi;
                if (deta<92 && deta>-92 && dphi<92 && dphi>-92) {
                    pfne_sel[t][i][j] = pfne[i][j];
                } else {
                    pfne_sel[t][i][j] = dummyn;
                }
            }
            for (int j = 0; j < NMU; j++) {
            #pragma HLS UNROLL
                etaphi_t deta = theeta - pfmu[i][j].hwEta;
                etaphi_t dphi = thephi - pfmu[i][j].hwPhi;
                if (deta<92 && deta>-92 && dphi<92 && dphi>-92) {
                    pfmu_sel[t][i][j] = pfmu[i][j];
                } else {
                    pfmu_sel[t][i][j] = dummyc;
                }
            }
        }
    }
}

void algo_unpacked(ap_uint<192> link_in[N_CH_IN], ap_uint<192> link_out[N_CH_OUT])
{

// !!! Retain these 4 #pragma directives below in your algo_unpacked implementation !!!
#pragma HLS ARRAY_PARTITION variable=link_in complete dim=0
#pragma HLS ARRAY_PARTITION variable=link_out complete dim=0
#pragma HLS PIPELINE II=3
#pragma HLS INTERFACE ap_ctrl_hs port=return

        ap_uint<6> decode[N_CH_IN] = {3,4,1,0,5,2,11,10,9,8,6,7,15,16,12,13,17,23,14,22,18,21,19,20,27,28,24,25,29,35,26,34,30,33,31,32,43,45,47,42,44,46,38,40,41,36,37,39};
        //ap_uint<6> decode[N_CH_IN] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47};

        MP7DataWord data_in[NREGIONS][NCHANN];
#pragma HLS ARRAY_PARTITION variable=data_in complete dim=0

        PFChargedObj pfch[NREGIONS][NTRACK];
        PFNeutralObj pfpho[NREGIONS][NPHOTON];
        PFNeutralObj pfne[NREGIONS][NSELCALO];
        PFChargedObj pfmu[NREGIONS][NMU];
#pragma HLS ARRAY_PARTITION variable=pfch  complete dim=0
#pragma HLS ARRAY_PARTITION variable=pfpho complete dim=0
#pragma HLS ARRAY_PARTITION variable=pfne  complete dim=0
#pragma HLS ARRAY_PARTITION variable=pfmu  complete dim=0

        for (int r = 0; r < NREGIONS; r++) {
#pragma HLS UNROLL
	    for (int lnk = 0; lnk < N_CH_L1; lnk++) {
#pragma HLS UNROLL
		data_in[r][lnk*5+0] = link_in[decode[lnk+(r*N_CH_L1)]].range( 63,  32);
		data_in[r][lnk*5+1] = link_in[decode[lnk+(r*N_CH_L1)]].range( 95,  64);
		data_in[r][lnk*5+2] = link_in[decode[lnk+(r*N_CH_L1)]].range(127,  96);
		data_in[r][lnk*5+3] = link_in[decode[lnk+(r*N_CH_L1)]].range(159, 128);
		data_in[r][lnk*5+4] = link_in[decode[lnk+(r*N_CH_L1)]].range(191, 160);
	    }
            mp7wrapped_unpack_out_necomb( data_in[r], pfch[r], pfpho[r], pfne[r], pfmu[r]);
        }

        input_t nn_data[NTAU][N_INPUTS];
        //#pragma HLS ARRAY_PARTITION variable=nn_data complete dim=1

        result_t taus[NTAU][N_OUTPUTS];
        //#pragma HLS ARRAY_PARTITION variable=taus complete dim=1

        pt_t taupt[NTAU];
        etaphi_t taueta[NTAU];
        etaphi_t tauphi[NTAU];
        #pragma HLS ARRAY_PARTITION variable=taupt  complete
        #pragma HLS ARRAY_PARTITION variable=taueta complete
        #pragma HLS ARRAY_PARTITION variable=tauphi complete

        PFChargedObj pfch_sel[NTAU][NREGIONS][NTRACK];
        PFNeutralObj pfpho_sel[NTAU][NREGIONS][NPHOTON];
        PFNeutralObj pfne_sel[NTAU][NREGIONS][NSELCALO];
        PFChargedObj pfmu_sel[NTAU][NREGIONS][NMU];
#pragma HLS ARRAY_PARTITION variable=pfch_sel  complete dim=0
#pragma HLS ARRAY_PARTITION variable=pfpho_sel complete dim=0
#pragma HLS ARRAY_PARTITION variable=pfne_sel  complete dim=0
#pragma HLS ARRAY_PARTITION variable=pfmu_sel  complete dim=0

        for (unsigned int itau = 0; itau < (NTAU < NREGIONS*NTRACK ? NTAU : NREGIONS*NTRACK); itau++) {
        #pragma HLS UNROLL
            for (unsigned int itmp = 0; itmp < N_INPUTS; itmp++) {
            #pragma HLS UNROLL
                nn_data[itau][itmp] = input_t(0);
            }
            taupt[itau] = pfch[0][itau].hwPt;//only seed from region 0
            taueta[itau] = pfch[0][itau].hwEta;
            tauphi[itau] = pfch[0][itau].hwPhi;
        }
        sel_data(pfch_sel, pfpho_sel, pfne_sel, pfmu_sel, pfch, pfpho, pfne, pfmu, taueta, tauphi);
        for (unsigned int itau = 0; itau < (NTAU < NREGIONS*NTRACK ? NTAU : NREGIONS*NTRACK); itau++) {
        #pragma HLS UNROLL
            prep_data(nn_data[itau], pfch_sel[itau], pfpho_sel[itau], pfne_sel[itau], pfmu_sel[itau]);
            tau_nn(nn_data[itau],taus[itau]);
        }

        //std::cout<<taus[0][0]<<std::endl;

	for (int lnk = 0; lnk < NTAU; lnk++) {
#pragma HLS UNROLL
		link_out[lnk].range(31, 0)  = 0;
		link_out[lnk].range(63, 32) = 0;
		link_out[lnk].range(81,  64) = taus[lnk][0].range(17,0);
		link_out[lnk].range(95,  82) = 0;
		link_out[lnk].range(101, 96) = ((ap_uint<6>)(lnk+1)).range(5,0);
		link_out[lnk].range(127, 102) = 0;
		link_out[lnk].range(143, 128) = taupt[lnk].range(15,0);
		link_out[lnk].range(159, 144) = 0;
		link_out[lnk].range(175, 160) = ((ap_int<16>)(taueta[lnk])).range(15,0);
		link_out[lnk].range(191, 176) = ((ap_int<16>)(tauphi[lnk])).range(15,0);
	}
	for (int lnk = NTAU; lnk < N_CH_IN; lnk++) {
#pragma HLS UNROLL
		link_out[lnk].range(31 ,   0) = 0;
		link_out[lnk].range(63 ,  32) = 0;
		link_out[lnk].range(95 ,  64) = 0;
		link_out[lnk].range(127,  96) = 0;
		link_out[lnk].range(159, 128) = 0;
		link_out[lnk].range(191, 160) = 0;
	}
}
