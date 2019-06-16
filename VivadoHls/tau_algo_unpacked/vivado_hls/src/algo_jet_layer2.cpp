#include "algo_jet_layer2.h"
#include "sorting_network.h"

//16+10+10+3+10 for pt,eta,phi,particleid,z0
//typedef ap_axis <64*NPART,1,1,1> axi_t;
//typedef hls::stream<axi_t> hls::stream<axi_t>;
//stream depth is TM6 and 240 MHz gives 24 clocks

// ===  FUNCTION  ============================================================
//         Name:  HexagonDR
//  Description:  /* cursor */
// ===========================================================================
inline bool HexagonDR(etaphi_t eta1, etaphi_t phi1,  etaphi_t eta2,  etaphi_t phi2)
{
  //hardcode for etaphi size
  int tmpe = eta1-eta2;
  etaphi_t deta = (tmpe > 0 ? tmpe : -tmpe);
  int tmpp = phi1-phi2;
  etaphi_t dphi = (tmpp > 0 ? tmpp : -tmpp);
  if (deta < DETAPHI && dphi < DETAPHI && (deta+dphi) < DETAPHI_HEX )
    return 1;
  else
    return 0;
}       // -----  end of function HexagonDR  -----


// ===  FUNCTION  ============================================================
//         Name:  RegrionCrossMap
//  Description:  
// ===========================================================================
void RegrionCrossMap(ap_uint<8> RegXMap[NREGIONS][6])
{
  for (int i = 0; i < NREGIONS; ++i)
  {
    RegXMap[i][0] = i;
    if (i % NPHIREGIONS == 0 && (i - 2 *NPHIREGIONS) + 1 > 0)
      RegXMap[i][1] = (i - 2 *NPHIREGIONS) + 1 ;
    if ((i - (NPHIREGIONS +1))> 0) 
      RegXMap[i][2] = i- (NPHIREGIONS+1);
    if ((i - NPHIREGIONS )> 0) 
      RegXMap[i][3] = i- NPHIREGIONS;
    if ((i - (NPHIREGIONS-1) )> 0) 
      RegXMap[i][4] = i- (NPHIREGIONS-1);
    if ((i - 1) > 0)
      RegXMap[i][5] = i- 1;
  }
}       // -----  end of function RegrionCrossMap  -----

// ===  FUNCTION  ============================================================
//         Name:  SelectPUPPIinRegion
//  Description:  
// ===========================================================================
template<unsigned int NParts, unsigned int NPUPPIREG>
inline void SelectPUPPIinRegion(PFChargedObj allparts[NParts], PFChargedObj allpuppis[NPUPPIREG], count_t regsize)
{
  #pragma HLS PIPELINE
  static ap_uint<1> IsPUPPI[NParts];
  static count_t IsPUPPIcnt[NParts];

  #pragma HLS ARRAY_PARTITION variable=IsPUPPI complete
  #pragma HLS ARRAY_PARTITION variable=IsPUPPIcnt complete
  #pragma HLS ARRAY_PARTITION variable=allparts complete

  for (int i = 0; i < NParts; ++i)
  {
    if(allparts[i].hwPt > 0)
      IsPUPPI[i] = 1;
    else
      IsPUPPI[i] = 0;
  }

  // Get accumulated PUPPI Count
  for (int i = 0; i < NParts; ++i)
  {
    IsPUPPIcnt[i] = 0;
    for (int j = 0; j < i; ++j)
    {
      IsPUPPIcnt[i] += IsPUPPI[j];
    }
  }

  for (int i = 0; i < NParts; ++i)
  {
    if (IsPUPPI[i])
    {
      allpuppis[IsPUPPIcnt[i]-1] = allparts[i];
    }
  }
  
  regsize = IsPUPPIcnt[NParts-1];
}       // -----  end of function SelectPUPPIinRegion  -----

template<unsigned int N> 
inline void mp7_unpack(axi_t data, PFChargedObj emcalo[N]) {
    #pragma HLS inline
    #pragma HLS PIPELINE
    for (unsigned int i = 0; i < N; ++i) {
      int pOffset = (i)*64;
      emcalo[i].hwPt       = data.data(15+pOffset, 0+pOffset);
      emcalo[i].hwEta      = data.data(24+pOffset,16+pOffset);
      emcalo[i].hwPhi      = data.data(34+pOffset,25+pOffset);
      emcalo[i].hwId       = data.data(38+pOffset,35+pOffset);
      emcalo[i].hwZ0       = data.data(48+pOffset,39+pOffset);
    }
}

template<unsigned int N> 
inline void mp7_unpack_seed(axi_t data, PFChargedObj &seed) {
    #pragma HLS inline
    #pragma HLS PIPELINE
    seed.hwPt          = data.data(15, 0);
    seed.hwEta         = data.data(24,16);
    seed.hwPhi         = data.data(34,25);
    seed.hwId          = data.data(38,35);
    seed.hwZ0          = data.data(48,39);
}

template<unsigned int N,unsigned int NR> 
inline void mp7_unpack(axi_t data, PFChargedObj emcalo[N*NR],int iBase) {
    #pragma HLS inline
    #pragma HLS PIPELINE
    for (unsigned int i = 0; i < N; ++i) {
      int pOffset = (i)*64;
      emcalo[iBase+i].hwPt       = data.data(15+pOffset, 0+pOffset);
      emcalo[iBase+i].hwEta      = data.data(24+pOffset,16+pOffset);
      emcalo[iBase+i].hwPhi      = data.data(34+pOffset,25+pOffset);
      emcalo[iBase+i].hwId       = data.data(38+pOffset,35+pOffset);
      emcalo[iBase+i].hwZ0       = data.data(48+pOffset,39+pOffset);
    }
}

template<unsigned int N> 
inline void mp7_unpack(axi_t data, PFNeutralObj emcalo[N]) { //,int iBase) {
    #pragma HLS inline
    #pragma HLS PIPELINE
    for (unsigned int i = 0; i < N; ++i) {
      int pOffset = (i)*64;
      emcalo[i].hwPt       = data.data(15+pOffset, 0+pOffset);
      emcalo[i].hwEta      = data.data(24+pOffset,16+pOffset);
      emcalo[i].hwPhi      = data.data(34+pOffset,25+pOffset);
      emcalo[i].hwId       = data.data(38+pOffset,35+pOffset);
      emcalo[i].hwPtPuppi  = data.data(55+pOffset,39+pOffset);
    }
}
template<unsigned int N, unsigned int OFFS> 
inline void mp7_pack(PFChargedObj obj[N], axi_t &data) {
    #pragma HLS inline
    for (unsigned int i = 0; i < N; ++i) {
      #pragma HLS UNROLL
      int pOffset = i*64;
      data.data.range(15+pOffset,0 +pOffset) = obj[i].hwPt;
      data.data.range(24+pOffset,16+pOffset) = obj[i].hwEta;
      data.data.range(34+pOffset,25+pOffset) = obj[i].hwPhi;
      data.data.range(38+pOffset,35+pOffset) = obj[i].hwId;
      data.data.range(48+pOffset,39+pOffset) = obj[i].hwZ0;
    }
}
template<unsigned int N, unsigned int OFFS> 
inline void mp7_pack(hls::stream<PFChargedObj> obj[], axi_t &data) {
    #pragma HLS inline
    for (unsigned int i = 0; i < N; ++i) {
      #pragma HLS UNROLL
      int pOffset = i*64;
      PFChargedObj tmpobj = obj[i].read();
      data.data.range(15+pOffset,0 +pOffset) = tmpobj.hwPt;
      data.data.range(24+pOffset,16+pOffset) = tmpobj.hwEta;
      data.data.range(34+pOffset,25+pOffset) = tmpobj.hwPhi;
      data.data.range(38+pOffset,35+pOffset) = tmpobj.hwId;
      data.data.range(48+pOffset,39+pOffset) = tmpobj.hwZ0;
    }
}
template<unsigned int N> void sumparts(pt_t &ipt,PFChargedObj iCol[4*N]) {
  #pragma HLS inline
  pt_t pt1 = 0; 
  pt_t pt2 = 0; 
  pt_t pt3 = 0; 
  pt_t pt4 = 0; 
  for(int i = 0; i < N; i++) { 
    pt1 = pt1 + iCol[i].hwPt;
  }
  for(int i = 0; i < N; i++) { 
    pt2 = pt2 + iCol[i+N].hwPt;
  }
  for(int i = 0; i < N; i++) { 
    pt3 = pt3 + iCol[i+N*2].hwPt;
  }
  for(int i = 0; i < N; i++) { 
    pt4 = pt4 + iCol[i+N*3].hwPt;
  }
  ipt = pt1 + pt2 + pt3 + pt4;
} 
void sumpt(pt_t &taupt, PFChargedObj pfch[NTRACK*4], PFChargedObj pfpho[NPHOTON*4], PFChargedObj pfne[NSELCALO*4], PFChargedObj pfmu[NMU*4]) {
  #pragma HLS inline
  #pragma HLS PIPELINE
  pt_t tauptch   = 0;
  pt_t tauptpho  = 0;
  pt_t tauptne   = 0;
  pt_t tauptmu   = 0;
  sumparts<NTRACK>  (tauptch ,pfch);
  sumparts<NPHOTON> (tauptpho,pfpho);
  sumparts<NSELCALO>(tauptne ,pfne);
  sumparts<NMU>     (tauptmu ,pfmu);
  taupt = tauptch + tauptpho + tauptne + tauptmu;
}
template<int N,int NMAX> 
void deltaR(int iOffSet, etaphi_t seedeta,etaphi_t seedphi,PFChargedObj pfch[N],PFChargedObj pfout[NMAX*4]) { 
  #pragma HLS inline
  #pragma HLS PIPELINE 
  const ap_int<16> eDR2MAX = DR2MAX;
  PFChargedObj dummyc; dummyc.hwPt = 0; dummyc.hwEta = 0; dummyc.hwPhi = 0; dummyc.hwId = 0; dummyc.hwZ0 = 0;
  for (int i = 0; i < NMAX; i++) {
   #pragma HLS UNROLL
    int drcheck = dr2_int_cap<12>(seedeta,seedphi,pfch[i].hwEta,pfch[i].hwPhi,eDR2MAX);
    if(drcheck < DRCONE) { 
     pfout[iOffSet+i] = pfch[i];
    } else { 
     pfout[iOffSet+i] = dummyc;
    }
  }
}

template<int N> 
void deltaR(int iOffSet, etaphi_t seedeta,etaphi_t seedphi,PFNeutralObj pfne[N],PFNeutralObj pfne_sel[N*4]) { 
  #pragma HLS inline 
  #pragma HLS PIPELINE 
  const ap_int<16> eDR2MAX = DR2MAX;
  PFNeutralObj dummyn; dummyn.hwPt = 0; dummyn.hwEta = 0; dummyn.hwPhi = 0; dummyn.hwId = 0; dummyn.hwPtPuppi = 0;
  for (int i = 0; i < N; i++) {
   #pragma HLS UNROLL
   int drcheck = dr2_int_cap<12>(seedeta,seedphi,pfne[i].hwEta,pfne[i].hwPhi,eDR2MAX);
   if(drcheck < DRCONE) { 
    pfne_sel[iOffSet+i] = pfne[i];
   } else { 
    pfne_sel[iOffSet+i] = dummyn;
   }
  }
}

void convert_out(PFChargedObj datas[NTAUPARTS], hls::stream<PFChargedObj> axis_out[DATA_SIZE]) {
	// input data from stream
	for (int i = 0; i < NTAUPARTS; ++i) {
	#pragma HLS UNROLL
	  datas[i] = axis_out[i].read();
	}
	/*  Below is needed to remove some errors when NTAUPARTS != DATA_SIZE
	for (int i = NTAUPARTS; i < DATA_SIZE; ++i) {
	  #pragma HLS unroll
	  axis_out[i].read();
	 }*/
}
void convert_in(PFChargedObj datas[DATA_SIZE], hls::stream<PFChargedObj > axis_in[DATA_SIZE]) {
	// output data to stream
	for (int i = 0; i < DATA_SIZE; ++i) {
	#pragma HLS unroll
		axis_in[i].write(datas[i]);
	}
}
template<unsigned int N0,unsigned int N1, unsigned int N2, unsigned int N3,unsigned int NOUT>
void MergeAllParts(PFChargedObj allparts[NOUT], PFChargedObj pfch[N0], PFNeutralObj pfpho[N1], PFNeutralObj pfne[N2], PFChargedObj pfmu[N3]) {
  #pragma HLS PIPELINE
  
  #pragma HLS ARRAY_PARTITION variable=allparts complete
  for(int i = 0; i < N0; i++) { 
    #pragma HLS UNROLL
    allparts[i] = pfch[i];
  }
  for(int i = 0; i < N1; i++) { 
    #pragma HLS UNROLL
    allparts[i+N0].hwPt = pfpho[i].hwPtPuppi;
    allparts[i+N0].hwEta = pfpho[i].hwEta;
    allparts[i+N0].hwPhi = pfpho[i].hwPhi;
    allparts[i+N0].hwId = pfpho[i].hwId;
    allparts[i+N0].hwZ0 = 0;
  }
  for(int i = 0; i < N2; i++) { 
    #pragma HLS UNROLL
    allparts[i+N0+N1].hwPt = pfne[i].hwPtPuppi;
    allparts[i+N0+N1].hwEta = pfne[i].hwEta;
    allparts[i+N0+N1].hwPhi = pfne[i].hwPhi;
    allparts[i+N0+N1].hwId = pfne[i].hwId;
    allparts[i+N0+N1].hwZ0 = 0;
  }
  for(int i = 0; i < N3; i++) { 
    #pragma HLS UNROLL
    allparts[i+N0+N1+N2] = pfmu[i];
  }
}
template<unsigned int N>
void make_inputs(input_t nn_data[N*8], hls::stream<PFChargedObj> pf[N]) {
  #pragma HLS inline
  #pragma HLS PIPELINE
  for (int i = 0; i < N; i++) {
        #pragma HLS PIPELINE II=1
    PFChargedObj tmpobj = pf[i].read();
    nn_data[i*8+0] = input_t(tmpobj.hwPt);
    nn_data[i*8+1] = input_t(tmpobj.hwEta);
    nn_data[i*8+2] = input_t(tmpobj.hwPhi);
    nn_data[i*8+3] = input_t(tmpobj.hwId == 2 ? 1 : 0);
    nn_data[i*8+4] = input_t(tmpobj.hwId == 3 ? 1 : 0);
    nn_data[i*8+5] = input_t(tmpobj.hwId == 4 ? 1 : 0);
    nn_data[i*8+6] = input_t(tmpobj.hwId == 1 ? 1 : 0);
    nn_data[i*8+7] = input_t(tmpobj.hwId == 0 ? 1 : 0);
  }
}
void algo_layer2(hls::stream<axi_t> &ch_link_in,hls::stream<axi_t> &ne_link_in, hls::stream<axi_t> &em_link_in, hls::stream<axi_t> &mu_link_in,
		 hls::stream<axi_t> &link_out) {
  #pragma HLS PIPELINE
  #pragma HLS INTERFACE axis port=ch_link_in
  #pragma HLS INTERFACE axis port=ne_link_in
  #pragma HLS INTERFACE axis port=em_link_in
  #pragma HLS INTERFACE axis port=mu_link_in
  #pragma HLS INTERFACE axis port=link_out
  #pragma HLS INTERFACE s_axilite port=return bundle=ctrl
  
  #pragma HLS stream variable=ch_link_in depth=36
  #pragma HLS stream variable=ne_link_in depth=36
  #pragma HLS stream variable=em_link_in depth=36
  #pragma HLS stream variable=mu_link_in depth=36
  #pragma HLS stream variable=link_out   depth=36

  PFChargedObj pfch[NREGIONS][NTRACK];
  PFChargedObj pfmu[NREGIONS][NMU];
  PFNeutralObj pfpho[NREGIONS][NPHOTON];
  PFNeutralObj pfne[NREGIONS][NCALO];
  #pragma HLS ARRAY_RESHAPE variable=pfch  block factor=25   dim=2
  #pragma HLS ARRAY_RESHAPE variable=pfpho block factor=20   dim=2
  #pragma HLS ARRAY_RESHAPE variable=pfne  block factor=20   dim=2
  #pragma HLS ARRAY_RESHAPE variable=pfmu  block factor=2    dim=2
  #pragma HLS RESOURCE      variable=pfch  core=RAM_2P_BRAM
  #pragma HLS RESOURCE      variable=pfpho core=RAM_2P_BRAM
  #pragma HLS RESOURCE      variable=pfne  core=RAM_2P_BRAM
  #pragma HLS RESOURCE      variable=pfmu  core=RAM_2P_BRAM
  
  //First take the seed from each region
  LoopA:
  for(int idepth =0; idepth < NREGIONS; idepth++) { 
    #pragma HLS PIPELINE II=1
    mp7_unpack<NTRACK> (ch_link_in.read(), pfch [idepth]);
    mp7_unpack<NEMCALO>(em_link_in.read(), pfpho[idepth]);
    mp7_unpack<NCALO>  (ne_link_in.read(), pfne [idepth]);
    mp7_unpack<NMU>    (mu_link_in.read(), pfmu [idepth]);
  }
  
//**************************************************************************//
//                         Start the Fancy jet algo                         //
//**************************************************************************//
  // Flatting out pt/eta/phi of ID for the next
  hls::stream<PFChargedObj > allparts_in [DATA_SIZE];
  #pragma HLS stream variable=allparts_in  depth=5
  const int NSUM = NTRACK+NCALO+NEMCALO+NMU;
  PFChargedObj allparts[NREGIONS][NSUM];
  FlatID:
  for(int idepth =0; idepth < NREGIONS; idepth++) { 
    #pragma HLS PIPELINE II=1
    MergeAllParts<NTRACK, NCALO, NEMCALO, NMU, NSUM> (allparts[idepth],pfch[idepth],pfne[idepth],pfne[idepth],pfmu[idepth]);  
  }


  // Reduce the size of the parts in the events
  // Region counting PUPPI Candidate
  PFChargedObj allpuppis[NREGIONS][NPUPPI_REG];
  count_t RegSize[NREGIONS];
  for(int idepth =0; idepth < NREGIONS; idepth++) { 
    #pragma HLS PIPELINE II=1
    // Get flag of PUPPI particle
    SelectPUPPIinRegion<NSUM, NPUPPI_REG>(allparts[idepth], allpuppis[idepth], RegSize[idepth]);
  }



  //ap_uint<8> IsPUPPIcnt[totParts] ;
  //ap_uint<1> IsPUPPIX[totParts] ;
  //PFChargedObj allpuppi[NPUPPI];

  //Next iterate through seeds and compute DR




  //
  //PFChargedObj seeds[NTAU];
  //#pragma HLS ARRAY_PARTITION variable=seeds complete
  //ptsort_hwopt_ind<PFChargedObj,NREGIONS,NTAU>(chseed, seeds);
  //const ap_int<16> eDR2MAX = DR2MAX;
  //ap_uint<6> seedsort2[NTAU];
  //#pragma HLS ARRAY_PARTITION variable=seedsort2 complete
  //LoopB:
  //for(int iseed0=0; iseed0 < NTAU; iseed0++) {
   //#pragma HLS UNROLL
    //seedsort2[iseed0]=iseed0;
  //}
  //LoopC:
  //for(int iseed0=0; iseed0 < NTAU; iseed0++) {
    //#pragma HLS UNROLL
    //for(int iseed1 = 1; iseed1 < NTAU; iseed1++) {
      //if(iseed1 < iseed0+1) continue;                                                                                                                                                                                                                                         
      //int drcheck = dr2_int_cap<16>(seeds[iseed0].hwEta,seeds[iseed0].hwPhi,seeds[iseed1].hwEta,seeds[iseed0].hwPhi,eDR2MAX);
      //if(drcheck < DRCONE && iseed1 > iseed0) {
        //seedsort2[iseed1] = 0;
      //}
    //}
  //}
  ////Now get all the candidates in each region => Note that a max of 4 regoins can fill array                                                                                                                                                                                  
  //PFChargedObj pfchreg [NTAU][NTAUPARTS*4];
  //PFChargedObj pfemreg [NTAU][NTAUPARTS*4];
  //PFChargedObj pfnereg [NTAU][NTAUPARTS*4];
  //PFChargedObj pfmureg [NTAU][NMU*4];
  //#pragma HLS ARRAY_RESHAPE variable=pfchreg  dim=0 complete
  //#pragma HLS ARRAY_RESHAPE variable=pfemreg  dim=0 complete
  //#pragma HLS ARRAY_RESHAPE variable=pfnereg  dim=0 complete
  //#pragma HLS ARRAY_RESHAPE variable=pfmureg  dim=0 complete
  //PFChargedObj dummyc; dummyc.hwPt = 0; dummyc.hwEta = 0; dummyc.hwPhi = 0; dummyc.hwId = 0; dummyc.hwZ0 = 0;
  //LoopD:
  //for(int itau = 0; itau < NTAU; itau++) {
    //#pragma HLS UNROLL
    //etaphi_t seedeta = seeds[itau].hwEta;
    //etaphi_t seedphi = seeds[itau].hwPhi;
    //ap_int<8> arr[4];
    //#pragma HLS ARRAY_PARTITION variable=arr  complete dim=0
    //arraymap(seedeta,seedphi,arr);
    //for(int iregion = 0; iregion < 4; iregion++) {
      //#pragma HLS UNROLL
      //int pRegion = arr[iregion];
      //int NTRACKOFFSET   = iregion*NTAUPARTS;//NTRACK;
      //int NPHOTONOFFSET  = iregion*NTAUPARTS;//NPHOTON;
      //int NSELCALOOFFSET = iregion*NTAUPARTS;//NSELCALO;
      //int NMUOFFSET      = iregion*NMU;
      //deltaR<NTRACK,NTAUPARTS>   (NTRACKOFFSET,  seedeta,seedphi,pfch [pRegion],pfchreg[itau]);
      //deltaR<NPHOTON,NTAUPARTS>  (NPHOTONOFFSET, seedeta,seedphi,pfpho[pRegion],pfemreg[itau]);
      //deltaR<NSELCALO,NTAUPARTS> (NSELCALOOFFSET,seedeta,seedphi,pfne [pRegion],pfnereg[itau]);
      //deltaR<NMU,NMU>            (NMUOFFSET,     seedeta,seedphi,pfmu [pRegion],pfmureg[itau]);
    //}
  //}
  
  //hls::stream<PFChargedObj > allparts_in [DATA_SIZE];
  //hls::stream<PFChargedObj > allparts_out[DATA_SIZE];
  //#pragma HLS stream variable=allparts_in  depth=5
  //#pragma HLS stream variable=allparts_out depth=5
  //for(int itau = 0; itau < NTAU; itau++) {
    //PFChargedObj allparts[DATA_SIZE];
    //#pragma HLS ARRAY_PARTITION variable=allparts complete
    //tausort<4*NTAUPARTS,4*NTAUPARTS,4*NTAUPARTS,4*NMU,DATA_SIZE>(allparts,pfchreg[itau],pfnereg[itau],pfemreg[itau],pfmureg[itau]);  
    //convert_in(allparts,allparts_in);
  //}
  //sorting_network_128(allparts_in,allparts_out);
  //
  
  //PFChargedObj taus[NTAU];
  //#pragma HLS ARRAY_PARTITION variable=taus complete
  //for(int itau = 0; itau < NTAU; itau++) {
    //input_t nn_data[NTAUPARTS*8];
    ////#pragma HLS ARRAY_PARTITION variable=nn_data complete
    //make_inputs<NTAUPARTS>(nn_data,allparts_out);            
    //result_t taunn[N_OUTPUTS];
    //tau_nn(nn_data,taunn);
    //PFChargedObj dummyc; 
    //dummyc.hwPt = taunn[0]*100;; dummyc.hwEta = nn_data[1]; dummyc.hwPhi = nn_data[2]; dummyc.hwId = 0; dummyc.hwZ0 = 0;
    //taus[itau] = dummyc;
  //}
  //ptsort_hwopt_ind<PFChargedObj,NTAU,NTAU>(taus, tausout);
  PFChargedObj jetsout[NJET]; 
  #pragma HLS ARRAY_PARTITION variable=jetsout complete
  for (int i = 0; i < NJET; ++i)
  {
    #pragma HLS PIPELINE
    jetsout[i].hwPt = allpuppis[0][i].hwPt;
    jetsout[i].hwEta = allpuppis[0][i].hwEta;
    jetsout[i].hwPhi = allpuppis[0][i].hwPhi;
    jetsout[i].hwId = allpuppis[0][i].hwId;
    jetsout[i].hwZ0 = 0;
  }
  axi_t data_out;
  mp7_pack<NJET,0>(jetsout,data_out);
  link_out.write(data_out);
}
