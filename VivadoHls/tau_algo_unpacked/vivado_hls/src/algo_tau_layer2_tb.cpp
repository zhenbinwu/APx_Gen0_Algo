#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

#include "ap_int.h"
#include "algo_tau_layer2.h"

using namespace std;

typedef ap_int<64>            tmpaxi_t;
typedef ap_axis <64*NPART,1,1,1> axi_t;
typedef hls::stream<axi_t>    stream_t;

int main(int argc, char ** argv) {
  hls::stream<PFChargedObj> link_in[DATA_SIZE];
  hls::stream<axi_t>        link_out;

  float dphi  = 0.02; 
  float deta  = 0.02; 
  float pt    = 1.;
  float pterr = 0.1;
  PFChargedObj parts[DEPTH][DATA_SIZE];  
  for(int idepth = 0; idepth < DEPTH; idepth++) { 
    int   pEta  = idepth % 4;
    int   pPhi  = idepth/3;
    for(int i0 = 0; i0 < DATA_SIZE; i0++) { 
      int shift   = (i0/20);
      float pTEta = 0. + pEta*0.7+dphi*(shift+int(i0/2));
      float pTPhi = 0. + pPhi*0.7+deta*(2*shift+1+int(i0/2));
      PFChargedObj pTmp;
      pTmp.hwPt  = pt*PT_SCALE;
      pTmp.hwEta = pTEta*ETAPHI_SCALE;
      pTmp.hwPhi = pTPhi*ETAPHI_SCALE;
      pTmp.hwId  = shift;
      pTmp.hwZ0  = 0;
      parts[idepth][i0] = pTmp;
    }
  } 
  for(int idepth = 0; idepth < DEPTH; idepth++) { 
    for(int i0 = 0; i0 < DATA_SIZE;  i0++)      link_in[i0].write(parts[idepth][i0]);
  }
  algo_tau_layer2(link_in, link_out);
  tmpaxi_t tau[DEPTH][NPART];  
  for(int idepth = 0; idepth < 1; idepth++) {
    axi_t tmp_out; 
    link_out.read(tmp_out);
    std::cout << " out " << tmp_out.data << std::endl;
    for(int i0 = 0; i0 < NPART; i0++) tau[idepth][i0] = tmp_out.data.range(63*(i0+1),64*(i0));
    for(int i0 = 0; i0 < NPART; i0++) { 
     tmpaxi_t pTmp = tau[idepth][i0];
     float pPt  = pTmp.range(15,0);pPt/=PT_SCALE;
     float pEta = pTmp.range(24,16);pEta/=ETAPHI_SCALE;
     float pPhi = pTmp.range(34,25);pPhi/=ETAPHI_SCALE;
     std::cout << "===> depth " << idepth << " -- part " << i0 << " vector " << pPt << "-- " << pEta << " -- " << pPhi << std::endl;
    }
  }
}

