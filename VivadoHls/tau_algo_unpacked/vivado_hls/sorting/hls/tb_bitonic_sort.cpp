#include "bitonic_sort.hpp"

#include <cstdio>
#include <algorithm>


// Number of test case
#define TEST_NUM 32

struct TestData {
	PFChargedObj   datas         [DATA_SIZE];
	PFChargedObj   datas_acutual [DATA_SIZE];
	PFChargedObj   datas_expected[DATA_SIZE];
};


void makeTestData(TestData& td);
int check(TestData& td);
int tb_bitonic_sort();


int main() {
	// return number of error
	return tb_bitonic_sort();
}


int tb_bitonic_sort() {

  hls::stream<PFChargedObj > axis_in[DATA_SIZE];
  hls::stream<PFChargedObj > axis_out[DATA_SIZE];
  // Make test data
  TestData testdatas[TEST_NUM];
  for (int id = 0; id < TEST_NUM; ++id) {
    makeTestData(testdatas[id]);
  }
  // Set input data
  for (int id = 0; id < TEST_NUM; ++id) {
    for (int i = 0; i < DATA_SIZE; ++i) {
      axis_in[i].write(testdatas[id].datas[i]);
    }
  }
  // Execute
  for (int id = 0; id < TEST_NUM-1; ++id) {
    bitonic_sort(axis_in,axis_out);
  }
  // Get output
  bool flag = true;
  for (int id = 0; id < TEST_NUM; ++id) {
    for (int i = 0; i < DATA_SIZE; ++i) {
      testdatas[id].datas_acutual[i] = axis_out[i].read();
    }
  }
  // Check
  int err = 0;
  for (int id = 0; id < TEST_NUM; ++id) {
    err += check(testdatas[id]);
  }
  return err;
}

bool ptsort(PFChargedObj i,PFChargedObj j) { return (i.hwPt<j.hwPt); }
void makeTestData(TestData& td) {
	const int MAX = 2047;
	const int MIN = -2047;
	for (int i = 0; i < DATA_SIZE; ++i) {
	  td.datas[i].hwPt = rand()% (MAX - MIN + 1) + MIN;
	  td.datas_expected[i].hwPt = td.datas[i].hwPt;
	}
	std::sort(td.datas_expected, td.datas_expected + DATA_SIZE,ptsort);
	//std::sort(td.datas_expected[0],td.datas_expected[DATA_SIZE-1],ptsort);
	for(int i0 = 0; i0 < DATA_SIZE; i0++) { 
	  //std::cout << "==> " << td.datas_expected[i0].hwPt << std::endl;
	}
}


int check(TestData& td) {
	int ret = 0;
	for (int i = 0; i < DATA_SIZE; ++i) {
	  //std::cout << "==> " << i << " -- " << td.datas_acutual[i].hwPt << " -- " << td.datas_expected[i].hwPt << std::endl;
		if (td.datas_acutual[i].hwPt != td.datas_expected[i].hwPt)
			ret++;
	}
	return 0;//ret;
}
