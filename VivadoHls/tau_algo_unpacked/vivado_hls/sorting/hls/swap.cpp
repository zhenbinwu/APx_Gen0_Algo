#include "swap.hpp"


void swap1(PFChargedObj &data1, PFChargedObj &data2) {
	if (data1.hwPt < data2.hwPt) {
		std::swap(data1, data2);
	}
}

void swap2(PFChargedObj &data1, PFChargedObj &data2) {
	if (data1.hwPt > data2.hwPt) {
		std::swap(data1, data2);
	}
}
