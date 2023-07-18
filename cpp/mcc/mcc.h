#pragma once
//
// mcc.h
//

#include "fingerprint_base.h"
#include "constants.h"

struct BitCylinderSet
{
	std::vector<LinearizedBitCylinder> cylinders;
	std::vector<LinearizedBitCylinder> validities;
	std::vector<unsigned char> minutiae_references;
	
	BitCylinderSet(std::vector<Minutia>& minutiae);
	size_t len();
};

std::vector<float> local_matching(BitCylinderSet& c1, std::vector<Minutia>& m1, BitCylinderSet& c2, std::vector<Minutia>& m2);
size_t calc_nP(BitCylinderSet& c1, BitCylinderSet& c2);
float global_matching(std::vector<float>& similarities, size_t nP);
float global_matching_LSA(std::vector<float>& similarities, BitCylinderSet& c1, BitCylinderSet& c2);

