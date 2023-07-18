//
// mcc.cpp
//

#include "mcc.h"
#include "utils.hpp"
#include "Munkres.h"

#include <array>
#include <numeric>
#include <algorithm>

using namespace std;

BitCylinderSet::BitCylinderSet(std::vector<Minutia>& minutiae)
{
	cylinders.reserve(minutiae.size());
	validities.reserve(minutiae.size());
	minutiae_references.reserve(minutiae.size());

	// convex hull of all minutiae
	auto chull = convex_hull(minutiae);

	// pre-calculate phi values as they repeat for all cells
	auto phi_vals = std::vector<float>();
	for (auto k = 1; k <= ND; k++)
		phi_vals.push_back(cell_angle(k));

	for (auto idx_probe = 0; idx_probe < minutiae.size(); idx_probe++)
	{
		// reference to ease readability
		auto m = &minutiae[idx_probe];

		// temp float_cylinder to hold the results
		std::array<float, CYLINDER_SIZE> tmp_cylinder = { 0.0f };

		// cylinders
		LinearizedBitCylinder cyl_values;
		LinearizedBitCylinder cyl_validity;

		// set all validities to 1
		cyl_validity = ~cyl_validity;

		// pre-calculate sin cos
		auto sin = sinf(m->arc());
		auto cos = cosf(m->arc());

		// invalidating checks
		auto contributing_minutiae = std::vector<int>(minutiae.size());
		auto invalid_cells = 0;

		// for each cell in a 'floor'
		for (size_t i = 0; i < NS; i++)
		{
			for (size_t j = 0; j < NS; j++)
			{
				// cell center
				auto center = cell_center(*m, i + 1, j + 1, sin, cos);

				// if cell center is further than radius or outside convexhull, invalidate cell
				if (ed_to_point(*m, center) > CYLINDER_RADIUS || !is_point_inside_chull(chull, center, OMEGA))
				{
					// invalidate on all 'floors'
					for (size_t k = 0; k < ND; k++)
						cyl_validity[linearize_idxs(i, j, k)] = false;

					invalid_cells++;
					continue;
				}

				// given a valid cell, we check the contribution
				// of each other minutiae on the template
				// filtering first those that are inside the radius 3*SIGMA_S
				for (auto idx_cand = 0; idx_cand < minutiae.size(); idx_cand++)
				{
					// skip itself
					if (idx_probe == idx_cand)
						continue;

					// reference to ease readability
					auto mt = &minutiae[idx_cand];

					auto dist = ed_to_point(*mt, center);
					if (dist > 3.0f * SIGMA_S)
						continue;

					auto spatial_contribution = gaussian(dist);
					contributing_minutiae[idx_cand] = 1;
					
					auto minutiae_ang_diff = ang_diff(m->arc(), mt->arc());
					for (size_t k = 0; k < ND; k++)
					{
						auto d_phi = phi_vals[k];
						auto ang_difference = ang_diff(d_phi, minutiae_ang_diff);
						auto directional_contribution = area_under_gaussian(ang_difference);

						tmp_cylinder[linearize_idxs(i, j, k)] += spatial_contribution * directional_contribution;
					} // k
				} // mt
			} // j
		} // i

		// check if the cylinder should be discarded
		auto sum_of_contributing = std::accumulate(contributing_minutiae.begin(), contributing_minutiae.end(), decltype(contributing_minutiae)::value_type(0));
		if ((1.0 - invalid_cells / (NS * NS) < MIN_VALID_CELLS) || (sum_of_contributing < MIN_CONTRIBUTING_MINUTIAE))
			continue;

		// linearize the results
		for (size_t i = 0; i < CYLINDER_SIZE; i++)
		{
			if (tmp_cylinder[i] > MICRO_PSI)
				cyl_values[i] = true;
		}

		cylinders.push_back(cyl_values);
		validities.push_back(cyl_validity);
		minutiae_references.push_back(idx_probe);
	} // m
}

size_t BitCylinderSet::len()
{
	return cylinders.size();
}

inline float bitwise_sum(LinearizedBitCylinder& b)
{
	return b.count();
}

inline float bitwise_norm(LinearizedBitCylinder& b)
{
	return sqrtf(b.count());
}

/// calculate np - np best similarities
size_t calc_nP(BitCylinderSet& c1, BitCylinderSet& c2)
{
	return size_t(NP_MIN + (NP_MAX - NP_MIN) * sigmoid(std::min(c1.len(), c2.len()), 20.0f, 2.0f / 5.0f) + 0.5f);
}

std::vector<float> local_matching(BitCylinderSet& c1, std::vector<Minutia>& m1, BitCylinderSet& c2, std::vector<Minutia>& m2)
{
	auto similarities = std::vector<float>();
	similarities.reserve(c1.len() * c2.len());

	// allocate memory once, afterall they all have fixed size
	LinearizedBitCylinder validitiy_ab, c_ab, c_ba, cxor;

	for (size_t i = 0; i < c1.len(); ++i) //ca = c1[i]
	{
		for (size_t j = 0; j < c2.len(); ++j) //cb = c2[j]
		{
			// CONDITION 1: angle difference
			if (ang_diff(m1[c1.minutiae_references[i]].arc(), m2[c2.minutiae_references[j]].arc()) > MAX_ANG_DIFF)
			{
				similarities.push_back(0.0f);
				continue;
			}

			validitiy_ab = c1.validities[i] & c2.validities[j];
			c_ab = c1.cylinders[i] & validitiy_ab;
			c_ba = c2.cylinders[j] & validitiy_ab;

			// CONDITION 2: number of matchable cells in cylinder
			if (bitwise_sum(validitiy_ab) / CYLINDER_SIZE < 0.6f)
			{
				similarities.push_back(0.0f);
				continue;
			}

			auto norm_c_ab = bitwise_norm(c_ab);
			auto norm_c_ba = bitwise_norm(c_ba);

			// CONDITION 3: norm sum != 0
			if (norm_c_ab + norm_c_ba == 0)
			{
				similarities.push_back(0.0f);
				continue;
			}

			cxor = c_ab ^ c_ba;
			auto lambda = 1.0f - (bitwise_norm(cxor) / (norm_c_ab + norm_c_ba));
			similarities.push_back(lambda);
		}
	}

	return similarities;
}

/// global score
/// (23) on paper
float global_matching(std::vector<float>& similarities, size_t nP)
{
	std::partial_sort(similarities.begin(), similarities.begin() + nP, similarities.end(), std::greater{});
	auto sum = std::accumulate(similarities.begin(), similarities.begin() + nP, 0.0f);

	return sum / static_cast<float>(nP);
}

/// global score using LSA
/// (23) on paper
float global_matching_LSA(std::vector<float>& similarities, BitCylinderSet& c1, BitCylinderSet& c2)
{
	auto nP = calc_nP(c1, c2);

	auto matrix = Matrix(similarities.data(), c1.len(), c2.len());
	auto macopy = Matrix(matrix);

	// inverte valores da matrix pois essa implementação de Munkres procura o menor path
	for (auto i = 0; i < c1.len(); ++i) {
		for (auto j = 0; j < c2.len(); ++j) {
			if (matrix(i, j) <= 0.0) {
				matrix(i, j) = 999;
			}
			else {
				matrix(i, j) = 1.0f / matrix(i, j);
			}
		}
	}

	Munkres m;
	m.solve(matrix); //TODO: alterar para poder calcular invertido sem copia

	std::vector<float> valid_scores;
	for (auto i = 0; i < c1.len(); ++i)
		for (auto j = 0; j < c2.len(); ++j)
			if (matrix(i, j) == 0)
				valid_scores.push_back(macopy(i, j));

	std::partial_sort(valid_scores.begin(), valid_scores.begin() + nP, valid_scores.end(), std::greater{});
	auto sum = std::accumulate(valid_scores.begin(), valid_scores.begin() + nP, 0.0f);

	return sum / static_cast<float>(nP);
}


