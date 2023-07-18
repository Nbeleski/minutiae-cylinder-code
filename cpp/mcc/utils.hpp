#pragma once
//
// utils.hpp
//

#include "fingerprint_base.h"
#include "constants.h"
#include "convexhull.hpp"

#include <math.h>

/// Calculate angle associated with a cell
/// (1) on paper
inline float cell_angle(int k)
{
	return -PI + (k - 0.5f) * DELTA_D;
}

/// Calculate cell-center based on minutiae and pre-calculated sin & cos
/// (2) on paper
inline Point cell_center(Minutia &m, size_t i, size_t j, float sin, float cos)
{
	auto i_term = i - (NS + 1.0f) / 2.0f;
	auto j_term = j - (NS + 1.0f) / 2.0f;

	auto x = 0.5f + m.position.x + DELTA_S * (cos * i_term + sin * j_term);
	auto y = 0.5f + m.position.y + DELTA_S * (-sin * i_term + cos * j_term);

	return Point(short(x), short(y));
}

/// Euclidean distance from minutia to point
inline float ed_to_point(Minutia& m, Point& p)
{
	auto dx = float(m.position.x - p.x);
	auto dy = float(m.position.y - p.y);

	return sqrtf(dx * dx + dy * dy);
}

/// Parametrized sigmoid
/// (5) on paper
inline float sigmoid(float v, float u, float t)
{
	return 1.0f / (1.0f + expf(-t * (v - u)));
}

/// Gaussian used in distance contribution calculation
/// (7) on paper
inline float gaussian(float t)
{
	return (1.0f / (SIGMA_S * SQRT_TWO_PI))* expf(-(t * t) / (2.0f * SIGMA_S * SIGMA_S));
}

/// Angle difference
/// (9) on paper

inline float ang_diff(float t1, float t2)
{
	auto diff = t1 - t2;
	if (diff >= -PI && diff < PI)
		return diff;
	else if (diff < -PI)
		return diff + TWO_PI;
	else
		return diff - TWO_PI;
}

/// Calculate de convex hull of a set of minutiae
std::vector<Point> convex_hull(std::vector<Minutia>& minutiae)
{
	std::vector<Point> coords;
	coords.reserve(minutiae.size());

	for (auto& m : minutiae)
		coords.push_back(m.position);

	return convex_hull_jarvis(coords);
}

bool point_inside_poly(std::vector<Point>& coords, Point& p)
{
	int i, j, c = 0;
	for (i = 0, j = coords.size() - 1; i < coords.size(); j = i++) {
		if (((coords[i].y > p.y) != (coords[j].y > p.y)) &&
			(p.x < (coords[j].x - coords[i].x) * (p.y - coords[i].y) / (coords[j].y - coords[i].y) + coords[i].x))
			c = ~c;
	}
	return c;
}

// TODO: find a faster impl combinint above
/// checks if a point is inside of a convex hull
/// with an offset omega
inline bool is_point_inside_chull(std::vector<Point>& chull, Point& p, float omega = OMEGA)
{
	if (point_inside_poly(chull, p))
		return true;

	// if the point is outside, we need the distance
	// to check if it is inside the offset omega

	auto min_dist = 99999.0f;

	for (auto i = 0; i < chull.size() - 1; ++i)
	{
		auto p1 = &chull[i];
		auto p2 = &chull[i + 1];

		auto vec_p1p0 = Point(p.x - p1->x, p.y - p1->y);
		auto vec_p1p2 = Point(p2->x - p1->x, p2->y - p1->y);
		auto vec_p2p0 = Point(p.x - p2->x, p.y - p2->y);

		auto norm = powf(vec_p1p2.norm(), 2);
		auto dot = vec_p1p0.x * vec_p1p2.x + vec_p1p0.y * vec_p1p2.y;
		auto r = dot / norm;

		auto dist = 0.0f;
		if (r > 1.0f)
			dist = vec_p2p0.norm();
		else if (r <= 0.0f)
			dist = vec_p1p0.norm();
		else { // between 0 and 1
			dist = sqrtf(powf(vec_p1p0.norm(), 2) - powf((r * vec_p1p2.norm()), 2));
		}

		min_dist = std::min(dist, min_dist);
	}

	if (min_dist < omega)
		return true;
	return false;
}

/// Area under curve, calculated using wolfram
/// https://www.wolframalpha.com/input?i2d=true&i=Divide%5B1%2Csigma+*+Sqrt%5B2+*+pi%5D%5D+*+Integrate%5BPower%5Be%2CDivide%5B-Power%5Bt%2C2%5D%2C2*Power%5Bsigma%2C2%5D%5D%5D%2C%7Bt%2Cx+-+Divide%5BD%2C2%5D%2Cx+%2B+Divide%5BD%2C2%5D%7D%5D
/// (11) on paper
inline float area_under_gaussian(float v)
{
	return 0.5f * (erff((0.5f * DELTA_D + v) / (sqrtf(2.0f) * SIGMA_D)) - erff((v - 0.5f * DELTA_D) / (sqrtf(2.0f) * SIGMA_D)));
}

/// Formula to linearize the idxs
/// slightly different from the paper as my arrays start at 0
/// (13) on paper
inline size_t linearize_idxs(size_t i, size_t j, size_t k)
{
	return static_cast<size_t>(k * NS * NS + j * NS + i);
}