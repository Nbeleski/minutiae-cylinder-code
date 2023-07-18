//
// fingerprint_base.h
//
#pragma once
#include <vector>

constexpr float PI = 3.14159265358f;
constexpr float TWO_PI = 2.0f * PI;
constexpr float RAD_TO_DEG = 180.0f / PI;
constexpr float DEG_TO_RAD = PI / 180.0f;

struct Point // 4 bytes
{
	short x = 0;
	short y = 0;

	Point(short _x, short _y);
	Point() = default;
	bool operator==(const Point& p) const;
	Point operator-(const Point& p) const;
	float norm() const;
};
struct Minutia
{
	Point position{ 0,0 };
	short direction = 0;
	unsigned char index = 0;

	//minutia();
	//minutia(short x, short y, short ang);
	bool operator==(const Minutia& m) const;
	float arc() const;
};

//struct pre_match
//{
//	point nucleo1{0,0};
//	point nucleo2{ 0,0 };
//	point delta1{ 0,0 };
//	point delta2{ 0,0 };
//	float grupos[6]{0,0,0,0,0,0};
//	unsigned char qual{0};
//};

std::vector<Minutia> decode_from_buffer(const unsigned char* buffer);
std::vector<unsigned char> read_file(const char* filename);


short compute_angle_diff_180(const short angle1, const short angle2);
int compute_squared_euclidean_distance(const Point& p1, const Point& p2);
float compute_angle_radians(const float dx, const float dy);
short compute_angle_degrees(const short dx, const short dy);
short compute_angle_degrees(const Point& p1, const Point& p2);
short compute_angle_diff_360(const short angle1, const short angle2);
std::vector<int> compute_adjacency_matrix(const std::vector<Minutia>& minutiae);