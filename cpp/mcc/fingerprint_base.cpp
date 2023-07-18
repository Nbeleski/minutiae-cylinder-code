//
// fingerprint_base.cpp
//

#include "fingerprint_base.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>


using namespace std;


short compute_angle_diff_180(const short angle1, const short angle2)
{
	auto diff = abs(angle1 - angle2);
	return (short)(diff < 180 ? diff : (360 - diff));
}
int compute_squared_euclidean_distance(const Point& p1, const Point& p2)
{
	int dx = p1.x - p2.x;
	int dy = p1.y - p2.y;
	return (int)sqrt((double)dx * (double)dx + (double)dy * (double)dy);
}
float compute_angle_radians(const float dx, const float dy)
{
	float angle = (float)atan2(dy, dx);
	if (angle < 0)
		angle += (float)(TWO_PI);
	return angle;
}
short compute_angle_degrees(const short dx, const short dy)
{
	return (short)(compute_angle_radians(dx, dy) * RAD_TO_DEG + 0.5f);
}
short compute_angle_degrees(const Point& p1, const Point& p2)
{
	return compute_angle_degrees(p2.x - p1.x, p1.y - p2.y);
}
short compute_angle_diff_360(const short angle1, const short angle2)
{
	short diff = angle1 - angle2;
	return  diff < 0 ? diff + 360 : diff;
}
vector<int> compute_adjacency_matrix(const vector<Minutia>& minutiae)
{
	int len = (int)minutiae.size();
	vector<int> adjacency(len * len, 0);

	for (int i = 0; i < len; i++)
	{
		for (int j = (i + 1); j < len; j++)
		{
			auto dist = (int)compute_squared_euclidean_distance(minutiae[i].position, minutiae[j].position);
			adjacency[i * len + j] = adjacency[j * len + i] = dist;
		}
	}
	return adjacency;
}

Point::Point(short _x, short _y) : x(_x), y(_y)
{

}

bool Point::operator==(const Point& p) const
{
	return (x == p.x) && (y == p.y);
}
Point Point::operator-(const Point& p) const
{
	Point ret;
	ret.x = x - p.x;
	ret.y = y - p.y;
	return ret;
}

float Point::norm() const
{
	return sqrtf(static_cast<float>(x * x + y * y));
}

bool Minutia::operator==(const Minutia& m) const
{
	return (position == m.position) && (direction == m.direction) && (index == m.index);
}

float Minutia::arc() const
{
	return 2.0f * direction * PI / 180.0f;
}


int ReadRecordReader(const unsigned char* buffer, int* currentCursorPosition)
{
	if (buffer == nullptr)
		return -1;


	int pos;

	//Format Identifier
	if (buffer[0] != 'F' || buffer[1] != 'M' || buffer[2] != 'R' || buffer[3] != 0x00)
		return -1;

	//Version of the standard
	if (buffer[4] != ' ' || buffer[5] != '2' || buffer[6] != '0' || buffer[7] != 0x00)
		return -2;

	//Length of total record in bytes
	long lengthOfRecord = 0;

	if (buffer[8] == 0x00 && buffer[9] == 0x00)
	{
		lengthOfRecord |= (buffer[10] << 24);
		lengthOfRecord |= (buffer[11] << 16);
		lengthOfRecord |= (buffer[12] << 8);
		lengthOfRecord |= buffer[13];
		pos = 14;
	}
	else
	{
		lengthOfRecord |= (buffer[8] << 8);
		lengthOfRecord |= buffer[9];
		pos = 10;
	}

	if (lengthOfRecord == 32)
		return -3; // TemplateAnsiIncitsErrorCodes::NULL_TEMPLATE;

	//Capture System: ignored to speed up
	pos += 6;

	//Image Size
	short imageWidth = buffer[pos++] << 8;
	imageWidth |= buffer[pos++];

	short imageHeight = buffer[pos++] << 8;
	imageHeight |= buffer[pos++];

	//Image Resolution
	short horizontalResolution = buffer[pos++] << 8;
	horizontalResolution |= buffer[pos++];

	horizontalResolution = static_cast<short>((horizontalResolution * 2.54f) + 0.5f);

	short verticalResolution = buffer[pos++] << 8;
	verticalResolution |= buffer[pos++];

	verticalResolution = static_cast<short>((verticalResolution * 2.54f) + 0.5f);

	//Number of Finger Views: here we're just considering one finger view
	pos++;

	//Reserved Byte
	pos++;

	*currentCursorPosition = pos;

	return 0;  //TemplateAnsiIncitsErrorCodes::SUCCESS;
}

int ReadFingerView(int actualFingerView, const unsigned char* buffer, int* currentCursorPosition, vector<Minutia>& minutiae)
{
	auto pos = *currentCursorPosition;

	//Finger Position
	unsigned char fgp = (unsigned char)(buffer[pos++]);

	//View Number: here we ignore because we're just considering the first view

	//Impression type
	unsigned char imp = (unsigned char)(buffer[pos++] & 0x0F);

	//Quality
	unsigned char fingerQuality = buffer[pos++];

	//NumMinutiae
	int numMinutiae = buffer[pos++];
	minutiae.resize(numMinutiae);

	for (int i = 0; i < numMinutiae; i++, pos += 6)
	{
		unsigned char type = (buffer[pos] >> 6) & 0x03;

		short posX = (buffer[pos] & 0x3F) << 8;
		posX |= buffer[pos + 1];

		short posY = (buffer[pos + 2] & 0x3F) << 8;
		posY |= buffer[pos + 3];

		short orientation = (short)(buffer[pos + 4]);

		unsigned char quality = buffer[pos + 5];

		//minutiae.push_back(MinutiaAnsiIncits(posX, posY, orientation, i, quality, type));

		minutiae[i].position.x = posX;
		minutiae[i].position.y = posY;
		minutiae[i].direction = orientation;
		minutiae[i].index = i;
		//minutiae[i].quality = quality;
		//minutiae[i].type = type;
	}

	//We do not check the extended are content
	*currentCursorPosition = pos;

	return 0;
}


vector<Minutia> decode_from_buffer(const unsigned char* buffer)
{
	vector<Minutia> minutiae;

	auto currentCursorPosition = 0;
	int status = ReadRecordReader(buffer, &currentCursorPosition);
	if (status == 0)
	{
		status = ReadFingerView(0, buffer, &currentCursorPosition, minutiae);
	}
	return minutiae;
}

vector<unsigned char> read_file(const char* filename)
{
	vector<unsigned char> buffer;
	std::ifstream ifs;
	ifs.open(filename, std::ifstream::in | std::ifstream::binary);

	if (ifs.is_open())
	{
		ifs.seekg(0, std::ios::end);
		const auto file_size = ifs.tellg();
		ifs.seekg(0, std::ios::beg);

		buffer = vector<unsigned char>(file_size);
		ifs.read((char*)buffer.data(), file_size);
	}
	return buffer;
}
