//
// main.cpp
//

#include <iostream>
#include <optional>
#include <assert.h>
#include <chrono>
#include <execution>

#include "..\mcc\fingerprint_base.h"
#include "..\mcc\mcc.h"

#include "base64.h"

using namespace std;


using us = std::chrono::duration<double, std::micro>;

// ┌──────────────────────────┐
// │   MCC16b usage example   │
// └──────────────────────────┘

std::optional<vector<Minutia>> incits378_from_buffer(unsigned char* buffer)
{
	std::vector<::Minutia> minucias = decode_from_buffer(buffer);
	if (minucias.size() < 2)
		return {};

	return minucias;
}

int main() try
{
	cout << "Teste mcc" << endl;

	auto b1 = base64_decode("Rk1SACAyMAAC8AAyAAAAAAIAAgAAxQDFAQAAAAB4QOYAt2sAQNwA4ZMAQOwA4TUAgPUA4zAAgO4AnFIAgQ0A1oUAQRAAs5cAQLwAxIMAgMEA7YwAQLMA5ocAQRAA+3oAQR0A/oAAQKoAsH0AgSYAqJsAgSoAwZMAQKEA2YIAQM4AeF8AgS4BA4cAgLcBCDIAQSYAcKMAQOIAb60AQToAf0gAQOsBF0sAQRwBGXgAgT8A8YwAQJgBGIkAgNgBHKAAQUEAY58AgOgAX68AgL0BIJAAQUwBGS4AgIcBJIcAgSkBJnkAgVAA1ZIAgVIAgkEAgVIArJUAgHwBIYMAgVIBGIsAgVUAyjoAgMABLJMAgToBMI0AgTMBLnoAgNQBL6MAgHYAongAgHQAu3oAQP0BNRUAgPABN10AQM4BOD0AQGkBKoIAgUYBQD8AgWkApJcAgTIAP0wAgToBQiMAQNoBQlEAgSMBRBkAgW4AR0AAgVEBRUIAQLcAOGMAgF4BM34AgW4A15MAgXMASjYAgW4BQEYAgFwBCXwAgXQAOJkAgJIBSDQAQXQAl5sAQFwAQRgAgWEBT1IAgXoAPDwAQSIALqIAgFQAj3YAgXwA6TsAQFEBNXYAQFsAIx4AQHsBVTcAgE4BLyIAgMcBWk4AgYABPkIAgYEBUlYAQFcAHnoAQQQBXwsAgYgAoz8AQYwAaJkAgWABZFwAgEQAPyoAQFgAEyUAQUcAGqgAQNYBaAUAQH4BaUwAgZgAFpEAgL0BalQAgDwAQh8AgZMA3TsAQZkAIDcAQTwAErIAQV0AEaIAQXkAEpEAgT0BbQEAQXMBcVgAQZMBWVQAgM8Bb10AgKEBb04AQDYAmXUAQDMA+nYAQNUBdgMAQDABWRQAQaIAGAcAQaIAkUAAgaQAezwAgCcBYQoAQaIBUEwAQRgBgQUAQCQAlRoAQIgBg1oAgbUBeFcAQaABmrAAgaUBoFoAgcgAM5IAQccBYVwAgc8ARDcAAAA=");
	auto r1 = incits378_from_buffer(b1.data());
	auto t1 = r1.value();
	auto c1 = BitCylinderSet(t1);

	auto b2 = base64_decode("Rk1SACAyMAAC8AAyAAAAAAJYAlgAxQDFAQAAAAB4QVUAtWwAgVUAx4cAQXQA3YIAgWIAnVYAQUMA4ZMAQVUA5TUAQX4AtJMAQSYAyYAAgSwA8YsAQR4A6IMAgXQA93kAQYMA/IAAgRYAsngAQUUAfWQAgZUArZYAgZQAxpAAQVoBASEAgZMBAYcAQV8AdLAAQQoA3oIAQaAAdZ0AgSABCzIAgaMA8owAgYgBDn0AQbEAiEEAgY4BIIAAgT0BHqEAgWoAW60AQbcA4Y4AQb0AbJkAQPsBIYcAQb0AuZMAgSYBJ5IAgRsBKosAQcMA1DgAgZcBL34AgVEBL1cAQaMBMYwAQa4ATEUAgOkBLYUAgOABLoAAgOAArXMAQTMARmMAgTwBMqgAgN8AwXkAgc0AljwAQWIBOBgAQTYBOzsAQdQA4ZMAQNcASxIAQSoBQJkAQacAOp0AgdoAqpUAQM4BN30AgYwBRBoAQccBRTgAQdwBLpgAgZ8BSSMAgawBSD8AgUsBSVwAgeUAVD8AgeMBRz4AgbYBTUIAQMIBEnsAQecA7TsAgMMAk3AAgeoAnZYAQPgBUi0AgLoBP3sAQegBKD0AQQUBU5MAQbsBYE8AQWwBYRAAgTIBYU4AQK4BQHQAQfwAqz4AgL4BZQMAQdsBaU8AQfQBblEAggMAcJcAQKQBNRwAQgMA4j0AQKwBbRAAgNYBbDwAQUIBcAgAgb4Bc1gAgScBd1EAQQIBd0kAgTsBd2AAgNQBfEwAQRcBfU4AQaQBgAEAQUEBgAUAQdkBgVYAQJABBXMAQPUBg0wAQYEBjAUAQIIBYhEAQIABcA0AQNkBlFQAgRcBlrAAgR0BnlcAQM0BvVwAQRIBybIAQSYBybIAgeUBy2UAgUIBygAAQJEBy2MAQVUBzV4AQNEBy2cAQPIBzlwAgSUBz1gAgQYB0q0AgM0B1AcAQKUB12EAgNkB1QcAgQ4B2lIAgb8B4F4AQV8B3VwAgUEB3lkAAAA=");
	auto r2 = incits378_from_buffer(b2.data());
	auto t2 = r2.value();
	auto c2 = BitCylinderSet(t2);

	cout << "num minutiae: " << t1.size() << endl;
	cout << "num cylinders: " << c1.cylinders.size() << endl;
	cout << "mem size: " << sizeof(c1.cylinders[0]) * c1.cylinders.size() * 2 + sizeof(c1.minutiae_references) << endl;

	cout << "num minutiae: " << t2.size() << endl;
	cout << "num cylinders: " << c2.cylinders.size() << endl;
	cout << "mem size: " << sizeof(c2.cylinders[0]) * c2.cylinders.size() * 2 + sizeof(c2.minutiae_references) << endl;

	auto similarities = local_matching(c1, t1, c2, t2);
	auto nP = calc_nP(c2, c1);
	auto score_LSS = global_matching(similarities, nP);
	auto score_LSA = global_matching_LSA(similarities, c1, c2);

	cout << "score_LSS: " << score_LSS << endl;
	cout << "score_LSA: " << score_LSA << endl;
	//assert(std::abs(score - 0.628780365) < 0.0000001); 

	std::vector<int> v;
	auto num_iters = 1000;
	for (auto i = 0; i < num_iters; ++i)
		v.push_back(i);

	const auto before = std::chrono::system_clock::now();

	std::for_each(std::execution::par, v.begin(), v.end(), [&c1, &c2, &t1, &t2](int i) {
		auto similarities = local_matching(c1, t1, c2, t2);
		auto nP = calc_nP(c1, c2);
		auto score = global_matching_LSA(similarities, c1, c2);
	});

	const us duration = std::chrono::system_clock::now() - before;

	std::cout << "average time: " << duration.count() / double(num_iters) << "us" << std::endl;
	cin.get();
}
catch (exception e)
{
	cout << e.what() << endl;
	cin.get();
}