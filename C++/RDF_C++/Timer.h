#pragma once
#include <chrono>

struct Timer
{
	std::chrono::high_resolution_clock::time_point start;
	std::chrono::high_resolution_clock::time_point stop;

	void begin()
	{
		start = std::chrono::high_resolution_clock::now();
	}

	long long end()
	{
		stop = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
		return duration.count();
	}
};