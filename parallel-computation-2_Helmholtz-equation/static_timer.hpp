#pragma once

#include <chrono>
#include <iostream>

///#include <iomanip> // for getting formated date
///#include <ctime>

// # StaticTimer #
// Small class used for time measurements, fully static aka does not require creating local instances
// - Thread safe due to 'chrono::steady_clock' guarantees
struct StaticTimer {
	using Clock = std::chrono::steady_clock;
	using Milliseconds = std::chrono::milliseconds;
	
	inline static void start() {
		_start_timepoint = Clock::now();
	}

	inline static long long end() {
		return std::chrono::duration_cast<Milliseconds>(Clock::now() - _start_timepoint).count();
			// time since last StaticTimer::start() call in ms 
	}

private:
	inline static Clock::time_point _start_timepoint = StaticTimer::Clock::now();
		// 'inline static' requires C++17
};

/*
std::string get_date_string() {
	auto start = std::chrono::system_clock::now();
	auto legacyStart = std::chrono::system_clock::to_time_t(start);
	char tmBuff[30];
	ctime_s(tmBuff, sizeof(tmBuff), &legacyStart);

	return std::string(tmBuff);
}
*/