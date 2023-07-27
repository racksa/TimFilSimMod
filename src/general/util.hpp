// util.hpp

// =============================================================================
// Include guard
#ifndef MY_UTIL_HEADER_INCLUDED
#define MY_UTIL_HEADER_INCLUDED

// =============================================================================
// Forward declared dependencies

// =============================================================================
// Included dependencies
#include <chrono>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include "../../config.hpp"

bool hasEnding (std::string const &fullString, std::string const &ending);


using clock_type    = std::chrono::high_resolution_clock;
using duration_type = std::chrono::duration<Real>;
static Real get_time_cilia() {
    static auto start_time = clock_type::now();
    return duration_type(clock_type::now()-start_time).count();
}

std::unordered_map<std::string, std::unordered_map<std::string, std::string>> parseINI(const std::string& filename);

std::string data_from_ini(std::string section, std::string variable);

#endif