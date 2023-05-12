// util.hpp

// =============================================================================
// Include guard
#ifndef MY_UTIL_HEADER_INCLUDED
#define MY_UTIL_HEADER_INCLUDED

#include <chrono>

#include "../../config.hpp"

// =============================================================================
// Forward declared dependencies

// =============================================================================
// Included dependencies
#include <iostream>
#include "../../config.hpp"

bool hasEnding (std::string const &fullString, std::string const &ending);


using clock_type    = std::chrono::high_resolution_clock;
using duration_type = std::chrono::duration<Real>;
static Real get_time() {
    static auto start_time = clock_type::now();
    return duration_type(clock_type::now()-start_time).count();
}

#endif