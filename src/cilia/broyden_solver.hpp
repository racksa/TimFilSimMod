// broyden_solver.hpp

// =============================================================================
// Include guard
#ifndef MY_BROYDEN_SOLVER_HEADER_INCLUDED
#define MY_BROYDEN_SOLVER_HEADER_INCLUDED

// =============================================================================
// Forward declared dependencies
class swimmer;

// =============================================================================
// Included dependencies
#include <vector>
#include "matrix.hpp"

class broyden_solver {

public:

  matrix C;
  matrix D;
  matrix error;
  matrix new_error;
  matrix update;
  int iter;
  int max_iter;
  int total_iter;
  double avg_iter;

  broyden_solver();
  ~broyden_solver();

  #if !(PRESCRIBED_CILIA || NO_CILIA_SQUIRMER)

    void find_update(const std::vector<swimmer>& swimmers, const int nt);
    void end_of_iter(const std::vector<swimmer>& swimmers, const int nt, const int nt_start, const bool error_is_too_large);

  #endif


}; // End of broyden_solver class definition.

#endif // MY_BROYDEN_SOLVER_HEADER_INCLUDED
