// stokes_drag_mobility_solver.hpp

// =============================================================================
// Include guard
#ifndef MY_STOKES_DRAG_MOBILITY_SOLVER_HEADER_INCLUDED
#define MY_STOKES_DRAG_MOBILITY_SOLVER_HEADER_INCLUDED

// =============================================================================
// Forward declared dependencies

// =============================================================================
// Included dependencies
#include "rpy_mobility_solver.hpp"

class stokes_drag_mobility_solver : public rpy_mobility_solver{

public:

  // =============================================================================
  // Everything we have to define for the base class:

  void evaluate_segment_segment_mobility();
  void evaluate_segment_blob_mobility();
  void evaluate_blob_blob_mobility();
  void evaluate_blob_segment_mobility();

  // =============================================================================
  // Everything unique to this derived class:

  ~stokes_drag_mobility_solver();
  stokes_drag_mobility_solver();

}; // End of stokes_drag_mobility_solver class definition

#endif // MY_STOKES_DRAG_MOBILITY_SOLVER_HEADER_INCLUDED
