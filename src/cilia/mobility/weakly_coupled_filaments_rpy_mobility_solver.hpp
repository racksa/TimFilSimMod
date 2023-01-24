// weakly_coupled_filaments_rpy_mobility_solver.hpp

// =============================================================================
// Include guard
#ifndef MY_WEAKLY_COUPLED_FILAMENTS_RPY_MOBILITY_SOLVER_HEADER_INCLUDED
#define MY_WEAKLY_COUPLED_FILAMENTS_RPY_MOBILITY_SOLVER_HEADER_INCLUDED

// =============================================================================
// Forward declared dependencies

// =============================================================================
// Included dependencies
#include "rpy_mobility_solver.hpp"

class weakly_coupled_filaments_rpy_mobility_solver : public rpy_mobility_solver{

public:

  // =============================================================================
  // Everything we have to define for the base class:

  void free_device_memory();
  void allocate_device_memory();

  void copy_segment_positions_to_device();
  void copy_segment_forces_to_device();

  void evaluate_segment_segment_mobility();

  // =============================================================================
  // Everything unique to this derived class:

  ~weakly_coupled_filaments_rpy_mobility_solver();
  weakly_coupled_filaments_rpy_mobility_solver();

  // GPU memory
  Real **f_fils_device; // Total force summed over all segments.
  Real **s_fils_device; // Total first-moment of force summed over all segments.
  Real **x_fils_device; // Average segment position for each filament.

}; // End of weakly_coupled_filaments_rpy_mobility_solver class definition

#endif // MY_WEAKLY_COUPLED_FILAMENTS_RPY_MOBILITY_SOLVER_HEADER_INCLUDED
