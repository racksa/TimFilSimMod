// fcm_mobility_solver.hpp

// =============================================================================
// Include guard
#ifndef MY_FCM_MOBILITY_SOLVER_HEADER_INCLUDED
#define MY_FCM_MOBILITY_SOLVER_HEADER_INCLUDED

// =============================================================================
// Forward declared dependencies

// =============================================================================
// Included dependencies
#include "mobility_solver.hpp"
#include "../../../../CUFCM/config.hpp"
#include "../../../../CUFCM/CUFCM_DATA.cuh"
#include "../../../../CUFCM/CUFCM_SOLVER.cuh"
#include "../../../../CUFCM/util/cuda_util.hpp"
#include "../../../../CUFCM/CUFCM_CELLLIST.cuh"

class fcm_mobility_solver : public mobility_solver{

public:

  // =============================================================================
  // Everything we have to define for the base class:

  void free_host_memory();
  void free_device_memory();
  void allocate_host_memory();
  void allocate_device_memory();

  void copy_segment_positions_to_device();
  void copy_segment_forces_to_device();
  void copy_blob_positions_to_device();
  void copy_blob_forces_to_device();

  void copy_interparticle_blob_forces_to_host();
  void copy_blob_velocities_to_host();
  void copy_segment_velocities_to_host();

  void apply_interparticle_forces();
  void wait_for_device();
  void evaluate_segment_segment_mobility();
  void evaluate_segment_blob_mobility();
  void evaluate_blob_blob_mobility();
  void evaluate_blob_segment_mobility();

  // =============================================================================
  // Everything unique to this derived class:

  // GPU info
  int num_gpus;
  int *num_segs;
  int *num_blobs;

  // GPU memory
  Real **v_segs_device;
  Real **v_blobs_device;
  Real **x_segs_device;
  Real **x_blobs_device;
  Real **f_segs_device;
  Real **f_blobs_device;
  Real **f_blobs_repulsion_device;

  Pars pars;

  FCM_solver *cufcm_solver;
  

  ~fcm_mobility_solver();
  fcm_mobility_solver();
  void write_repulsion();

}; // End of fcm_mobility_solver class definition

#endif // MY_FCM_MOBILITY_SOLVER_HEADER_INCLUDED
