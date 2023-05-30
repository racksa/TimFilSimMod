// mobility_solver.hpp

// =============================================================================
// Include guard
#ifndef MY_GENERIC_MOBILITY_SOLVER_HEADER_INCLUDED
#define MY_GENERIC_MOBILITY_SOLVER_HEADER_INCLUDED

// =============================================================================
// Forward declared dependencies
class swimmer;

// =============================================================================
// Included dependencies
#include <vector>
#include <chrono>
#include "../../general/matrix.hpp"
#include "../../../config.hpp"

class mobility_solver{

public:

  // =============================================================================
  // These pure virtual methods MUST be defined by the derived classes:

  virtual void free_host_memory() = 0;
  virtual void free_device_memory() = 0;
  virtual void allocate_host_memory() = 0;
  virtual void allocate_device_memory() = 0;

  virtual void copy_segment_positions_to_device() = 0;
  virtual void copy_segment_forces_to_device() = 0;
  virtual void copy_blob_positions_to_device() = 0;
  virtual void copy_blob_forces_to_device() = 0;

  virtual void copy_interparticle_blob_forces_to_host() = 0;
  virtual void copy_blob_velocities_to_host() = 0;
  virtual void copy_segment_velocities_to_host() = 0;

  virtual void apply_interparticle_forces() = 0;
  virtual void wait_for_device() = 0;
  virtual void evaluate_segment_segment_mobility() = 0;
  virtual void evaluate_segment_blob_mobility() = 0;
  virtual void evaluate_blob_blob_mobility() = 0;
  virtual void evaluate_blob_segment_mobility() = 0;
  virtual void evaluate_full_mobility() = 0;

        

  // =============================================================================
  // Everything below here is defined at the generic level:

  // CPU-side/host storage is used by the filaments and segments to store forces and positions directly to avoid copying.
  // As such, it must exist for any mobility solver. Moreover, these external classes expect to receive
  // raw pointers to locations for storing these data, and so that is the form we will require here, even
  // if it just ends up pointing to the storage inside a std::vector or something similar inside the derived class.
  Real *v_segs_host;
  Real *v_blobs_host;
  Real *x_segs_host;
  Real *x_blobs_host;
  Real *f_segs_host;
  Real *f_blobs_host;
  Real *f_blobs_repulsion_host;

  #if DISPLAYTIME
    float time_start;
    float read_position_time;
    float assemble_rhs_time;
    float solve_linear_system_time;
      float preconditioned_time;
      float system_mult_time;
        float mobility_mult_time;
        float k_mult_time;
      float gmres_solve_time;
      float found_solution_finishing_time;
  #endif

  ~mobility_solver();
  mobility_solver();

  void initialise();
  void finalise();
  void read_positions_and_forces(std::vector<swimmer>& swimmers);
  void compute_velocities(std::vector<swimmer>& swimmers, int& num_gmres_iterations, const int nt);
  bool compute_errors(matrix& error, const std::vector<swimmer>& swimmers, const int nt);
  void write_data(const int nt, const std::vector<swimmer>& swimmers);
  

  #if !USE_BROYDEN_FOR_EVERYTHING

    // linear system storage
    matrix rhs;
    matrix v_bodies;
    matrix Q;
    matrix H;
    matrix beta;
    matrix SN;
    matrix CS;

    // Reference matrices for preconditioning
    matrix body_mobility_reference;
    matrix KTKinv_reference;

    void assemble_rhs(const std::vector<swimmer>& swimmers, const int nt);
    matrix apply_preconditioner(const matrix& in, const std::vector<swimmer>& swimmers);
    matrix apply_preconditioner_new(const matrix& in, const std::vector<swimmer>& swimmers);
    matrix system_matrix_mult(const matrix& in, const std::vector<swimmer>& swimmers);
    matrix system_matrix_mult_new(const matrix& in, const std::vector<swimmer>& swimmers);
    int solve_linear_system(std::vector<swimmer>& swimmers);
    void make_body_reference_matrices(std::vector<swimmer>& swimmers);

  #endif

  #if DYNAMIC_PHASE_EVOLUTION

    std::vector<Real> gen_phase_force_refs;

  #endif

  #if DYNAMIC_SHAPE_ROTATION

    std::vector<Real> gen_angle_force_refs;

  #endif


}; // End of mobility_solver class definition

#endif // MY_GENERIC_MOBILITY_SOLVER_HEADER_INCLUDED
