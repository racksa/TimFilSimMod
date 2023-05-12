// swimmer.hpp

// =============================================================================
// Include guard
#ifndef MY_SWIMMER_HEADER_INCLUDED
#define MY_SWIMMER_HEADER_INCLUDED

// =============================================================================
// Forward declared dependencies

// =============================================================================
// Included dependencies
#include <iostream>
#include <vector>
#include "../general/matrix.hpp"
#include "rigid_body.hpp"
#include "filament.hpp"
#include "config.hpp"

class swimmer{

public:

  rigid_body body;

  std::vector<Real> filament_references;
  std::vector<Real> polar_dir_refs;
  std::vector<Real> azi_dir_refs;
  std::vector<Real> normal_refs;
  
  std::vector<filament> filaments;

  matrix f;

  ~swimmer();
  swimmer();

  void initial_setup(const int id, const Real *const data_from_file, Real *const x_segs_address, Real *const f_segs_address, Real *const f_blobs_address);
  void initial_guess(const int nt);
  void forces_and_torques(const int nt, int id);
  void end_of_step(const int nt);
  void write_reference_positions() const;
  void write_data(std::ofstream& seg_state_file, std::ofstream& body_state_file, std::ofstream& tether_force_file) const;
  void write_backup(std::ofstream& backup_file) const;

  #if PRESCRIBED_CILIA

    matrix KTMinvK_inv;

    void make_precon_mat();

  #else

    #if !INFINITE_PLANE_WALL

      matrix schur_mat_inv;
      std::vector<matrix> jacobian_B_blocks;
      matrix body_mobility_reference;

    #endif

    void prepare_jacobian_inv(const int nt);
    matrix jacobian_inv_mult(const matrix& in, const int nt) const;
    void update(const Real *const swimmer_update);

  #endif

}; // End of swimmer class definition

#endif // MY_SWIMMER_HEADER_INCLUDED
