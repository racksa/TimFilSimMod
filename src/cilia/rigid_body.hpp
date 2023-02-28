// rigid_body.hpp

// =============================================================================
// Include guard
#ifndef MY_RIGID_BODY_HEADER_INCLUDED
#define MY_RIGID_BODY_HEADER_INCLUDED

// =============================================================================
// Forward declared dependencies

// =============================================================================
// Included dependencies
#include <vector>
#include <iostream>
#include "matrix.hpp"
#include "quaternion.hpp"
#include <random>
#include "../../config.hpp"

class rigid_body{

public:

  Real x[3];
  Real xm1[3];
  Real xm2[3];

  Real u[3];
  Real um1[3];
  quaternion q;
  quaternion qm1;
  matrix Q_init; // Rotation matrix associated with the initial guess for the body's quaternion.
  std::vector<Real> blob_references;
  std::vector<Real> polar_dir_refs;
  std::vector<Real> azi_dir_refs;
  std::vector<Real> normal_refs;

  ~rigid_body();
  rigid_body();

  void initial_setup(const int id, Real *const f_address, const Real *const data_from_file);
  void initial_guess(const int nt);
  void blob_positions(Real *const x_array) const;
  void update(const Real *const body_update);
  void write_reference_positions() const;
  void write_data(std::ofstream& body_state_file) const;
  void write_backup(std::ofstream& backup_file) const;

  #if PRESCRIBED_BODY_VELOCITIES

    void prescribed_translational_velocity(Real *const V, const Real t, const int id) const;
    void prescribed_rotational_velocity(Real *const W, const Real t, const int id) const;

  #endif

  #if USE_BROYDEN_FOR_EVERYTHING

    Real *blob_forces; // Stored globally in the mobility solver.
    std::vector<Real> blob_forces_m1;
    std::vector<Real> blob_forces_m2;

  #endif

  #if NO_CILIA_SQUIRMER

    Real max_cylindrical_radius;

  #endif

}; // End of rigid body class definition

#endif // MY_RIGID_BODY_HEADER_INCLUDED
