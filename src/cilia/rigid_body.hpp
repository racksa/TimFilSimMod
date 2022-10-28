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
#include "config.hpp"

class rigid_body{

public:

  double x[3];
  double xm1[3];
  double xm2[3];

  double u[3];
  double um1[3];
  quaternion q;
  quaternion qm1;
  matrix Q_init; // Rotation matrix associated with the initial guess for the body's quaternion.

  std::vector<double> blob_references;
  std::vector<double> polar_dir_refs;
  std::vector<double> azi_dir_refs;
  std::vector<double> normal_refs;

  ~rigid_body();
  rigid_body();

  void initial_setup(const int id, double *const f_address, const double *const data_from_file);
  void initial_guess(const int nt);
  void blob_positions(double *const x_array) const;
  void update(const double *const body_update);
  void write_reference_positions() const;
  void write_data(std::ofstream& body_state_file) const;
  void write_backup(std::ofstream& backup_file) const;

  #if PRESCRIBED_BODY_VELOCITIES

    void prescribed_translational_velocity(double *const V, const double t, const int id) const;
    void prescribed_rotational_velocity(double *const W, const double t, const int id) const;

  #endif

  #if USE_BROYDEN_FOR_EVERYTHING

    double *blob_forces; // Stored globally in the mobility solver.
    std::vector<double> blob_forces_m1;
    std::vector<double> blob_forces_m2;

  #endif

  #if NO_CILIA_SQUIRMER

    double max_cylindrical_radius;

  #endif

}; // End of rigid body class definition

#endif // MY_RIGID_BODY_HEADER_INCLUDED
