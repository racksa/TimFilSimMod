// filament.hpp

// =============================================================================
// Include guard
#ifndef MY_FILAMENT_HEADER_INCLUDED
#define MY_FILAMENT_HEADER_INCLUDED

// =============================================================================
// Forward declared dependencies

// =============================================================================
// Included dependencies
#include <vector>
#include <iostream>
#include "segment.hpp" // The filament inherits knowledge of quaternions through here
#include "matrix.hpp"
#include "../../config.hpp"

matrix body_frame_moment_lie_derivative(const quaternion& q1, const quaternion& q2, const bool upper);

class filament{

public:

  std::vector<segment> segments;
  Real strain_twist[3];

  // Jacobian-related members
  matrix inverse_jacobian;
  matrix elastic_clamping_block1;
  matrix elastic_clamping_block2;

  // Lagrange multipliers
  std::vector<Real> lambda;
  std::vector<Real> lambdam1;
  std::vector<Real> lambdam2;
  Real tether_lambda[3];
  Real tether_lambdam1[3];
  Real tether_lambdam2[3];

  // The forces and torques on the segments. They appear in the mobility solve, so we store them globally to avoid copying etc.
  Real *f;

  ~filament();
  filament();

  void initial_setup(const Real *const base_pos,
                      const Real *const dir,
                      const Real *const strain_twist_in,
                      const Real *const data_from_file,
                      Real *const x_address,
                      Real *const f_address,
                      const int fil_id);
  void robot_arm();
  void accept_state_from_rigid_body(const Real *const x_in, const Real *const u_in);
  void initial_guess(const int nt, const Real *const x_in, const Real *const u_in);
  void end_of_step(const int nt);
  matrix body_frame_moment(const int link_id);
  void internal_forces_and_torques(const int nt);
  matrix jacobian_lie_algebra_block(const int nt);
  void invert_approx_jacobian(const int nt);
  void update(const Real *const u);
  void write_data(std::ofstream& data_file, std::ofstream& tether_force_file) const;
  void write_backup(std::ofstream& data_file) const;

  #if INSTABILITY_CILIA

    Real clamp_lambda[3];
    Real clamp_lambdam1[3];
    Real clamp_lambdam2[3];
    Real perturbation1[3];
    Real perturbation2[3];

  #elif GEOMETRIC_CILIA

    vec3 my_vertical;
    vec3 fast_beating_direction;
    int step_at_start_of_transition;
    bool just_switched_from_fast_to_slow;
    bool just_switched_from_slow_to_fast;
    Real base_torque_magnitude_factor;
    quaternion body_q;
    quaternion body_qm1;

    Real find_my_angle();

  #elif PRESCRIBED_CILIA

    Real phase;
    Real phase_dot;
    Real omega0;
    quaternion body_q;
    quaternion body_qm1;
    std::vector<Real> vel_dir_phase;

    #if (DYNAMIC_SHAPE_ROTATION || WRITE_GENERALISED_FORCES)

      // If it isn't solved for dynamically, there is no rotation.
      // So we don't need to even store the angle in such cases; everything to do with it reduces to the identity.
      Real shape_rotation_angle;
      Real shape_rotation_angle_dot;
      std::vector<Real> vel_dir_angle;

    #endif

    #if FIT_TO_DATA_BEAT

      std::vector<Real> s_to_use;
      matrix Ax;
      matrix Ay;
      matrix Bx;
      matrix By;

      #if (FIT_TO_DATA_BEAT && !WRITE_GENERALISED_FORCES)

        matrix *s_to_use_ref_ptr;

      #endif

      void fill_fourier_coeff_mats();
      void fitted_shape_tangent(Real& tx, Real& ty, const Real s) const;
      matrix fitted_shape(const Real s) const;
      matrix fitted_shape_velocity_direction(const Real s) const;
      Real fitted_curve_length(const Real s) const;
      void find_fitted_shape_s();

    #elif BUILD_A_BEAT

      Real beat_switch_theta;

      Real build_a_beat_tangent_angle(const Real s) const;
      void build_a_beat_tangent(matrix& t, const Real s) const;
      void build_a_beat_tangent_phase_deriv(matrix& k, const Real s) const;

    #endif

  #endif

}; // End of filament class definition

  #if PRESCRIBED_CILIA

    #include <string>

    std::string reference_phase_generalised_force_file_name();
    std::string reference_angle_generalised_force_file_name();
    std::string reference_s_values_file_name();

  #endif

#endif // MY_FILAMENT_HEADER_INCLUDED
