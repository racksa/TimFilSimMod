// mobility_solver.cpp

#include <iomanip>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <string>
#include "mobility_solver.hpp"
#include "swimmer.hpp"
#include "omp.h"

mobility_solver::~mobility_solver(){}

mobility_solver::mobility_solver(){}

void mobility_solver::finalise(){

  free_host_memory();

  free_device_memory();

}

void mobility_solver::initialise(){

  allocate_host_memory();

  allocate_device_memory();

  // Initialise these forces as zeros as they will be used for GMRES initial guesses.
  for (int n = 0; n < 3*NSWIM*NBLOB; n++){

    f_blobs_host[n] = 0.0;

  }

  #if !USE_BROYDEN_FOR_EVERYTHING

    #if DYNAMIC_PHASE_EVOLUTION

      std::ifstream generalised_phase_force_file(reference_phase_generalised_force_file_name());

      if (generalised_phase_force_file.good()){

        int num_saves;
        generalised_phase_force_file >> num_saves;
        gen_phase_force_refs = std::vector<Real>(num_saves);
        for (int n = 0; n < num_saves; n++){

          generalised_phase_force_file >> gen_phase_force_refs[n];

        }

        generalised_phase_force_file.close();

      } else {

        std::cout << std::endl;
        std::cout << std::endl;
        std::cout << "Required reference file for generalised phase forces was not found." << std::endl;
        std::cout << "Run an appropriate simulation with the WRITE_GENERALISED_FORCES option enabled." << std::endl;
        std::cout << std::endl;
        std::cout << std::endl;

        std::exit(-1);

      }

    #endif

    #if DYNAMIC_SHAPE_ROTATION

      std::ifstream generalised_angle_force_file(reference_angle_generalised_force_file_name());

      if (generalised_angle_force_file.good()){

        int num_saves;
        generalised_angle_force_file >> num_saves;
        gen_angle_force_refs = std::vector<Real>(num_saves);
        for (int n = 0; n < num_saves; n++){

          generalised_angle_force_file >> gen_angle_force_refs[n];

        }

        generalised_angle_force_file.close();

      } else {

        std::cout << std::endl;
        std::cout << std::endl;
        std::cout << "Required reference file for generalised angle forces was not found." << std::endl;
        std::cout << "Run an appropriate simulation with the WRITE_GENERALISED_FORCES option enabled." << std::endl;
        std::cout << std::endl;
        std::cout << std::endl;

        std::exit(-1);

      }

    #endif

    int system_size = 0;

    #if !INFINITE_PLANE_WALL

      system_size += 3*NSWIM*NBLOB;

    #endif

    #if !PRESCRIBED_BODY_VELOCITIES

      system_size += 6*NSWIM;

      body_mobility_reference = identity(6);
      KTKinv_reference = identity(6);

    #endif

    #if PRESCRIBED_CILIA

      system_size += 3*NSWIM*NFIL*NSEG;

      #if DYNAMIC_PHASE_EVOLUTION

        system_size += NSWIM*NFIL;

      #endif

      #if DYNAMIC_SHAPE_ROTATION

        system_size += NSWIM*NFIL;

      #endif

    #endif

    if (system_size > 0){

      rhs = matrix(system_size, 1);

      v_bodies = matrix(6*NSWIM, 1);
      v_bodies.zero();

      #if USE_GMRES_FOR_LINEAR_SYSTEM

        const int max_iter = std::min<int>(system_size, MAX_LINEAR_SYSTEM_ITER);

        Q = matrix(system_size, max_iter+1);
        Q.zero();

        beta = matrix(max_iter+1, 1);
        beta.zero();

        H = matrix(max_iter+1, max_iter+1);
        H.zero();

        SN = matrix(max_iter, 1);
        SN.zero();

        CS = matrix(max_iter, 1);
        CS.zero();

      #endif

    }

  #endif

}

void mobility_solver::read_positions_and_forces(std::vector<swimmer>& swimmers){

  copy_segment_positions_to_device();
  copy_segment_forces_to_device();

  #if (USE_BROYDEN_FOR_EVERYTHING && !INFINITE_PLANE_WALL)

      copy_blob_forces_to_device();

  #endif

  // Calculate the blob positions, storing them in the appropriate host memory
  for (int n = 0; n < NSWIM; n++){

    swimmers[n].body.blob_positions(&x_blobs_host[3*n*NBLOB]);

  }

  copy_blob_positions_to_device();

  #if !PRESCRIBED_CILIA

      apply_interparticle_forces();

      #if !PRESCRIBED_BODY_VELOCITIES

        copy_interparticle_blob_forces_to_host();

      #endif

  #endif

  #if !(PRESCRIBED_BODY_VELOCITIES || PRESCRIBED_CILIA)

    wait_for_device();
    
    // #pragma omp parallel
    // {
    //   int swimmer_per_thread = NSWIM/omp_get_num_threads();
    //   for(int i = 0; i < swimmer_per_thread; i++){
    //     int n = swimmer_per_thread*omp_get_thread_num() + i;
    //     matrix& f = swimmers[n].f;
    //     const matrix R = swimmers[n].body.q.rot_mat();

    //     for (int m = 0; m < NBLOB; m++){

    //       const int id = 3*(n*NBLOB + m);
    //       const matrix r = R*matrix(3, 1, &swimmers[n].body.blob_references[3*m]);

    //       f(0) += f_blobs_repulsion_host[id];
    //       f(1) += f_blobs_repulsion_host[id + 1];
    //       f(2) += f_blobs_repulsion_host[id + 2];

    //       f(3) += r(1)*f_blobs_repulsion_host[id + 2] - r(2)*f_blobs_repulsion_host[id + 1];
    //       f(4) += r(2)*f_blobs_repulsion_host[id] - r(0)*f_blobs_repulsion_host[id + 2];
    //       f(5) += r(0)*f_blobs_repulsion_host[id + 1] - r(1)*f_blobs_repulsion_host[id];

    //     }
    //   }
    // }

    for (int n = 0; n < NSWIM; n++){

      matrix& f = swimmers[n].f;
      const matrix R = swimmers[n].body.q.rot_mat();

      for (int m = 0; m < NBLOB; m++){

        const int id = 3*(n*NBLOB + m);
        const matrix r = R*matrix(3, 1, &swimmers[n].body.blob_references[3*m]);

        f(0) += f_blobs_repulsion_host[id];
        f(1) += f_blobs_repulsion_host[id + 1];
        f(2) += f_blobs_repulsion_host[id + 2];

        f(3) += r(1)*f_blobs_repulsion_host[id + 2] - r(2)*f_blobs_repulsion_host[id + 1];
        f(4) += r(2)*f_blobs_repulsion_host[id] - r(0)*f_blobs_repulsion_host[id + 2];
        f(5) += r(0)*f_blobs_repulsion_host[id + 1] - r(1)*f_blobs_repulsion_host[id];

      }

    }

  #endif

}

#if !USE_BROYDEN_FOR_EVERYTHING

  Real cyclic_cubic_spline_interpolation(const Real phase, const std::vector<Real>& in){

    Real phi_index = 0.5*phase/PI; // = phase/(2*pi)
    phi_index -= std::floor(phi_index); // Map this ratio into [0,1]
    phi_index *= in.size();

    Real q;

    // Cubic Hermite spline interpolation
    int phi_index_int_lower_bound = int(phi_index); // Rounds towards 0.

    if (phi_index_int_lower_bound == in.size()){

      // This can only be true if phi_index == in.size()
      // I don't think we can ever actually enter here unless rounding error messes things up,
      // but better safe than sorry.
      q = in[0];

    } else {

      // Otherwise, we're safely in the interior
      const int phi_index_int_upper_bound = (phi_index_int_lower_bound == in.size() - 1) ? 0 : (phi_index_int_lower_bound + 1);

      const Real q_lower = in[phi_index_int_lower_bound];
      const Real q_prime_lower = (phi_index_int_lower_bound == 0) ? (in[1] - in[0]) : 0.5*(in[phi_index_int_upper_bound] - in[phi_index_int_lower_bound-1]);

      const Real q_upper = in[phi_index_int_upper_bound];
      const Real q_prime_upper = (phi_index_int_upper_bound == 0) ? (in[0] - in[in.size()-1])
                                                                    : 0.5*(in[(phi_index_int_upper_bound == in.size()-1) ? 0 : phi_index_int_upper_bound+1] - in[phi_index_int_lower_bound]);

      const Real t = phi_index - phi_index_int_lower_bound;
      const Real t2 = t*t;
      const Real t3 = t2*t;

      q = q_lower*(2.0*t3 - 3.0*t2 + 1.0) + q_prime_lower*(t3 - 2.0*t2 + t) + q_upper*(3.0*t2 - 2.0*t3) + q_prime_upper*(t3 - t2);

    }

    return q;

  };

  void mobility_solver::assemble_rhs(const std::vector<swimmer>& swimmers, const int nt){

    #if PRESCRIBED_BODY_VELOCITIES

      const Real t = DT*(nt+1.0);

      for (int n = 0; n < NSWIM; n++){

        swimmers[n].body.prescribed_translational_velocity(&v_bodies.data[6*n], t, n);
        swimmers[n].body.prescribed_rotational_velocity(&v_bodies.data[6*n + 3], t, n);

      }

    #endif

    #if PRESCRIBED_CILIA

      rhs.zero();

      for (int n = 0; n < NSWIM; n++){

        for (int i = 0; i < NFIL; i++){

          #if DYNAMIC_SHAPE_ROTATION

            // Interpolate to find the generalised driving force
            Real q_angle = cyclic_cubic_spline_interpolation(swimmers[n].filaments[i].phase, gen_angle_force_refs);

            // Scale if the natural frequency of this cilium differs from the reference case
            q_angle *= 0.5*swimmers[n].filaments[i].omega0/PI;

            // Apply the torsional spring force
            Real Qbar = 0.0;
            for (int ii = 0; ii < gen_angle_force_refs.size(); ii++){
              Qbar += std::abs(gen_angle_force_refs[ii]);
            }
            Qbar /= gen_angle_force_refs.size();
            q_angle -= TORSIONAL_SPRING_MAGNITUDE_FACTOR*Qbar*swimmers[n].filaments[i].shape_rotation_angle;

          #endif

          #if DYNAMIC_PHASE_EVOLUTION

            // Interpolate to find the generalised driving force
            Real q_phase = cyclic_cubic_spline_interpolation(swimmers[n].filaments[i].phase, gen_phase_force_refs);

            // Scale if the natural frequency of this cilium differs from the reference case
            q_phase *= 0.5*swimmers[n].filaments[i].omega0/PI;

            // Store minus the generalised force in the RHS
            #if PRESCRIBED_BODY_VELOCITIES

              rhs(3*NSWIM*(NFIL*NSEG + NBLOB) + n*NFIL + i) = -q_phase;

            #else

              rhs(3*NSWIM*(NFIL*NSEG + NBLOB + 2) + n*NFIL + i) = -q_phase;

            #endif

            #if DYNAMIC_SHAPE_ROTATION

              #if PRESCRIBED_BODY_VELOCITIES

                rhs(3*NSWIM*(NFIL*NSEG + NBLOB) + n*NFIL + i + NSWIM*NFIL) = -q_angle;

              #else

                rhs(3*NSWIM*(NFIL*NSEG + NBLOB + 2) + n*NFIL + i + NSWIM*NFIL) = -q_angle;

              #endif

            #endif

          #else

            const int index = 3*NSEG*(i + NFIL*n);

            for (int m = 0; m < 3*NSEG; m++){

              rhs(index + m) = swimmers[n].filaments[i].phase_dot*swimmers[n].filaments[i].vel_dir_phase[m];

            }

            #if DYNAMIC_SHAPE_ROTATION

              #if PRESCRIBED_BODY_VELOCITIES

                rhs(3*NSWIM*(NFIL*NSEG + NBLOB) + n*NFIL + i) = -q_angle;

              #else

                rhs(3*NSWIM*(NFIL*NSEG + NBLOB + 2) + n*NFIL + i) = -q_angle;

              #endif

            #endif

          #endif

        }

        #if (PRESCRIBED_BODY_VELOCITIES && !INFINITE_PLANE_WALL) // All contributions are 0 in the half-space special case.

          const matrix R = swimmers[n].body.q.rot_mat();

          for (int m = 0; m < NBLOB; m++){

            const matrix r = R*matrix(3, 1, &swimmers[n].body.blob_references[3*m]);

            const int id = 3*(NSWIM*NFIL*NSEG + n*NBLOB + m);

            // K mult
            rhs(id) += v_bodies(6*n) + v_bodies(6*n + 4)*r(2) - v_bodies(6*n + 5)*r(1);
            rhs(id + 1) += v_bodies(6*n + 1) + v_bodies(6*n + 5)*r(0) - v_bodies(6*n + 3)*r(2);
            rhs(id + 2) += v_bodies(6*n + 2) + v_bodies(6*n + 3)*r(1) - v_bodies(6*n + 4)*r(0);

          }

          for (int m = 0; m < NFIL; m++){

            for (int k = 0; k < NSEG; k++){

              const int id = 3*(n*NFIL*NSEG + m*NSEG + k);

              matrix r(3,1);

              r(0) = swimmers[n].filaments[m].segments[k].x[0] - swimmers[n].body.x[0];
              r(1) = swimmers[n].filaments[m].segments[k].x[1] - swimmers[n].body.x[1];
              r(2) = swimmers[n].filaments[m].segments[k].x[2] - swimmers[n].body.x[2];

              // K mult
              rhs(id) += v_bodies(6*n) + v_bodies(6*n + 4)*r(2) - v_bodies(6*n + 5)*r(1);
              rhs(id + 1) += v_bodies(6*n + 1) + v_bodies(6*n + 5)*r(0) - v_bodies(6*n + 3)*r(2);
              rhs(id + 2) += v_bodies(6*n + 2) + v_bodies(6*n + 3)*r(1) - v_bodies(6*n + 4)*r(0);

            }

          }

        #endif

        // N.B. If there is a net external force or torque on the swimmer, then -F
        // should appear in the RHS instead of zeros in the rows corresponding to the force-
        // and torque-balances. Note that this will need to be the total -F -- i.e. including
        // any external forces and torques on the segments -- and thus is NOT the same as swimmer[n].f.
        // Indeed, if it is important that only blobs are seen as having external forces, for example,
        // then we may even have to change how we setup the system.

      }

    #else

      wait_for_device();

      for (int n = 0; n < 3*NSWIM*NBLOB; n++){

        rhs(n) = -v_blobs_host[n];

      }

      #if NO_CILIA_SQUIRMER

        // This is all but copied verbatim from the filament object.
        // I should probably come up with a better system.
        matrix Ap(3,3);
        Ap(0,0) = 3.091205e-01; Ap(0,1) = -1.993712e+00; Ap(0,2) = 1.091847e+00;
        Ap(1,0) = -2.952120e-01; Ap(1,1) = 7.189021e-01; Ap(1,2) = -6.328589e-01;
        Ap(2,0) = -9.868018e-02; Ap(2,1) = 3.312812e-01; Ap(2,2) = -2.578206e-01;

        matrix An(3,3);
        An(0,0) = -2.468656e-01; An(0,1) = 1.215869e+00; An(0,2) = -5.867236e-01;
        An(1,0) = 3.075836e-01; An(1,1) = -6.966322e-01; An(1,2) = 4.554729e-01;
        An(2,0) = 2.378489e-02; An(2,1) = -1.355582e-01; An(2,2) = 7.315831e-02;

        matrix Bp(3,3);
        Bp(0,0) = 9.165546e-01; Bp(0,1) = 8.848315e-02; Bp(0,2) = -4.838832e-01;
        Bp(1,0) = -2.946211e-01; Bp(1,1) = 1.057964e+00; Bp(1,2) = -6.335721e-01;
        Bp(2,0) = 6.122744e-02; Bp(2,1) = -2.319094e-01; Bp(2,2) = 1.906292e-01;

        matrix Bn(3,3);
        Bn(0,0) = 1.336807e-01; Bn(0,1) = -6.807852e-01; Bn(0,2) = 6.008592e-01;
        Bn(1,0) = 8.146824e-02; Bn(1,1) = 3.472676e-01; Bn(1,2) = -2.744220e-01;
        Bn(2,0) = 3.615272e-02; Bn(2,1) = 8.619119e-02; Bn(2,2) = -6.122992e-02;

        // Absorb the omega = 2pi
        Ap *= 2.0*PI;
        An *= 2.0*PI;
        Bp *= 2.0*PI;
        Bn *= 2.0*PI;

        // Match the ratio of cilium length to body length.
        const Real cilia_length_scale = 18.0*(2.0*FIL_LENGTH)/785.0;
        Ap *= cilia_length_scale;
        An *= cilia_length_scale;
        Bp *= cilia_length_scale;
        Bn *= cilia_length_scale;

        const Real MCW_length = cilia_length_scale; // TODO: This should be in the config, but will require other things being added to the config too.

        // Slip velocities for pure squirmer:
        for (int n = 0; n < NSWIM; n++){

          const std::vector<Real>& polar = swimmers[n].body.polar_dir_refs;
          const std::vector<Real>& normal = swimmers[n].body.normal_refs;
          const std::vector<Real>& refs = swimmers[n].body.blob_references;

          // Some other things we will need to evaluate velocities:
          const int num_fourier_modes = Ap.num_rows;

          matrix s_vec(3,1);
          s_vec(0) = 1.0;
          s_vec(1) = 1.0;
          s_vec(2) = 1.0;

          const matrix Q = swimmers[n].body.q.rot_mat();

          const Real MCW_length_in_radians = MCW_length/swimmers[n].body.max_cylindrical_radius;

          for (int m = 0; m < NBLOB; m++){

            const Real phi = std::atan2(refs[3*m + 1], refs[3*m]); // Azimuthal angle coordinate of the point on the surface in a body-fixed spherical cordinate system.

            // What would be the phase of a filament attached to this point on the surface?
            Real phase = 2.0*PI*nt*DT; // = 2.0*PI*t/T, since the period is T = 1.
            phase += 2.0*PI*phi/MCW_length_in_radians; // Shift due to the MCW.

         /*   // What would the velocity of the tip of this filament be?
            matrix cos_vec(1, num_fourier_modes);
            matrix sin_vec(1, num_fourier_modes);
            for (int n = 1; n <= num_fourier_modes; n++){
              cos_vec(n-1) = n*std::cos(n*phase);
              sin_vec(n-1) = -n*std::sin(n*phase);
            }*/
            const Real p_coeff = (0.1019 + 0.0064*std::sin(phase))*2.0*FIL_LENGTH; //*AXIS_DIR_BODY_LENGTH; //(cos_vec*Bp + sin_vec*Ap)*s_vec;
            const Real n_coeff = 0.0; //(cos_vec*Bn + sin_vec*An)*s_vec;
            matrix p_vec(3,1);
            p_vec(0) = polar[3*m];
            p_vec(1) = polar[3*m + 1];
            p_vec(2) = polar[3*m + 2];
            matrix n_vec(3,1);
            n_vec(0) = normal[3*m];
            n_vec(1) = normal[3*m + 1];
            n_vec(2) = normal[3*m + 2];
            const matrix v = Q*(p_coeff*p_vec + n_coeff*n_vec);

            // Store the slip velocity.
            rhs.add_to_block(3*(n*NBLOB + m), 3, v);
          }

        }

      #endif

      for (int n = 0; n < NSWIM; n++){

        #if PRESCRIBED_BODY_VELOCITIES

          const matrix R = swimmers[n].body.q.rot_mat();

          for (int m = 0; m < NBLOB; m++){

            const matrix r = R*matrix(3, 1, &swimmers[n].body.blob_references[3*m]);

            // K mult
            rhs(3*(n*NBLOB + m)) += v_bodies(6*n) + v_bodies(6*n + 4)*r(2) - v_bodies(6*n + 5)*r(1);
            rhs(3*(n*NBLOB + m) + 1) += v_bodies(6*n + 1) + v_bodies(6*n + 5)*r(0) - v_bodies(6*n + 3)*r(2);
            rhs(3*(n*NBLOB + m) + 2) += v_bodies(6*n + 2) + v_bodies(6*n + 3)*r(1) - v_bodies(6*n + 4)*r(0);

          }

        #else

          // This part of the vector should contain minus the total force and torque on the rigid body,
          // which is equal to the total force and torque across all filament segments minus the net
          // external force and torque on the swimmer by the force- and torque-free swimming conditions.
          rhs.set_block(3*NSWIM*NBLOB + 6*n, 6, -swimmers[n].f);

        #endif

      }

    #endif

  }

  matrix mobility_solver::apply_preconditioner(const matrix& in, const std::vector<swimmer>& swimmers){

    #if PRESCRIBED_CILIA

      matrix out = in;

      const Real seg_mob_inv = 6.0*PI*MU*RSEG;
      const Real blob_mob_inv = 6.0*PI*MU*RBLOB;

      out.multiply_block(0, 3*NSWIM*NFIL*NSEG, seg_mob_inv);
      out.multiply_block(3*NSWIM*NFIL*NSEG, 3*NSWIM*NBLOB, blob_mob_inv);

      #if PRESCRIBED_BODY_VELOCITIES

        #if DYNAMIC_PHASE_EVOLUTION

          out.multiply_block(3*NSWIM*(NBLOB + NFIL*NSEG), NSWIM*NFIL, -1.0);

          #if DYNAMIC_SHAPE_ROTATION

            out.multiply_block(3*NSWIM*(NBLOB + NFIL*NSEG) + NSWIM*NFIL, NSWIM*NFIL, -1.0);

          #endif

          for (int n = 0; n < NSWIM; n++){

            for (int m = 0; m < NFIL; m++){

              Real dot_phase = 0.0;
              Real Knorm_phase = 0.0;

              #if DYNAMIC_SHAPE_ROTATION

                Real dot_angle = 0.0;
                Real Knorm_angle = 0.0;

              #endif

              for (int k = 0; k < 3*NSEG; k++){

                dot_phase += swimmers[n].filaments[m].vel_dir_phase[k]*out(3*NSEG*(n*NFIL + m) + k);
                Knorm_phase += swimmers[n].filaments[m].vel_dir_phase[k]*swimmers[n].filaments[m].vel_dir_phase[k];

                #if DYNAMIC_SHAPE_ROTATION

                  dot_angle += swimmers[n].filaments[m].vel_dir_angle[k]*out(3*NSEG*(n*NFIL + m) + k);
                  Knorm_angle += swimmers[n].filaments[m].vel_dir_angle[k]*swimmers[n].filaments[m].vel_dir_angle[k];

                #endif

              }

              out(3*NSWIM*(NBLOB + NFIL*NSEG) + n*NFIL + m) -= dot_phase;
              out(3*NSWIM*(NBLOB + NFIL*NSEG) + n*NFIL + m) /= Knorm_phase*seg_mob_inv;

              #if DYNAMIC_SHAPE_ROTATION

                out(3*NSWIM*(NBLOB + NFIL*NSEG) + n*NFIL + m + NSWIM*NFIL) -= dot_angle;
                out(3*NSWIM*(NBLOB + NFIL*NSEG) + n*NFIL + m + NSWIM*NFIL) /= Knorm_angle*seg_mob_inv;

              #endif

              for (int k = 0; k < 3*NSEG; k++){

                out(3*NSEG*(n*NFIL + m) + k) += seg_mob_inv*out(3*NSWIM*(NBLOB + NFIL*NSEG) + n*NFIL + m)*swimmers[n].filaments[m].vel_dir_phase[k];

                #if DYNAMIC_SHAPE_ROTATION

                  out(3*NSEG*(n*NFIL + m) + k) += seg_mob_inv*out(3*NSWIM*(NBLOB + NFIL*NSEG) + n*NFIL + m + NSWIM*NFIL)*swimmers[n].filaments[m].vel_dir_angle[k];

                #endif

              }

            }

          }

        #elif DYNAMIC_SHAPE_ROTATION

          out.multiply_block(3*NSWIM*(NBLOB + NFIL*NSEG), NSWIM*NFIL, -1.0);

          for (int n = 0; n < NSWIM; n++){

            for (int m = 0; m < NFIL; m++){

              Real dot_angle = 0.0;
              Real Knorm_angle = 0.0;

              for (int k = 0; k < 3*NSEG; k++){

                dot_angle += swimmers[n].filaments[m].vel_dir_angle[k]*out(3*NSEG*(n*NFIL + m) + k);
                Knorm_angle += swimmers[n].filaments[m].vel_dir_angle[k]*swimmers[n].filaments[m].vel_dir_angle[k];

              }

              out(3*NSWIM*(NBLOB + NFIL*NSEG) + n*NFIL + m) -= dot_angle;
              out(3*NSWIM*(NBLOB + NFIL*NSEG) + n*NFIL + m) /= Knorm_angle*seg_mob_inv;

              for (int k = 0; k < 3*NSEG; k++){

                out(3*NSEG*(n*NFIL + m) + k) += seg_mob_inv*out(3*NSWIM*(NBLOB + NFIL*NSEG) + n*NFIL + m)*swimmers[n].filaments[m].vel_dir_angle[k];

              }

            }

          }

        #endif

      #else

        out.multiply_block(3*NSWIM*(NBLOB + NFIL*NSEG), 6*NSWIM, -1.0);

        for (int n = 0; n < NSWIM; n++){

          const int out_id = 3*NSWIM*(NBLOB + NFIL*NSEG) + 6*n;

          for (int m = 0; m < NFIL; m++){

            const std::vector<segment>& segments = swimmers[n].filaments[m].segments;

            for (int k = 0; k < NSEG; k++){

              const int seg_force_id = 3*(n*NFIL*NSEG + m*NSEG + k);

              out(out_id) -= out(seg_force_id);
              out(out_id + 1) -= out(seg_force_id + 1);
              out(out_id + 2) -= out(seg_force_id + 2);

              const Real diff[3] = {segments[k].x[0] - swimmers[n].body.x[0], segments[k].x[1] - swimmers[n].body.x[1], segments[k].x[2] - swimmers[n].body.x[2]};
              out(out_id + 3) -= diff[1]*out(seg_force_id + 2) - diff[2]*out(seg_force_id + 1);
              out(out_id + 4) -= diff[2]*out(seg_force_id) - diff[0]*out(seg_force_id + 2);
              out(out_id + 5) -= diff[0]*out(seg_force_id + 1) - diff[1]*out(seg_force_id);

            }

          }

          const matrix Q = swimmers[n].body.q.rot_mat();

          for (int m = 0; m < NBLOB; m++){

            const int blob_force_id = 3*(NSWIM*NFIL*NSEG + n*NBLOB + m);

            out(out_id) -= out(blob_force_id);
            out(out_id + 1) -= out(blob_force_id + 1);
            out(out_id + 2) -= out(blob_force_id + 2);

            const matrix diff = Q*matrix(3, 1, &swimmers[n].body.blob_references[3*m]);
            out(out_id + 3) -= diff(1)*out(blob_force_id + 2) - diff(2)*out(blob_force_id + 1);
            out(out_id + 4) -= diff(2)*out(blob_force_id) - diff(0)*out(blob_force_id + 2);
            out(out_id + 5) -= diff(0)*out(blob_force_id + 1) - diff(1)*out(blob_force_id + 0);

          }

          #if DYNAMIC_PHASE_EVOLUTION

            out.multiply_block(3*NSWIM*(NBLOB + NFIL*NSEG + 2), NFIL*NSWIM, -1.0);

            for (int m = 0; m < NFIL; m++){

              const int phase_id = 3*NSWIM*(NBLOB + NFIL*NSEG + 2) + NFIL*n + m;

              Real Knormsq = 0.0;

              matrix v1(3,1), v2(3,1);
              v1.zero();
              v2.zero();

              for (int k = 0; k < NSEG; k++){

                const int seg_force_id = 3*(n*NFIL*NSEG + m*NSEG + k);

                matrix K(3,1);
                K(0) = swimmers[n].filaments[m].vel_dir_phase[3*k];
                K(1) = swimmers[n].filaments[m].vel_dir_phase[3*k + 1];
                K(2) = swimmers[n].filaments[m].vel_dir_phase[3*k + 2];

                // Additional -K^T part
                out(phase_id) -= dot(K, out.get_block(seg_force_id, 3));

                Knormsq += dot(K,K);
                v1 += K;

                matrix diff(3,1);
                diff(0) = swimmers[n].filaments[m].segments[k].x[0] - swimmers[n].body.x[0];
                diff(1) = swimmers[n].filaments[m].segments[k].x[1] - swimmers[n].body.x[1];
                diff(2) = swimmers[n].filaments[m].segments[k].x[2] - swimmers[n].body.x[2];
                v2 += cross(diff, K);

              }

              // For the (K^T * M^(-1) * K)^(-1) part
              out(phase_id) /= Knormsq;
              out.subtract_from_block(out_id, 3, out(phase_id)*v1);
              out.subtract_from_block(out_id + 3, 3, out(phase_id)*v2);

            }

            // Apply inverse of (K^T * M^(-1) * K)
            out.set_block(out_id, 6, swimmers[n].KTMinvK_inv*out.get_block(out_id, 6));
            out.divide_block(3*NSWIM*(NBLOB + NFIL*NSEG + 2), NSWIM*NFIL, seg_mob_inv);

            for (int m = 0; m < NFIL; m++){

              const int phase_id = 3*NSWIM*(NBLOB + NFIL*NSEG + 2) + NFIL*n + m;

              Real Knormsq = 0.0;

              matrix v1(3,1), v2(3,1);
              v1.zero();
              v2.zero();

              for (int k = 0; k < NSEG; k++){

                const int seg_force_id = 3*(n*NFIL*NSEG + m*NSEG + k);

                matrix K(3,1);
                K(0) = swimmers[n].filaments[m].vel_dir_phase[3*k];
                K(1) = swimmers[n].filaments[m].vel_dir_phase[3*k + 1];
                K(2) = swimmers[n].filaments[m].vel_dir_phase[3*k + 2];

                Knormsq += dot(K,K);
                v1 += K;

                matrix diff(3,1);
                diff(0) = swimmers[n].filaments[m].segments[k].x[0] - swimmers[n].body.x[0];
                diff(1) = swimmers[n].filaments[m].segments[k].x[1] - swimmers[n].body.x[1];
                diff(2) = swimmers[n].filaments[m].segments[k].x[2] - swimmers[n].body.x[2];
                v2 += cross(diff, K);

              }

              out(phase_id) -= (dot(v1, out.get_block(out_id,3)) + dot(v2, out.get_block(out_id+3,3)))/Knormsq;

            }

          #else

            out.set_block(out_id, 6, swimmers[n].KTMinvK_inv*out.get_block(out_id, 6));

          #endif

          for (int m = 0; m < NFIL; m++){

            const std::vector<segment>& segments = swimmers[n].filaments[m].segments;

            #if DYNAMIC_PHASE_EVOLUTION

              const int phase_id = 3*NSWIM*(NBLOB + NFIL*NSEG + 2) + NFIL*n + m;

            #endif

            for (int k = 0; k < NSEG; k++){

              const int seg_force_id = 3*(n*NFIL*NSEG + m*NSEG + k);
              const Real diff[3] = {segments[k].x[0] - swimmers[n].body.x[0], segments[k].x[1] - swimmers[n].body.x[1], segments[k].x[2] - swimmers[n].body.x[2]};
              out(seg_force_id) += seg_mob_inv*(out(out_id) + out(out_id + 4)*diff[2] - out(out_id + 5)*diff[1]);
              out(seg_force_id + 1) += seg_mob_inv*(out(out_id + 1) + out(out_id + 5)*diff[0] - out(out_id + 3)*diff[2]);
              out(seg_force_id + 2) += seg_mob_inv*(out(out_id + 2) + out(out_id + 3)*diff[1] - out(out_id + 4)*diff[0]);

              #if DYNAMIC_PHASE_EVOLUTION

                // Additional K mult part
                out(seg_force_id) += seg_mob_inv*out(phase_id)*swimmers[n].filaments[m].vel_dir_phase[3*k];
                out(seg_force_id + 1) += seg_mob_inv*out(phase_id)*swimmers[n].filaments[m].vel_dir_phase[3*k + 1];
                out(seg_force_id + 2) += seg_mob_inv*out(phase_id)*swimmers[n].filaments[m].vel_dir_phase[3*k + 2];

              #endif

            }

          }

          for (int m = 0; m < NBLOB; m++){

            const int blob_force_id = 3*(NSWIM*NFIL*NSEG + n*NBLOB + m);
            const matrix diff = Q*matrix(3, 1, &swimmers[n].body.blob_references[3*m]);
            out(blob_force_id) += blob_mob_inv*(out(out_id) + out(out_id + 4)*diff(2) - out(out_id + 5)*diff(1));
            out(blob_force_id + 1) += blob_mob_inv*(out(out_id + 1) + out(out_id + 5)*diff(0) - out(out_id + 3)*diff(2));
            out(blob_force_id + 2) += blob_mob_inv*(out(out_id + 2) + out(out_id + 3)*diff(1) - out(out_id + 4)*diff(0));

          }

        }

      #endif

    #else

      const Real Minv = 6.0*PI*MU*RBLOB; // Diagonal approx.

      #if PRESCRIBED_BODY_VELOCITIES

        matrix out = Minv*in;

      #else

        matrix out(in.num_rows, 1);

        for (int n = 0; n < NSWIM; n++){

          const matrix Q = swimmers[n].body.q.rot_mat();
          const matrix Qinv = transpose(Q);

          matrix body_mobility = body_mobility_reference;
          body_mobility.set_block(0, 0, 3, 3, Q*body_mobility.get_block(0, 0, 3, 3)*Qinv);
          body_mobility.set_block(0, 3, 3, 3, Q*body_mobility.get_block(0, 3, 3, 3)*Qinv);
          body_mobility.set_block(3, 0, 3, 3, Q*body_mobility.get_block(3, 0, 3, 3)*Qinv);
          body_mobility.set_block(3, 3, 3, 3, Q*body_mobility.get_block(3, 3, 3, 3)*Qinv);

          matrix KTKinv = KTKinv_reference;
          KTKinv.set_block(0, 0, 3, 3, Q*KTKinv.get_block(0, 0, 3, 3)*Qinv);
          KTKinv.set_block(0, 3, 3, 3, Q*KTKinv.get_block(0, 3, 3, 3)*Qinv);
          KTKinv.set_block(3, 0, 3, 3, Q*KTKinv.get_block(3, 0, 3, 3)*Qinv);
          KTKinv.set_block(3, 3, 3, 3, Q*KTKinv.get_block(3, 3, 3, 3)*Qinv);

          matrix b(6,1);
          b.zero();

          for (int m = 0; m < NBLOB; m++){

            const int id = 3*(n*NBLOB + m);

            b(0) -= in(id);
            b(1) -= in(id + 1);
            b(2) -= in(id + 2);

            const matrix diff = Q*matrix(3, 1, &swimmers[n].body.blob_references[3*m]);

            b(3) -= diff(1)*in(id + 2) - diff(2)*in(id + 1);
            b(4) -= diff(2)*in(id) - diff(0)*in(id + 2);
            b(5) -= diff(0)*in(id + 1) - diff(1)*in(id);

          }

          const int vel_id = 3*NSWIM*NBLOB + 6*n;

          out.set_block(vel_id, 6, KTKinv*b - body_mobility*in.get_block(vel_id, 6));

          for (int m = 0; m < NBLOB; m++){

            const int id = 3*(n*NBLOB + m);
            const matrix diff = Q*matrix(3, 1, &swimmers[n].body.blob_references[3*m]);

            matrix block = in.get_block(id, 3) + out.get_block(vel_id, 3);
            block(0) += out(vel_id + 4)*diff(2) - out(vel_id + 5)*diff(1);
            block(1) += out(vel_id + 5)*diff(0) - out(vel_id + 3)*diff(2);
            block(2) += out(vel_id + 3)*diff(1) - out(vel_id + 4)*diff(0);

            out.set_block(id, 3, block);

          }

          out.multiply_block(0, 3*NSWIM*NBLOB, Minv);

        }

      #endif

    #endif

    return out;

  }

  matrix mobility_solver::system_matrix_mult(const matrix& in, const std::vector<swimmer>& swimmers){

    matrix out(rhs.num_rows, 1);

    #if PRESCRIBED_CILIA

      for (int n = 0; n < NSWIM*NFIL*NSEG; n++){

        // We don't need to be scared that we're writing here because there are no
        // conventional forces and torques on the segments when the motion is prescribed.
        f_segs_host[6*n] = in(3*n);
        f_segs_host[6*n + 1] = in(3*n + 1);
        f_segs_host[6*n + 2] = in(3*n + 2);

      }

      for (int n = 0; n < 3*NSWIM*NBLOB; n++){

        f_blobs_host[n] = in(n + 3*NSWIM*NFIL*NSEG);

      }

      // Get everything running on the GPUs that we can at this point.
      copy_segment_forces_to_device();
      evaluate_segment_segment_mobility();

    #if !INFINITE_PLANE_WALL

        copy_blob_forces_to_device();
        evaluate_segment_blob_mobility();
        evaluate_blob_blob_mobility();
        copy_blob_velocities_to_host();

      #endif

      copy_segment_velocities_to_host();

      // Do CPU work (i.e. everything involving the geometric matrices) while the GPU kernels run.
      out.zero();

      #if !PRESCRIBED_BODY_VELOCITIES

        for (int n = 0; n < NSWIM; n++){

          const int swim_id = 3*NSWIM*(NBLOB + NFIL*NSEG) + 6*n;

          for (int m = 0; m < NFIL*NSEG; m++){

            int out_id = 3*(n*NFIL*NSEG + m);

            // -K mult
            out(out_id) -= in(swim_id);
            out(out_id + 1) -= in(swim_id + 1);
            out(out_id + 2) -= in(swim_id + 2);

            const Real diff[3] = {x_segs_host[out_id] - swimmers[n].body.x[0], x_segs_host[out_id + 1] - swimmers[n].body.x[1], x_segs_host[out_id + 2] - swimmers[n].body.x[2]};
            out(out_id) -= in(swim_id + 4)*diff[2] - in(swim_id + 5)*diff[1];
            out(out_id + 1) -= in(swim_id + 5)*diff[0] - in(swim_id + 3)*diff[2];
            out(out_id + 2) -= in(swim_id + 3)*diff[1] - in(swim_id + 4)*diff[0];

            // -K^T mult
            const Real seg_force[3] = {f_segs_host[2*out_id], f_segs_host[2*out_id + 1], f_segs_host[2*out_id + 2]};
            out(swim_id) -= seg_force[0];
            out(swim_id + 1) -= seg_force[1];
            out(swim_id + 2) -= seg_force[2];
            out(swim_id + 3) -= diff[1]*seg_force[2] - diff[2]*seg_force[1];
            out(swim_id + 4) -= diff[2]*seg_force[0] - diff[0]*seg_force[2];
            out(swim_id + 5) -= diff[0]*seg_force[1] - diff[1]*seg_force[0];

          }

          const matrix R = swimmers[n].body.q.rot_mat();

          for (int m = 0; m < NBLOB; m++){

            const int out_id = 3*(n*NBLOB + NSWIM*NFIL*NSEG + m);

            // -K mult
            out(out_id) -= in(swim_id);
            out(out_id + 1) -= in(swim_id + 1);
            out(out_id + 2) -= in(swim_id + 2);

            const matrix diff = R*matrix(3, 1, &swimmers[n].body.blob_references[3*m]);
            out(out_id) -= in(swim_id + 4)*diff(2) - in(swim_id + 5)*diff(1);
            out(out_id + 1) -= in(swim_id + 5)*diff(0) - in(swim_id + 3)*diff(2);
            out(out_id + 2) -= in(swim_id + 3)*diff(1) - in(swim_id + 4)*diff(0);

            // -K^T mult
            const Real blob_force[3] = {f_blobs_host[3*(n*NBLOB + m)], f_blobs_host[3*(n*NBLOB + m) + 1], f_blobs_host[3*(n*NBLOB + m) + 2]};
            out(swim_id) -= blob_force[0];
            out(swim_id + 1) -= blob_force[1];
            out(swim_id + 2) -= blob_force[2];
            out(swim_id + 3) -= diff(1)*blob_force[2] - diff(2)*blob_force[1];
            out(swim_id + 4) -= diff(2)*blob_force[0] - diff(0)*blob_force[2];
            out(swim_id + 5) -= diff(0)*blob_force[1] - diff(1)*blob_force[0];

          }

        }

      #endif

      #if DYNAMIC_PHASE_EVOLUTION

        for (int n = 0; n < NSWIM; n++){

          for (int m = 0; m < NFIL; m++){

            const int out_pos = 3*NSEG*(n*NFIL + m);

            #if PRESCRIBED_BODY_VELOCITIES

              const int in_pos = 3*NSWIM*(NFIL*NSEG + NBLOB) + n*NFIL + m;

            #else

              const int in_pos = 3*NSWIM*(NFIL*NSEG + NBLOB + 2) + n*NFIL + m;

            #endif

            const Real phi_dot = in(in_pos);

            #if DYNAMIC_SHAPE_ROTATION

              const Real shape_rotation_angle_dot = in(in_pos + NSWIM*NFIL);

            #endif

            for (int k = 0; k < 3*NSEG; k++){

              // "-K" mult
              out(out_pos + k) -= phi_dot*swimmers[n].filaments[m].vel_dir_phase[k];

              // "-K^T" mult
              out(in_pos) -= in(out_pos + k)*swimmers[n].filaments[m].vel_dir_phase[k];

              #if DYNAMIC_SHAPE_ROTATION

                // "-K" mult
                out(out_pos + k) -= shape_rotation_angle_dot*swimmers[n].filaments[m].vel_dir_angle[k];

                // "-K^T" mult
                out(in_pos + NSWIM*NFIL) -= in(out_pos + k)*swimmers[n].filaments[m].vel_dir_angle[k];

              #endif

            }

          }

        }

      #elif DYNAMIC_SHAPE_ROTATION

        for (int n = 0; n < NSWIM; n++){

          for (int m = 0; m < NFIL; m++){

            const int out_pos = 3*NSEG*(n*NFIL + m);

            #if PRESCRIBED_BODY_VELOCITIES

              const int in_pos = 3*NSWIM*(NFIL*NSEG + NBLOB) + n*NFIL + m;

            #else

              const int in_pos = 3*NSWIM*(NFIL*NSEG + NBLOB + 2) + n*NFIL + m;

            #endif

            const Real shape_rotation_angle_dot = in(in_pos);

            for (int k = 0; k < 3*NSEG; k++){

              // "-K" mult
              out(out_pos + k) -= shape_rotation_angle_dot*swimmers[n].filaments[m].vel_dir_angle[k];

              // "-K^T" mult
              out(in_pos) -= in(out_pos + k)*swimmers[n].filaments[m].vel_dir_angle[k];

            }

          }

        }

      #endif

      // Now we need the results from the GPU. Get the last bit of GPU work running while we use them.
      wait_for_device();

      #if !INFINITE_PLANE_WALL

        // Do the bit using v_blobs here so we can make next set of kernel calls
        for (int n = 0; n < 3*NSWIM*NBLOB; n++){

          out(3*NSWIM*NFIL*NSEG + n) += v_blobs_host[n];

        }

        evaluate_blob_segment_mobility();
        copy_blob_velocities_to_host();

      #endif

      for (int n = 0; n < NSWIM*NFIL*NSEG; n++){

        out(3*n) += v_segs_host[6*n];
        out(3*n + 1) += v_segs_host[6*n + 1];
        out(3*n + 2) += v_segs_host[6*n + 2];

      }

      #if !INFINITE_PLANE_WALL

        // Retrieve and use the final GPU results.
        wait_for_device();

        for (int n = 0; n < 3*NSWIM*NBLOB; n++){

          out(3*NSWIM*NFIL*NSEG + n) += v_blobs_host[n];

        }

      #endif

    #else

      #if !INFINITE_PLANE_WALL

        for (int n = 0; n < 3*NSWIM*NBLOB; n++){

          f_blobs_host[n] = in(n);

        }

        copy_blob_forces_to_device();
        evaluate_blob_blob_mobility();
        copy_blob_velocities_to_host();

        out.zero();

        #if !PRESCRIBED_BODY_VELOCITIES

          for (int n = 0; n < NSWIM; n++){

            const matrix Q = swimmers[n].body.q.rot_mat();

            const int swim_id = 3*NSWIM*NBLOB + 6*n;

            for (int m = 0; m < NBLOB; m++){

              const matrix diff = Q*matrix(3, 1, &swimmers[n].body.blob_references[3*m]);

              const int blob_id = 3*(n*NBLOB + m);

              // -K^T mult
              out(swim_id) -= in(blob_id);
              out(swim_id + 1) -= in(blob_id + 1);
              out(swim_id + 2) -= in(blob_id + 2);
              out(swim_id + 3) -= diff(1)*in(blob_id + 2) - diff(2)*in(blob_id + 1);
              out(swim_id + 4) -= diff(2)*in(blob_id) - diff(0)*in(blob_id + 2);
              out(swim_id + 5) -= diff(0)*in(blob_id + 1) - diff(1)*in(blob_id);

              // -K mult
              out(blob_id) -= in(swim_id) - diff(1)*in(swim_id + 5) + diff(2)*in(swim_id + 4);
              out(blob_id + 1) -= in(swim_id + 1) - diff(2)*in(swim_id + 3) + diff(0)*in(swim_id + 5);
              out(blob_id + 2) -= in(swim_id + 2) - diff(0)*in(swim_id + 4) + diff(1)*in(swim_id + 3);

            }

          }

        #endif

        wait_for_device();

        for (int n = 0; n < 3*NSWIM*NBLOB; n++){

          out(n) += v_blobs_host[n];

        }

      #endif

    #endif

    return out;

  }

  int mobility_solver::solve_linear_system(std::vector<swimmer>& swimmers){

    /*

    In this function, we solve the system A * x = RHS, where A is the saddle-point matrix or a preconditioned version of it.

    We search iteratively for a solution x_n in the n-th Krylov subspace; i.e. the space spanned by RHS, A * RHS, ..., A^(n-1) * RHS.

    Since these spanning vectors may be close to linearly dependent, we use the stabilised Gram-Schmidt process to construct an orthonormal basis
    q_1, ..., q_n, which we store as the columns of the (system_size)-by-n matrix Q_n. Hence, there exists some vector y_n of coefficients such that
    x_n = Q_n * y_n.

    As we produce this orthonormal basis, we can also construct the (n+1)-by-n upper Hessenburg matrix H_n which satisfies A * Q_n = Q_(n+1) * H_n;
    the first n elements in column n are the projections of the corresponding Gram-Schmidt iterate onto q_n, and the (n+1)-th element is the norm of
    q_(n+1) before it is projected to unit length. The previous columns are unchanged -- we form H_n from H_(n-1) by adding a row of zeros and the
    column just described.

    Since the columns of Q_n are orthonormal,
      error = ||A * x_n - RHS|| = ||H_n * y_n - Q^T_(n+1) * RHS || = ||H_n * y_n - ||RHS|| * e_1 ||.

    Thus we can find x_n by solving a least-squares problem for y_n. Although H_n is not a square matrix, we can convert the problem to a square one
    as we proceed by performing successive Givens rotations to both the H matrices and the vector which starts as ||RHS|| * e_1. This process will
    also give us the error in the least-squares problem, and hence the error in the overall linear system.

    N.B. We actually solve one of the pre-conditioned systems A*P^(-1)*P*x = RHS or P^(-1)*A*x = P^(-1)*RHS, where P is a readily invertible approximation to A.

    */

    #if !USE_RIGHT_PRECON

      rhs = apply_preconditioner(rhs, swimmers);

    #endif

    const Real norm_of_rhs = norm(rhs);

    // If the right hand side is zero, then the solution is zero.
    if (norm_of_rhs == 0.0){

      #if PRESCRIBED_CILIA

        for (int n = 0; n < NSWIM*NFIL*NSEG; n++){

          f_segs_host[6*n] = 0.0;
          f_segs_host[6*n + 1] = 0.0;
          f_segs_host[6*n + 2] = 0.0;

        }

        for (int n = 0; n < 3*NSWIM*NBLOB; n++){

          f_blobs_host[n] = 0.0;

        }

        #if !PRESCRIBED_BODY_VELOCITIES

          v_bodies.zero();

        #endif

        #if DYNAMIC_PHASE_EVOLUTION

          for (int n = 0; n < NSWIM; n++){

            for (int m = 0; m < NFIL; m++){

              swimmers[n].filaments[m].phase_dot = 0.0;

            }

          }

        #endif

         #if DYNAMIC_SHAPE_ROTATION

          for (int n = 0; n < NSWIM; n++){

            for (int m = 0; m < NFIL; m++){

              swimmers[n].filaments[m].shape_rotation_angle_dot = 0.0;

            }

          }

        #endif

      #else

        for (int n = 0; n < 3*NSWIM*NBLOB; n++){

          f_blobs_host[n] = 0.0;

        }

        #if !PRESCRIBED_BODY_VELOCITIES

          v_bodies.zero();

        #endif

      #endif

      return 0;

    }

    // If 0 isn't the solution, try the previous solution as the initial guess. This should get progressively better as the time-step goes on.

    matrix soln(rhs.num_rows, 1);

    #if PRESCRIBED_CILIA

      for (int n = 0; n < NSWIM*NFIL*NSEG; n++){

        soln(3*n) = f_segs_host[6*n];
        soln(3*n + 1) = f_segs_host[6*n + 1];
        soln(3*n + 2) = f_segs_host[6*n + 2];

      }

      for (int n = 0; n < 3*NSWIM*NBLOB; n++){

        soln(n + 3*NSWIM*NFIL*NSEG) = f_blobs_host[n];

      }

      #if !PRESCRIBED_BODY_VELOCITIES

        soln.set_block(3*NSWIM*(NBLOB + NFIL*NSEG), 6*NSWIM, v_bodies);

      #endif

      #if DYNAMIC_PHASE_EVOLUTION

        for (int n = 0; n < NSWIM; n++){

          for (int m = 0; m < NFIL; m++){

            #if PRESCRIBED_BODY_VELOCITIES

              const int in_pos = 3*NSWIM*(NFIL*NSEG + NBLOB) + n*NFIL + m;

            #else

              const int in_pos = 3*NSWIM*(NFIL*NSEG + NBLOB + 2) + n*NFIL + m;

            #endif

            soln(in_pos) = swimmers[n].filaments[m].phase_dot;

            #if DYNAMIC_SHAPE_ROTATION

              soln(in_pos + NSWIM*NFIL) = swimmers[n].filaments[m].shape_rotation_angle_dot;

            #endif

          }

        }

      #elif DYNAMIC_SHAPE_ROTATION

        for (int n = 0; n < NSWIM; n++){

          for (int m = 0; m < NFIL; m++){

            #if PRESCRIBED_BODY_VELOCITIES

              const int in_pos = 3*NSWIM*(NFIL*NSEG + NBLOB) + n*NFIL + m;

            #else

              const int in_pos = 3*NSWIM*(NFIL*NSEG + NBLOB + 2) + n*NFIL + m;

            #endif

            soln(in_pos) = swimmers[n].filaments[m].shape_rotation_angle_dot;

          }

        }

      #endif

    #else

      for (int n = 0; n < 3*NSWIM*NBLOB; n++){

        soln(n) = f_blobs_host[n];

      }

      #if !PRESCRIBED_BODY_VELOCITIES

        soln.set_block(3*NSWIM*NBLOB, 6*NSWIM, v_bodies);

      #endif

    #endif

    #if USE_RIGHT_PRECON

      Q.set_col(0, rhs - system_matrix_mult(soln, swimmers));

    #else

      Q.set_col(0, rhs - apply_preconditioner(system_matrix_mult(soln, swimmers), swimmers));

    #endif

    beta(0) = norm(Q.get_col(0));

    if (std::abs(beta(0)) <= LINEAR_SYSTEM_TOL*norm_of_rhs){

      return 0; // Initial guess was good enough so do nothing. We don't even have to write it back into the relevant arrays because that's where we read it from in the first place.

    }

    if (std::abs(beta(0)) > norm_of_rhs){

      // Don't use the initial guess if it happens to be rubbish.
      soln.zero();
      beta(0) = norm_of_rhs;
      Q.set_col(0, rhs);

    }

    Q.divide_col(0, beta(0));

    const int max_iter = std::min<int>(rhs.num_rows, MAX_LINEAR_SYSTEM_ITER);

    for (int iter = 1; iter <= max_iter; iter++){

      // Produce the new orthonormal vector, using the appropriate values to update H as we do.
      #if USE_RIGHT_PRECON

        Q.set_col(iter, system_matrix_mult(apply_preconditioner(Q.get_col(iter-1), swimmers), swimmers));

      #else

        Q.set_col(iter, apply_preconditioner(system_matrix_mult(Q.get_col(iter-1), swimmers), swimmers));

      #endif

      for (int i = 0; i < iter; i++){

        H(i,iter-1) = dot(Q.get_col(iter), Q.get_col(i));

        Q.subtract_from_col(iter, H(i,iter-1)*Q.get_col(i));

      }

      const Real q_norm = norm(Q.get_col(iter));

      H(iter,iter-1) = q_norm;

      Q.divide_col(iter, q_norm);

      // Apply Givens rotations to transform the least-squares problem to a square one.
      // Apply any previous rotations IN ORDER. These have all already been applied to the RHS vector beta.
      for (int i = 1; i < iter; i++){

        const Real temp = CS(i-1)*H(i-1,iter-1) + SN(i-1)*H(i,iter-1);
        H(i,iter-1) = CS(i-1)*H(i,iter-1) - SN(i-1)*H(i-1,iter-1);
        H(i-1,iter-1) = temp;

      }

      // Then calculate and apply the new Givens rotation. This one is also applied to the RHS of the least-squares problem, beta.
      if (H(iter,iter-1)==0.0){

        // It's already lower diagonal.
        CS(iter-1) = 1.0;
        SN(iter-1) = 0.0;

      } else {

        if (std::abs(H(iter,iter-1)) > std::abs(H(iter-1,iter-1))){

          const Real temp = H(iter-1,iter-1)/H(iter,iter-1);
          SN(iter-1) = (2.0*Real(H(iter,iter-1) > 0) - 1.0)/sqrt(1.0 + temp*temp);
          CS(iter-1) = temp*SN(iter-1);

        } else {

          const Real temp = H(iter,iter-1)/H(iter-1,iter-1);
          CS(iter-1) = (2.0*Real(H(iter-1,iter-1) > 0) - 1.0)/sqrt(1.0 + temp*temp);
          SN(iter-1) = temp*CS(iter-1);

        }

      }

      H(iter-1,iter-1) = CS(iter-1)*H(iter-1,iter-1) + SN(iter-1)*H(iter,iter-1);
      H(iter,iter-1) = 0.0;

      beta(iter) = -SN(iter-1)*beta(iter-1);
      beta(iter-1) = CS(iter-1)*beta(iter-1);

      // beta(iter) now contains the (signed) error in the least-squares system, and hence the error in the linear system.
      // If it is small enough, or we've reached the maximum number of iterations, we generate the solution and return.
      const Real relative_error = std::abs(beta(iter))/norm_of_rhs;

      if ((relative_error <= LINEAR_SYSTEM_TOL) || (iter == max_iter)){

        matrix y(iter, 1);
        y.zero();

        y(iter-1) = beta(iter-1)/H(iter-1,iter-1);

        for (int i = iter-2; i >= 0; i--){

          y(i) = beta(i);

          for (int j = iter-1; j > i; j--){

            y(i) -= H(i,j)*y(j);

          }

          y(i) /= H(i,i);

        }

        matrix temp_soln(rhs.num_rows, 1);
        temp_soln.zero();

        for (int i = 0; i < iter; i++){

          temp_soln += y(i)*Q.get_col(i);

        }

        #if USE_RIGHT_PRECON

          soln += apply_preconditioner(temp_soln, swimmers);

          // When we use a right preconditioner, the velocity arrays are storing M*M_approx_inv*F
          // rather than M*F (and thus the force arrays store M_approx_inv*F rather than F). We need to
          // have them store the latter so that what we save in the segment velocity and force files is
          // meaningful, which we can achieve by simply calling system_matrix_mult(...) on the solution.
          // We don't want to do that here because then we'd be doing many more matrix multiplications
          // than are required, but rather as and when we want to save data. To do this, we save the
          // solution into the rhs vector, which is persistent storage across function calls.
          rhs = soln;

        #else

          soln += temp_soln;

        #endif

        // Store the solution back in the separate arrays and vectors
        #if PRESCRIBED_CILIA

          for (int n = 0; n < NSWIM*NFIL*NSEG; n++){

            f_segs_host[6*n] = soln(3*n);
            f_segs_host[6*n + 1] = soln(3*n + 1);
            f_segs_host[6*n + 2] = soln(3*n + 2);

          }

          for (int n = 0; n < 3*NSWIM*NBLOB; n++){

            f_blobs_host[n] = soln(n + 3*NSWIM*NFIL*NSEG);

          }

          #if !PRESCRIBED_BODY_VELOCITIES

            v_bodies = soln.get_block(3*NSWIM*(NBLOB + NFIL*NSEG), 6*NSWIM);

          #endif

          #if DYNAMIC_PHASE_EVOLUTION

            for (int n = 0; n < NSWIM; n++){

              for (int m = 0; m < NFIL; m++){

                #if PRESCRIBED_BODY_VELOCITIES

                  const int in_pos = 3*NSWIM*(NFIL*NSEG + NBLOB) + n*NFIL + m;

                #else

                  const int in_pos = 3*NSWIM*(NFIL*NSEG + NBLOB + 2) + n*NFIL + m;

                #endif

                swimmers[n].filaments[m].phase_dot = soln(in_pos);

                #if DYNAMIC_SHAPE_ROTATION

                  swimmers[n].filaments[m].shape_rotation_angle_dot = soln(in_pos + NSWIM*NFIL);

                #endif

              }

            }

          #elif DYNAMIC_SHAPE_ROTATION

            for (int n = 0; n < NSWIM; n++){

              for (int m = 0; m < NFIL; m++){

                #if PRESCRIBED_BODY_VELOCITIES

                  const int in_pos = 3*NSWIM*(NFIL*NSEG + NBLOB) + n*NFIL + m;

                #else

                  const int in_pos = 3*NSWIM*(NFIL*NSEG + NBLOB + 2) + n*NFIL + m;

                #endif

                swimmers[n].filaments[m].shape_rotation_angle_dot = soln(in_pos);

              }

            }

          #endif

        #else

          for (int n = 0; n < 3*NSWIM*NBLOB; n++){

            f_blobs_host[n] = soln(n);

          }

          #if !PRESCRIBED_BODY_VELOCITIES

            v_bodies = soln.get_block(3*NSWIM*NBLOB, 6*NSWIM);

          #endif

        #endif

        return iter;

      }

    }

  return 0;

  }

#endif // End of section to compile if not using Broyden's method for everything.

void mobility_solver::compute_velocities(std::vector<swimmer>& swimmers, int& num_gmres_iterations, const int nt){
  
  read_positions_and_forces(swimmers);  

  #if PRESCRIBED_CILIA

    assemble_rhs(swimmers, nt);

    num_gmres_iterations = solve_linear_system(swimmers);

  #else
    #if !ROD
      evaluate_segment_segment_mobility();
    #endif

    #if !INFINITE_PLANE_WALL

      #if USE_BROYDEN_FOR_EVERYTHING

          evaluate_blob_blob_mobility();

      #endif

      #if !ROD
        evaluate_blob_segment_mobility();
      #endif
      copy_blob_velocities_to_host();
      
      #if !USE_BROYDEN_FOR_EVERYTHING

        assemble_rhs(swimmers, nt);

        num_gmres_iterations = solve_linear_system(swimmers);

      #endif

      #if !ROD
        copy_blob_forces_to_device();
        evaluate_segment_blob_mobility();
      #endif

    #endif

    copy_segment_velocities_to_host();

    wait_for_device();

  #endif

}

bool mobility_solver::compute_errors(matrix& error, const std::vector<swimmer>& swimmers, const int nt){

  // Storage order is one body at a time, and within each body it is all filaments first, then
  // the body state/force- and torque-balance and then finally the blob velocities (if applicable).

  bool error_is_too_large = false;

  const Real c1 = -4.0/3.0;
  const Real c2 = 1.0/3.0;
  const Real c3 = -2.0/3.0;

  for (int n = 0; n < NSWIM; n++){

    #if USE_BROYDEN_FOR_EVERYTHING

      const int body_id = n*(6*NFIL*NSEG + 3*NBLOB + 6) + 6*NFIL*NSEG;

      #if !PRESCRIBED_BODY_VELOCITIES

        error.set_block(body_id, 6, -swimmers[n].f);

      #endif

    #endif

    const std::vector<filament>& filaments = swimmers[n].filaments;

    for (int i = 0; i < NFIL; i++){

      const std::vector<segment>& segments = filaments[i].segments;

      for (int j = 0; j < NSEG; j++){

        int id = 6*(j + NSEG*(i + n*NFIL));

        const Real *const v = &v_segs_host[id];
        const Real *const omega = &v_segs_host[id + 3];

        #if INFINITE_PLANE_WALL

          id = 6*i*NSEG + 3*j;

        #elif USE_BROYDEN_FOR_EVERYTHING

          id = n*(6*NFIL*NSEG + 3*NBLOB + 6) + 6*i*NSEG + 3*j;

        #else

          id = n*6*(NFIL*NSEG + 1) + 6*i*NSEG + 3*j;

        #endif

        if (nt < NUM_EULER_STEPS){

          error(id) = segments[j].x[0] - segments[j].xm1[0] - DT*v[0];
          error(id + 1) = segments[j].x[1] - segments[j].xm1[1] - DT*v[1];
          error(id + 2) = segments[j].x[2] - segments[j].xm1[2] - DT*v[2];

        } else {

          error(id) = segments[j].x[0] + c1*segments[j].xm1[0] + c2*segments[j].xm2[0]+ c3*DT*v[0];
          error(id + 1) = segments[j].x[1] + c1*segments[j].xm1[1] + c2*segments[j].xm2[1]+ c3*DT*v[1];
          error(id + 2) = segments[j].x[2] + c1*segments[j].xm1[2] + c2*segments[j].xm2[2]+ c3*DT*v[2];

        }

        if (!error_is_too_large){

          // Scale TOL by DL for position equations
          error_is_too_large = (std::abs(error(id)) > DL*TOL) || (std::abs(error(id + 1)) > DL*TOL) || (std::abs(error(id + 2)) > DL*TOL);

        }

        id += 3*NSEG;

        Real W[3];
        dexpinv(W, segments[j].u, omega);

        if (nt < NUM_EULER_STEPS){

          error(id) = segments[j].u[0] - DT*W[0];
          error(id + 1) = segments[j].u[1] - DT*W[1];
          error(id + 2) = segments[j].u[2] - DT*W[2];

        } else {

          error(id) = segments[j].u[0]- c2*segments[j].um1[0] + c3*DT*W[0];
          error(id + 1) = segments[j].u[1]- c2*segments[j].um1[1] + c3*DT*W[1];
          error(id + 2) = segments[j].u[2]- c2*segments[j].um1[2] + c3*DT*W[2];

        }

        if (!error_is_too_large){

          error_is_too_large = (std::abs(error(id)) > TOL) || (std::abs(error(id + 1)) > TOL) || (std::abs(error(id + 2)) > TOL);

        }

      }

    }

    #if !INFINITE_PLANE_WALL

    // We still let Broyden's method find the updated body state when it's velocities are prescribed.
    // This is obviously overkill for the positions because we can re-arrange both backwards-Euler and BDF2
    // to get the new position explicitly in terms of the previous ones and the imposed velocity, but this
    // is not true for the Lie algebra element (technically we can solve the Euler steps but I don't really
    // want the number of equations to change over time...). It might be worth having the positions handled
    // differently and effectively solve the position update equations to machine precision? As it stands,
    // I only ever use this to hold the body fixed, meaning that the error is always zero, Broyden's method
    // never updates the position or orientation and they're both satisified exactly anyway.

      #if USE_BROYDEN_FOR_EVERYTHING

        Real v[3], omega[3];

        #if PRESCRIBED_BODY_VELOCITIES

          swimmers[n].body.prescribed_translational_velocity(v, DT*(nt + 1.0), n);
          swimmers[n].body.prescribed_rotational_velocity(omega, DT*(nt + 1.0), n);

        #else

          if (nt < NUM_EULER_STEPS){

            v[0] = (swimmers[n].body.x[0] - swimmers[n].body.xm1[0])/DT;
            v[1] = (swimmers[n].body.x[1] - swimmers[n].body.xm1[1])/DT;
            v[2] = (swimmers[n].body.x[2] - swimmers[n].body.xm1[2])/DT;

            omega[0] = swimmers[n].body.u[0]/DT; // = dexp(swimmers[n].u, swimmers[n].u)/DT
            omega[1] = swimmers[n].body.u[1]/DT;
            omega[2] = swimmers[n].body.u[2]/DT;

          } else {

            v[0] = -(swimmers[n].body.x[0] + c1*swimmers[n].body.xm1[0] + c2*swimmers[n].body.xm2[0])/(c3*DT);
            v[1] = -(swimmers[n].body.x[1] + c1*swimmers[n].body.xm1[1] + c2*swimmers[n].body.xm2[1])/(c3*DT);
            v[2] = -(swimmers[n].body.x[2] + c1*swimmers[n].body.xm1[2] + c2*swimmers[n].body.xm2[2])/(c3*DT);

            const Real w[3] = {(c2*swimmers[n].body.um1[0] - swimmers[n].body.u[0])/(c3*DT),
                                  (c2*swimmers[n].body.um1[1] - swimmers[n].body.u[1])/(c3*DT),
                                  (c2*swimmers[n].body.um1[2] - swimmers[n].body.u[2])/(c3*DT)};

            dexp(omega, swimmers[n].body.u, w);

          }

        #endif

        const matrix Q = swimmers[n].body.q.rot_mat();

        // Blob equations
        for (int m = 0; m < NBLOB; m++){

          const int id = n*(6*NFIL*NSEG + 3*NBLOB + 6) + 6*(NFIL*NSEG + 1) + 3*m;

          const matrix diff = Q*matrix(3, 1, &swimmers[n].body.blob_references[3*m]);

          matrix w = cross(omega, diff);

          error(id) = v[0] + w(0) - v_blobs_host[3*(n*NBLOB + m)];
          error(id + 1) = v[1] + w(1) - v_blobs_host[3*(n*NBLOB + m) + 1];
          error(id + 2) = v[2] + w(2) - v_blobs_host[3*(n*NBLOB + m) + 2];

          if (!error_is_too_large){

            error_is_too_large = (std::abs(error(id)) > TOL) || (std::abs(error(id + 1)) > TOL) || (std::abs(error(id + 2)) > TOL);

          }

          #if !PRESCRIBED_BODY_VELOCITIES

            error(body_id) += f_blobs_host[3*(n*NBLOB + m)];
            error(body_id + 1) += f_blobs_host[3*(n*NBLOB + m) + 1];
            error(body_id + 2) += f_blobs_host[3*(n*NBLOB + m) + 2];

            w = cross(diff, &f_blobs_host[3*(n*NBLOB + m)]);

            error(body_id + 3) += w(0);
            error(body_id + 4) += w(1);
            error(body_id + 5) += w(2);

          #endif

        }

        #if PRESCRIBED_BODY_VELOCITIES

          Real w[3];
          dexpinv(w, swimmers[n].body.u, omega);

          if (nt < NUM_EULER_STEPS){

            error(body_id) = swimmers[n].body.x[0] - swimmers[n].body.xm1[0] - DT*v[0];
            error(body_id + 1) = swimmers[n].body.x[1] - swimmers[n].body.xm1[1] - DT*v[1];
            error(body_id + 2) = swimmers[n].body.x[2] - swimmers[n].body.xm1[2] - DT*v[2];
            error(body_id + 3) = swimmers[n].body.u[0] - DT*w[0];
            error(body_id + 4) = swimmers[n].body.u[1] - DT*w[1];
            error(body_id + 5) = swimmers[n].body.u[2] - DT*w[2];

          } else {

            error(body_id) = swimmers[n].body.x[0] + c1*swimmers[n].body.xm1[0] + c2*swimmers[n].body.xm2[0] + c3*DT*v[0];
            error(body_id + 1) = swimmers[n].body.x[1] + c1*swimmers[n].body.xm1[1] + c2*swimmers[n].body.xm2[1] + c3*DT*v[1];
            error(body_id + 2) = swimmers[n].body.x[2] + c1*swimmers[n].body.xm1[2] + c2*swimmers[n].body.xm2[2] + c3*DT*v[2];
            error(body_id + 3) = swimmers[n].body.u[0] - c2*swimmers[n].body.um1[0] + c3*DT*w[0];
            error(body_id + 4) = swimmers[n].body.u[1] - c2*swimmers[n].body.um1[1] + c3*DT*w[1];
            error(body_id + 5) = swimmers[n].body.u[2] - c2*swimmers[n].body.um1[2] + c3*DT*w[2];

          }

          if (!error_is_too_large){

            error_is_too_large = (std::abs(error(body_id)) > TOL*DL) || (std::abs(error(body_id + 1)) > TOL*DL) || (std::abs(error(body_id + 2)) > TOL*DL);

          }

          if (!error_is_too_large){

            error_is_too_large = (std::abs(error(body_id + 3)) > TOL) || (std::abs(error(body_id + 4)) > TOL) || (std::abs(error(body_id + 5)) > TOL);

          }

        #else

          // Force and torque balances have been formed as we've gone along. We just need to check the error.
          if (!error_is_too_large){

            error_is_too_large = (std::abs(error(body_id)) > TOL) || (std::abs(error(body_id + 1)) > TOL) ||
                                  (std::abs(error(body_id + 2)) > TOL) || (std::abs(error(body_id + 3)) > TOL) ||
                                  (std::abs(error(body_id + 4)) > TOL) || (std::abs(error(body_id + 5)) > TOL);

          }

        #endif

      #else

        const int body_id = 6*n;
        const int id = n*6*(NFIL*NSEG + 1) + 6*NFIL*NSEG;

        Real w[3];
        dexpinv(w, swimmers[n].body.u, &v_bodies.data[body_id+3]);

        if (nt < NUM_EULER_STEPS){

          error(id) = swimmers[n].body.x[0] - swimmers[n].body.xm1[0] - DT*v_bodies(body_id);
          error(id + 1) = swimmers[n].body.x[1] - swimmers[n].body.xm1[1] - DT*v_bodies(body_id + 1);
          error(id + 2) = swimmers[n].body.x[2] - swimmers[n].body.xm1[2] - DT*v_bodies(body_id + 2);
          error(id + 3) = swimmers[n].body.u[0] - DT*w[0];
          error(id + 4) = swimmers[n].body.u[1] - DT*w[1];
          error(id + 5) = swimmers[n].body.u[2] - DT*w[2];

        } else {

          error(id) = swimmers[n].body.x[0] + c1*swimmers[n].body.xm1[0] + c2*swimmers[n].body.xm2[0] + c3*DT*v_bodies(body_id);
          error(id + 1) = swimmers[n].body.x[1] + c1*swimmers[n].body.xm1[1] + c2*swimmers[n].body.xm2[1] + c3*DT*v_bodies(body_id + 1);
          error(id + 2) = swimmers[n].body.x[2] + c1*swimmers[n].body.xm1[2] + c2*swimmers[n].body.xm2[2] + c3*DT*v_bodies(body_id + 2);
          error(id + 3) = swimmers[n].body.u[0] - c2*swimmers[n].body.um1[0] + c3*DT*w[0];
          error(id + 4) = swimmers[n].body.u[1] - c2*swimmers[n].body.um1[1] + c3*DT*w[1];
          error(id + 5) = swimmers[n].body.u[2] - c2*swimmers[n].body.um1[2] + c3*DT*w[2];

        }

        if (!error_is_too_large){

          error_is_too_large = (std::abs(error(id)) > TOL*DL) || (std::abs(error(id + 1)) > TOL*DL) || (std::abs(error(id + 2)) > TOL*DL);

        }

        if (!error_is_too_large){

          error_is_too_large = (std::abs(error(id + 3)) > TOL) || (std::abs(error(id + 4)) > TOL) || (std::abs(error(id + 5)) > TOL);

        }

      #endif

    #endif

  }

  return error_is_too_large;

}

#if !USE_BROYDEN_FOR_EVERYTHING

  void mobility_solver::make_body_reference_matrices(std::vector<swimmer>& swimmers){

    // Reference (K^T * K)^(-1)
    // Do this one first so that it can be used in the preconditioner for the next part.
    matrix KTK(6,6);
    KTK.zero();

    KTK(0,0) = NBLOB;
    KTK(1,1) = NBLOB;
    KTK(2,2) = NBLOB;

    for (int m = 0; m < NBLOB; m++){

      const matrix rcross_mat = rcross(&swimmers[0].body.blob_references[3*m]);

      KTK.add_to_block(0, 3, 3, 3, rcross_mat);
      KTK.subtract_from_block(3, 0, 3, 3, rcross_mat);
      KTK.subtract_from_block(3, 3, 3, 3, rcross_mat*rcross_mat);

    }

    KTKinv_reference = inverse(KTK);

    // Reference mobility matrix
    Real* original_posns = new Real[3*(NSWIM-1)];

    for (int n = 0; n < NSWIM-1; n++){

      original_posns[3*n] = swimmers[n].body.x[0];
      original_posns[3*n + 1] = swimmers[n].body.x[1];
      original_posns[3*n + 2] = swimmers[n].body.x[2];

      swimmers[n].body.x[0] += 1e10*DL*NSEG;
      swimmers[n].body.x[1] += 1e10*DL*NSEG;
      swimmers[n].body.x[2] += 1e10*DL*NSEG;

    }

    read_positions_and_forces(swimmers);

    #if SURFACE_OF_REVOLUTION_BODIES

      std::vector<std::string> names = {"WZ", "WY", "WX", "VZ", "VY", "VX"};

    #endif

    for (int n = 0; n < 6; n++){

      rhs.zero();
      rhs(rhs.num_rows - 1 - n) = -1.0;

      const int temp = solve_linear_system(swimmers); // Not interested in the number of required iterations here.

      body_mobility_reference.set_col(5-n, v_bodies.get_block(6*(NSWIM-1), 6));

      #if SURFACE_OF_REVOLUTION_BODIES

        std::ofstream file(std::string(GENERATRIX_FILE_NAME)+std::to_string(NBLOB)+"_"+names[n]+"_blob_force_distribution.dat");
        for (int ii = 0; ii < 3*NBLOB; ii++){
          file << rhs(ii) << " ";
        }
        file.close();

      #endif

    }

    #if SURFACE_OF_REVOLUTION_BODIES

      if (std::string(GENERATRIX_FILE_NAME) == std::string("sphere")){

        const Real R = 0.5*AXIS_DIR_BODY_LENGTH;
        matrix Ntrue(6, 6);
        Ntrue.zero();
        Ntrue(0,0) = 1.0/(6.0*PI*MU*R);
        Ntrue(1,1) = Ntrue(0,0);
        Ntrue(2,2) = Ntrue(0,0);
        Ntrue(3,3) = 1.0/(8.0*PI*MU*R*R*R);
        Ntrue(4,4) = Ntrue(3,3);
        Ntrue(5,5) = Ntrue(3,3);
        std::cout << std::endl << "Percentage error in mobility matrix compared to Stokes drag is " << 100.0*norm(body_mobility_reference - Ntrue)/norm(Ntrue) << "%" << std::endl;

      }

    #endif

    for (int n = 0; n < NSWIM-1; n++){

      swimmers[n].body.x[0] = original_posns[3*n];
      swimmers[n].body.x[1] = original_posns[3*n + 1];
      swimmers[n].body.x[2] = original_posns[3*n + 2];

    }

    delete[] original_posns;

  }

#endif

void mobility_solver::write_data(const int nt, const std::vector<swimmer>& swimmers){

  #if (!USE_BROYDEN_FOR_EVERYTHING && USE_RIGHT_PRECON && PRESCRIBED_CILIA)

    // See solve_linear_system(...)
    const matrix junk = system_matrix_mult(rhs, swimmers);

    // N.B. The use of rhs here means write_data(...) cannot be a const method.

  #endif

  #if !INFINITE_PLANE_WALL

    std::ofstream body_vel_file(SIMULATION_BODY_VEL_NAME, std::ios::app);
    body_vel_file << nt << " ";
    body_vel_file << std::scientific << std::setprecision(10);

    std::ofstream blob_forces_file(SIMULATION_BLOB_FORCES_NAME, std::ios::app);
    blob_forces_file << nt << " ";
    blob_forces_file << std::scientific << std::setprecision(10);

  #endif

  std::ofstream seg_vel_file(SIMULATION_SEG_VEL_NAME, std::ios::app);
  seg_vel_file << nt << " ";
  seg_vel_file << std::scientific << std::setprecision(10);

  std::ofstream seg_forces_file(SIMULATION_SEG_FORCES_NAME, std::ios::app);
  seg_forces_file << nt << " ";
  seg_forces_file << std::scientific << std::setprecision(10);

  for (int n = 0; n < NSWIM; n++){

    #if !INFINITE_PLANE_WALL

      #if USE_BROYDEN_FOR_EVERYTHING

        Real v[3], omega[3];

        #if PRESCRIBED_BODY_VELOCITIES

          swimmers[n].body.prescribed_translational_velocity(v, DT*(nt + 1.0), n);
          swimmers[n].body.prescribed_rotational_velocity(omega, DT*(nt + 1.0), n);

        #else

          if (nt < NUM_EULER_STEPS){

            v[0] = (swimmers[n].body.x[0] - swimmers[n].body.xm1[0])/DT;
            v[1] = (swimmers[n].body.x[1] - swimmers[n].body.xm1[1])/DT;
            v[2] = (swimmers[n].body.x[2] - swimmers[n].body.xm1[2])/DT;

            omega[0] = swimmers[n].body.u[0]/DT; // = dexp(swimmers[n].u, swimmers[n].u)/DT
            omega[1] = swimmers[n].body.u[1]/DT;
            omega[2] = swimmers[n].body.u[2]/DT;

          } else {

            v[0] = 0.5*(3.0*swimmers[n].body.x[0] - 4.0*swimmers[n].body.xm1[0] + swimmers[n].body.xm2[0])/DT;
            v[1] = 0.5*(3.0*swimmers[n].body.x[1] - 4.0*swimmers[n].body.xm1[1] + swimmers[n].body.xm2[1])/DT;
            v[2] = 0.5*(3.0*swimmers[n].body.x[2] - 4.0*swimmers[n].body.xm1[2] + swimmers[n].body.xm2[2])/DT;

            const Real w[3] = {(1.5*swimmers[n].body.u[0] - 0.5*swimmers[n].body.um1[0])/DT,
                                  (1.5*swimmers[n].body.u[1] - 0.5*swimmers[n].body.um1[1])/DT,
                                   (1.5*swimmers[n].body.u[2] - 0.5*swimmers[n].body.um1[2])/DT};

            dexp(omega, swimmers[n].body.u, w);

          }

        #endif

        body_vel_file << v[0] << " " << v[1] << " " << v[2] << " " << omega[0] << " " << omega[1] << " " << omega[2] << " ";

      #else

        body_vel_file << v_bodies(6*n) << " " << v_bodies(6*n + 1) << " " << v_bodies(6*n + 2) << " " << v_bodies(6*n + 3) << " " << v_bodies(6*n + 4) << " " << v_bodies(6*n + 5) << " ";

      #endif

      for (int i = 0; i < NBLOB; i++){

        const int id = 3*(i + n*NBLOB);

        blob_forces_file << f_blobs_host[id] << " " << f_blobs_host[id + 1] << " " << f_blobs_host[id + 2] << " ";

      }

    #endif

    for (int i = 0; i < NFIL; i++){
      for (int j = 0; j < NSEG; j++){

        const int id = 6*(j + NSEG*(i + n*NFIL));

        seg_vel_file << v_segs_host[id] << " " << v_segs_host[id + 1] << " " << v_segs_host[id + 2] << " " << v_segs_host[id + 3] << " " << v_segs_host[id + 4] << " " << v_segs_host[id + 5] << " ";

        seg_forces_file << f_segs_host[id] << " " << f_segs_host[id + 1] << " " << f_segs_host[id + 2] << " " << f_segs_host[id + 3] << " " << f_segs_host[id + 4] << " " << f_segs_host[id + 5] << " ";

      }
    }

  }

  #if !INFINITE_PLANE_WALL

    body_vel_file << std::endl;
    body_vel_file.close();

    blob_forces_file << std::endl;
    blob_forces_file.close();

  #endif

  seg_vel_file << std::endl;
  seg_vel_file.close();

  seg_forces_file << std::endl;
  seg_forces_file.close();

}