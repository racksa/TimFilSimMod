// rigid_body.cpp

#include <cmath>
#include <string>
#include <algorithm>
#include "rigid_body.hpp"
#include "seeding.hpp"
#include "../../config.hpp"
#include "../general/matrix.hpp"
#include "../general/quaternion.hpp"

rigid_body::~rigid_body(){}

rigid_body::rigid_body(){}

void rigid_body::initial_setup(const int id, Real *const f_address, const Real *const data_from_file){

  blob_references = std::vector<Real>(3*NBLOB);
  polar_dir_refs = std::vector<Real>(3*NBLOB);
  azi_dir_refs = std::vector<Real>(3*NBLOB);
  normal_refs = std::vector<Real>(3*NBLOB);

  Q_init = matrix(3,3);

  #if (USE_BROYDEN_FOR_EVERYTHING && !INFINITE_PLANE_WALL)

    blob_forces = f_address;

    blob_forces_m1 = std::vector<Real>(3*NBLOB);
    blob_forces_m2 = std::vector<Real>(3*NBLOB);

  #endif

  #if READ_INITIAL_CONDITIONS_FROM_BACKUP

    int p = 0;

    q.scalar_part = data_from_file[p++];
    q.vector_part[0] = data_from_file[p++];
    q.vector_part[1] = data_from_file[p++];
    q.vector_part[2] = data_from_file[p++];

    qm1 = q;

    x[0] = data_from_file[p++];
    x[1] = data_from_file[p++];
    x[2] = data_from_file[p++];
    xm1[0] = data_from_file[p++];
    xm1[1] = data_from_file[p++];
    xm1[2] = data_from_file[p++];
    xm2[0] = data_from_file[p++];
    xm2[1] = data_from_file[p++];
    xm2[2] = data_from_file[p++];

    u[0] = data_from_file[p++];
    u[1] = data_from_file[p++];
    u[2] = data_from_file[p++];
    um1[0] = data_from_file[p++];
    um1[1] = data_from_file[p++];
    um1[2] = data_from_file[p++];

    #if USE_BROYDEN_FOR_EVERYTHING

      // Do this last so we can start GMRES sims from broyden-only backups.
      for (int n = 0; n < 3*NBLOB; n++){

        blob_forces[n] = data_from_file[p++];
        blob_forces_m1[n] = data_from_file[p++];
        blob_forces_m2[n] = data_from_file[p++];

      }

    #endif

  #endif

  #if INFINITE_PLANE_WALL

    #if !READ_INITIAL_CONDITIONS_FROM_BACKUP

      x[0] = 0.0;
      x[1] = 0.0;
      x[2] = 0.0;
      xm1[0] = x[0];
      xm1[1] = x[1];
      xm1[2] = x[2];
      xm2[0] = x[0];
      xm2[1] = x[1];
      xm2[2] = x[2];

      u[0] = 0.0;
      u[1] = 0.0;
      u[2] = 0.0;
      um1[0] = 0.0;
      um1[1] = 0.0;
      um1[2] = 0.0;

      q = quaternion(1.0, 0.0, 0.0, 0.0);
      qm1 = q;

    #endif

  #elif SADDLE_BODIES

    // Setting both to INFINITY will form a plane.
    const Real radius_of_curvature_x = INFINITY;
    const Real radius_of_curvature_y = INFINITY;

    const Real grid_length_x = DL*NSEG;
    const Real grid_length_y = DL*NSEG;
    const int surf_grid_dim_x = round(sqrt(grid_length_x*NBLOB/grid_length_y));
    const int surf_grid_dim_y = ceil(NBLOB/Real(surf_grid_dim_x));
    const Real surf_grid_step_x = (surf_grid_dim_x > 1) ? grid_length_x/(Real(surf_grid_dim_x) - 1.0) : 0.0;
    const Real surf_grid_step_y = (surf_grid_dim_y > 1) ? grid_length_y/(Real(surf_grid_dim_y) - 1.0) : 0.0;

    #if !READ_INITIAL_CONDITIONS_FROM_BACKUP

      const Real body_spacing = std::max<Real>(grid_length_x, grid_length_y) + 3.0*DL*NSEG;

      x[0] = id*body_spacing;
      x[1] = 0.0;
      x[2] = 0.0;
      xm1[0] = x[0];
      xm1[1] = x[1];
      xm1[2] = x[2];
      xm2[0] = x[0];
      xm2[1] = x[1];
      xm2[2] = x[2];

      u[0] = 0.0;
      u[1] = 0.0;
      u[2] = 0.0;
      um1[0] = 0.0;
      um1[1] = 0.0;
      um1[2] = 0.0;

      q = quaternion(1.0, 0.0, 0.0, 0.0);
      qm1 = q;

    #endif

    const Real im = 0.5*(surf_grid_dim_x - 1.0);
    const Real jm = 0.5*(surf_grid_dim_y - 1.0);

    const Real theta_x = surf_grid_step_x/radius_of_curvature_x;
    const Real theta_y = surf_grid_step_y/radius_of_curvature_y;

    // The multiblob method won't work if the point we designate as the position of the body
    // coincides with one of the blob positions, so we displace it from the body.
    const Real height_offset = RBLOB;

    for (int i = 0; i < surf_grid_dim_x; i++){
      for (int j = 0; j < surf_grid_dim_y; j++){

        const int id = 3*(j + i*surf_grid_dim_y);

        if (blob_id < NBLOB){

          if (std::isinf(radius_of_curvature_x) && std::isinf(radius_of_curvature_y)){

            blob_references[id] = (i-im)*surf_grid_step_x;
            blob_references[id + 1] = (j-jm)*surf_grid_step_y;
            blob_references[id + 2] = height_offset;

          } else if (std::isinf(radius_of_curvature_x)) {

            blob_references[id] = (i-im)*surf_grid_step_x;
            blob_references[id + 1] = radius_of_curvature_y*sin(theta_y*(j - jm));
            blob_references[id + 2] = height_offset + radius_of_curvature_y*(1.0 - cos(theta_y*(j - jm)));

          } else if (std::isinf(radius_of_curvature_y)) {

            blob_references[id] = radius_of_curvature_x*sin(theta_x*(i - im));
            blob_references[id + 1] = (j-jm)*surf_grid_step_y;
            blob_references[id + 2] = height_offset + radius_of_curvature_x*(1.0 - cos(theta_x*(i - im)));

          } else {

            blob_references[id] = radius_of_curvature_x*sin(theta_x*(i - im));
            blob_references[id + 1] = radius_of_curvature_y*sin(theta_y*(j - jm));
            blob_references[id + 2] = height_offset + radius_of_curvature_x*(1.0 - cos(theta_x*(i - im))) + radius_of_curvature_y*(1.0 - cos(theta_y*(j - jm)));

          }

        }

      }
    }

  #elif SURFACE_OF_REVOLUTION_BODIES or ROD or RIGIDWALL

    #if !READ_INITIAL_CONDITIONS_FROM_BACKUP

      #if SURFACE_OF_REVOLUTION_BODIES
        const Real body_spacing = 3.0*(AXIS_DIR_BODY_LENGTH + FIL_LENGTH);

        Real sep;
        std::ifstream in("separation.dat"); // input
        in >> sep;

        x[0] = id*sep;
        x[1] = 0.0;
        x[2] = 0.0;
        xm1[0] = x[0];
        xm1[1] = x[1];
        xm1[2] = x[2];
        xm2[0] = x[0];
        xm2[1] = x[1];
        xm2[2] = x[2];

        u[0] = 0.0;
        u[1] = 0.0;
        u[2] = 0.0;
        um1[0] = 0.0;
        um1[1] = 0.0;
        um1[2] = 0.0;

        q = quaternion(1.0, 0.0, 0.0, 0.0);
        qm1 = q;
        
      #elif ROD
        Real body_spacing = 1.8*(AXIS_DIR_BODY_LENGTH);

        int relative_Lx = 8;
        int relative_Ly = 1;
        int relative_Lz = 1;
        int b = relative_Lx/relative_Ly;
        int c = relative_Lx/relative_Lz;

        // initialise plane 
        // int Nx = ceil(sqrt(NSWIM));
        // const int k = 0;
        // const int j = id/Nx;
        // const int i = id - j*Nx;

        // initialise lattice        
        int Nx = ceil(cbrt(NSWIM*b*c));
        int Ny = Nx/b;
        int k = id/(Nx*Ny);
        int j = (id - k*Nx*Ny)/Nx;
        int i = id - k*Nx*Ny - j*Nx;
        
        std::ifstream in("separation.dat"); // input
        in >> body_spacing;

        x[0] = i*body_spacing + body_spacing;
        x[1] = j*body_spacing + body_spacing;
        x[2] = k*body_spacing + body_spacing;

        // if(id==0){
        //   x[0] = body_spacing;
        //   x[1] = body_spacing;
        // }
        // if(id==1){
        //   x[0] = -body_spacing;
        //   x[1] = body_spacing;
        // }
        // if(id==2){
        //   x[0] = body_spacing;
        //   x[1] = -body_spacing;
        // }
        // if(id==3){
        //   x[0] = -body_spacing;
        //   x[1] = -body_spacing;
        // }

        xm1[0] = x[0];
        xm1[1] = x[1];
        xm1[2] = x[2];
        xm2[0] = x[0];
        xm2[1] = x[1];
        xm2[2] = x[2];

        u[0] = 0.0;
        u[1] = 0.0;
        u[2] = 0.0;
        um1[0] = 0.0;
        um1[1] = 0.0;
        um1[2] = 0.0;

        // Randomise orientation in 2d
        // q = quaternion(1.0, 0.0, 0.0, 0.0);
        // std::random_device rd{};
        // std::mt19937 gen = std::mt19937{rd()};
        // std::normal_distribution<Real> d(0.0, 1.0);
        // q(3) = d(gen);
        // q.normalise_in_place();
        
        // Randomise orientation in 3d
        q.randomise();

        // if(id==0){
        //   q = quaternion(1.0, 0.0, 1.0, 0.0);
        //   q.normalise_in_place();
        // }
        // if(id==1){
        //   q = quaternion(1.0, 0.0, 0.0, 1.0);
        //   q.normalise_in_place();
        // }
        // if(id==2){
        //   q = quaternion(1.0, 0.0, 0.0, 1.0);
        //   q.normalise_in_place();
        // }
        // if(id==3){
        //   q = quaternion(1.0, 0.0, 1.0, 0.0);
        //   q.normalise_in_place();
        // }

        qm1 = q;
      #elif RIGIDWALL

        const Real body_spacing = 1.6*(AXIS_DIR_BODY_LENGTH);

        // initialise plane 

        x[0] = 0;
        x[0] = 0;
        x[0] = 0;

        xm1[0] = x[0];
        xm1[1] = x[1];
        xm1[2] = x[2];
        xm2[0] = x[0];
        xm2[1] = x[1];
        xm2[2] = x[2];

        u[0] = 0.0;
        u[1] = 0.0;
        u[2] = 0.0;
        um1[0] = 0.0;
        um1[1] = 0.0;
        um1[2] = 0.0;

        q = quaternion(1.0, 0.0, 0.0, 0.0);
        
      #endif

    #endif

    const std::string file_name_trunk = GENERATRIX_FILE_NAME+std::to_string(NBLOB);

    std::ifstream pos_file(file_name_trunk + ".seed");
    std::ifstream polar_file(file_name_trunk + ".polar_dir");
    std::ifstream azi_file(file_name_trunk + ".azi_dir");
    std::ifstream normal_file(file_name_trunk + ".normal");

    #if ROD or RIGIDWALL

      seed_rod_blobs(&blob_references[0], &polar_dir_refs[0], &azi_dir_refs[0], &normal_refs[0]);

    #endif

    #if SURFACE_OF_REVOLUTION_BODIES

      if (pos_file.good() && polar_file.good() && azi_file.good() && normal_file.good()){

        for (int i = 0; i < 3*NBLOB; i++){

          pos_file >> blob_references[i];
          polar_file >> polar_dir_refs[i];
          azi_file >> azi_dir_refs[i];
          normal_file >> normal_refs[i];

        }

      } else { // If any are missing, we'll have to re-make them all...

        seed_blobs(&blob_references[0], &polar_dir_refs[0], &azi_dir_refs[0], &normal_refs[0]);

      }
   

    // The seeding functions work on unit-length surfaces, so the scaling must be done after we read or calculate.
    for (int i = 0; i < 3*NBLOB; i++){

      blob_references[i] *= AXIS_DIR_BODY_LENGTH;

    }

    pos_file.close();
    polar_file.close();
    azi_file.close();
    normal_file.close();
    
    #endif

    #if NO_CILIA_SQUIRMER

      // We want to know this because when we have experimental observations for the wavelength of an MCW,
      // we assume these values are for the 'widest' part of the swimmer.
      max_cylindrical_radius = 0.0;

      for (int n = 0; n < NBLOB; n++){

        const Real my_cyl_rad = std::sqrt(blob_references[3*n]*blob_references[3*n] + blob_references[3*n + 1]*blob_references[3*n + 1]);
        max_cylindrical_radius = std::max<Real>(max_cylindrical_radius, my_cyl_rad);

      }

    #endif

  #elif TORUS_BODIES

    // NOT YET IMPLEMENTED
  #endif

}

void rigid_body::initial_guess(const int nt){

  Real initial_guess_x[3], initial_guess_u[3];

  if (nt <= NUM_EULER_STEPS){

    initial_guess_x[0] = 2.0*x[0] - xm1[0];
    initial_guess_x[1] = 2.0*x[1] - xm1[1];
    initial_guess_x[2] = 2.0*x[2] - xm1[2];

    initial_guess_u[0] = u[0];
    initial_guess_u[1] = u[1];
    initial_guess_u[2] = u[2];

  } else {

    initial_guess_x[0] = 3.0*(x[0] - xm1[0]) + xm2[0];
    initial_guess_x[1] = 3.0*(x[1] - xm1[1]) + xm2[1];
    initial_guess_x[2] = 3.0*(x[2] - xm1[2]) + xm2[2];

    initial_guess_u[0] = 2.0*u[0] - um1[0];
    initial_guess_u[1] = 2.0*u[1] - um1[1];
    initial_guess_u[2] = 2.0*u[2] - um1[2];

  }

  xm2[0] = xm1[0];
  xm2[1] = xm1[1];
  xm2[2] = xm1[2];
  xm1[0] = x[0];
  xm1[1] = x[1];
  xm1[2] = x[2];
  x[0] = initial_guess_x[0];
  x[1] = initial_guess_x[1];
  x[2] = initial_guess_x[2];

  um1[0] = u[0];
  um1[1] = u[1];
  um1[2] = u[2];
  u[0] = initial_guess_u[0];
  u[1] = initial_guess_u[1];
  u[2] = initial_guess_u[2];

  qm1 = q;
  lie_exp(q, u);
  q *= qm1;

  q.rot_mat(Q_init);

  #if USE_BROYDEN_FOR_EVERYTHING

    if (nt == 1){

      for (int n = 0; n < 3*NBLOB; n++){

        blob_forces_m1[n] = blob_forces[n];

      }

    } else if (nt == 2){

      for (int n = 0; n < 3*NBLOB; n++){

        blob_forces_m2[n] = blob_forces_m1[n];
        blob_forces_m1[n] = blob_forces[n];
        blob_forces[n] = 2.0*blob_forces_m1[n] - blob_forces_m2[n];

      }

    } else if (nt > 2){

      for (int n = 0; n < 3*NBLOB; n++){

        const Real temp = 3.0*(blob_forces[n] - blob_forces_m1[n]) + blob_forces_m2[n];

        blob_forces_m2[n] = blob_forces_m1[n];
        blob_forces_m1[n] = blob_forces[n];
        blob_forces[n] = temp;

      }

    }

  #endif

}

void rigid_body::blob_positions(Real *const x_array) const {

  const matrix R = q.rot_mat();

  const matrix X(3, 1, x);

  for (int i = 0; i < NBLOB; i++){

    const matrix ref(3, 1, &blob_references[3*i]);

    const matrix pos = X + R*ref;

    x_array[3*i] = pos(0);
    x_array[3*i + 1] = pos(1);
    x_array[3*i + 2] = pos(2);

  }

}

void rigid_body::update(const Real *const body_update){

  #if !INFINITE_PLANE_WALL

    #if USE_BROYDEN_FOR_EVERYTHING

      for (int n = 0; n < 3*NBLOB; n++){

        blob_forces[n] += body_update[6 + n];

      }

    #endif

    x[0] += body_update[0];
    x[1] += body_update[1];
    x[2] += body_update[2];

    u[0] += body_update[3];
    u[1] += body_update[4];
    u[2] += body_update[5];

    lie_exp(q, u);
    q *= qm1;

  #endif

}

void rigid_body::write_reference_positions() const {

  #if !INFINITE_PLANE_WALL

    std::ofstream blob_ref_file(SIMULATION_NAME+std::string("_blob_references.dat"));

    for (int n = 0; n < 3*NBLOB; n++){

      blob_ref_file << blob_references[n] << " " ;

    }

    blob_ref_file << "\n";
    blob_ref_file.close();

  #endif

}

void rigid_body::write_data(std::ofstream& body_state_file) const {

  body_state_file << x[0] << " " << x[1] << " " << x[2] << " " << q.scalar_part << " " << q.vector_part[0] << " " << q.vector_part[1] << " " << q.vector_part[2] << " ";

}

void rigid_body::write_backup(std::ofstream& backup_file) const {

  backup_file << q.scalar_part << " " << q.vector_part[0] << " " << q.vector_part[1] << " " << q.vector_part[2] << " ";

  backup_file << x[0] << " " << x[1] << " " << x[2] << " ";
  backup_file << xm1[0] << " " << xm1[1] << " " << xm1[2] << " ";
  backup_file << xm2[0] << " " << xm2[1] << " " << xm2[2] << " ";

  backup_file << u[0] << " " << u[1] << " " << u[2] << " ";
  backup_file << um1[0] << " " << um1[1] << " " << um1[2] << " ";

  #if USE_BROYDEN_FOR_EVERYTHING

    // Do this after everything else so we can start GMRES from Broyden-only backups.
    for (int n = 0; n < 3*NBLOB; n++){

      backup_file << blob_forces[n] << " ";
      backup_file << blob_forces_m1[n] << " ";
      backup_file << blob_forces_m2[n] << " ";

    }

  #endif

}

#if PRESCRIBED_BODY_VELOCITIES

  void rigid_body::prescribed_translational_velocity(Real *const V, const Real t, const int id) const {

    V[0] = 0.0;
    V[1] = 0.0;
    V[2] = 0.0;

  }

  void rigid_body::prescribed_rotational_velocity(Real *const W, const Real t, const int id) const {

    W[0] = 0.0;
    W[1] = 0.0;
    W[2] = 0.0;

  }

#endif
