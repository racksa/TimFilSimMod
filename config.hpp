// config.hpp
#include <string>

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Include guard
#ifndef MY_CONFIG_HEADER_INCLUDED
#define MY_CONFIG_HEADER_INCLUDED

// Only define it if it is not defined in the makefile
// #ifndef SIMULATION_NAME
//   #define SIMULATION_DIR "data/expr_sims/global/"
//   #define SIMULATION_FILE "cilia"
//   #define SIMULATION_NAME SIMULATION_DIR SIMULATION_FILE
// #endif

extern std::string SIMULATION_DIR;
extern std::string SIMULATION_FILE;
extern std::string SIMULATION_NAME;

extern std::string SIMULATION_CONFIG_NAME;
extern std::string SIMULATION_BACKUP_NAME;
extern std::string SIMULATION_BODY_STATE_NAME; // Blob states are recoverable from body states.
extern std::string SIMULATION_SEG_STATE_NAME;
extern std::string SIMULATION_BODY_VEL_NAME; // Blob velocities are recoverable from body velocities.
extern std::string SIMULATION_SEG_VEL_NAME;
extern std::string SIMULATION_BLOB_FORCES_NAME; // Body forces are recoverable from blob forces.
extern std::string SIMULATION_SEG_FORCES_NAME;
extern std::string SIMULATION_TIME_NAME;
extern std::string SIMULATION_TETHERLAM_NAME;


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Simulation type
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#define CILIA_TYPE 3
// Valid options:
// 0 = Instability-driven cilia. This choice has some sub-types (see below).
// 1 = Geometrically-switching cilia (partially implemented)
// 2 = Constant base rotation (partially implemented)
// 3 = Cilia follow a prescribed sequence of shapes. This choice has some sub-types (see below).
// 4 = Squirmer-type simulation; i.e. there aren't actually any filaments/cilia. The slip velocity can be set in the mobility solver.

#if CILIA_TYPE==0

  #define CILIA_IC_TYPE 2
  // Valid options:
  // 0 = All cilia have identical planar perturbations.
  // 1 = All cilia have identical out-of-plane perturbations.
  // 2 = Cilia have random planar perturbations.
  // 3 = Cilia have random out-of-plane perturbations.

#elif CILIA_TYPE==3

  #define SHAPE_SEQUENCE 1
  // Valid options:
  // 0 = 'Build-a-beat'. This choice has some parameters to set (see below).
  // 1 = The 'Fulford and Blake' beat pattern for mammalian airway cilia. See the data-fitting description in  "A model for the micro-structure in ciliated organisms", Blake (1972).
  // 2 = Coral larvae beat pattern. Data fitting done by me, in the same way as Blake (1972).

  #if SHAPE_SEQUENCE==0

    #define SCALED_BEAT_AMPLITUDE 1.99 // A value in (0,2), giving the beat amplitude in units of filament length.
    #define RECOVERY_STROKE_WINDOW_LENGTH (1.0/12.0) // A value in (0,1), giving the fraction of the beat cycle over which a given point on the filament completes its recovery-stroke tangent angle change.
    #define EFFECTIVE_STROKE_LENGTH 0.75 // A value in (0,1), giving the fraction of the cycle spent in the effective stroke.
    // N.B. We must have RECOVERY_STROKE_WINDOW_LENGTH + EFFECTIVE_STROKE_LENGTH < 1
    #define ZERO_VELOCITY_AVOIDANCE_LENGTH 0.05 // A value in (0,1), giving the maximum fraction of the cycle by which we shift the tangent angle curve to ensure the velocity cannot be zero everywhere along the filament at once.

  #endif

  #define DYNAMIC_PHASE_EVOLUTION true
  // If true, cilia phase speeds are solved for as part of the dynamics. Note that this requires having run a reference simulation with WRITE_GENERALISED_FORCES=true previously.
  // If false, phase_dot = omega0 is constant for each cilium.

  #define DYNAMIC_SHAPE_ROTATION true
  // If true, the vertical in the cilia reference configuration can rotate with respect to the surface normal.
  // Essentially, the cilia can 'tip backwards or forwards' in their beat planes.
  // If false, no such rotation ever occurs.

  #define INITIAL_PHASE 3
  // 0 = Random
  // 1 = All zeros
  // 2 = Ishikawa
  // 3 = Diaplectic

  #define WRITE_GENERALISED_FORCES false
  // If true, this simulation will save its generalised forces to file for use as the reference values.
  // It will also generate reference s-values for shape sequences which don't result in inextensible filaments.
  // NOTE: This will overwrite any existing reference files unless their names have been changed.

  #define CILIA_IC_TYPE 1
  // Valid options:
  // 0 = All cilia start in-phase with phase 0.
  // 1 = Cilia start with a (uniformly) random initial phase.
  // 2 = A metachronal wave (MCW). Its wavelength and direction are defined below.

  #if CILIA_IC_TYPE==2

    #define CILIA_MCW_WAVELENGTH FIL_LENGTH

    // This angle is measured anti-clockwise from the positive beat-wise reference-coordinate direction. Exactly which stroke this equates to will depend on the specific beat used.
    // Note: Only the sign of this angle is considered on spherical surfaces, where the nature of the seeding means we do not control the beat-wise displacements
    // and hence we only allow the wave to travel azimuthally. Positive values will have the wave propagate in the direction of increasing azimuthal angle.
    #define CILIA_MCW_ANGLE (0.5*PI)

  #endif

#endif

#define BODY_OR_SURFACE_TYPE 2
// Valid options:
// 0 = An infinite plane wall at z = 0. This choice has some sub-types (see below).
// 1 = Deformed planes with 2 principal curvatures (partially implemented)
// 2 = Surface-of-revolution bodies. This choice has some sub-types (see below).
// 3 = Toroidal bodies (partially implemented)
// 4 = Rigid Rod
// 5 = Filament on rigid wall

#define BASE_HEIGHT_ABOVE_SURFACE (0.5*DL)

#if BODY_OR_SURFACE_TYPE==0

  #define SEEDING_TYPE 2
  // Valid options:
  // 0 = Filaments are placed on a rectangular grid.
  // 1 = Filaments are placed on a hexagonal grid.
  // 2 = Filaments are placed on a rectangular grid compatable with cuFCM.
  // 3 = Filaments are placed on a lattice compatable with cuFCM.

  // Define one lattice size and leave the other blank to have it automatically calculated to fit the number of filaments.
  // Leave both blank to have the code attempt to make a regular lattice.
  // If values are given for both, the code will use FIL_LATTICE_X_NUM and re-calculate FIL_LATTICE_Y_NUM.
  #define FIL_LATTICE_X_NUM 1
  #define FIL_LATTICE_Y_NUM
  #define FIL_LATTICE_Y_SPACING (1.8*FIL_LENGTH) // For prescribed-shape cilia, this is the beat-wise separation.
  #define FIL_LATTICE_X_SPACING (2.0*DL) // Leave blank for a regular grid (i.e. squares or regular hexagons (as appropriate) based on FIL_LATTICE_Y_SPACING).

#elif BODY_OR_SURFACE_TYPE==2 or BODY_OR_SURFACE_TYPE==4 or BODY_OR_SURFACE_TYPE==5

  #define SEEDING_TYPE 7
  // Valid options:
  // 0 = Filaments are evenly distributed over the surface.
  // 1 = Filaments are seeded in an equatorial band.
  // 2 = Platynereis-style seeding. Most filament are in an equatorial band but some form a small ring at the rear of the swimmer.
  // 3 = Hexagonal grid seeding. (rigidbody plane)
  // 4 = Meridian seeding
  // 5 = Icosa seeding
  // 6 = Mismatched seeding
  // 7 = Filaments are evenly distributed over the surface but with potential at poles

  #define FOURIER_DIR "data/fourier_modes/"
  #define GENERATRIX_FILE_NAME FOURIER_DIR "sphere"
  // #define GENERATRIX_FILE_NAME FOURIER_DIR "avg_shape_plus_1.5_first_mode"
  // #define GENERATRIX_FILE_NAME FOURIER_DIR "rod_mode"

  // The code will search for a file called GENERATRIX_FILE_NAME.fourier_modes
  // The first entry in this file should be a positive integer N, which is the number of Fourier modes used to describe the surface's generatrix.
  // The remaining 2N entries are the Fourier coefficients, in (cos_coeff, sin_coeff) pairs for each mode in turn.
  // Note that this will only work for convex surfaces which contain the body centre in their interior.
  // Note that it also assumes the surface is centred at the origin and that it has unit length parallel to the rotation axis for the revolution
  // of the surface. The script format_shape_fourier_modes.m can be used to ensure this is true.

#endif

// Define whether the motion of the rigid bodies is imposed or allowed to evolve dynamically.
#define PRESCRIBED_BODY_VELOCITIES false

#define MOBILITY_TYPE 4
// Valid options:
// 0 = Basic Stokes drag. No hydrodynamic interactions between particles.
// 1 = Rotne-Prager-Yamakawa (RPY) mobility matrices (with the corrections due to Swan and Brady if an infinite plane wall is selected).
// 2 = Weakly-coupled-filaments RPY mobility. Intra-filament hydrodynamics are as in option 1, but we use an approximation for inter-filament interactions.
// 3 = The Force Coupling Method (FCM). Implemented using the UAMMD library (https://github.com/RaulPPelaez/UAMMD).
// 4 = cuFCM

#define INITIAL_CONDITIONS_TYPE 0
// Valid options:
// 0 = Default
// 1 = Resume from backup file (only really valid if we ensure DT is the same etc., but this is NOT checked automatically)
// 2 = Fresh start from backup file (i.e. uses saved state but resets t=0. Use for hold-then-release style simulations)
#define INITIAL_CONDITIONS_FILE_NAME SIMULATION_NAME // SIMULATION_NAME or "a_different_sim_name"
// N.B. Simulations using GMRES can resume/start using backup files from Broyden-only simulations, but the opposite is not true.
// N.B. For options 0 and 2, whilst the simulation state will be fresh, all saved data will still be appended to any data from a previous simulation of the same name.

#if MOBILITY_TYPE==4
  #define CUFCM_CONFIG_FILE_NAME "../CUFCM/simulation_info_cilia"
#endif

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Physical parameters
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
extern int NSWIM;
extern int NSEG;
extern int NFIL;
extern int NBLOB;
extern float AR;
extern float AXIS_DIR_BODY_LENGTH;
extern float TORSIONAL_SPRING_MAGNITUDE_FACTOR; // Pre-multiplies the mean generalised driving force magnitude to give the spring constant that resists rigid-body rotations of the shape.
extern int NTOTAL;
extern float END_FORCE_MAGNITUDE;

#define MU 1.0 // Fluid viscosity.

#define RSEG 1.0 // Segment radius.
#define DL (2.2*RSEG) // Inter-segment distance.
#define RBLOB 1.0 // Surface blob radius.

#define KB 1800.0 // Bending modulus.
#define KT 1800.0 // Twist modulus.

// #if BODY_OR_SURFACE_TYPE==2
//   #define AXIS_DIR_BODY_LENGTH (AR*FIL_LENGTH)
//   // #define AXIS_DIR_BODY_LENGTH (1.6496*2.0*FIL_LENGTH) // The length of the body parallel to the axis of rotation for the surface of revolution.
// #elif BODY_OR_SURFACE_TYPE==4
//   #define AXIS_DIR_BODY_LENGTH (0.5*NBLOB*RBLOB) // The length of the body parallel to the axis of rotation for the surface of revolution.
// #elif BODY_OR_SURFACE_TYPE==5
//   #define AXIS_DIR_BODY_LENGTH (NBLOB*RBLOB)
// #endif

#if CILIA_TYPE==0 or CILIA_TYPE==3

  #define DIMENSIONLESS_FORCE 220.0 

#elif CILIA_TYPE==1

  #define DRIVING_TORQUE_MAGNITUDE_RATIO 3.0 // How much stronger is the driving torque in the fast stroke than in the slow
  #define DEFAULT_BASE_TORQUE_MAGNITUDE (0.1*KB) // Driving torque magnitude in the fast stroke
  #define CRITICAL_ANGLE (0.4*PI) // The angle at which the cilia changes stroke

#elif CILIA_TYPE==2

  #define BASE_ROTATION_RATE 0.1

#endif

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Computational parameters
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// Threads per block for kernel execution. Should be a multiple of 32, the warp size.
#define THREADS_PER_BLOCK 64

// This factor avoids over-adjusting during the Broyden's iterations.
// Setting it to 1 will give standard behaviour, and values smaller than 1 should help with convergence problems by having Broyden's method do more of the heavy lifting.
// Choosing 0.4 seems to work well when using GMRES, whilst 0.1 appears to be a good choice for a "Broyden's only" simulation.
#define JACOBIAN_CONFIDENCE_FACTOR 0.1

#define MAX_BROYDEN_ITER 400 // Maximum number of Broyden's method iterations per time-step.
#define TOL 1e-4 // Tolerance to be reached by the Broyden's method solve.

#define SOLVER_TYPE 1
// Valid options:
// 0: Use Broyden's method for absolutely everything. When there is a rigid body with forces (and velocities if they're not prescribed) to solve for,
//    the associated linear system is embedded in the wider Broyden's solve, rather than being solved for the current iterate at each iteration.
// 1: Use GMRES to solve the linear system at each iteration of Broyden's method.

#if SOLVER_TYPE==1

  #define MAX_LINEAR_SYSTEM_ITER 350 // Maximum number of iterations used to solve the linear system in each mobility solve.
  #define LINEAR_SYSTEM_TOL 1e-4 // Relative tolerance in the linear system solves.

  // GMRES preconditioner type.
  // Uses left preconditioning if set to false; if you don't want a preconditioner,
  // simply write "return in;" at the start of rpy_mobility_solver::apply_preconditioner(...) and it should be optimised away.
  // Right preconditioning seems to result in fewer GMRES iterations and also means that the error in GMRES
  // is the same as the error in the original system we want to solve, so this is the default option.
  // This ought to be checked whenever the preconditioner is changed though.
  #define USE_RIGHT_PRECON true

#endif

#define TOTAL_TIME_STEPS (30*STEPS_PER_PERIOD) // Total number of time-steps in the simulation.
#define NUM_EULER_STEPS 1 // Number of time-steps to use backwards-Euler before switching to BDF2.

#if CILIA_TYPE==1

  #define TRANSITION_STEPS 21 // The number of time-steps over which to change the driving torque between strokes, must be >= 1
  #define DT 2.0
  #define PLOT_FREQUENCY_IN_STEPS 10

#else

  #define STEPS_PER_PERIOD 300
  #define SAVES_PER_PERIOD 30

#endif


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Derived/redefined parameters (these should be left alone)
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#define GEOMETRIC_CILIA (CILIA_TYPE==1)
#define INSTABILITY_CILIA (CILIA_TYPE==0)
#define CONSTANT_BASE_ROTATION (CILIA_TYPE==2)
#define PRESCRIBED_CILIA (CILIA_TYPE==3)
#define NO_CILIA_SQUIRMER (CILIA_TYPE==4)

#define INFINITE_PLANE_WALL (BODY_OR_SURFACE_TYPE==0)
#define SADDLE_BODIES (BODY_OR_SURFACE_TYPE==1)
#define SURFACE_OF_REVOLUTION_BODIES (BODY_OR_SURFACE_TYPE==2)
#define TORUS_BODIES (BODY_OR_SURFACE_TYPE==3)
#define ROD (BODY_OR_SURFACE_TYPE==4)
#define RIGIDWALL (BODY_OR_SURFACE_TYPE==5)

#define STOKES_DRAG_MOBILITY (MOBILITY_TYPE==0)
#define RPY_MOBILITY (MOBILITY_TYPE==1)
#define WEAKLY_COUPLED_FILAMENTS_RPY_MOBILITY (MOBILITY_TYPE==2)
#define UAMMD_FCM (MOBILITY_TYPE==3)
#define CUFCM (MOBILITY_TYPE==4)

#define USE_BROYDEN_FOR_EVERYTHING (SOLVER_TYPE==0)
#define USE_GMRES_FOR_LINEAR_SYSTEM (SOLVER_TYPE==1)

#define READ_INITIAL_CONDITIONS_FROM_BACKUP (INITIAL_CONDITIONS_TYPE != 0)

#define PI 3.14159265358979323846264338327950288

// #define SIMULATION_CONFIG_NAME SIMULATION_NAME ".par"
// #define SIMULATION_BACKUP_NAME SIMULATION_NAME ".backup"
// #define SIMULATION_BODY_STATE_NAME SIMULATION_NAME "_body_states.dat" // Blob states are recoverable from body states.
// #define SIMULATION_SEG_STATE_NAME SIMULATION_NAME "_seg_states.dat"
// #define SIMULATION_BODY_VEL_NAME SIMULATION_NAME "_body_vels.dat" // Blob velocities are recoverable from body velocities.
// #define SIMULATION_SEG_VEL_NAME SIMULATION_NAME "_seg_vels.dat"
// #define SIMULATION_BLOB_FORCES_NAME SIMULATION_NAME "_blob_forces.dat" // Body forces are recoverable from blob forces.
// #define SIMULATION_SEG_FORCES_NAME SIMULATION_NAME "_seg_forces.dat"
// #define SIMULATION_TIME_NAME SIMULATION_NAME "_time.dat"
// #define SIMULATION_TETHERLAM_NAME SIMULATION_NAME "_tether_force.dat"

#define DELETE_CURRENT_LINE "                                                                                                               " << "\r"

#if USE_BROYDEN_FOR_EVERYTHING

  #define NBROY (3*NSWIM*(NBLOB + 2*(NFIL*NSEG + 1)))

#else

  #define NBROY (6*(NSWIM*NFIL*NSEG + NSWIM))

#endif

#if !GEOMETRIC_CILIA

  #define PLOT_FREQUENCY_IN_STEPS (STEPS_PER_PERIOD/SAVES_PER_PERIOD)

#endif

#if INSTABILITY_CILIA
  
  // #define END_FORCE_MAGNITUDE (DIMENSIONLESS_FORCE*KB/(DL*DL*NSEG*NSEG))
  #define REPULSIVE_FORCE_FACTOR 2.0 // How much stronger is the barrier force than the driving force.
  #define DT (36.3833/STEPS_PER_PERIOD) // Based on the period of a single DIMENSIONLESS_FORCE = 220.0 filament above a no-slip wall.

#elif CONSTANT_BASE_ROTATION

  #define DT (2.0*PI/(BASE_ROTATION_RATE*STEPS_PER_PERIOD))
  #define REPULSIVE_FORCE_FACTOR 1.0
  #define END_FORCE_MAGNITUDE (0.1*KB) // There isn't really an end force, it's just used to define the repulsion.

#elif PRESCRIBED_CILIA

  #define DT (1.0/STEPS_PER_PERIOD) // Pick T = 1.

  #if USE_BROYDEN_FOR_EVERYTHING

    #error "The fully-prescribed cilia motion option currently doesn't support using Broyden's method to solve the system."

  #endif

  #if WRITE_GENERALISED_FORCES

    #define PRESCRIBED_BODY_VELOCITIES true
    #define DYNAMIC_PHASE_EVOLUTION false
    #define DYNAMIC_SHAPE_ROTATION false
    #define TOTAL_TIME_STEPS STEPS_PER_PERIOD
    #define INITIAL_CONDITIONS_TYPE 0

  #endif

  #define BUILD_A_BEAT (SHAPE_SEQUENCE==0)

  #define FIT_TO_DATA_BEAT (SHAPE_SEQUENCE != 0)
  #define FULFORD_AND_BLAKE_BEAT (SHAPE_SEQUENCE==1)
  #define CORAL_LARVAE_BEAT (SHAPE_SEQUENCE==2)

#endif

#if INFINITE_PLANE_WALL
  #define NBROY (6*NFIL*NSEG)
  #define PRESCRIBED_BODY_VELOCITIES true
  #if RPY_MOBILITY
    #define RSEG 1.0
  #endif
  #define MU 1.0

  #define RECTANGULAR_SEEDING (SEEDING_TYPE==0)
  #define HEXAGONAL_SEEDING (SEEDING_TYPE==1)
  #define FCM_RECTANGULAR_SEEDING (SEEDING_TYPE==2)
  #define FCM_LATTICE_SEEDING (SEEDING_TYPE==3)

#endif

#if ROD
  #define DT (0.1) 
#endif

#if NO_CILIA_SQUIRMER

  #define NFIL 0
  #define NTOTAL (NSWIM*NBLOB)
  #define DT (1.0/STEPS_PER_PERIOD) // Pick T = 1.

  #if !SURFACE_OF_REVOLUTION_BODIES

    // #error "Squirmer-type simulations are only compatible with surface-of-revolution bodies."

  #endif

#endif

#if SURFACE_OF_REVOLUTION_BODIES or RIGIDWALL

  #define UNIFORM_SEEDING (SEEDING_TYPE==0)
  #define EQUATORIAL_SEEDING (SEEDING_TYPE==1)
  #define PLATY_SEEDING (SEEDING_TYPE==2)
  #define HEXAGONAL_WALL_SEEDING (SEEDING_TYPE==3)
  #define MERIDIAN_SEEDING (SEEDING_TYPE==4)
  #define ICOSA_SEEDING (SEEDING_TYPE==5)
  #define MISMATCH_SEEDING (SEEDING_TYPE==6)
  #define UNIFORM_SEEDING_POLE (SEEDING_TYPE==7)

#endif

#if WEAKLY_COUPLED_FILAMENTS_RPY_MOBILITY

  #if !PRESCRIBED_CILIA

    // This is because I've only worked it out for forces along the filament, not torques as well.
    // It also means we don't have to be concerned about calculating the total force and total first moment on each filament prior to
    // the barrier_forces(...) call, because this does nothing for prescribed-shape filaments.
    #error "The weakly-coupled-filaments mobility currently only supports filaments which follow a prescribed sequence of shapes."

  #endif

#endif

#if UAMMD_FCM

  #define RBLOB RSEG // FCM requires that all radii are the same.

#endif

#define DEFINED_BUT_EMPTY(VAR) ((~(~VAR + 0) == 0) && (~(~VAR + 1) == 1)) // Macro magic...

// This is just how we've ended up defining the length in these different cases.
// The length is always a bit ambiguous up to DL anyway; e.g. when we try to reproduce
// the results in Landau and Lifshitz for an elastic filament under a load at the free end,
// we actually need to use the PRESCRIBED_CILIA definition instead
// (because nothing elastic happens outside of that part of the filament).
#if PRESCRIBED_CILIA

  #define FIL_LENGTH (DL*(NSEG-1))

#else

  #define FIL_LENGTH (DL*NSEG)

#endif

// #if (DYNAMIC_SHAPE_ROTATION && !PRESCRIBED_BODY_VELOCITIES)

//   #error "I haven't implemented the preconditioner for free bodies with DYNAMIC_SHAPE_ROTATION enabled yet..."

// #endif

#define DISPLAYTIME true

#define FIL_USE_DOUBLE_PRECISION false

#if FIL_USE_DOUBLE_PRECISION
    typedef double Real;
    typedef long Integer;
    #define myfil_rint rint
    #define myfil_exp exp
    #define myfil_floor floor
    #define myfil_fmod fmod
    #define myfil_getrf_ dgetrf_
    #define myfil_getri_ dgetri_
    #define myfil_gemm_ dgemm_
#else
    typedef float Real;
    typedef int Integer;
    #define myfil_rint rintf
    #define myfil_exp expf
    #define myfil_floor floorf
    #define myfil_fmod fmodf
    #define myfil_getrf_ sgetrf_
    #define myfil_getri_ sgetri_
    #define myfil_gemm_ sgemm_
    
#endif

#endif // MY_CONFIG_HEADER_INCLUDED
