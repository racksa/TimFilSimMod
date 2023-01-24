// uammd_fcm_mobility_solver.cuh

// =============================================================================
// Include guard
#ifndef MY_UAMMD_FCM_MOBILITY_SOLVER_HEADER_INCLUDED
#define MY_UAMMD_FCM_MOBILITY_SOLVER_HEADER_INCLUDED

// =============================================================================
// Forward declared dependencies

// =============================================================================
// Included dependencies
#include "mobility_solver.hpp"

//
// An interface for UAMMD FCM, provided by Raul Pelaez.
//
#include "Integrator/BDHI/FCM/FCM_kernels.cuh"
#include "uammd.cuh"
#include "Integrator/BDHI/FCM/FCM_impl.cuh"
using Kernel = uammd::BDHI::FCM_ns::Kernels::Gaussian;
using KernelTorque = uammd::BDHI::FCM_ns::Kernels::GaussianTorque;
using FCM = uammd::BDHI::FCM_impl<Kernel, KernelTorque>;
namespace uammd_fcm{

    template<class T> using cached_vector = uammd::BDHI::cached_vector<T>;
    using uammd::real3;
    using uammd::real4;
    using uammd::real;
  
    auto initializeFCM(real hydrodynamicRadius, real viscosity,
               real Lx, real Ly, real Lz,
               real tolerance){
  
      FCM::Parameters par;
      par.viscosity = viscosity;
      par.box = uammd::Box({Lx, Ly, Lz});
      par.hydrodynamicRadius = hydrodynamicRadius;
      par.tolerance = tolerance;
      par.adaptBoxSize = true;
      auto fcm = std::make_shared<FCM>(par);
      return fcm;
    }
  
    namespace detail{
  
      struct InterleavedPermute{
        int instance;
        InterleavedPermute(int instance):instance(instance){}
  
        __host__ __device__ auto operator()(int i){
      return 2*i+instance;
        }
      };
    }
  
    template<class Iter>
    auto make_interleaved_iterator(Iter it, int instance){
      auto cit = thrust::make_transform_iterator(thrust::make_counting_iterator(0),
                             detail::InterleavedPermute(instance));
      return thrust::make_permutation_iterator(it, cit);
    }
  
    __global__ void interleaveResults(real* f, real* t, real* result, int n){
      int id = blockIdx.x*blockDim.x + threadIdx.x;
      if(id>=n) return;
      for(int i=0; i<3; i++){
        result[6*id+i] = f[3*id+i];
        result[6*id+3+i] = t[3*id+i];
      }
    }
    auto computeHydrodynamicDisplacements(std::shared_ptr<FCM> fcm,
                      real* i_pos, real* i_forceTorque,
                      int i_numberParticles,
                      real* o_pos, int o_numberParticles,
                      cudaStream_t st = 0){
      real temperature = 0;
      real prefactor = 0;
      auto i_ft = (real3*)i_forceTorque;
      auto i_force = make_interleaved_iterator(i_ft, 0);
      auto i_torque = make_interleaved_iterator(i_ft, 1);
  
      auto res = fcm->computeHydrodynamicDisplacements((real3*)i_pos, i_force, i_torque, i_numberParticles,
                               (real3*)o_pos, o_numberParticles,
                               temperature, prefactor, st);
  
      cached_vector<real> result(6*o_numberParticles);
      int nthreads = THREADS_PER_BLOCK;
      int nblocks = o_numberParticles/nthreads +1;
      interleaveResults<<<nblocks, nthreads>>>((real*)thrust::raw_pointer_cast(res.first.data()),
                           (real*)thrust::raw_pointer_cast(res.second.data()),
                           (real*)thrust::raw_pointer_cast(result.data()),
                           o_numberParticles);
      return result;
    }
  
}

// Because UAMMD is a header library, we cannot include it in multiple source files or the compiler
// will throw "multiple definition" errors. The workaround is to define the UAMMD-based mobility class
// in its header too, and not have a separate .cu file.

class uammd_fcm_mobility_solver : public mobility_solver {

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

    // These vectors will be aliased by the raw pointers required by the base interface class.
    thrust::host_vector<Real> v_segs_host_thrust;
    thrust::host_vector<Real> v_blobs_host_thrust;
    thrust::host_vector<Real> x_segs_host_thrust;
    thrust::host_vector<Real> x_blobs_host_thrust;
    thrust::host_vector<Real> f_segs_host_thrust;
    thrust::host_vector<Real> f_blobs_host_thrust;
    thrust::host_vector<Real> f_blobs_repulsion_host_thrust;

    thrust::device_vector<Real> v_segs_device_thrust;
    thrust::device_vector<Real> v_blobs_device_thrust;
    thrust::device_vector<Real> x_segs_device_thrust;
    thrust::device_vector<Real> x_blobs_device_thrust;
    thrust::device_vector<Real> f_segs_device_thrust;
    thrust::device_vector<Real> f_blobs_device_thrust;
    thrust::device_vector<Real> f_blobs_repulsion_device_thrust;

    std::shared_ptr<FCM> fcm;

    ~uammd_fcm_mobility_solver();
    uammd_fcm_mobility_solver();

};

#include "config.hpp"

uammd_fcm_mobility_solver::~uammd_fcm_mobility_solver(){}
uammd_fcm_mobility_solver::uammd_fcm_mobility_solver(){}

// All memory is freed for us with UAMMD (seemingly by virtue of its use of the thrust library)
void uammd_fcm_mobility_solver::free_host_memory(){}
void uammd_fcm_mobility_solver::free_device_memory(){}

void uammd_fcm_mobility_solver::allocate_host_memory(){

    v_segs_host_thrust.resize(6*NSWIM*NFIL*NSEG);
    v_segs_host = (Real *) thrust::raw_pointer_cast(v_segs_host_thrust.data());

    v_blobs_host_thrust.resize(3*NSWIM*NBLOB);
    v_blobs_host = (Real *) thrust::raw_pointer_cast(v_blobs_host_thrust.data());

    x_segs_host_thrust.resize(3*NSWIM*NFIL*NSEG);
    x_segs_host = (Real *) thrust::raw_pointer_cast(x_segs_host_thrust.data());

    x_blobs_host_thrust.resize(3*NSWIM*NBLOB);
    x_blobs_host = (Real *) thrust::raw_pointer_cast(x_blobs_host_thrust.data());

    f_segs_host_thrust.resize(6*NSWIM*NFIL*NSEG);
    f_segs_host = (Real *) thrust::raw_pointer_cast(f_segs_host_thrust.data());

    f_blobs_host_thrust.resize(3*NSWIM*NBLOB);
    f_blobs_host = (Real *) thrust::raw_pointer_cast(f_blobs_host_thrust.data());

    f_blobs_repulsion_host_thrust.resize(3*NSWIM*NBLOB);
    f_blobs_repulsion_host = (Real *) thrust::raw_pointer_cast(f_blobs_repulsion_host_thrust.data());

    // These should be provided in the config file eventually.
    int domain_length = 128; // Measured in grid points. The domain must be a cube because of the way Raul forced it to keep the radius as RSEG ( = RBLOB).
    Real tolerance = 1e-4;

    fcm = uammd_fcm::initializeFCM(RSEG, MU, domain_length, domain_length, domain_length, tolerance);

}

void uammd_fcm_mobility_solver::allocate_device_memory(){

    v_segs_device_thrust.resize(6*NSWIM*NFIL*NSEG);
    v_blobs_device_thrust.resize(3*NSWIM*NBLOB);
    x_segs_device_thrust.resize(3*NSWIM*NFIL*NSEG);
    x_blobs_device_thrust.resize(3*NSWIM*NBLOB);
    f_segs_device_thrust.resize(6*NSWIM*NFIL*NSEG);
    f_blobs_device_thrust.resize(3*NSWIM*NBLOB);
    f_blobs_repulsion_device_thrust.resize(3*NSWIM*NBLOB);

}

void uammd_fcm_mobility_solver::copy_segment_positions_to_device(){

    x_segs_device_thrust = x_segs_host_thrust;

}

void uammd_fcm_mobility_solver::copy_segment_forces_to_device(){

    f_segs_device_thrust = f_segs_host_thrust;

}

void uammd_fcm_mobility_solver::copy_blob_positions_to_device(){

    x_blobs_device_thrust = x_blobs_host_thrust;

}

void uammd_fcm_mobility_solver::copy_blob_forces_to_device(){

    f_blobs_device_thrust = f_blobs_host_thrust;

}

void uammd_fcm_mobility_solver::copy_interparticle_blob_forces_to_host(){

    f_blobs_repulsion_host_thrust = f_blobs_repulsion_device_thrust;

}

void uammd_fcm_mobility_solver::copy_blob_velocities_to_host(){

    v_blobs_host_thrust = v_blobs_device_thrust;

}

void uammd_fcm_mobility_solver::copy_segment_velocities_to_host(){

    v_segs_host_thrust = v_segs_device_thrust;

}

void uammd_fcm_mobility_solver::apply_interparticle_forces(){

    // I'll worry about this later, because it's going to require me learning about how UAMMD implements its linked list,
    // or what built-in force options it provides.

}

// All copying is done for me by UAMMD (read: thrust).
void uammd_fcm_mobility_solver::wait_for_device(){}

void uammd_fcm_mobility_solver::evaluate_segment_segment_mobility(){

    auto i_pos_ptr = thrust::raw_pointer_cast(x_segs_device_thrust.data());
    auto o_pos_ptr = thrust::raw_pointer_cast(x_segs_device_thrust.data());
    auto ft_ptr = thrust::raw_pointer_cast(f_segs_device_thrust.data());
    auto result = uammd_fcm::computeHydrodynamicDisplacements(fcm, i_pos_ptr, ft_ptr, NSWIM*NFIL*NSEG, o_pos_ptr, NSWIM*NFIL*NSEG);
    thrust::copy(result.begin(), result.end(), v_segs_device_thrust.begin()); // There must be a way of putting the input straight into my own thrust vector...

}

void uammd_fcm_mobility_solver::evaluate_segment_blob_mobility(){

    auto i_pos_ptr = thrust::raw_pointer_cast(x_blobs_device_thrust.data());
    auto o_pos_ptr = thrust::raw_pointer_cast(x_segs_device_thrust.data());
    auto ft_ptr = thrust::raw_pointer_cast(f_blobs_device_thrust.data());
    auto result = uammd_fcm::computeHydrodynamicDisplacements(fcm, i_pos_ptr, ft_ptr, NSWIM*NFIL*NSEG, o_pos_ptr, NSWIM*NBLOB);
    thrust::copy(result.begin(), result.end(), v_segs_device_thrust.begin()); // There must be a way of putting the input straight into my own thrust vector...

}

void uammd_fcm_mobility_solver::evaluate_blob_blob_mobility(){

    auto i_pos_ptr = thrust::raw_pointer_cast(x_blobs_device_thrust.data());
    auto o_pos_ptr = thrust::raw_pointer_cast(x_blobs_device_thrust.data());
    auto ft_ptr = thrust::raw_pointer_cast(f_blobs_device_thrust.data());
    auto result = uammd_fcm::computeHydrodynamicDisplacements(fcm, i_pos_ptr, ft_ptr, NSWIM*NBLOB, o_pos_ptr, NSWIM*NBLOB);
    thrust::copy(result.begin(), result.end(), v_blobs_device_thrust.begin()); // There must be a way of putting the input straight into my own thrust vector...

}

void uammd_fcm_mobility_solver::evaluate_blob_segment_mobility(){

    auto i_pos_ptr = thrust::raw_pointer_cast(x_segs_device_thrust.data());
    auto o_pos_ptr = thrust::raw_pointer_cast(x_blobs_device_thrust.data());
    auto ft_ptr = thrust::raw_pointer_cast(f_segs_device_thrust.data());
    auto result = uammd_fcm::computeHydrodynamicDisplacements(fcm, i_pos_ptr, ft_ptr, NSWIM*NBLOB, o_pos_ptr, NSWIM*NFIL*NSEG);
    thrust::copy(result.begin(), result.end(), v_blobs_device_thrust.begin()); // There must be a way of putting the input straight into my own thrust vector...

}

#endif // MY_UAMMD_FCM_MOBILITY_SOLVER_HEADER_INCLUDED
