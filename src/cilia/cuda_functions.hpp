// cuda_functions.hpp
#include "../../config.hpp"
// =============================================================================
// Include guard
#ifndef MY_CUDA_FUNCTIONS_HEADER_INCLUDED
#define MY_CUDA_FUNCTIONS_HEADER_INCLUDED

// RPY kernels
__global__ void Mss_mult(Real *V, const Real *const F, const Real *const X, const int start_seg, const int num_segs);
__global__ void Mbb_mult(Real *V, const Real *const F, const Real *const X, const int start_blob, const int num_blobs);
__global__ void Msb_mult(Real *V, const Real *const F, const Real *const Xs, const Real *const Xb, const int start_seg, const int num_segs);
__global__ void Mbs_mult(Real *V, const Real *const F, const Real *const Xb, const Real *const Xs, const int start_blob, const int num_blobs);
__global__ void Mbs_mult_add(Real *V, const Real *const F, const Real *const Xb, const Real *const Xs, const int start_blob, const int num_blobs);


// Stokes drag kernels
__global__ void Ms_mult(Real *V, const Real *const F, const int start_seg, const int num_segs);
__global__ void Mb_mult(Real *V, const Real *const F, const int start_blob, const int num_blobs);
__global__ void Mb_fill_zero(Real *V, const int start_blob, const int num_blobs);

// Weakly-coupled-filaments kernels
__global__ void reset_mean_fil_posns(Real *Z);
__global__ void write_mean_fil_posns(Real *Z, const Real *const X);
__global__ void reset_mean_fil_force_and_moment(Real * __restrict__ F, Real * __restrict__ S);
__global__ void find_total_fil_force_and_first_moment(Real *F, Real *S, const Real *const Fs, const Real *const Xs, const Real *const Z);
__global__ void Msf_mult(Real *V, const Real *const F, const Real *const Ftot, const Real *const Stot, const Real *const X, const Real *const Z, const int start_seg, const int num_segs);

// Generic interaction kernels
__global__ void barrier_forces(Real *f_segs, Real *f_blobs_repulsion, const Real *const x_segs, const Real *const x_blobs, const int start_seg, const int num_segs, const int start_blob, const int num_blobs);
__global__ void periodic_barrier_forces(Real *f_segs, Real *f_blobs_repulsion, const Real *const x_segs, const Real *const x_blobs, const int num_segs, const int num_blobs, const Real boxsize);

__host__ __device__
void box_images(Real &x, Real box_size);

__global__
void sync_var(int nswim, int nseg, int nfil, int nblob, int end_force_magnitude, Real seg_sep);

#endif // MY_CUDA_FUNCTIONS_HEADER_INCLUDED
