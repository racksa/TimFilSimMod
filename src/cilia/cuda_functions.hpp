// cuda_functions.hpp

// =============================================================================
// Include guard
#ifndef MY_CUDA_FUNCTIONS_HEADER_INCLUDED
#define MY_CUDA_FUNCTIONS_HEADER_INCLUDED

// RPY kernels
__global__ void Mss_mult(double *V, const double *const F, const double *const X, const int start_seg, const int num_segs);
__global__ void Mbb_mult(double *V, const double *const F, const double *const X, const int start_blob, const int num_blobs);
__global__ void Msb_mult(double *V, const double *const F, const double *const Xs, const double *const Xb, const int start_seg, const int num_segs);
__global__ void Mbs_mult(double *V, const double *const F, const double *const Xb, const double *const Xs, const int start_blob, const int num_blobs);

// Stokes drag kernels
__global__ void Ms_mult(double *V, const double *const F, const int start_seg, const int num_segs);
__global__ void Mb_mult(double *V, const double *const F, const int start_blob, const int num_blobs);
__global__ void Mb_fill_zero(double *V, const int start_blob, const int num_blobs);

// Weakly-coupled-filaments kernels
__global__ void reset_mean_fil_posns(double *Z);
__global__ void write_mean_fil_posns(double *Z, const double *const X);
__global__ void reset_mean_fil_force_and_moment(double * __restrict__ F, double * __restrict__ S);
__global__ void find_total_fil_force_and_first_moment(double *F, double *S, const double *const Fs, const double *const Xs, const double *const Z);
__global__ void Msf_mult(double *V, const double *const F, const double *const Ftot, const double *const Stot, const double *const X, const double *const Z, const int start_seg, const int num_segs);

// Generic interaction kernels
__global__ void barrier_forces(double *f_segs, double *f_blobs_repulsion, const double *const x_segs, const double *const x_blobs, const int start_seg, const int num_segs, const int start_blob, const int num_blobs);
__global__ void periodic_barrier_forces(double *f_segs, double *f_blobs_repulsion, const double *const x_segs, const double *const x_blobs, const int num_segs, const int num_blobs, const double boxsize);

__host__ __device__
void box_images(double &x, double box_size);


#endif // MY_CUDA_FUNCTIONS_HEADER_INCLUDED
