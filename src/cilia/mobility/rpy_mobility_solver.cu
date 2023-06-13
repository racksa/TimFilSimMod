// rpy_mobility_solver.cu

#include <iomanip>
#include <iostream>
#include <cmath>
#include <algorithm>
#include "rpy_mobility_solver.hpp"
#include "../cuda_functions.hpp"

rpy_mobility_solver::~rpy_mobility_solver(){}

rpy_mobility_solver::rpy_mobility_solver(){}

void rpy_mobility_solver::free_host_memory(){
  
  delete[] num_segs;
  delete[] num_blobs;

  cudaFreeHost(v_segs_host);
  cudaFreeHost(v_blobs_host);
  cudaFreeHost(x_segs_host);
  cudaFreeHost(x_blobs_host);
  cudaFreeHost(f_segs_host);
  cudaFreeHost(f_blobs_host);
  cudaFreeHost(f_blobs_repulsion_host);

}

void rpy_mobility_solver::free_device_memory(){

  for (int n = 0; n < num_gpus; n++){

    cudaSetDevice(n);
    cudaFree(v_segs_device[n]);
    cudaFree(v_blobs_device[n]);
    cudaFree(x_segs_device[n]);
    cudaFree(x_blobs_device[n]);
    cudaFree(f_segs_device[n]);
    cudaFree(f_blobs_device[n]);
    cudaFree(f_blobs_repulsion_device[n]);

  }

  delete[] v_segs_device;
  delete[] v_blobs_device;
  delete[] x_segs_device;
  delete[] x_blobs_device;
  delete[] f_segs_device;
  delete[] f_blobs_device;
  delete[] f_blobs_repulsion_device;

}

void rpy_mobility_solver::allocate_host_memory(){
  
  std::cout << std::endl << std::endl << "Running on all GPUs visible to this shell environment, as defined by the environment variable CUDA_VISIBLE_DEVICES." << std::endl;

  cudaGetDeviceCount(&num_gpus);

  std::cout <<  "Found " << num_gpus << " GPU(s)." << std::endl;

  cudaSetDevice(0);

  // Allocate pinned host memory to allow async copying
  cudaHostAlloc(&v_segs_host, 6*NSWIM*NFIL*NSEG*sizeof(Real), cudaHostAllocPortable);
  cudaHostAlloc(&v_blobs_host, 3*NSWIM*NBLOB*sizeof(Real), cudaHostAllocPortable);
  cudaHostAlloc(&x_segs_host, 3*NSWIM*NFIL*NSEG*sizeof(Real), cudaHostAllocPortable);
  cudaHostAlloc(&x_blobs_host, 3*NSWIM*NBLOB*sizeof(Real), cudaHostAllocPortable);
  cudaHostAlloc(&f_segs_host, 6*NSWIM*NFIL*NSEG*sizeof(Real), cudaHostAllocPortable);
  cudaHostAlloc(&f_blobs_host, 3*NSWIM*NBLOB*sizeof(Real), cudaHostAllocPortable);
  cudaHostAlloc(&f_blobs_repulsion_host, 3*NSWIM*NBLOB*sizeof(Real), cudaHostAllocPortable);

  num_segs = new int[num_gpus];
  num_blobs = new int[num_gpus];

  num_segs[0] = NSWIM*NFIL*NSEG;
  num_blobs[0] = NSWIM*NBLOB;

  for (int n = 1; n < num_gpus; n++){

    num_segs[n] = (NSWIM*NFIL*NSEG)/num_gpus;
    num_segs[0] -= num_segs[n];

    num_blobs[n] = (NSWIM*NBLOB)/num_gpus;
    num_blobs[0] -= num_blobs[n];

  }

}

void rpy_mobility_solver::allocate_device_memory(){
  
  v_segs_device = new Real*[num_gpus];
  v_blobs_device = new Real*[num_gpus];
  x_segs_device = new Real*[num_gpus];
  x_blobs_device = new Real*[num_gpus];
  f_segs_device = new Real*[num_gpus];
  f_blobs_device = new Real*[num_gpus];
  f_blobs_repulsion_device = new Real*[num_gpus];

  for (int n = 0; n < num_gpus; n++){

    cudaSetDevice(n);

    cudaMalloc(&v_segs_device[n], 6*num_segs[n]*sizeof(Real));
    cudaMalloc(&v_blobs_device[n], 3*num_blobs[n]*sizeof(Real));

    cudaMalloc(&x_segs_device[n], 3*NSWIM*NFIL*NSEG*sizeof(Real));
    cudaMalloc(&x_blobs_device[n], 3*NSWIM*NBLOB*sizeof(Real));

    cudaMalloc(&f_segs_device[n], 6*NSWIM*NFIL*NSEG*sizeof(Real));
    cudaMalloc(&f_blobs_device[n], 3*NSWIM*NBLOB*sizeof(Real));
    cudaMalloc(&f_blobs_repulsion_device[n], 3*num_blobs[n]*sizeof(Real));

  }

}

void rpy_mobility_solver::copy_segment_positions_to_device(){

  for (int n = 0; n < num_gpus; n++){

    cudaSetDevice(n);
    cudaMemcpyAsync(x_segs_device[n], x_segs_host, 3*NSWIM*NFIL*NSEG*sizeof(Real), cudaMemcpyHostToDevice);

  }

}

void rpy_mobility_solver::copy_segment_forces_to_device(){

  for (int n = 0; n < num_gpus; n++){

    cudaSetDevice(n);
    cudaMemcpyAsync(f_segs_device[n], f_segs_host, 6*NSWIM*NFIL*NSEG*sizeof(Real), cudaMemcpyHostToDevice);

  }

}

void rpy_mobility_solver::copy_blob_positions_to_device(){

  for (int n = 0; n < num_gpus; n++){

    cudaSetDevice(n);
    cudaMemcpyAsync(x_blobs_device[n], x_blobs_host, 3*NSWIM*NBLOB*sizeof(Real), cudaMemcpyHostToDevice);

  }

}

void rpy_mobility_solver::copy_blob_forces_to_device(){

  for (int n = 0; n < num_gpus; n++){

    cudaSetDevice(n);
    cudaMemcpyAsync(f_blobs_device[n], f_blobs_host, 3*NSWIM*NBLOB*sizeof(Real), cudaMemcpyHostToDevice);

  }

}

void rpy_mobility_solver::copy_interparticle_blob_forces_to_host(){

  int start_blob = 0;

  for (int n = 0; n < num_gpus; n++){

    cudaSetDevice(n);
    cudaMemcpyAsync(&f_blobs_repulsion_host[3*start_blob], f_blobs_repulsion_device[n], 3*num_blobs[n]*sizeof(Real), cudaMemcpyDeviceToHost);
    start_blob += num_blobs[n];

  }

}

void rpy_mobility_solver::copy_blob_velocities_to_host(){

  int start_blob = 0;

  for (int n = 0; n < num_gpus; n++){

    cudaSetDevice(n);
    cudaMemcpyAsync(&v_blobs_host[3*start_blob], v_blobs_device[n], 3*num_blobs[n]*sizeof(Real), cudaMemcpyDeviceToHost);
    start_blob += num_blobs[n];

  }

}

void rpy_mobility_solver::copy_segment_velocities_to_host(){

  int start_seg = 0;

  for (int n = 0; n < num_gpus; n++){

    cudaSetDevice(n);

    cudaMemcpyAsync(&v_segs_host[6*start_seg], v_segs_device[n], 6*num_segs[n]*sizeof(Real), cudaMemcpyDeviceToHost);

    start_seg += num_segs[n];

  }

}

void rpy_mobility_solver::apply_interparticle_forces(){

  #if !PRESCRIBED_CILIA

    int start_seg = 0;
    int start_blob = 0;

    for (int n = 0; n < num_gpus; n++){

      cudaSetDevice(n);

      const int num_thread_blocks = (std::max<int>(num_segs[n], num_blobs[n]) + THREADS_PER_BLOCK - 1)/THREADS_PER_BLOCK;
      
      barrier_forces<<<num_thread_blocks, THREADS_PER_BLOCK>>>(f_segs_device[n], f_blobs_repulsion_device[n], x_segs_device[n], x_blobs_device[n], start_seg, num_segs[n], start_blob, num_blobs[n]);

      if (num_gpus > 1){

        // In the barrier_forces(...) kernel, each GPU only evaluates barrier forces for the same subset of segments it evaluates velocities for.
        // So if we have multiple GPUs but don't copy back the forces to share, HIs between segments will only include any barrier force contributions
        // if their velocities happen to be computed on the same GPU.
        cudaMemcpyAsync(&f_segs_host[6*start_seg], &f_segs_device[n][6*start_seg], 6*num_segs[n]*sizeof(Real), cudaMemcpyDeviceToHost);

      }

      start_seg += num_segs[n];
      start_blob += num_blobs[n];

    }

    if (num_gpus > 1){

      rpy_mobility_solver::wait_for_device();
      rpy_mobility_solver::copy_segment_forces_to_device();

    }

  #endif

}

void rpy_mobility_solver::wait_for_device(){

  for (int n = 0; n < num_gpus; n++){

      cudaSetDevice(n);
      cudaDeviceSynchronize();

    }

}

void rpy_mobility_solver::evaluate_segment_segment_mobility(){

  int start_seg = 0;

    for (int n = 0; n < num_gpus; n++){

      cudaSetDevice(n);

      int num_thread_blocks = (num_segs[n] + THREADS_PER_BLOCK - 1)/THREADS_PER_BLOCK;

      Mss_mult<<<num_thread_blocks, THREADS_PER_BLOCK>>>(v_segs_device[n], f_segs_device[n], x_segs_device[n], start_seg, num_segs[n]);

      start_seg += num_segs[n];

    }

}

void rpy_mobility_solver::evaluate_segment_blob_mobility(){

  int start_seg = 0;

  for (int n = 0; n < num_gpus; n++){

    cudaSetDevice(n);

    int num_thread_blocks = (num_segs[n] + THREADS_PER_BLOCK - 1)/THREADS_PER_BLOCK;

    Msb_mult<<<num_thread_blocks, THREADS_PER_BLOCK>>>(v_segs_device[n], f_blobs_device[n], x_segs_device[n], x_blobs_device[n], start_seg, num_segs[n]);

    start_seg += num_segs[n];

  }

}

void rpy_mobility_solver::evaluate_blob_blob_mobility(){

  int start_blob = 0;

  for (int n = 0; n < num_gpus; n++){

    cudaSetDevice(n);

    int num_thread_blocks = (num_blobs[n] + THREADS_PER_BLOCK - 1)/THREADS_PER_BLOCK;

    Mbb_mult<<<num_thread_blocks, THREADS_PER_BLOCK>>>(v_blobs_device[n], f_blobs_device[n], x_blobs_device[n], start_blob, num_blobs[n]);

    start_blob += num_blobs[n];

  }

}

void rpy_mobility_solver::evaluate_blob_segment_mobility(){

  int start_blob = 0;

  for (int n = 0; n < num_gpus; n++){

    cudaSetDevice(n);

    const int num_thread_blocks = (num_blobs[n] + THREADS_PER_BLOCK - 1)/THREADS_PER_BLOCK;

    Mbs_mult<<<num_thread_blocks, THREADS_PER_BLOCK>>>(v_blobs_device[n], f_segs_device[n], x_blobs_device[n], x_segs_device[n], start_blob, num_blobs[n]);

    start_blob += num_blobs[n];

  }

}

void rpy_mobility_solver::evaluate_full_mobility(){
  evaluate_segment_segment_mobility();
  evaluate_segment_blob_mobility();
  copy_segment_velocities_to_host();

  evaluate_blob_blob_mobility();
  evaluate_blob_segment_mobility();
  copy_blob_velocities_to_host();
}