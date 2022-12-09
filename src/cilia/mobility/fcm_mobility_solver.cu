// fcm_mobility_solver.cu

#include <iomanip>
#include <iostream>
#include <cmath>
#include <algorithm>
#include "fcm_mobility_solver.hpp"
#include "../cuda_functions.hpp"

fcm_mobility_solver::~fcm_mobility_solver(){}

fcm_mobility_solver::fcm_mobility_solver(){
  Real values[100];
  std::string fcm_folder = "../CUFCM/";
  std::vector<std::string> datafile_names{3};
  read_config(values, datafile_names, "../CUFCM/simulation_info_long");
  for(int i = 0; i < 3; i++){
    datafile_names[i] = fcm_folder + datafile_names[i];
  }

  pars.N = NSWIM*NFIL*NSEG + NSWIM*NBLOB;
  pars.rh = RSEG;
  pars.alpha = values[2];
  pars.beta = values[3];
  pars.eta = values[4];
  int npts = values[5];
  pars.nx = npts;
  pars.ny = npts;
  pars.nz = npts;
  pars.repeat = values[8];
  pars.prompt = values[9];
  pars.boxsize = values[13];
  // pars.boxsize = pars.rh/1.7724538509055159 * pars.nx / pars.alpha;

  cufcm_solver = new FCM_solver(pars);
  cufcm_solver->init_aux_for_filament();
}

void fcm_mobility_solver::free_host_memory(){
  
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

void fcm_mobility_solver::free_device_memory(){

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

void fcm_mobility_solver::allocate_host_memory(){
  
  std::cout << std::endl << std::endl << "Running on all GPUs visible to this shell environment, as defined by the environment variable CUDA_VISIBLE_DEVICES." << std::endl;

  cudaGetDeviceCount(&num_gpus);

  std::cout <<  "Found " << num_gpus << " GPU(s)." << std::endl;

  cudaSetDevice(0);

  // Allocate pinned host memory to allow async copying
  cudaHostAlloc(&v_segs_host, 6*NSWIM*NFIL*NSEG*sizeof(double), cudaHostAllocPortable);
  cudaHostAlloc(&v_blobs_host, 3*NSWIM*NBLOB*sizeof(double), cudaHostAllocPortable);
  cudaHostAlloc(&x_segs_host, 3*NSWIM*NFIL*NSEG*sizeof(double), cudaHostAllocPortable);
  cudaHostAlloc(&x_blobs_host, 3*NSWIM*NBLOB*sizeof(double), cudaHostAllocPortable);
  cudaHostAlloc(&f_segs_host, 6*NSWIM*NFIL*NSEG*sizeof(double), cudaHostAllocPortable);
  cudaHostAlloc(&f_blobs_host, 3*NSWIM*NBLOB*sizeof(double), cudaHostAllocPortable);
  cudaHostAlloc(&f_blobs_repulsion_host, 3*NSWIM*NBLOB*sizeof(double), cudaHostAllocPortable);

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

void fcm_mobility_solver::allocate_device_memory(){
  
  v_segs_device = new double*[num_gpus];
  v_blobs_device = new double*[num_gpus];
  x_segs_device = new double*[num_gpus];
  x_blobs_device = new double*[num_gpus];
  f_segs_device = new double*[num_gpus];
  f_blobs_device = new double*[num_gpus];
  f_blobs_repulsion_device = new double*[num_gpus];

  for (int n = 0; n < num_gpus; n++){

    cudaSetDevice(n);

    cudaMalloc(&v_segs_device[n], 6*num_segs[n]*sizeof(double));
    cudaMalloc(&v_blobs_device[n], 3*num_blobs[n]*sizeof(double));

    cudaMalloc(&x_segs_device[n], 3*NSWIM*NFIL*NSEG*sizeof(double));
    cudaMalloc(&x_blobs_device[n], 3*NSWIM*NBLOB*sizeof(double));

    cudaMalloc(&f_segs_device[n], 6*NSWIM*NFIL*NSEG*sizeof(double));
    cudaMalloc(&f_blobs_device[n], 3*NSWIM*NBLOB*sizeof(double));
    cudaMalloc(&f_blobs_repulsion_device[n], 3*num_blobs[n]*sizeof(double));

  }

}

void fcm_mobility_solver::copy_segment_positions_to_device(){

  for (int n = 0; n < num_gpus; n++){

    cudaSetDevice(n);
    cudaMemcpyAsync(x_segs_device[n], x_segs_host, 3*NSWIM*NFIL*NSEG*sizeof(double), cudaMemcpyHostToDevice);

  }

}

void fcm_mobility_solver::copy_segment_forces_to_device(){

  for (int n = 0; n < num_gpus; n++){

    cudaSetDevice(n);
    cudaMemcpyAsync(f_segs_device[n], f_segs_host, 6*NSWIM*NFIL*NSEG*sizeof(double), cudaMemcpyHostToDevice);

  }

}

void fcm_mobility_solver::copy_blob_positions_to_device(){

  for (int n = 0; n < num_gpus; n++){

    cudaSetDevice(n);
    cudaMemcpyAsync(x_blobs_device[n], x_blobs_host, 3*NSWIM*NBLOB*sizeof(double), cudaMemcpyHostToDevice);

  }

}

void fcm_mobility_solver::copy_blob_forces_to_device(){

  for (int n = 0; n < num_gpus; n++){

    cudaSetDevice(n);
    cudaMemcpyAsync(f_blobs_device[n], f_blobs_host, 3*NSWIM*NBLOB*sizeof(double), cudaMemcpyHostToDevice);

  }

}

void fcm_mobility_solver::copy_interparticle_blob_forces_to_host(){

  int start_blob = 0;

  for (int n = 0; n < num_gpus; n++){

    cudaSetDevice(n);
    cudaMemcpyAsync(&f_blobs_repulsion_host[3*start_blob], f_blobs_repulsion_device[n], 3*num_blobs[n]*sizeof(double), cudaMemcpyDeviceToHost);
    start_blob += num_blobs[n];

  }

}

void fcm_mobility_solver::copy_blob_velocities_to_host(){

  int start_blob = 0;

  for (int n = 0; n < num_gpus; n++){

    cudaSetDevice(n);
    cudaMemcpyAsync(&v_blobs_host[3*start_blob], v_blobs_device[n], 3*num_blobs[n]*sizeof(double), cudaMemcpyDeviceToHost);
    start_blob += num_blobs[n];

  }

}

void fcm_mobility_solver::copy_segment_velocities_to_host(){

  int start_seg = 0;

  for (int n = 0; n < num_gpus; n++){

    cudaSetDevice(n);

    cudaMemcpyAsync(&v_segs_host[6*start_seg], v_segs_device[n], 6*num_segs[n]*sizeof(double), cudaMemcpyDeviceToHost);

    start_seg += num_segs[n];

  }

}

void fcm_mobility_solver::apply_interparticle_forces(){

  // #if !PRESCRIBED_CILIA

  //   int start_seg = 0;
  //   int start_blob = 0;

  //   for (int n = 0; n < num_gpus; n++){

  //     cudaSetDevice(n);

  //     const int num_thread_blocks = (std::max<int>(num_segs[n], num_blobs[n]) + THREADS_PER_BLOCK - 1)/THREADS_PER_BLOCK;
  //     barrier_forces<<<num_thread_blocks, THREADS_PER_BLOCK>>>(f_segs_device[n], f_blobs_repulsion_device[n], x_segs_device[n], x_blobs_device[n], start_seg, num_segs[n], start_blob, num_blobs[n]);

  //     start_seg += num_segs[n];
  //     start_blob += num_blobs[n];

  //   }

  // #endif

}

void fcm_mobility_solver::wait_for_device(){

  for (int n = 0; n < num_gpus; n++){

      cudaSetDevice(n);
      cudaDeviceSynchronize();

    }

}

void fcm_mobility_solver::evaluate_segment_segment_mobility(){

  int start_seg = 0;

  for (int n = 0; n < num_gpus; n++){

    cudaSetDevice(n);

    int num_thread_blocks = (num_segs[n] + THREADS_PER_BLOCK - 1)/THREADS_PER_BLOCK;

    // Mss_mult<<<num_thread_blocks, THREADS_PER_BLOCK>>>(v_segs_device[n], f_segs_device[n], x_segs_device[n], start_seg, num_segs[n]);

    cufcm_solver->reform_data(x_segs_device[n], f_segs_device[n], v_segs_device[n],
                              x_blobs_device[n], f_blobs_device[n], v_blobs_device[n],
                              num_segs[n], num_blobs[n]);
                    
    cufcm_solver->Mss();

    cufcm_solver->reform_data_back(x_segs_device[n], f_segs_device[n], v_segs_device[n],
                                   x_blobs_device[n], f_blobs_device[n], v_blobs_device[n],
                                   num_segs[n], num_blobs[n]);

    start_seg += num_segs[n];

  }

}

void fcm_mobility_solver::evaluate_segment_blob_mobility(){

  printf("NOT IMPLEMENTED ERROR\n");

  int start_seg = 0;

  for (int n = 0; n < num_gpus; n++){

    cudaSetDevice(n);

    int num_thread_blocks = (num_segs[n] + THREADS_PER_BLOCK - 1)/THREADS_PER_BLOCK;

    Msb_mult<<<num_thread_blocks, THREADS_PER_BLOCK>>>(v_segs_device[n], f_blobs_device[n], x_segs_device[n], x_blobs_device[n], start_seg, num_segs[n]);

    start_seg += num_segs[n];

  }

}

void fcm_mobility_solver::evaluate_blob_blob_mobility(){

  int start_blob = 0;

  for (int n = 0; n < num_gpus; n++){

    cudaSetDevice(n);

    int num_thread_blocks = (num_segs[n] + THREADS_PER_BLOCK - 1)/THREADS_PER_BLOCK;

    cufcm_solver->reform_data(x_segs_device[n], f_segs_device[n], v_segs_device[n],
                              x_blobs_device[n], f_blobs_device[n], v_blobs_device[n],
                              num_segs[n], num_blobs[n]);
                    
    cufcm_solver->Mss();

    cufcm_solver->reform_data_back(x_segs_device[n], f_segs_device[n], v_segs_device[n],
                                   x_blobs_device[n], f_blobs_device[n], v_blobs_device[n],
                                   num_segs[n], num_blobs[n]);

    start_blob += num_blobs[n];

  }

}

void fcm_mobility_solver::evaluate_blob_segment_mobility(){

  printf("NOT IMPLEMENTED ERROR\n");

  int start_blob = 0;

  for (int n = 0; n < num_gpus; n++){

    cudaSetDevice(n);

    const int num_thread_blocks = (num_blobs[n] + THREADS_PER_BLOCK - 1)/THREADS_PER_BLOCK;

    Mbs_mult<<<num_thread_blocks, THREADS_PER_BLOCK>>>(v_blobs_device[n], f_segs_device[n], x_blobs_device[n], x_segs_device[n], start_blob, num_blobs[n]);

    start_blob += num_blobs[n];

  }

}
