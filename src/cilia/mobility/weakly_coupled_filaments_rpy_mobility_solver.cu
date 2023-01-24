// weakly_coupled_filaments_rpy_mobility_solver.cu

#include "weakly_coupled_filaments_rpy_mobility_solver.hpp"
#include "cuda_functions.hpp"

weakly_coupled_filaments_rpy_mobility_solver::~weakly_coupled_filaments_rpy_mobility_solver(){}

weakly_coupled_filaments_rpy_mobility_solver::weakly_coupled_filaments_rpy_mobility_solver(){}

void weakly_coupled_filaments_rpy_mobility_solver::free_device_memory(){

  rpy_mobility_solver::free_device_memory();

  for (int n = 0; n < num_gpus; n++){

    cudaSetDevice(n);

    cudaFree(&f_fils_device[n]);
    cudaFree(&s_fils_device[n]);
    cudaFree(&x_fils_device[n]);

  }

  delete[] f_fils_device;
  delete[] s_fils_device;
  delete[] x_fils_device;

}

void weakly_coupled_filaments_rpy_mobility_solver::allocate_device_memory(){
  
  rpy_mobility_solver::allocate_device_memory();

  f_fils_device = new Real*[num_gpus];
  s_fils_device = new Real*[num_gpus];
  x_fils_device = new Real*[num_gpus];

  for (int n = 0; n < num_gpus; n++){

    cudaSetDevice(n);

    cudaMalloc(&f_fils_device[n], 3*NSWIM*NFIL*sizeof(Real));
    cudaMalloc(&s_fils_device[n], 1*sizeof(Real)); // cudaMalloc(&s_fils_device[n], 9*NSWIM*NFIL*sizeof(Real));
    cudaMalloc(&x_fils_device[n], 3*NSWIM*NFIL*sizeof(Real));

  }

}

void weakly_coupled_filaments_rpy_mobility_solver::copy_segment_positions_to_device(){

  for (int n = 0; n < num_gpus; n++){

    cudaSetDevice(n);
    cudaMemcpyAsync(x_segs_device[n], x_segs_host, 3*NSWIM*NFIL*NSEG*sizeof(Real), cudaMemcpyHostToDevice);

    int num_thread_blocks = (3*NSWIM*NFIL*NSEG + THREADS_PER_BLOCK - 1)/THREADS_PER_BLOCK;

    reset_mean_fil_posns<<<num_thread_blocks, THREADS_PER_BLOCK>>>(x_fils_device[n]);

    num_thread_blocks = (NSWIM*NFIL*NSEG + THREADS_PER_BLOCK - 1)/THREADS_PER_BLOCK;

    write_mean_fil_posns<<<num_thread_blocks, THREADS_PER_BLOCK>>>(x_fils_device[n], x_segs_device[n]);

  }

}

void weakly_coupled_filaments_rpy_mobility_solver::copy_segment_forces_to_device(){

  for (int n = 0; n < num_gpus; n++){

    cudaSetDevice(n);
    cudaMemcpyAsync(f_segs_device[n], f_segs_host, 6*NSWIM*NFIL*NSEG*sizeof(Real), cudaMemcpyHostToDevice);

    int num_thread_blocks = (3*NSWIM*NFIL*NSEG + THREADS_PER_BLOCK - 1)/THREADS_PER_BLOCK;

    reset_mean_fil_force_and_moment<<<num_thread_blocks, THREADS_PER_BLOCK>>>(f_fils_device[n], s_fils_device[n]);

    num_thread_blocks = (NSWIM*NFIL*NSEG + THREADS_PER_BLOCK - 1)/THREADS_PER_BLOCK;

    find_total_fil_force_and_first_moment<<<num_thread_blocks, THREADS_PER_BLOCK>>>(f_fils_device[n], s_fils_device[n], f_segs_device[n], x_segs_device[n], x_fils_device[n]);

  }

}

void weakly_coupled_filaments_rpy_mobility_solver::evaluate_segment_segment_mobility(){

  int start_seg = 0;

    for (int n = 0; n < num_gpus; n++){

      cudaSetDevice(n);

      int num_thread_blocks = (num_segs[n] + THREADS_PER_BLOCK - 1)/THREADS_PER_BLOCK;

      Msf_mult<<<num_thread_blocks, THREADS_PER_BLOCK>>>(v_segs_device[n], f_segs_device[n], f_fils_device[n], s_fils_device[n], x_segs_device[n], x_fils_device[n], start_seg, num_segs[n]);

      start_seg += num_segs[n];

    }

}
