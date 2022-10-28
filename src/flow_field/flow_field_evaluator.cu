// flow_field_evaluator.cu

#include <iostream>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include "quaternion.hpp"
#include "flow_field_evaluator.hpp"

#define PI 3.14159265358979323846264338327950288
#define THREADS_PER_BLOCK 64

__global__ void flow_velocities(double *const __restrict__ Vf, const double *const __restrict__ Xf, const double *const __restrict__ Fs,
                                  const double *const __restrict__ Fb, const double *const __restrict__ Xs, const double *const __restrict__ Xbody, const double *const __restrict__ Xb,
                                  const int Nf, const int Ns, const int Nbods, const int Nb, const double prefac, const int remove_net_force){

  const int index = threadIdx.x + blockIdx.x*blockDim.x;
  const int stride = blockDim.x*gridDim.x;

  // Declare the shared memory for this thread block
  __shared__ double x_shared[THREADS_PER_BLOCK];
  __shared__ double y_shared[THREADS_PER_BLOCK];
  __shared__ double z_shared[THREADS_PER_BLOCK];
  __shared__ double fx_shared[THREADS_PER_BLOCK];
  __shared__ double fy_shared[THREADS_PER_BLOCK];
  __shared__ double fz_shared[THREADS_PER_BLOCK];
  __shared__ double taux_shared[THREADS_PER_BLOCK];
  __shared__ double tauy_shared[THREADS_PER_BLOCK];
  __shared__ double tauz_shared[THREADS_PER_BLOCK];

  double vx, vy, vz;
  double xi, yi, zi;

  // Stay in the loop as long as any thread in the block still needs to compute velocities.
  for (int i = index; (i-threadIdx.x) < Nf; i += stride){

    if (i < Nf){

      vx = 0.0; vy = 0.0; vz = 0.0;
      xi = Xf[3*i]; yi = Xf[3*i + 1]; zi = Xf[3*i + 2];

    }

    // Loop over the segments
    for (int j_start = 0; j_start < Ns*Nbods; j_start += THREADS_PER_BLOCK){

      const int j_to_read = j_start + threadIdx.x;

      if (j_to_read < Ns*Nbods){

        x_shared[threadIdx.x] = Xs[3*j_to_read];
        y_shared[threadIdx.x] = Xs[3*j_to_read + 1];
        z_shared[threadIdx.x] = Xs[3*j_to_read + 2];
        fx_shared[threadIdx.x] = Fs[6*j_to_read];
        fy_shared[threadIdx.x] = Fs[6*j_to_read + 1];
        fz_shared[threadIdx.x] = Fs[6*j_to_read + 2];
        taux_shared[threadIdx.x] = Fs[6*j_to_read + 3];
        tauy_shared[threadIdx.x] = Fs[6*j_to_read + 4];
        tauz_shared[threadIdx.x] = Fs[6*j_to_read + 5];

      }

      __syncthreads();

      if (i < Nf){

        for (int j = 0; (j < THREADS_PER_BLOCK) && (j_start + j < Ns*Nbods); j++){

          // Calculate the flow velocity
          double xdiff = xi - x_shared[j];
          double ydiff = yi - y_shared[j];
          double zdiff = zi - z_shared[j];

          // Stokeslet part
          double rm2 = 1.0/(xdiff*xdiff + ydiff*ydiff + zdiff*zdiff);
          double rm1 = sqrt(rm2);
          double rm2_times_r_dot_f = rm2*(xdiff*fx_shared[j] + ydiff*fy_shared[j] + zdiff*fz_shared[j]);

          vx += rm1*(fx_shared[j] + xdiff*rm2_times_r_dot_f);
          vy += rm1*(fy_shared[j] + ydiff*rm2_times_r_dot_f);
          vz += rm1*(fz_shared[j] + zdiff*rm2_times_r_dot_f);

          if (Nb == 0){

            // If there are no blobs, we assume the simulation is in a half-space and hence
            // the Stokeslet must be corrected to the Blakelet.

            // First, remove a Stokeslet contribution at the image location
            const double Rz = zi + z_shared[j];
            const double R2 = xdiff*xdiff + ydiff*ydiff + Rz*Rz;
            const double Rm1 = sqrt(1.0/R2);
            const double Rm2_times_R_dot_f = (xdiff*fx_shared[j] + ydiff*fy_shared[j] + Rz*fz_shared[j])*Rm1*Rm1;

            vx -= Rm1*(fx_shared[j] + xdiff*Rm2_times_R_dot_f);
            vy -= Rm1*(fy_shared[j] + ydiff*Rm2_times_R_dot_f);
            vz -= Rm1*(fz_shared[j] + Rz*Rm2_times_R_dot_f);

            // Then add the higher-order correction, which doesn't depend on the z-component of the force
            const double Dxx = (3.0*xdiff*xdiff - R2)*zi;
            const double Dxy = 3.0*xdiff*ydiff*zi;
            const double Dxz = (R2 - 3.0*Rz*zi)*xdiff;

            // Dyx = Dxy so no new term
            const double Dyy = (3.0*ydiff*ydiff - R2)*zi;
            const double Dyz = (R2 - 3.0*Rz*zi)*ydiff;

            const double Dzx = (3.0*Rz*zi + R2)*xdiff;
            const double Dzy = (3.0*Rz*zi + R2)*ydiff;
            const double Dzz = (R2 - 3.0*Rz*Rz)*zi;

            const double Dfac = 2.0*z_shared[j]*Rm1*Rm1*Rm1*Rm1*Rm1;
            vx += Dfac*(Dxx*fx_shared[j] + Dxy*fy_shared[j] + Dxz*fz_shared[j]);
            vy += Dfac*(Dxy*fx_shared[j] + Dyy*fy_shared[j] + Dyz*fz_shared[j]);
            vz += Dfac*(Dzx*fx_shared[j] + Dzy*fy_shared[j] + Dzz*fz_shared[j]);

          }

          // Rotlet part
          const double rm3 = rm1*rm2;

          vx += rm3*(tauy_shared[j]*zdiff - tauz_shared[j]*ydiff);
          vy += rm3*(tauz_shared[j]*xdiff - taux_shared[j]*zdiff);
          vz += rm3*(taux_shared[j]*ydiff - tauy_shared[j]*xdiff);

          if (remove_net_force){

            // Which body does this segment belong to?
            const int body_id = (j_start + j)/Ns;

            double xdiff = xi - Xbody[3*body_id];
            double ydiff = yi - Xbody[3*body_id + 1];
            double zdiff = zi - Xbody[3*body_id + 2];

            // Stokeslet part
            double rm2 = 1.0/(xdiff*xdiff + ydiff*ydiff + zdiff*zdiff);
            double rm1 = sqrt(rm2);
            double rm2_times_r_dot_f = rm2*(xdiff*fx_shared[j] + ydiff*fy_shared[j] + zdiff*fz_shared[j]);

            vx -= rm1*(fx_shared[j] + xdiff*rm2_times_r_dot_f);
            vy -= rm1*(fy_shared[j] + ydiff*rm2_times_r_dot_f);
            vz -= rm1*(fz_shared[j] + zdiff*rm2_times_r_dot_f);

          }

        }

      }

      __syncthreads();

    }

    // Loop over the blobs
    for (int j_start = 0; j_start < Nb*Nbods; j_start += THREADS_PER_BLOCK){

      const int j_to_read = j_start + threadIdx.x;

      if (j_to_read < Nb*Nbods){

        x_shared[threadIdx.x] = Xb[3*j_to_read];
        y_shared[threadIdx.x] = Xb[3*j_to_read + 1];
        z_shared[threadIdx.x] = Xb[3*j_to_read + 2];
        fx_shared[threadIdx.x] = Fb[3*j_to_read];
        fy_shared[threadIdx.x] = Fb[3*j_to_read + 1];
        fz_shared[threadIdx.x] = Fb[3*j_to_read + 2];

      }

      __syncthreads();

      if (i < Nf){

        for (int j = 0; (j < THREADS_PER_BLOCK) && (j_start + j < Nb*Nbods); j++){

          // Calculate the flow velocity
          double xdiff = xi - x_shared[j];
          double ydiff = yi - y_shared[j];
          double zdiff = zi - z_shared[j];

          // Only have the Stokeslet part for blobs
          double rm2 = 1.0/(xdiff*xdiff + ydiff*ydiff + zdiff*zdiff);
          double rm1 = sqrt(rm2);
          double rm2_times_r_dot_f = rm2*(xdiff*fx_shared[j] + ydiff*fy_shared[j] + zdiff*fz_shared[j]);

          vx += rm1*(fx_shared[j] + xdiff*rm2_times_r_dot_f);
          vy += rm1*(fy_shared[j] + ydiff*rm2_times_r_dot_f);
          vz += rm1*(fz_shared[j] + zdiff*rm2_times_r_dot_f);

          if (remove_net_force){

            // Which body does this blob belong to?
            const int body_id = (j_start + j)/Nb;

            double xdiff = xi - Xbody[3*body_id];
            double ydiff = yi - Xbody[3*body_id + 1];
            double zdiff = zi - Xbody[3*body_id + 2];

            // Stokeslet part
            double rm2 = 1.0/(xdiff*xdiff + ydiff*ydiff + zdiff*zdiff);
            double rm1 = sqrt(rm2);
            double rm2_times_r_dot_f = rm2*(xdiff*fx_shared[j] + ydiff*fy_shared[j] + zdiff*fz_shared[j]);

            vx -= rm1*(fx_shared[j] + xdiff*rm2_times_r_dot_f);
            vy -= rm1*(fy_shared[j] + ydiff*rm2_times_r_dot_f);
            vz -= rm1*(fz_shared[j] + zdiff*rm2_times_r_dot_f);

          }

        }

      }

      __syncthreads();

    }

    // Write the results and finish
    if (i < Nf){

      Vf[3*i] = prefac*vx;
      Vf[3*i + 1] = prefac*vy;
      Vf[3*i + 2] = prefac*vz;

    }

  }

};

flow_field_evaluator::~flow_field_evaluator(){

  // Free all memory
  delete[] v_flow_device;
  delete[] x_flow_device;
  delete[] x_segs_device;
  delete[] x_blobs_device;
  delete[] f_segs_device;
  delete[] f_blobs_device;
  delete[] num_flow_points;

  cudaFreeHost(v_flow_host);
  cudaFreeHost(x_flow_host);
  cudaFreeHost(x_segs_host);
  cudaFreeHost(x_blobs_host);
  cudaFreeHost(f_segs_host);
  cudaFreeHost(f_blobs_host);

  for (int n = 0; n < num_gpus; n++){

    cudaSetDevice(n);

    cudaFree(v_flow_device[n]);
    cudaFree(x_flow_device[n]);
    cudaFree(x_segs_device[n]);
    cudaFree(x_blobs_device[n]);
    cudaFree(f_segs_device[n]);
    cudaFree(f_blobs_device[n]);

  }

};

flow_field_evaluator::flow_field_evaluator(){};

void flow_field_evaluator::setup(FILE *fil_ref_file, FILE *blob_ref_file){

  std::cout << std::endl << std::endl << "Running on all GPUs visible to this shell environment, as defined by the environment variable CUDA_VISIBLE_DEVICES." << std::endl;
  cudaGetDeviceCount(&num_gpus);
  std::cout <<  "Found " << num_gpus << " GPU(s)." << std::endl;

  cudaSetDevice(0);

  // Allocate pinned host memory to allow async copying
  cudaHostAlloc(&v_flow_host, 3*NFLOW*sizeof(double), cudaHostAllocPortable);
  cudaHostAlloc(&x_flow_host, 3*NFLOW*sizeof(double), cudaHostAllocPortable);
  cudaHostAlloc(&x_segs_host, 3*NBOD*NFIL*NSEG*sizeof(double), cudaHostAllocPortable);
  cudaHostAlloc(&x_bods_host, 3*NBOD*sizeof(double), cudaHostAllocPortable);
  cudaHostAlloc(&x_blobs_host, 3*NBOD*NBLOB*sizeof(double), cudaHostAllocPortable);
  cudaHostAlloc(&f_segs_host, 6*NBOD*NFIL*NSEG*sizeof(double), cudaHostAllocPortable);
  cudaHostAlloc(&f_blobs_host, 3*NBOD*NBLOB*sizeof(double), cudaHostAllocPortable);

  v_flow_device = new double*[num_gpus];
  x_flow_device = new double*[num_gpus];
  x_segs_device = new double*[num_gpus];
  x_bods_device = new double*[num_gpus];
  x_blobs_device = new double*[num_gpus];
  f_segs_device = new double*[num_gpus];
  f_blobs_device = new double*[num_gpus];

  num_flow_points = new int[num_gpus];
  num_flow_points[0] = NFLOW;

  for (int n = 1; n < num_gpus; n++){

    num_flow_points[n] = NFLOW/num_gpus;
    num_flow_points[0] -= num_flow_points[n];

  }

  // Allocate device memory
  for (int n = 0; n < num_gpus; n++){

    cudaSetDevice(n);

    cudaMalloc(&v_flow_device[n], 3*num_flow_points[n]*sizeof(double));
    cudaMalloc(&x_flow_device[n], 3*num_flow_points[n]*sizeof(double));

    cudaMalloc(&x_segs_device[n], 3*NBOD*NFIL*NSEG*sizeof(double));
    cudaMalloc(&x_bods_device[n], 3*NBOD*sizeof(double));
    cudaMalloc(&x_blobs_device[n], 3*NBOD*NBLOB*sizeof(double));

    cudaMalloc(&f_segs_device[n], 6*NBOD*NFIL*NSEG*sizeof(double));
    cudaMalloc(&f_blobs_device[n], 3*NBOD*NBLOB*sizeof(double));

  }

  fil_refs = matrix(3, NFIL);
  for (int n = 0; n < NFIL; n++){

    std::fscanf(fil_ref_file, "%lf", &fil_refs(0,n));
    std::fscanf(fil_ref_file, "%lf", &fil_refs(1,n));
    std::fscanf(fil_ref_file, "%lf", &fil_refs(2,n));

  }

  if (NBLOB != 0){

    blob_refs = matrix(3, NBLOB);
    for (int n = 0; n < NBLOB; n++){

      std::fscanf(blob_ref_file, "%lf", &blob_refs(0,n));
      std::fscanf(blob_ref_file, "%lf", &blob_refs(1,n));
      std::fscanf(blob_ref_file, "%lf", &blob_refs(2,n));

    }

  }

  // Provide initial flow positions
  const double DL = 2.2*RSEG;
  const double r = 2.3*DL*NSEG;
  double phi = 0.0;

  const double b = DL; // Only for planar spiral seeding

  const double x_length = 36.0*DL*NSEG; // These are for rectangular planar seeding
  const double y_length = x_length;
  const int num_flow_x = int(std::sqrt(double(NFLOW)));
  const int num_flow_y = std::max<int>(1, int(ceil(NFLOW/double(num_flow_x))));

  for (int n = 0; n < NFLOW; n++){

    /* // Spiral seeding on a sphere of radius r:
    x_flow_host[3*n + 2] = (NFLOW == 1) ? -1.0 : 2.0*n/double(NFLOW-1) - 1.0;

    double r_local = std::sqrt(1.0 - x_flow_host[3*n + 2]*x_flow_host[3*n + 2]);

    if ((n == 0) || (n == NFLOW-1)){

      phi = 0.0;

    } else {

      phi += 3.6/(r_local*std::sqrt(NFLOW));

    }

    x_flow_host[3*n] = r*r_local*std::cos(phi);
    x_flow_host[3*n + 1] = r*r_local*std::sin(phi);
    x_flow_host[3*n + 2] *= r;

    */

    /* // Spiral seeding in a plane, starting at distance r from the swimmer's centre:
    const double r_n = r + b*phi;

    x_flow_host[3*n] = r_n*std::cos(phi);
    x_flow_host[3*n + 1] = 0.0;
    x_flow_host[3*n + 2] = r_n*std::sin(phi);

    phi += 2.0*PI*b/r_n;

    */

    // Rectangular seeding in a plane. This is primarily for producing time-averaged flow speed images
    // and so will not worry about tracers being inside the swimmer -- speeds at these locations can simply be ignored.
    x_flow_host[3*n] = 0.0; //(double(n - num_flow_x*std::floor(double(n)/double(num_flow_x)))/double(num_flow_x) - 0.5)*x_length;
    x_flow_host[3*n + 1] = -2.2032e+03 + 8.0*(n - 551.0*std::floor(n/551.0)); //0.0;
    x_flow_host[3*n + 2] = 8.0*(1.0 + std::floor(n/551.0)); //(std::floor(double(n)/double(num_flow_x))/double(num_flow_y) - 0.5)*y_length;


    // Fill velocity array with zeros initially so the update does nothing at the start of the first step.
    v_flow_host[3*n] = 0.0;
    v_flow_host[3*n + 1] = 0.0;
    v_flow_host[3*n + 2] = 0.0;

  }

};

void flow_field_evaluator::read_positions_and_forces(FILE *body_state_file, FILE *seg_state_file, FILE *seg_force_file, FILE *blob_force_file){

  matrix X(3,1); // Body position
  quaternion q; // Body orientation

  int i = 0; // Entry in segment position arrays
  int j = 0; // Entry in segment force arrays
  int k = 0; // Entry for blob arrays

  const double DL = 2.2*RSEG;

  // Write the segment and blob data to the host arrays
  for (int b = 0; b < NBOD; b++){

    std::fscanf(body_state_file, "%lf", &X(0));
    std::fscanf(body_state_file, "%lf", &X(1));
    std::fscanf(body_state_file, "%lf", &X(2));

    x_bods_host[3*b] = X(0);
    x_bods_host[3*b + 1] = X(1);
    x_bods_host[3*b + 2] = X(2);

    std::fscanf(body_state_file, "%lf", &q.scalar_part);
    std::fscanf(body_state_file, "%lf", &q.vector_part[0]);
    std::fscanf(body_state_file, "%lf", &q.vector_part[1]);
    std::fscanf(body_state_file, "%lf", &q.vector_part[2]);

    const matrix R = q.rot_mat();

    const matrix fil_base_diffs = R*fil_refs;

    for (int f = 0; f < NFIL; f++){

      // Read the 6 force and torque values per segment straight into the CUDA array
      for (int s = 0; s < 6*NSEG; s++){
        std::fscanf(seg_force_file, "%lf", &f_segs_host[j++]);
      }

      if (QUATERNION_BASED_DATA){

        matrix seg_pos = X + fil_base_diffs.get_col(f);

        x_segs_host[i++] = seg_pos(0);
        x_segs_host[i++] = seg_pos(1);
        x_segs_host[i++] = seg_pos(2);

        matrix t1(3,1), t2(3,1);

        quaternion qseg;

        std::fscanf(seg_state_file, "%lf", &qseg.scalar_part);
        std::fscanf(seg_state_file, "%lf", &qseg.vector_part[0]);
        std::fscanf(seg_state_file, "%lf", &qseg.vector_part[1]);
        std::fscanf(seg_state_file, "%lf", &qseg.vector_part[2]);

        qseg.tangent(t1);

        for (int s = 1; s < NSEG; s++){

          std::fscanf(seg_state_file, "%lf", &qseg.scalar_part);
          std::fscanf(seg_state_file, "%lf", &qseg.vector_part[0]);
          std::fscanf(seg_state_file, "%lf", &qseg.vector_part[1]);
          std::fscanf(seg_state_file, "%lf", &qseg.vector_part[2]);

          qseg.tangent(t2);

          seg_pos += 0.5*DL*(t1 + t2);

          x_segs_host[i++] = seg_pos(0);
          x_segs_host[i++] = seg_pos(1);
          x_segs_host[i++] = seg_pos(2);

          t1 = t2;

        }

      } else {

        for (int s = 0; s < 3*NSEG; s++){
          std::fscanf(seg_state_file, "%lf", &x_segs_host[i++]);
        }

      }

    }

    if (NBLOB != 0){

      const matrix blob_diffs = R*blob_refs;

      for (int bl = 0; bl < NBLOB; bl++){

        const matrix blob_pos = X + blob_diffs.get_col(bl);

        x_blobs_host[k] = blob_pos(0);
        std::fscanf(blob_force_file, "%lf", &f_blobs_host[k++]);
        x_blobs_host[k] = blob_pos(1);
        std::fscanf(blob_force_file, "%lf", &f_blobs_host[k++]);
        x_blobs_host[k] = blob_pos(2);
        std::fscanf(blob_force_file, "%lf", &f_blobs_host[k++]);

      }

    }

  }

  // Send the segment and blob data to the GPUs
  for (int n = 0; n < num_gpus; n++){

    cudaSetDevice(n);
    cudaMemcpyAsync(x_segs_device[n], x_segs_host, 3*NBOD*NFIL*NSEG*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpyAsync(x_bods_device[n], x_bods_host, 3*NBOD*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpyAsync(f_segs_device[n], f_segs_host, 6*NBOD*NFIL*NSEG*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpyAsync(x_blobs_device[n], x_blobs_host, 3*NBOD*NBLOB*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpyAsync(f_blobs_device[n], f_blobs_host, 3*NBOD*NBLOB*sizeof(double), cudaMemcpyHostToDevice);

  }

  if (ADVANCE_TRACER_POSITIONS){

    // While that data transfers, determine the new positions to evaluate the flow at
    for (int n = 0; n < NFLOW; n++){

      // Euler's method
      x_flow_host[3*n] += DT*v_flow_host[3*n];
      x_flow_host[3*n + 1] += DT*v_flow_host[3*n + 1];
      x_flow_host[3*n + 2] += DT*v_flow_host[3*n + 2];

    }

  }

  // Transfer the flow evaluation points
  int flow_points_offset = 0;
  for (int n = 0; n < num_gpus; n++){

    cudaSetDevice(n);
    cudaMemcpyAsync(x_flow_device[n], &x_flow_host[3*flow_points_offset], 3*num_flow_points[n]*sizeof(double), cudaMemcpyHostToDevice);
    flow_points_offset += num_flow_points[n];

  }

};

void flow_field_evaluator::velocities(){

  // Launch the kernel
  const double prefac = 1.0/(8.0*PI*MU);

  for (int n = 0; n < num_gpus; n++){

    cudaSetDevice(n);
    const int num_thread_blocks = (num_flow_points[n] + THREADS_PER_BLOCK - 1)/THREADS_PER_BLOCK;
    flow_velocities<<<num_thread_blocks, THREADS_PER_BLOCK>>>(v_flow_device[n], x_flow_device[n], f_segs_device[n],
                                      f_blobs_device[n], x_segs_device[n], x_bods_device[n], x_blobs_device[n],
                                      num_flow_points[n], NFIL*NSEG, NBOD, NBLOB, prefac, REMOVE_NET_FORCE_STOKESLET);

  }

  // Copy the flow velocities back to the host
  int flow_points_offset = 0;
  for (int n = 0; n < num_gpus; n++){

    cudaSetDevice(n);
    cudaMemcpyAsync(&v_flow_host[3*flow_points_offset], v_flow_device[n], 3*num_flow_points[n]*sizeof(double), cudaMemcpyDeviceToHost);
    flow_points_offset += num_flow_points[n];

  }

  // Make sure we've got them before we move on
  for (int n = 0; n < num_gpus; n++){

    cudaSetDevice(n);
    cudaDeviceSynchronize();

  }

};
