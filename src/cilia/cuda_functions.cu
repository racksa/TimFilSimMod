// cuda_functions.cu

#include "cuda_functions.hpp"
#include <stdio.h>

// RPY kernels

template<int vel_dof, int force_dof>
__device__ void rpy_interaction(Real *const V, const Real *const F, const int i, const Real xi, const Real yi, const Real zi, const Real ai, const int j, const Real xj, const Real yj, const Real zj, const Real aj){

  #if INFINITE_PLANE_WALL

    // Calculates the velocities in a domain bounded by an infinite no-slip wall at z = 0.
    // It does so using the expressions given by Swan and Brady (2007).
    // NOTE: There seems to be a typo/sign error/unusual convention in the alternating
    // tensors that appear in the appendices of Swan and Brady's work. To get the correct
    // behaviour, it suffices to change the sign of all rot-trans (and hence trans-rot) terms.
    // This is highlighted in the code by evaluating these terms as they appear in the paper
    // but using -= rather than += to evaluate their contributions.

    // N.B. This function assumes that MU = 1 and ai = aj = 1.

    if (i == j){

      // Self-mobility
      const Real hm1 = 1.0/zi;
      const Real hm3 = hm1*hm1*hm1;
      const Real hm4 = hm3*hm1;
      const Real hm5 = hm4*hm1;

      Real a0 = 0.053051647697298 - (9.0*hm1 - 2.0*hm3 + hm5)*0.003315727981081; // 0.053051647697298 = 1/(6*PI), 0.003315727981081 = 1/(6*PI*16)
      Real a1 = 0.053051647697298 - (9.0*hm1 - 4.0*hm3 + hm5)*0.006631455962162; // 0.006631455962162 = 1/(6*PI*8)
      const Real a2 = 0.004973591971622*hm4; // 0.004973591971622 = 3/(6*PI*32)

      V[0] += a0*F[0];
      V[1] += a0*F[1];
      V[2] += a1*F[2];

      if (force_dof > 3){

        V[0] -= -a2*F[4];
        V[1] -= a2*F[3];

      }

      if (vel_dof > 3){

        V[3] -= a2*F[1];
        V[4] -= -a2*F[0];

        if (force_dof > 3){

          a0 = 0.039788735772974 - 0.012433979929054*hm3; // 0.039788735772974 = 1/(8*PI), 0.012433979929054 = 15/(6*PI*64)
          a1 = 0.039788735772974 - 0.004973591971622*hm3;

          V[3] += a0*F[3];
          V[4] += a0*F[4];
          V[5] += a1*F[5];

        }

      }

    } else {

      Real rm1 = rsqrt((xi - xj)*(xi - xj) + (yi - yj)*(yi - yj) + (zi - zj)*(zi - zj));

      // Begin with the unbounded terms
      Real rhatx = (xi - xj)*rm1;
      Real rhaty = (yi - yj)*rm1;
      Real rhatz = (zi - zj)*rm1;
      Real rhat_dot_force = rhatx*F[0] + rhaty*F[1] + rhatz*F[2];
      Real rhat_dot_torque;

      Real a0 = 0.039788735772974*rm1; // 0.039788735772974 = 1/(8*PI)
      Real a1 = 1.0 + 0.666666666666667*rm1*rm1;
      Real a2 = 1.0 - 2.0*rm1*rm1;

      V[0] += a0*(a1*F[0] + a2*rhat_dot_force*rhatx);
      V[1] += a0*(a1*F[1] + a2*rhat_dot_force*rhaty);
      V[2] += a0*(a1*F[2] + a2*rhat_dot_force*rhatz);

      if ((force_dof > 3) || (vel_dof > 3)){

        a0 *= rm1; // = 0.039788735772974*rm1*rm1

      }

      if (force_dof > 3){

        rhat_dot_torque = rhatx*F[3] + rhaty*F[4] + rhatz*F[5];

        V[0] += a0*(F[4]*rhatz - F[5]*rhaty);
        V[1] += a0*(F[5]*rhatx - F[3]*rhatz);
        V[2] += a0*(F[3]*rhaty - F[4]*rhatx);

      }

      if (vel_dof > 3){

        V[3] += a0*(F[1]*rhatz - F[2]*rhaty);
        V[4] += a0*(F[2]*rhatx - F[0]*rhatz);
        V[5] += a0*(F[0]*rhaty - F[1]*rhatx);

        if (force_dof > 3){

          a0 *= 0.5*rm1; // = 0.019894367886487*rm1*rm1*rm1, where 0.019894367886487 = 1/(16*PI)

          V[3] += a0*(3.0*rhat_dot_torque*rhatx - F[3]);
          V[4] += a0*(3.0*rhat_dot_torque*rhaty - F[4]);
          V[5] += a0*(3.0*rhat_dot_torque*rhatz - F[5]);

        }

      }

      // Now for the wall-induced corrections
      rm1 = rsqrt((xi - xj)*(xi - xj) + (yi - yj)*(yi - yj) + (zi + zj)*(zi + zj));
      Real rm2 = rm1*rm1;
      Real rm3 = rm2*rm1;
      Real rm4 = rm3*rm1;
      Real rm5 = rm4*rm1;

      Real ex = (xi - xj)*rm1;
      Real ey = (yi - yj)*rm1;
      Real ez = (zi + zj)*rm1;

      Real h = zj/(zi + zj);

      Real e_dot_force = ex*F[0] + ey*F[1] + ez*F[2];
      Real e_dot_torque;

      if (force_dof > 3){

        e_dot_torque = ex*F[3] + ey*F[4] + ez*F[5];

      }

      Real ez2 = ez*ez;

      a0 = -(0.75*rm1*(1.0 + 2.0*h*(1.0-h)*ez2) + 0.5*(rm3*(1.0 - 3.0*ez2) - rm5*(1.0 - 5.0*ez2)));
      a1 = -(0.75*rm1*(1.0 - 6.0*h*(1.0-h)*ez2) - 1.5*rm3*(1.0 - 5.0*ez2) + 2.5*rm5*(1.0 - 7.0*ez2))*e_dot_force;
      a2 = ez*F[2]*(1.5*h*rm1*(1.0 - 6.0*(1.0-h)*ez2) - 3.0*rm3*(1.0 - 5.0*ez2) + 5.0*rm5*(2.0 - 7.0*ez2));
      Real a3 = (1.5*h*rm1 - 5.0*rm5)*e_dot_force;
      Real a4 = -((3.0*h*h*rm1 + 3.0*rm3)*ez2 + (2.0 - 15.0*ez2)*rm5)*F[2];

      V[0] += 0.053051647697298*(a0*F[0] + (a1 + a2)*ex);
      V[1] += 0.053051647697298*(a0*F[1] + (a1 + a2)*ey);
      V[2] += 0.053051647697298*(a0*F[2] + (a1 + a2 + a3)*ez + a4);

      if (vel_dof > 3){

        a0 = 0.75*rm2;
        a1 = F[2]*(9.0*h*ez2*rm2 + (1.5 - 15.0*ez2)*rm4);
        a2 = -ez*e_dot_force*(4.5*h*rm2 - 7.5*rm4);
        a3 = -1.5*ez*(h*rm2 - rm4);

        V[3] -= 0.053051647697298*(a0*(F[1]*ez - F[2]*ey) - (a1+a2)*ey + a3*F[1]);
        V[4] -= 0.053051647697298*(a0*(F[2]*ex - F[0]*ez) + (a1+a2)*ex - a3*F[0]);
        V[5] -= 0.053051647697298*a0*(F[0]*ey - F[1]*ex);

        if (force_dof > 3){

          a0 = rm3*(0.375 - 2.25*ez2);
          a1 = -1.125*rm3*e_dot_torque;
          a2 = 2.25*ez*rm3*e_dot_torque;
          a3 = 2.25*rm3*(ex*F[4] - ey*F[3]);

          V[3] += 0.053051647697298*(a0*F[3]+ a1*ex - a3*ey);
          V[4] += 0.053051647697298*(a0*F[4] + a1*ey + a3*ex);
          V[5] += 0.053051647697298*(a0*F[5] + a1*ez + a2);

        }

      }

      if (force_dof > 3){

        // We use symmetry of the grand mobility matrix to calculate this contribution as
        // the torque multiplied by the transpose of the sub-block relating the angular velocity
        // of particle j to the force on particle i.

        // We do this interaction last because all others use the same h, ex and ey.

        h = 1.0 - h;
        ex *= -1.0;
        ey *= -1.0;

        const Real e_cross_T_3 = ex*F[4] - ey*F[3];
        a0 = 0.75*rm2;
        a1 = (9.0*h*ez2*rm2 + (1.5 - 15.0*ez2)*rm4)*e_cross_T_3;
        a2 = -ez*(4.5*h*rm2 - 7.5*rm4)*e_cross_T_3;
        a3 = -1.5*ez*(h*rm2 - rm4);

        V[0] -= 0.053051647697298*(a0*(F[5]*ey - F[4]*ez) + a2*ex - a3*F[4]);
        V[1] -= 0.053051647697298*(a0*(F[3]*ez - F[5]*ex) + a2*ey + a3*F[3]);
        V[2] -= 0.053051647697298*(a0*(F[4]*ex - F[3]*ey) + a2*ez + a1);

      }


    }

  #else

    // Calculates the velocities in an unbounded domain.
    Real rx = xi - xj;
    Real ry = yi - yj;
    Real rz = zi - zj;

    const Real r = sqrt(rx*rx + ry*ry + rz*rz);
    const Real rm1 = 1.0/(r + 1e-20); // We don't treat r = 0 separately.

    const Real amax_minus_amin = abs(ai - aj);

    rx *= rm1; ry *= rm1; rz *= rm1;

    Real A, B;

    // translation-translation
    Real r_dot = rx*F[0] + ry*F[1] + rz*F[2];

    if (r > ai + aj){

      A = rm1*rm1*(ai*ai + aj*aj)/3.0;
      B = 1.0 - 3.0*A;
      A += 1.0;

      Real temp = rm1/(8.0*PI*MU);

      A *= temp;
      B *= temp*r_dot;

    } else if (r > amax_minus_amin){

      Real temp = 32.0*r*r*r;

      A = amax_minus_amin*amax_minus_amin + 3.0*r*r;
      A *= -A/temp;
      A += 0.5*(ai + aj);

      B = amax_minus_amin*amax_minus_amin - r*r;
      B *= 3.0*B/temp;

      temp = 1.0/(6.0*PI*MU*ai*aj);

      A *= temp;
      B *= temp*r_dot;

    } else {

      A = 1.0/(6.0*PI*MU*((ai > aj) ? ai : aj)); // ternary operator gives largest radius
      B = 0.0;

    }

    V[0] += A*F[0] + B*rx;
    V[1] += A*F[1] + B*ry;
    V[2] += A*F[2] + B*rz;

    if (force_dof > 3){

      // translation-rotation
      if (r > ai + aj){

        A = rm1*rm1/(8.0*PI*MU);

      } else if (r > amax_minus_amin){

        A = aj - ai + r;
        A *= rm1*rm1*A*(ai*ai + 2.0*ai*(r + aj) - 3.0*(aj - r)*(aj - r))/(128.0*PI*MU*ai*aj*aj*aj);

      } else {

        A = (aj > ai)*r/(8.0*PI*MU*aj*aj*aj);

      }

      V[0] += A*(F[4]*rz - F[5]*ry);
      V[1] += A*(F[5]*rx - F[3]*rz);
      V[2] += A*(F[3]*ry - F[4]*rx);

    }

    if (vel_dof > 3){

      // rotation-translation ( = translation-rotation with ai and aj swapped; see Zuk et al. (2014))
      if (r > ai + aj){

        A = rm1*rm1/(8.0*PI*MU);

      } else if (r > amax_minus_amin){

        A = ai - aj + r;
        A *= rm1*rm1*A*(aj*aj + 2.0*aj*(r + ai) - 3.0*(ai - r)*(ai - r))/(128.0*PI*MU*aj*ai*ai*ai);

      } else {

        A = (ai > aj)*r/(8.0*PI*MU*ai*ai*ai);

      }

      V[3] += A*(F[1]*rz - F[2]*ry);
      V[4] += A*(F[2]*rx - F[0]*rz);
      V[5] += A*(F[0]*ry - F[1]*rx);

      if (force_dof > 3){

        // rotation-rotation
        r_dot = rx*F[3] + ry*F[4] + rz*F[5];

        if (r > ai + aj){

          A = -rm1*rm1*rm1/(16.0*PI*MU);
          B = -3.0*r_dot*A;

        } else if (r > amax_minus_amin){

          // Don't worry about optimising this part, we shouldn't ever actually use it anyway.

          A = 5.0*(r*r*r)*(r*r*r);
          B = (ai - aj)*(ai - aj) - r*r;

          A -= 27.0*r*r*r*r*(ai*ai + aj*aj);
          B *= 3.0*B;

          A += 32.0*(r*r*r)*((ai*ai*ai) + (aj*aj*aj));
          B *= ai*ai + aj*aj + 4.0*ai*aj - r*r;

          A -= 9.0*r*r*(ai*ai - aj*aj)*(ai*ai - aj*aj);
          A -= amax_minus_amin*amax_minus_amin*amax_minus_amin*amax_minus_amin*(ai*ai + aj*aj + 4.0*ai*aj);

          const Real temp = 1.0/(512.0*PI*MU*(r*r*r)*(ai*ai*ai)*(aj*aj*aj));

          A *= temp;
          B *= temp*r_dot;

        } else {

          const Real amax = ((ai > aj) ? ai : aj);

          A = 1.0/(8.0*PI*MU*amax*amax*amax);
          B = 0.0;

        }

        V[3] += A*F[3] + B*rx;
        V[4] += A*F[4] + B*ry;
        V[5] += A*F[5] + B*rz;

      }

    }

  #endif

}



__global__ void Mss_mult(Real * __restrict__ V, const Real *const __restrict__ F, const Real *const __restrict__ X, const int start_seg, const int num_segs){

  // Calculates the velocities of filament segments given the forces and torques
  // on the segments.

  const int index = threadIdx.x + blockIdx.x*blockDim.x;
  const int stride = blockDim.x*gridDim.x;

  // Declare the shared memory for this thread block
  __shared__ Real x_shared[THREADS_PER_BLOCK];
  __shared__ Real y_shared[THREADS_PER_BLOCK];
  __shared__ Real z_shared[THREADS_PER_BLOCK];
  __shared__ Real fx_shared[THREADS_PER_BLOCK];
  __shared__ Real fy_shared[THREADS_PER_BLOCK];
  __shared__ Real fz_shared[THREADS_PER_BLOCK];

  #if !PRESCRIBED_CILIA

    __shared__ Real taux_shared[THREADS_PER_BLOCK];
    __shared__ Real tauy_shared[THREADS_PER_BLOCK];
    __shared__ Real tauz_shared[THREADS_PER_BLOCK];

  #endif

  Real v[6];
  Real xi, yi, zi;

  // Stay in the loop as long as any thread in the block still needs to compute velocities.
  for (int i = (start_seg + index); (i-threadIdx.x) < (start_seg + num_segs); i+=stride){

    if (i < (start_seg + num_segs)){

      v[0] = 0.0; v[1] = 0.0; v[2] = 0.0; v[3] = 0.0; v[4] = 0.0; v[5] = 0.0;
      xi = X[3*i]; yi = X[3*i + 1]; zi = X[3*i + 2];

    }

    for (int j_start = 0; j_start < NSWIM*NFIL*NSEG; j_start += THREADS_PER_BLOCK){

      const int j_to_read = j_start + threadIdx.x;

      if (j_to_read < NSWIM*NFIL*NSEG){

        x_shared[threadIdx.x] = X[3*j_to_read];
        y_shared[threadIdx.x] = X[3*j_to_read + 1];
        z_shared[threadIdx.x] = X[3*j_to_read + 2];
        fx_shared[threadIdx.x] = F[6*j_to_read];
        fy_shared[threadIdx.x] = F[6*j_to_read + 1];
        fz_shared[threadIdx.x] = F[6*j_to_read + 2];

        #if !PRESCRIBED_CILIA

          taux_shared[threadIdx.x] = F[6*j_to_read + 3];
          tauy_shared[threadIdx.x] = F[6*j_to_read + 4];
          tauz_shared[threadIdx.x] = F[6*j_to_read + 5];

        #endif

      }

      __syncthreads();

      if (i < (start_seg + num_segs)){

        for (int j=0; (j < THREADS_PER_BLOCK) && (j_start + j < NSWIM*NFIL*NSEG); j++){

          Real f[6];
          f[0] = fx_shared[j];
          f[1] = fy_shared[j];
          f[2] = fz_shared[j];

          #if !PRESCRIBED_CILIA

            f[3] = taux_shared[j];
            f[4] = tauy_shared[j];
            f[5] = tauz_shared[j];

          #endif

          #if PRESCRIBED_CILIA

            rpy_interaction<3,3>(v, f, i, xi, yi, zi, RSEG, j + j_start, x_shared[j], y_shared[j], z_shared[j], RSEG);

          #else

            rpy_interaction<6,6>(v, f, i, xi, yi, zi, RSEG, j + j_start, x_shared[j], y_shared[j], z_shared[j], RSEG);

          #endif

        }

      }

      __syncthreads();

    } // End of loop over filament segment forces and torques.

    if (i < (start_seg + num_segs)){

      const int p = 6*(i - start_seg);

      V[p] = v[0];
      V[p + 1] = v[1];
      V[p + 2] = v[2];
      V[p + 3] = v[3];
      V[p + 4] = v[4];
      V[p + 5] = v[5];

    }

  } // End of striding loop over filament segment velocities.

} // End of Mss_mult kernel.

__global__ void Mbb_mult(Real * __restrict__ V, const Real *const __restrict__ F, const Real *const __restrict__ X, const int start_blob, const int num_blobs){

  // Calculates the velocities of rigid-body blobs given the forces they experience.

  const int index = threadIdx.x + blockIdx.x*blockDim.x;
  const int stride = blockDim.x*gridDim.x;

  // Declare the shared memory for this thread block
  __shared__ Real x_shared[THREADS_PER_BLOCK];
  __shared__ Real y_shared[THREADS_PER_BLOCK];
  __shared__ Real z_shared[THREADS_PER_BLOCK];
  __shared__ Real fx_shared[THREADS_PER_BLOCK];
  __shared__ Real fy_shared[THREADS_PER_BLOCK];
  __shared__ Real fz_shared[THREADS_PER_BLOCK];

  Real v[3];
  Real xi, yi, zi;

  // Stay in the loop as long as any thread in the block still needs to compute velocities.
  for (int i = (start_blob + index); (i-threadIdx.x) < (start_blob + num_blobs); i+=stride){

    if (i < (start_blob + num_blobs)){

      v[0] = 0.0; v[1] = 0.0; v[2] = 0.0;
      xi = X[3*i]; yi = X[3*i + 1]; zi = X[3*i + 2];

    }

    for (int j_start = 0; j_start < NSWIM*NBLOB; j_start += THREADS_PER_BLOCK){

      const int j_to_read = j_start + threadIdx.x;

      if (j_to_read < NSWIM*NBLOB){

        x_shared[threadIdx.x] = X[3*j_to_read];
        y_shared[threadIdx.x] = X[3*j_to_read + 1];
        z_shared[threadIdx.x] = X[3*j_to_read + 2];
        fx_shared[threadIdx.x] = F[3*j_to_read];
        fy_shared[threadIdx.x] = F[3*j_to_read + 1];
        fz_shared[threadIdx.x] = F[3*j_to_read + 2];

      }

      __syncthreads();

      if (i < (start_blob + num_blobs)){

        for (int j=0; (j < THREADS_PER_BLOCK) && (j_start + j < NSWIM*NBLOB); j++){

          const Real f[3] = {fx_shared[j], fy_shared[j], fz_shared[j]};

          rpy_interaction<3,3>(v, f, i, xi, yi, zi, RBLOB, j_start + j, x_shared[j], y_shared[j], z_shared[j], RBLOB);

        }

      }

      __syncthreads();

    } // End of loop over blob forces.

    if (i < (start_blob + num_blobs)){

      const int p = 3*(i - start_blob);

      V[p] = v[0];
      V[p + 1] = v[1];
      V[p + 2] = v[2];

    }

  } // End of striding loop over blob velocities.

} // End of Mbb_mult kernel.

__global__ void Msb_mult(Real * __restrict__ V, const Real *const __restrict__ F, const Real *const __restrict__ Xs, const Real *const __restrict__ Xb, const int start_seg, const int num_segs){

  // Calculates the velocities of filament segments given the forces on
  // rigid-body blobs.

  const int index = threadIdx.x + blockIdx.x*blockDim.x;
  const int stride = blockDim.x*gridDim.x;

  // Declare the shared memory for this thread block
  __shared__ Real x_shared[THREADS_PER_BLOCK];
  __shared__ Real y_shared[THREADS_PER_BLOCK];
  __shared__ Real z_shared[THREADS_PER_BLOCK];
  __shared__ Real fx_shared[THREADS_PER_BLOCK];
  __shared__ Real fy_shared[THREADS_PER_BLOCK];
  __shared__ Real fz_shared[THREADS_PER_BLOCK];

  Real v[6];
  Real xi, yi, zi;

  // Stay in the loop as long as any thread in the block still needs to compute velocities.
  for (int i = (start_seg + index); (i-threadIdx.x) < (start_seg + num_segs); i+=stride){

    if (i < (start_seg + num_segs)){

      v[0] = 0.0; v[1] = 0.0; v[2] = 0.0; v[3] = 0.0; v[4] = 0.0; v[5] = 0.0;
      xi = Xs[3*i]; yi = Xs[3*i + 1]; zi = Xs[3*i + 2];

    }

    for (int j_start = 0; j_start < NSWIM*NBLOB; j_start += THREADS_PER_BLOCK){

      const int j_to_read = j_start + threadIdx.x;

      if (j_to_read < NSWIM*NBLOB){

        x_shared[threadIdx.x] = Xb[3*j_to_read];
        y_shared[threadIdx.x] = Xb[3*j_to_read + 1];
        z_shared[threadIdx.x] = Xb[3*j_to_read + 2];
        fx_shared[threadIdx.x] = F[3*j_to_read];
        fy_shared[threadIdx.x] = F[3*j_to_read + 1];
        fz_shared[threadIdx.x] = F[3*j_to_read + 2];

      }

      __syncthreads();

      if (i < (start_seg + num_segs)){

        for (int j=0; (j < THREADS_PER_BLOCK) && (j_start + j < NSWIM*NBLOB); j++){

          const Real f[3] = {fx_shared[j], fy_shared[j], fz_shared[j]};

          #if PRESCRIBED_CILIA

            rpy_interaction<3,3>(v, f, i, xi, yi, zi, RSEG, j_start + j, x_shared[j], y_shared[j], z_shared[j], RBLOB);

          #else

            rpy_interaction<6,3>(v, f, i, xi, yi, zi, RSEG, j_start + j, x_shared[j], y_shared[j], z_shared[j], RBLOB);

          #endif

        }

      }

      __syncthreads();

    } // End of loop over blob forces.

    if (i < (start_seg + num_segs)){

      const int p = 6*(i - start_seg);

      V[p] += v[0];
      V[p + 1] += v[1];
      V[p + 2] += v[2];
      V[p + 3] += v[3];
      V[p + 4] += v[4];
      V[p + 5] += v[5];

    }

  } // End of striding loop over filament segment velocities.

} // End of Msb_mult kernel.

__global__ void Mbs_mult(Real * __restrict__ V, const Real *const __restrict__ F, const Real *const __restrict__ Xb, const Real *const __restrict__ Xs, const int start_blob, const int num_blobs){

  // Calculates the velocities of blobs given the forces and torques on
  // filament segments.

  const int index = threadIdx.x + blockIdx.x*blockDim.x;
  const int stride = blockDim.x*gridDim.x;

  // Declare the shared memory for this thread block
  __shared__ Real x_shared[THREADS_PER_BLOCK];
  __shared__ Real y_shared[THREADS_PER_BLOCK];
  __shared__ Real z_shared[THREADS_PER_BLOCK];
  __shared__ Real fx_shared[THREADS_PER_BLOCK];
  __shared__ Real fy_shared[THREADS_PER_BLOCK];
  __shared__ Real fz_shared[THREADS_PER_BLOCK];

  #if !PRESCRIBED_CILIA

    __shared__ Real taux_shared[THREADS_PER_BLOCK];
    __shared__ Real tauy_shared[THREADS_PER_BLOCK];
    __shared__ Real tauz_shared[THREADS_PER_BLOCK];

  #endif

  Real v[3];
  Real xi, yi, zi;

  // Stay in the loop as long as any thread in the block still needs to compute velocities.
  for (int i = (start_blob + index); (i-threadIdx.x) < (start_blob + num_blobs); i+=stride){

    if (i < (start_blob + num_blobs)){

      v[0] = 0.0; v[1] = 0.0; v[2] = 0.0;
      xi = Xb[3*i]; yi = Xb[3*i + 1]; zi = Xb[3*i + 2];

    }

    for (int j_start = 0; j_start < NSWIM*NFIL*NSEG; j_start += THREADS_PER_BLOCK){

      const int j_to_read = j_start + threadIdx.x;

      if (j_to_read < NSWIM*NFIL*NSEG){

        x_shared[threadIdx.x] = Xs[3*j_to_read];
        y_shared[threadIdx.x] = Xs[3*j_to_read + 1];
        z_shared[threadIdx.x] = Xs[3*j_to_read + 2];
        fx_shared[threadIdx.x] = F[6*j_to_read];
        fy_shared[threadIdx.x] = F[6*j_to_read + 1];
        fz_shared[threadIdx.x] = F[6*j_to_read + 2];

        #if !PRESCRIBED_CILIA

          taux_shared[threadIdx.x] = F[6*j_to_read + 3];
          tauy_shared[threadIdx.x] = F[6*j_to_read + 4];
          tauz_shared[threadIdx.x] = F[6*j_to_read + 5];

        #endif

      }

      __syncthreads();

      if (i < (start_blob + num_blobs)){

        for (int j=0; (j < THREADS_PER_BLOCK) && (j_start + j < NSWIM*NFIL*NSEG); j++){

          Real f[6];
          f[0] = fx_shared[j];
          f[1] = fy_shared[j];
          f[2] = fz_shared[j];

          #if !PRESCRIBED_CILIA

            f[3] = taux_shared[j];
            f[4] = tauy_shared[j];
            f[5] = tauz_shared[j];

          #endif

          #if PRESCRIBED_CILIA

            rpy_interaction<3,3>(v, f, i, xi, yi, zi, RBLOB, j_start + j, x_shared[j], y_shared[j], z_shared[j], RSEG);

          #else

            rpy_interaction<3,6>(v, f, i, xi, yi, zi, RBLOB, j_start + j, x_shared[j], y_shared[j], z_shared[j], RSEG);

          #endif

        }

      }

      __syncthreads();

    } // End of loop over segment forces and torques.

    if (i < (start_blob + num_blobs)){

      const int p = 3*(i - start_blob);

      #if USE_BROYDEN_FOR_EVERYTHING || PRESCRIBED_CILIA

        V[p] += v[0];
        V[p + 1] += v[1];
        V[p + 2] += v[2];

      #else

        V[p] = v[0];
        V[p + 1] = v[1];
        V[p + 2] = v[2];

      #endif

    }

  } // End of striding loop over blob velocities.

} // End of Mbs_mult kernel.


__global__ void Mbs_mult_add(Real * __restrict__ V, const Real *const __restrict__ F, const Real *const __restrict__ Xb, const Real *const __restrict__ Xs, const int start_blob, const int num_blobs){

  // Calculates the velocities of blobs given the forces and torques on
  // filament segments.

  const int index = threadIdx.x + blockIdx.x*blockDim.x;
  const int stride = blockDim.x*gridDim.x;

  // Declare the shared memory for this thread block
  __shared__ Real x_shared[THREADS_PER_BLOCK];
  __shared__ Real y_shared[THREADS_PER_BLOCK];
  __shared__ Real z_shared[THREADS_PER_BLOCK];
  __shared__ Real fx_shared[THREADS_PER_BLOCK];
  __shared__ Real fy_shared[THREADS_PER_BLOCK];
  __shared__ Real fz_shared[THREADS_PER_BLOCK];

  #if !PRESCRIBED_CILIA

    __shared__ Real taux_shared[THREADS_PER_BLOCK];
    __shared__ Real tauy_shared[THREADS_PER_BLOCK];
    __shared__ Real tauz_shared[THREADS_PER_BLOCK];

  #endif

  Real v[3];
  Real xi, yi, zi;

  // Stay in the loop as long as any thread in the block still needs to compute velocities.
  for (int i = (start_blob + index); (i-threadIdx.x) < (start_blob + num_blobs); i+=stride){

    if (i < (start_blob + num_blobs)){

      v[0] = 0.0; v[1] = 0.0; v[2] = 0.0;
      xi = Xb[3*i]; yi = Xb[3*i + 1]; zi = Xb[3*i + 2];

    }

    for (int j_start = 0; j_start < NSWIM*NFIL*NSEG; j_start += THREADS_PER_BLOCK){

      const int j_to_read = j_start + threadIdx.x;

      if (j_to_read < NSWIM*NFIL*NSEG){

        x_shared[threadIdx.x] = Xs[3*j_to_read];
        y_shared[threadIdx.x] = Xs[3*j_to_read + 1];
        z_shared[threadIdx.x] = Xs[3*j_to_read + 2];
        fx_shared[threadIdx.x] = F[6*j_to_read];
        fy_shared[threadIdx.x] = F[6*j_to_read + 1];
        fz_shared[threadIdx.x] = F[6*j_to_read + 2];

        #if !PRESCRIBED_CILIA

          taux_shared[threadIdx.x] = F[6*j_to_read + 3];
          tauy_shared[threadIdx.x] = F[6*j_to_read + 4];
          tauz_shared[threadIdx.x] = F[6*j_to_read + 5];

        #endif

      }

      __syncthreads();

      if (i < (start_blob + num_blobs)){

        for (int j=0; (j < THREADS_PER_BLOCK) && (j_start + j < NSWIM*NFIL*NSEG); j++){

          Real f[6];
          f[0] = fx_shared[j];
          f[1] = fy_shared[j];
          f[2] = fz_shared[j];

          #if !PRESCRIBED_CILIA

            f[3] = taux_shared[j];
            f[4] = tauy_shared[j];
            f[5] = tauz_shared[j];

          #endif

          #if PRESCRIBED_CILIA

            rpy_interaction<3,3>(v, f, i, xi, yi, zi, RBLOB, j_start + j, x_shared[j], y_shared[j], z_shared[j], RSEG);

          #else

            rpy_interaction<3,6>(v, f, i, xi, yi, zi, RBLOB, j_start + j, x_shared[j], y_shared[j], z_shared[j], RSEG);

          #endif

        }

      }

      __syncthreads();

    } // End of loop over segment forces and torques.

    if (i < (start_blob + num_blobs)){

      const int p = 3*(i - start_blob);

      #if USE_BROYDEN_FOR_EVERYTHING || PRESCRIBED_CILIA

        V[p] += v[0];
        V[p + 1] += v[1];
        V[p + 2] += v[2];

      #else

        V[p] = v[0];
        V[p + 1] = v[1];
        V[p + 2] = v[2];

      #endif

    }

  } // End of striding loop over blob velocities.

} // End of Mbs_mult_add kernel.



// Stokes drag kernels

__global__ void Ms_mult(Real * __restrict__ V, const Real *const __restrict__ F, const int start_seg, const int num_segs){

    // Calculates the velocities of filament segments given the forces and torques
    // on the segments.

    const int index = threadIdx.x + blockIdx.x*blockDim.x;
    const int stride = blockDim.x*gridDim.x;

    const Real temp = 1.0/(6.0*PI*MU*RSEG);
    const Real temp2 = 1.0/(8.0*PI*MU*RSEG*RSEG*RSEG);

    // Stay in the loop as long as any thread in the block still needs to compute velocities.
    for (int i = (start_seg + index); i < (start_seg + num_segs); i+=stride){

      const int p = 6*(i - start_seg);

      V[p] = temp*F[6*i];
      V[p + 1] = temp*F[6*i + 1];
      V[p + 2] = temp*F[6*i + 2];
      V[p + 3] = temp2*F[6*i + 3];
      V[p + 4] = temp2*F[6*i + 4];
      V[p + 5] = temp2*F[6*i + 5];

    } // End of striding loop over filament segment velocities.

  } // End of Ms_mult kernel.

__global__ void Mb_mult(Real * __restrict__ V, const Real *const __restrict__ F, const int start_blob, const int num_blobs){

  // Calculates the velocities of rigid-body blobs given the forces they experience.

  const int index = threadIdx.x + blockIdx.x*blockDim.x;
  const int stride = blockDim.x*gridDim.x;

  const Real temp = 1.0/(6.0*PI*MU*RBLOB);

  // Stay in the loop as long as any thread in the block still needs to compute velocities.
  for (int i = (start_blob + index); i < (start_blob + num_blobs); i+=stride){

    const int p = 3*(i - start_blob);

    V[p] = temp*F[3*i];
    V[p + 1] = temp*F[3*i + 1];
    V[p + 2] = temp*F[3*i + 2];

  } // End of striding loop over blob velocities.

} // End of Mb_mult kernel.

__global__ void Mb_fill_zero(Real * __restrict__ V, const int start_blob, const int num_blobs){

  // Fill zero velocity arrays

  const int index = threadIdx.x + blockIdx.x*blockDim.x;
  const int stride = blockDim.x*gridDim.x;

  // Stay in the loop as long as any thread in the block still needs to fill zeros.
  for (int i = (start_blob + index); i < (start_blob + num_blobs); i+=stride){

    const int p = 3*(i - start_blob);

    #if USE_BROYDEN_FOR_EVERYTHING

      V[p] += 0.0;
      V[p + 1] += 0.0;
      V[p + 2] += 0.0;

    #else

      V[p] = 0.0;
      V[p + 1] = 0.0;
      V[p + 2] = 0.0;

    #endif

  } // End of striding loop over blob velocities.

} // End of Mb_fill_zero kernel.








// Weakly-coupled-filaments kernels

// See https://docs.nvidia.com/cuda/cuda-c-programming-guide/index.html#atomic-functions
__device__ Real atomicAdd_arch_indep(Real* address, Real val){
    unsigned long long int* address_as_ull =
                              (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;

    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,
                        __double_as_longlong(val +
                               __longlong_as_double(assumed)));

    // Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
    } while (assumed != old);

    return __longlong_as_double(old);

}

__global__ void reset_mean_fil_posns(Real *Z){

  // Sets the mean position of each filament to zero, allowing us to assign-via-addition in write_mean_fil_posns(...).

  const int index = threadIdx.x + blockIdx.x*blockDim.x;
  const int stride = blockDim.x*gridDim.x;

  for (int i = index; i < 3*NSWIM*NFIL; i += stride){

      Z[i] = 0.0;

  }

}

__global__ void write_mean_fil_posns(Real * __restrict__ Z, const Real *const __restrict__ X){

  // Calculates the mean position of each filament.

  const int index = threadIdx.x + blockIdx.x*blockDim.x;
  const int stride = blockDim.x*gridDim.x;

  for (int i = index; i < NSWIM*NFIL*NSEG; i += stride){

      const int fil_id = i/NSEG;

      atomicAdd_arch_indep(&Z[3*fil_id], X[3*i]/Real(NSEG));
      atomicAdd_arch_indep(&Z[3*fil_id + 1], X[3*i + 1]/Real(NSEG));
      atomicAdd_arch_indep(&Z[3*fil_id + 2], X[3*i + 2]/Real(NSEG));

  }

}

__global__ void reset_mean_fil_force_and_moment(Real * __restrict__ F, Real * __restrict__ S){

  // Sets the total force and moment on each filament to zero, allowing us to assign-via-addition in find_total_fil_force_and_first_moment(...).

  const int index = threadIdx.x + blockIdx.x*blockDim.x;
  const int stride = blockDim.x*gridDim.x;

  for (int i = index; i < 3*NSWIM*NFIL; i += stride){

      F[i] = 0.0;

      /*S[3*i] = 0.0;
      S[3*i + 1] = 0.0;
      S[3*i + 2] = 0.0;*/

  }

}

__global__ void find_total_fil_force_and_first_moment(Real * __restrict__ F, Real * __restrict__ S, const Real *const __restrict__ Fs, const Real *const __restrict__ Xs, const Real *const __restrict__ Z){


  // Calculates the total force and total first moment of force on each filament.

  const int index = threadIdx.x + blockIdx.x*blockDim.x;
  const int stride = blockDim.x*gridDim.x;

  for (int i = index; i < NSWIM*NFIL*NSEG; i += stride){

    const int fil_id = i/NSEG;

    Real f[3], r[3];

    f[0] = Fs[6*i];
    f[1] = Fs[6*i + 1];
    f[2] = Fs[6*i + 2];

    /*r[0] = Xs[3*i] - Z[3*fil_id];
    r[0] = Xs[3*i + 1] - Z[3*fil_id + 1];
    r[0] = Xs[3*i + 2] - Z[3*fil_id + 2];*/

    for (int n = 0; n < 3; n++){

      atomicAdd_arch_indep(&F[3*fil_id + n], f[n]);

      /*for (int m = 0; m < 3; m++){

        atomicAdd_arch_indep(&S[9*fil_id + 3*n + m], f[m]*r[n]);

      }*/

    }

  }

}

__device__ void grad_oseen_tensor_times_force_moment(Real& vx, Real& vy, Real& vz, const Real ForceFactor, const Real rx, const Real ry, const Real rz, const Real rm3, const Real rm5, const Real *const S){

  // x-derivative
  Real Hxx = ForceFactor*(rm3*rx - 3.0*rm5*rx*rx*rx);
  Real Hxy = ForceFactor*(rm3*ry - 3.0*rm5*rx*rx*ry);
  Real Hxz = ForceFactor*(rm3*rz - 3.0*rm5*rx*rx*rz);

  Real Hyy = ForceFactor*(-rm3*rx - 3.0*rm5*rx*ry*ry);
  Real Hyz = ForceFactor*(-3.0*rm5*rx*ry*rz);

  Real Hzz = ForceFactor*(-rm3*rx - 3.0*rm5*rx*rz*rz);

  vx += Hxx*S[0] + Hxy*S[1] + Hxz*S[2];
  vy += Hxy*S[0] + Hyy*S[1] + Hyz*S[2];
  vz += Hxz*S[0] + Hyz*S[1] + Hzz*S[2];

  // y-derivative
  Hxx = ForceFactor*(-rm3*ry - 3.0*rm5*ry*rx*rx);
  Hxy = ForceFactor*(rm3*rx - 3.0*rm5*ry*rx*ry);
  Hxz = Hyz; // We haven't re-defined Hyz yet

  Hyy = ForceFactor*(rm3*ry - 3.0*rm5*ry*ry*ry);
  Hyz = ForceFactor*(rm3*rz - 3.0*rm5*ry*ry*rz);

  Hzz = ForceFactor*(-rm3*ry - 3.0*rm5*ry*rz*rz);

  vx += Hxx*S[3] + Hxy*S[4] + Hxz*S[5];
  vy += Hxy*S[3] + Hyy*S[4] + Hyz*S[5];
  vz += Hxz*S[3] + Hyz*S[4] + Hzz*S[5];

  // z-derivative
  Hxx = ForceFactor*(-rm3*rz - 3.0*rm5*rz*rx*rx);
  Hxy = Hxz;
  Hxz = ForceFactor*(rm3*rx - 3.0*rm5*rz*rx*rz);

  Hyy = ForceFactor*(-rm3*rz - 3.0*rm5*rz*ry*ry);
  Hyz = ForceFactor*(rm3*ry - 3.0*rm5*rz*ry*rz);

  Hzz = ForceFactor*(rm3*rz - 3.0*rm5*rz*rz*rz);

  vx += Hxx*S[6] + Hxy*S[7] + Hxz*S[8];
  vy += Hxy*S[6] + Hyy*S[7] + Hyz*S[8];
  vz += Hxz*S[6] + Hyz*S[7] + Hzz*S[8];

}

__global__ void Msf_mult(Real *V, const Real *const F, const Real *const Ftot, const Real *const Stot, const Real *const X, const Real *const Z, const int start_seg, const int num_segs){

  // Calculates the segment velocities.
  // This kernel ASSUMES that we're using prescribed-motion filaments, rather than checking for the macro etc.

  const int index = threadIdx.x + blockIdx.x*blockDim.x;
  const int stride = blockDim.x*gridDim.x;

  // This kernel won't make use of shared memory etc. to start with and I'll optimise it later.
  for (int i = (start_seg + index); i < (start_seg + num_segs); i += stride){

    const int fil_id = i/NSEG;

    Real v[3] = {0.0, 0.0, 0.0};
    Real xi = X[3*i], yi = X[3*i + 1], zi = X[3*i + 2];

    for (int j = 0; j < NSWIM*NFIL; j++){

      if (fil_id == j){

        // In the same filament, we do normal RPY interactions.
        for (int k = 0; k < NSEG; k++){

          const int jj = k + j*NSEG; // Global ID of other segment.

          const Real f[3] = {F[6*jj], F[6*jj + 1], F[6*jj + 2]};

          rpy_interaction<3,3>(v, f, i, xi, yi, zi, RSEG, jj, X[3*jj], X[3*jj + 1], X[3*jj + 2], RSEG);

        }

      } else {

        // Between different filaments, we have weak coupling.
        const Real rx = X[3*i] - Z[3*j];
        const Real ry = X[3*i + 1] - Z[3*j + 1];
        const Real rz = X[3*i + 2] - Z[3*j + 2];

        const Real rm1 = rsqrt(rx*rx + ry*ry + rz*rz); // norm(r) cannot possibly be 0 so we don't need to add a small number to the denominator.
        const Real rm3 = rm1*rm1*rm1;

        const Real ForceFactor = 0.03978873577297383394/MU; // = 1/(8*PI*MU)

        // Stokeslet term
        Real r_dot_Ftot = rx*Ftot[3*j] + ry*Ftot[3*j + 1] + rz*Ftot[3*j + 2];
        v[0] += ForceFactor*(rm1*Ftot[3*j] + rm3*rx*r_dot_Ftot);
        v[1] += ForceFactor*(rm1*Ftot[3*j + 1] + rm3*ry*r_dot_Ftot);
        v[2] += ForceFactor*(rm1*Ftot[3*j + 2] + rm3*rz*r_dot_Ftot);

        #if INFINITE_PLANE_WALL

          const Real Rz = X[3*i + 2] + Z[3*j + 2];
          const Real Rm1 = rsqrt(rx*rx + ry*ry + Rz*Rz);
          const Real Rm3 = Rm1*Rm1*Rm1;
          const Real Rm5 = Rm3*Rm1*Rm1;

          // Blakelet term (this just involves additional terms to the Stokeslet above)
          r_dot_Ftot = rx*Ftot[3*j] + ry*Ftot[3*j + 1] + Rz*Ftot[3*j + 2];
          v[0] -= ForceFactor*(Rm1*Ftot[3*j] + Rm3*rx*r_dot_Ftot);
          v[1] -= ForceFactor*(Rm1*Ftot[3*j + 1] + Rm3*ry*r_dot_Ftot);
          v[2] -= ForceFactor*(Rm1*Ftot[3*j + 2] + Rm3*Rz*r_dot_Ftot);

          const Real Dxx = Rm5*3.0*rx*rx*X[3*i + 2] + Rm3*X[3*i + 2];
          const Real Dxy = Rm5*3.0*rx*ry*X[3*i + 2];
          const Real Dxz = -(Rm5*3.0*rx*Rz*X[3*i + 2] - Rm3*rx);
          const Real Dyy = Rm5*3.0*ry*ry*X[3*i + 2] + Rm3*X[3*i + 2];
          const Real Dyz = -(Rm5*3.0*ry*Rz*X[3*i + 2] - Rm3*ry);
          const Real Dzx = Rm5*3.0*Rz*rx*X[3*i + 2] + Rm3*rx;
          const Real Dzy = Rm5*3.0*Rz*ry*X[3*i + 2] + Rm3*ry;
          const Real Dzz = -(Rm5*3.0*Rz*Rz*X[3*i + 2] + Rm3*X[3*i + 2]);

          v[0] += 2.0*ForceFactor*Z[3*j + 2]*(Dxx*Ftot[3*j] + Dxy*Ftot[3*j + 1] + Dxz*Ftot[3*j + 2]);
          v[1] += 2.0*ForceFactor*Z[3*j + 2]*(Dxy*Ftot[3*j] + Dyy*Ftot[3*j + 1] + Dyz*Ftot[3*j + 2]);
          v[2] += 2.0*ForceFactor*Z[3*j + 2]*(Dzx*Ftot[3*j] + Dzy*Ftot[3*j + 1] + Dzz*Ftot[3*j + 2]);

          // Gradient term

        #else

         /* // Gradient term
          const Real rm5 = rm3*rm1*rm1;
          grad_oseen_tensor_times_force_moment(v[0], v[1], v[2], ForceFactor, rx, ry, rz, rm3, rm5, &Stot[9*j]); */

        #endif

      }

    }

    const int p = 6*(i - start_seg);

    V[p] = v[0];
    V[p + 1] = v[1];
    V[p + 2] = v[2];
    V[p + 3] = 0.0;
    V[p + 4] = 0.0;
    V[p + 5] = 0.0;

  }

}








// Generic interaction kernels

__global__ void barrier_forces(Real * __restrict__ f_segs, Real * __restrict__ f_blobs_repulsion, const Real *const __restrict__ x_segs, const Real *const __restrict__ x_blobs, const int start_seg, const int num_segs, const int start_blob, const int num_blobs){

  #if !(PRESCRIBED_CILIA || NO_CILIA_SQUIRMER)

    // Work out which particle(s) this thread will compute the force for
    const int index = threadIdx.x + blockIdx.x*blockDim.x;
    const int stride = blockDim.x*gridDim.x;

    // Declare the shared memory for this thread block
    __shared__ Real x_shared[THREADS_PER_BLOCK];
    __shared__ Real y_shared[THREADS_PER_BLOCK];
    __shared__ Real z_shared[THREADS_PER_BLOCK];

    Real fx, fy, fz;
    Real xi, yi, zi;
    int fili;

    for (int i = (start_seg + index); (i-threadIdx.x) < (start_seg + num_segs); i+=stride){

      if (i < (start_seg + num_segs)){

        fx = 0.0;
        fy = 0.0;
        fz = 0.0;

        xi = x_segs[3*i];
        yi = x_segs[3*i + 1];
        zi = x_segs[3*i + 2];

        fili = i/NSEG;

        #if INFINITE_PLANE_WALL

          if (zi < BASE_HEIGHT_ABOVE_SURFACE){

            fz = fmin(1.0, 1.0 - (zi - RSEG)/(0.5*DL - RSEG)); // max magnitude one radius from wall
            fz *= REPULSIVE_FORCE_FACTOR*END_FORCE_MAGNITUDE*fz*fz*fz; // 4th power

          }

        #endif

      }

      for (int j_start = 0; j_start < NTOTAL; j_start += THREADS_PER_BLOCK){

        const int j_to_read = j_start + threadIdx.x;

        if (j_to_read < NSWIM*NFIL*NSEG){

          x_shared[threadIdx.x] = x_segs[3*j_to_read];
          y_shared[threadIdx.x] = x_segs[3*j_to_read + 1];
          z_shared[threadIdx.x] = x_segs[3*j_to_read + 2];

        } else if (j_to_read < NTOTAL){

          x_shared[threadIdx.x] = x_blobs[3*(j_to_read - NSWIM*NFIL*NSEG)];
          y_shared[threadIdx.x] = x_blobs[3*(j_to_read - NSWIM*NFIL*NSEG) + 1];
          z_shared[threadIdx.x] = x_blobs[3*(j_to_read - NSWIM*NFIL*NSEG) + 2];

        }

        __syncthreads();

        if (i < (start_seg + num_segs)){

          for (int j=0; (j < THREADS_PER_BLOCK) && (j_start + j < NTOTAL); j++){

            const Real a_sum = (j_start + j < NSWIM*NFIL*NSEG) ? 2.0*RSEG : RSEG + RBLOB;
            const Real chi_fac = 10.0/a_sum; // 1.0/(1.1*a_sum - a_sum)

            const Real dx = xi - x_shared[j];
            const Real dy = yi - y_shared[j];
            const Real dz = zi - z_shared[j];

            const Real dist = sqrt(dx*dx + dy*dy + dz*dz);

            int filj = (j_start + j)/NSEG;

            if (!(fili==filj && abs(i -(j_start + j))<=1) && (dist < 1.1*a_sum)){

              Real fac = fmin(1.0, 1.0 - chi_fac*(dist - a_sum));
              fac *= REPULSIVE_FORCE_FACTOR*END_FORCE_MAGNITUDE*fac*fac*fac;

              const Real dm1 = 1.0/dist;

              fx += fac*dx*dm1;
              fy += fac*dy*dm1;
              fz += fac*dz*dm1;

            }

          }

        }

        __syncthreads();

      }

      if (i < (start_seg + num_segs)){

        f_segs[6*i] += fx; // Don't shift index as f_segs has global size.
        f_segs[6*i + 1] += fy;
        f_segs[6*i + 2] += fz;

      }

    }

    #if !PRESCRIBED_BODY_VELOCITIES

      for (int i = (start_blob + index); (i-threadIdx.x) < (start_blob + num_blobs); i+=stride){

        if (i < (start_blob + num_blobs)){

          fx = 0.0;
          fy = 0.0;
          fz = 0.0;

          xi = x_blobs[3*i];
          yi = x_blobs[3*i + 1];
          zi = x_blobs[3*i + 2];

        }

        for (int j_start = 0; j_start < NTOTAL; j_start += THREADS_PER_BLOCK){

          const int j_to_read = j_start + threadIdx.x;

          if (j_to_read < NSWIM*NFIL*NSEG){

            x_shared[threadIdx.x] = x_segs[3*j_to_read];
            y_shared[threadIdx.x] = x_segs[3*j_to_read + 1];
            z_shared[threadIdx.x] = x_segs[3*j_to_read + 2];

          } else if (j_to_read < NTOTAL){

            x_shared[threadIdx.x] = x_blobs[3*(j_to_read - NSWIM*NFIL*NSEG)];
            y_shared[threadIdx.x] = x_blobs[3*(j_to_read - NSWIM*NFIL*NSEG) + 1];
            z_shared[threadIdx.x] = x_blobs[3*(j_to_read - NSWIM*NFIL*NSEG) + 2];

          }

          __syncthreads();

          if (i < (start_blob + num_blobs)){

            for (int j=0; (j < THREADS_PER_BLOCK) && (j_start + j < NTOTAL); j++){

              const Real a_sum = (j_start + j < NSWIM*NFIL*NSEG) ? RBLOB + RSEG : 2.0*RBLOB;
              const Real chi_fac = 10.0/a_sum; // 1.0/(1.1*a_sum - a_sum)

              const Real dx = xi - x_shared[j];
              const Real dy = yi - y_shared[j];
              const Real dz = zi - z_shared[j];

              const Real dist = sqrt(dx*dx + dy*dy + dz*dz);

              if (((i + NSWIM*NFIL*NSEG) != (j_start + j)) && (dist < 1.1*a_sum)){

                Real fac = fmin(1.0, 1.0 - chi_fac*(dist - a_sum));
                fac *= REPULSIVE_FORCE_FACTOR*END_FORCE_MAGNITUDE*fac*fac*fac;

                const Real dm1 = 1.0/dist;

                fx += fac*dx*dm1;
                fy += fac*dy*dm1;
                fz += fac*dz*dm1;

              }

            }

          }

          __syncthreads();

        }

        if (i < (start_blob + num_blobs)){

          const int p = 3*(i - start_blob);

          f_blobs_repulsion[p] = fx;
          f_blobs_repulsion[p + 1] = fy;
          f_blobs_repulsion[p + 2] = fz;

        }

      }

    #endif

  #endif

} // End of barrier forces kernel.


// Periodic interaction kernels

__global__ void periodic_barrier_forces(Real * __restrict__ f_segs,
                                        Real * __restrict__ f_blobs_repulsion,
                                        const Real *const __restrict__ x_segs,
                                        const Real *const __restrict__ x_blobs,
                                        const int num_segs,
                                        const int num_blobs,
                                        const Real boxsize){

  #if !(PRESCRIBED_CILIA || NO_CILIA_SQUIRMER)

    // Work out which particle(s) this thread will compute the force for
    const int index = threadIdx.x + blockIdx.x*blockDim.x;
    const int stride = blockDim.x*gridDim.x;

    Real fx, fy, fz;
    Real xi, yi, zi;
    int fili;

    #if !PRESCRIBED_BODY_VELOCITIES

      for (int i = index; i < num_blobs; i+=stride){

        fx = 0.0;
        fy = 0.0;
        fz = 0.0;

        xi = x_blobs[3*i];
        yi = x_blobs[3*i + 1];
        zi = x_blobs[3*i + 2];

        box_images(xi, boxsize);
        box_images(yi, boxsize);
        box_images(zi, boxsize);

        for (int j=0; j < NTOTAL; j++){

          Real xj = x_blobs[3*j];
          Real yj = x_blobs[3*j+1];
          Real zj = x_blobs[3*j+2];

          box_images(xj, boxsize);
          box_images(yj, boxsize);
          box_images(zj, boxsize);

          const Real a_sum = (j < NSWIM*NFIL*NSEG) ? RBLOB + RSEG : 2.0*RBLOB;
          const Real chi_fac = 10.0/a_sum; // 1.0/(1.1*a_sum - a_sum)

          Real dx = xi - xj;
          Real dy = yi - yj;
          Real dz = zi - zj;

          dx -= boxsize * Real(int(dx/(0.5*boxsize)));
          dy -= boxsize * Real(int(dy/(0.5*boxsize)));
          dz -= boxsize * Real(int(dz/(0.5*boxsize)));

          // printf("**********xi=%.4f xj=%.4f dx=%.4f ratio=%.4f reset = %.4f\n", 
          //       xi, xj, dx, dx/(0.5*boxsize), Real(int(dx/(0.5*boxsize))));
          
          // if(Real(int(dx/(0.5*boxsize))) != 0){
          //   printf("**********xi=%.4f xj=%.4f dx=%.4f\n", xi, xj, dx);
          // }

          const Real dist = sqrt(dx*dx + dy*dy + dz*dz);

          if (((i + NSWIM*NFIL*NSEG) != j) && (dist < 1.1*a_sum)){

            Real fac = fmin(1.0, 1.0 - chi_fac*(dist - a_sum));
            fac *= REPULSIVE_FORCE_FACTOR*END_FORCE_MAGNITUDE*fac*fac*fac;

            const Real dm1 = 1.0/dist;

            fx += fac*dx*dm1;
            fy += fac*dy*dm1;
            fz += fac*dz*dm1;

          }

        }

        __syncthreads();

        f_blobs_repulsion[3*i] = fx;
        f_blobs_repulsion[3*i + 1] = fy;
        f_blobs_repulsion[3*i + 2] = fz;

      }

    #endif
    
  #endif

} // End of barrier forces kernel.


__host__ __device__
void box_images(Real &x, Real box_size){
    x -= floor(x/box_size)*box_size;
}
