// broyden_solver.cpp

#include <algorithm>
#include "broyden_solver.hpp"
#include "swimmer.hpp"

broyden_solver::~broyden_solver(){}

broyden_solver::broyden_solver(){

  C = matrix(NBROY, MAX_BROYDEN_ITER);
  D = matrix(NBROY, MAX_BROYDEN_ITER);
  error = matrix(NBROY, 1);
  new_error = matrix(NBROY, 1);
  update = matrix(NBROY, 1);

  iter = 0;
  max_iter = 0;
  total_iter = 0;
  avg_iter = 0.0;

}

#if !(PRESCRIBED_CILIA || NO_CILIA_SQUIRMER)

  void broyden_solver::find_update(const std::vector<swimmer>& swimmers, const int nt){

    // Note all of the minus signs because we want to end up with update = -Jinv * error;

    #if INFINITE_PLANE_WALL

      const int per_body = 6*NFIL*NSEG;

    #elif USE_BROYDEN_FOR_EVERYTHING

      const int per_body = 6*NFIL*NSEG + 3*NBLOB + 6;

    #else

      const int per_body = 6*NFIL*NSEG + 6;

    #endif

    // Start-of-step part
    for (int n = 0; n < NSWIM; n++){

      update.set_block(n*per_body, per_body, swimmers[n].jacobian_inv_mult(error.get_block(n*per_body, per_body), nt));

    }

    update *= -JACOBIAN_CONFIDENCE_FACTOR;

    // Update part
    for (int i = 0; i < iter; i++){

      double error_dot_D_col_i = 0.0;

      for (int n = 0; n < NBROY; n++){

        error_dot_D_col_i += error(n)*D(n,i);

      }

      update -= error_dot_D_col_i*C.get_col(i);

    }

  }

  void broyden_solver::end_of_iter(const std::vector<swimmer>& swimmers, const int nt, const int nt_start, const bool error_is_too_large){

    // Note all of the minus signs because we want to end up with C.col(iter) = -Jinv * new_error;

    if (error_is_too_large && (iter < MAX_BROYDEN_ITER-1)){

      #if INFINITE_PLANE_WALL

        const int per_body = 6*NFIL*NSEG;

      #elif USE_BROYDEN_FOR_EVERYTHING

        const int per_body = 6*NFIL*NSEG + 3*NBLOB + 6;

      #else

        const int per_body = 6*NFIL*NSEG + 6;

      #endif

      // Start-of-step part
      for (int n = 0; n < NSWIM; n++){

        update.set_block(n*per_body, per_body, swimmers[n].jacobian_inv_mult(new_error.get_block(n*per_body, per_body), nt));

      }

      C.set_col(iter, -JACOBIAN_CONFIDENCE_FACTOR*update);

      // Update part
      for (int i = 0; i < iter; i++){

        double error_dot_D_col_i = 0.0;

        for (int n = 0; n < NBROY; n++){

          error_dot_D_col_i += new_error(n)*D(n,i);

        }

        C.subtract_from_col(iter, error_dot_D_col_i*C.get_col(i));

      }

      D.set_col(iter, new_error - error);

      const double delta_error = norm(D.get_col(iter));

      C.divide_col(iter, delta_error);
      D.divide_col(iter, delta_error);

      error.swap(new_error);

    }

    iter++;
    total_iter++;
    max_iter = std::max<int>(max_iter, iter);
    avg_iter = double(total_iter)/(nt + 1.0 - nt_start);

  }

#endif
