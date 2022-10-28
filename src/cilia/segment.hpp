// segment.hpp

// =============================================================================
// Include guard
#ifndef MY_SEGMENT_HEADER_INCLUDED
#define MY_SEGMENT_HEADER_INCLUDED

// =============================================================================
// Forward declared dependencies
class matrix;

// =============================================================================
// Included dependencies
#include "quaternion.hpp"

class segment{

public:

  // current positions are used in the mobility solve so we store them in global arrays to avoid copying and needlessly storing duplicates
  double *x;

  // orientations and old positions don't appear in the mobility solve and can be stored locally
  double xm1[3];
  double xm2[3];
  quaternion q;
  quaternion qm1;
  double u[3];
  double um1[3];

  ~segment();
  segment();

  void initial_setup(const double *const x_in, const quaternion& q_in, const double *const u_in);
  void initial_guess(const int nt);
  void update(const double *const u_update);
  void tangent(double *const t) const;
  void normal(double *const n) const;
  void binormal(double *const b) const;
  void tangent(matrix& t) const;
  void normal(matrix& n) const;
  void binormal(matrix& b) const;

};

#endif // MY_SEGMENT_HEADER_INCLUDED
