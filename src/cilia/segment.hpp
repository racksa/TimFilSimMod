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
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>


class segment{

public:

  // current positions are used in the mobility solve so we store them in global arrays to avoid copying and needlessly storing duplicates
  Real *x;

  // orientations and old positions don't appear in the mobility solve and can be stored locally
  Real xm1[3];
  Real xm2[3];
  quaternion q;
  quaternion qm1;
  Real u[3];
  Real um1[3];

  ~segment();
  segment();

  void initial_setup(const Real *const x_in, const quaternion& q_in, const Real *const u_in);
  void initial_guess(const int nt);
  void update(const Real *const u_update);
  void tangent(Real *const t) const;
  void normal(Real *const n) const;
  void binormal(Real *const b) const;
  void tangent(matrix& t) const;
  void normal(matrix& n) const;
  void binormal(matrix& b) const;

};

#endif // MY_SEGMENT_HEADER_INCLUDED
