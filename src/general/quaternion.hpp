// quaternion.hpp

// =============================================================================
// Include guard
#ifndef MY_QUATERNION_HEADER_INCLUDED
#define MY_QUATERNION_HEADER_INCLUDED

// =============================================================================
// Forward declared dependencies
class matrix;

// =============================================================================
// Included dependencies
#include <iostream>
#include <fstream>
#include "../../config.hpp"

class quaternion {

public:

  Real scalar_part;
  Real vector_part[3];

  ~quaternion();
  quaternion();
  quaternion(const Real q0, const Real q1, const Real q2, const Real q3);
  quaternion(const Real q0, const Real *const q);
  quaternion(const quaternion& q);

  //
  // OPERATOR OVERLOADING
  //

  // Element access
  Real& operator()(const int elem);
  Real operator()(const int elem) const;

  // Assignment
  quaternion& operator =(const quaternion& q); // q1 = q

  // In-place scalar multiplication
  quaternion& operator *=(const Real s); // q *= s

  // In-place scalar division
  quaternion& operator /=(const Real s); // q /= s

  // Sign swap
  quaternion operator -() const; // For calls of the form q1 = -q2;

  // In-place right quaternion multiplication
  quaternion& operator *=(const quaternion& q); // q1 *= q2, resulting in q1 = q1*q2

  // In-place addition
  quaternion& operator +=(const quaternion& q); // q1 += q2

  // In-place subtraction
  quaternion& operator -=(const quaternion& q); // q1 -= q2

  //
  // OTHER METHODS
  //

  void randomise();
  Real norm() const;
  void normalise_in_place();
  void conj_in_place();
  void sqrt_in_place();
  void sqrt_in_place_ortho(Real *polar_dir);
  void write_data(std::ofstream& data_file) const;

  // Array versions for if we just want to see the values and maybe do basic addition etc.
  void tangent(Real *const t) const;
  void normal(Real *const n) const;
  void binormal(Real *const b) const;

  // Matrix versions for if we want to rotate them etc.
  void tangent(matrix& t) const;
  void normal(matrix& n) const;
  void binormal(matrix& b) const;

  // No array versions for these; I don't know why we would form these matrices if we
  // weren't planning on multiplying something by them. We can have in-place and
  // out-of-place versions though, so we can choose based on what's already been allocated
  // in a given function call.
  void rot_mat(matrix& R) const;
  matrix rot_mat() const;
  void psi_mat(matrix& Psi) const;
  matrix psi_mat() const;
  void left_mult_mat(matrix& qdot) const;
  matrix left_mult_mat() const;
  void right_mult_mat(matrix& dotq) const;
  matrix right_mult_mat() const;

}; // End of quaternion class definition

//
// BINARY OPERATOR OVERLOADING
//

quaternion operator /(quaternion q, const Real s);
quaternion operator *(quaternion q, const Real s);
quaternion operator *(const Real s, quaternion q);
quaternion operator *(quaternion p, const quaternion& q);
quaternion operator +(quaternion p, const quaternion& q);
quaternion operator -(quaternion p, const quaternion& q);
std::ostream& operator <<(std::ostream& stream, const quaternion& q);

//
// OTHER FUNCTIONS ASSOCIATED WITH QUATERNIONS
//

void lie_exp(quaternion& q, const Real *const u);
quaternion lie_exp(const Real *const u);
void midpoint_quaternion(quaternion& qmid, const quaternion& q1, const quaternion& q2);
quaternion midpoint_quaternion(const quaternion& q1, const quaternion& q2);
void dexp(Real *const out, const Real *const u, const Real *const v);
void dexpinv(Real *const out, const Real *const u, const Real *const v);
void dexpinv_transpose(Real *const out, const Real *const u, const Real *const v);
void bch(Real *const out, const Real *const u, const Real *const v);

#endif // MY_QUATERNION_HEADER_INCLUDED
