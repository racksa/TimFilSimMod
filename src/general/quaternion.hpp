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

class quaternion {

public:

  double scalar_part;
  double vector_part[3];

  ~quaternion();
  quaternion();
  quaternion(const double q0, const double q1, const double q2, const double q3);
  quaternion(const double q0, const double *const q);
  quaternion(const quaternion& q);

  //
  // OPERATOR OVERLOADING
  //

  // Element access
  double& operator()(const int elem);
  double operator()(const int elem) const;

  // Assignment
  quaternion& operator =(const quaternion& q); // q1 = q

  // In-place scalar multiplication
  quaternion& operator *=(const double s); // q *= s

  // In-place scalar division
  quaternion& operator /=(const double s); // q /= s

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

  double norm() const;
  void normalise_in_place();
  void conj_in_place();
  void sqrt_in_place();
  void write_data(std::ofstream& data_file) const;

  // Array versions for if we just want to see the values and maybe do basic addition etc.
  void tangent(double *const t) const;
  void normal(double *const n) const;
  void binormal(double *const b) const;

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

quaternion operator /(quaternion q, const double s);
quaternion operator *(quaternion q, const double s);
quaternion operator *(const double s, quaternion q);
quaternion operator *(quaternion p, const quaternion& q);
quaternion operator +(quaternion p, const quaternion& q);
quaternion operator -(quaternion p, const quaternion& q);
std::ostream& operator <<(std::ostream& stream, const quaternion& q);

//
// OTHER FUNCTIONS ASSOCIATED WITH QUATERNIONS
//

void lie_exp(quaternion& q, const double *const u);
quaternion lie_exp(const double *const u);
void midpoint_quaternion(quaternion& qmid, const quaternion& q1, const quaternion& q2);
quaternion midpoint_quaternion(const quaternion& q1, const quaternion& q2);
void dexp(double *const out, const double *const u, const double *const v);
void dexpinv(double *const out, const double *const u, const double *const v);
void dexpinv_transpose(double *const out, const double *const u, const double *const v);
void bch(double *const out, const double *const u, const double *const v);

#endif // MY_QUATERNION_HEADER_INCLUDED
