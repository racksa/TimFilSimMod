// quaternion.cpp

#include <cmath>
#include "matrix.hpp"
#include "quaternion.hpp"

#define SMALL_ANGLE 1e-5 // Angle below which we use special versions of dexpinv etc.

using namespace std;

// Constructors, destructors etc.
quaternion::~quaternion(){}

quaternion::quaternion(){

  scalar_part = 1.0;
  vector_part[0] = 0.0;
  vector_part[1] = 0.0;
  vector_part[2] = 0.0;

}

quaternion::quaternion(const double q0, const double q1, const double q2, const double q3){

  scalar_part = q0;
  vector_part[0] = q1;
  vector_part[1] = q2;
  vector_part[2] = q3;

}

quaternion::quaternion(const double q0, const double *const q){

  scalar_part = q0;
  vector_part[0] = q[0];
  vector_part[1] = q[1];
  vector_part[2] = q[2];

}

quaternion::quaternion(const quaternion& q){

  scalar_part = q.scalar_part;
  vector_part[0] = q.vector_part[0];
  vector_part[1] = q.vector_part[1];
  vector_part[2] = q.vector_part[2];

}

//
// OPERATOR OVERLOADING
//

// Element access
double& quaternion::operator()(const int elem){

  return (elem == 0) ? scalar_part : vector_part[elem-1];

}

double quaternion::operator()(const int elem) const {

  return (elem == 0) ? scalar_part : vector_part[elem-1];

}

// Assignment
quaternion& quaternion::operator =(const quaternion& q){

  scalar_part = q.scalar_part;
  vector_part[0] = q.vector_part[0];
  vector_part[1] = q.vector_part[1];
  vector_part[2] = q.vector_part[2];

  return *this;

}

// In-place scalar multiplication
quaternion& quaternion::operator *=(const double s){

  scalar_part *= s;
  vector_part[0] *= s;
  vector_part[1] *= s;
  vector_part[2] *= s;

  return *this;

}

// In-place scalar division
quaternion& quaternion::operator /=(const double s){

  scalar_part /= s;
  vector_part[0] /= s;
  vector_part[1] /= s;
  vector_part[2] /= s;

  return *this;

}

// Sign swap
quaternion quaternion::operator -() const {

  quaternion q = *this;
  q *= -1.0;
  return q;

}

// In-place quaternion multiplication
quaternion& quaternion::operator *=(const quaternion& q){

  // It's possible there's a more efficient way to write this, but this should
  // be "q *= q safe", which is nice.

  const double q0 = scalar_part*q.scalar_part - (vector_part[0]*q.vector_part[0] + vector_part[1]*q.vector_part[1] + vector_part[2]*q.vector_part[2]);
  const double q1 = scalar_part*q.vector_part[0] + q.scalar_part*vector_part[0] + vector_part[1]*q.vector_part[2] - vector_part[2]*q.vector_part[1];
  const double q2 = scalar_part*q.vector_part[1] + q.scalar_part*vector_part[1] + vector_part[2]*q.vector_part[0] - vector_part[0]*q.vector_part[2];
  const double q3 = scalar_part*q.vector_part[2] + q.scalar_part*vector_part[2] + vector_part[0]*q.vector_part[1] - vector_part[1]*q.vector_part[0];

  scalar_part = q0;
  vector_part[0] = q1;
  vector_part[1] = q2;
  vector_part[2] = q3;

  return *this;

}

// In-place addition
quaternion& quaternion::operator +=(const quaternion& q){

  scalar_part += q.scalar_part;
  vector_part[0] += q.vector_part[0];
  vector_part[1] += q.vector_part[1];
  vector_part[2] += q.vector_part[2];

  return *this;

}

// In-place subtraction
quaternion& quaternion::operator -=(const quaternion& q){

  scalar_part -= q.scalar_part;
  vector_part[0] -= q.vector_part[0];
  vector_part[1] -= q.vector_part[1];
  vector_part[2] -= q.vector_part[2];

  return *this;

}

//
// OTHER METHODS
//

  double quaternion::norm() const {

    return sqrt(scalar_part*scalar_part + vector_part[0]*vector_part[0] + vector_part[1]*vector_part[1] + vector_part[2]*vector_part[2]);

  }

  void quaternion::normalise_in_place(){

    (*this) /= this->norm();

  }

  void quaternion::conj_in_place(){

    vector_part[0] *= -1.0;
    vector_part[1] *= -1.0;
    vector_part[2] *= -1.0;

  }

  void quaternion::sqrt_in_place(){

    if (scalar_part < -0.999999999) {

      scalar_part = 0.0;
      vector_part[0] = 0.0;
      vector_part[1] = 0.0;
      vector_part[2] = 1.0;
      // Any unit vector would be technically correct here.
      // If we're taking the sqrt to find the rotation mapping two vectors to one another
      // and we enter this case -- i.e. we're mapping v to -v -- then this unit vector should be
      // orthogonal to v. This function isn't supposed to handle this -- it's primary purpose
      // is to do the sqrts appearing in the elastic moment expression where the rotation should be small --
      // and care should be taken around this case when seeding filements etc. This default choice of the
      // z-direction is so that it is perpendicular to the x-direction, which is mapped to the filament
      // tangent and hence is the v --> -v scenario most likely to occur, hopefully avoiding most issues
      // in practice.

    } else {

      scalar_part = sqrt(0.5*(1.0 + scalar_part));

      const double temp = 2.0*scalar_part;

      vector_part[0] /= temp;
      vector_part[1] /= temp;
      vector_part[2] /= temp;

    }

  }

  void quaternion::write_data(ofstream& data_file) const {

    data_file << scalar_part << " " << vector_part[0] << " " << vector_part[1] << " " << vector_part[2] << " ";

  }

  void quaternion::tangent(double *const t) const {

    t[0] = 1.0 - 2.0*(vector_part[1]*vector_part[1] + vector_part[2]*vector_part[2]);
    t[1] = 2.0*(vector_part[0]*vector_part[1] + scalar_part*vector_part[2]);
    t[2] = 2.0*(vector_part[0]*vector_part[2] - scalar_part*vector_part[1]);

  }

  void quaternion::normal(double *const n) const {

    n[0] = 2.0*(vector_part[0]*vector_part[1] - scalar_part*vector_part[2]);
    n[1] = 1.0 - 2.0*(vector_part[0]*vector_part[0] + vector_part[2]*vector_part[2]);
    n[2] = 2.0*(vector_part[1]*vector_part[2] + scalar_part*vector_part[0]);

  }

  void quaternion::binormal(double *const b) const {

    b[0] = 2.0*(vector_part[0]*vector_part[2] + scalar_part*vector_part[1]);
    b[1] = 2.0*(vector_part[1]*vector_part[2] - scalar_part*vector_part[0]);
    b[2] = 1.0 - 2.0*(vector_part[0]*vector_part[0] + vector_part[1]*vector_part[1]);

  }

  void quaternion::tangent(matrix& t) const {

    this->quaternion::tangent(&t.data[0]);

  }

  void quaternion::normal(matrix& n) const {

    this->quaternion::normal(&n.data[0]);

  }

  void quaternion::binormal(matrix& b) const {

    this->quaternion::binormal(&b.data[0]);

  }

  void quaternion::rot_mat(matrix& R) const {

    R(0,0) = 1.0;
    R(1,1) = 1.0;
    R(2,2) = 1.0;

    double temp = 2.0*vector_part[0]*vector_part[0];
    R(1,1) -= temp;
    R(2,2) -= temp;

    temp = 2.0*vector_part[1]*vector_part[1];
    R(0,0) -= temp;
    R(2,2) -= temp;

    temp = 2.0*vector_part[2]*vector_part[2];
    R(0,0) -= temp;
    R(1,1) -= temp;

    temp = 2.0*vector_part[0]*vector_part[1];
    R(1,0) = temp;
    R(0,1) = temp;

    temp = 2.0*vector_part[0]*vector_part[2];
    R(2,0) = temp;
    R(0,2) = temp;

    temp = 2.0*vector_part[1]*vector_part[2];
    R(1,2) = temp;
    R(2,1) = temp;

    temp = 2.0*scalar_part*vector_part[2];
    R(1,0) += temp;
    R(0,1) -= temp;

    temp = 2.0*scalar_part*vector_part[1];
    R(2,0) -= temp;
    R(0,2) += temp;

    temp = 2.0*scalar_part*vector_part[0];
    R(2,1) += temp;
    R(1,2) -= temp;

  }

  matrix quaternion::rot_mat() const {

    matrix R(3,3);
    this->quaternion::rot_mat(R);
    return R;

  }

  // Psi is the matrix satisfying dq/dt = Psi*Omega when the quaternion is viewed as a 4-vector.
  // It's used as an approximation to dq/du.

  void quaternion::psi_mat(matrix& Psi) const {

    Psi(0,0) = -vector_part[0];
    Psi(0,1) = -vector_part[1];
    Psi(0,2) = -vector_part[2];

    Psi.set_block(1, 0, 3, 3, rcross(vector_part));

    Psi(1,0) += scalar_part;
    Psi(2,1) += scalar_part;
    Psi(3,2) += scalar_part;

    Psi *= 0.5;

  }

  matrix quaternion::psi_mat() const {

    matrix Psi(4,3);
    this->quaternion::psi_mat(Psi);
    return Psi;

  }

  void quaternion::left_mult_mat(matrix& qdot) const {

    qdot(0,0) = scalar_part;

    qdot(0,1) = -vector_part[0];
    qdot(1,0) = vector_part[0];

    qdot(0,2) = -vector_part[1];
    qdot(2,0) = vector_part[1];

    qdot(0,3) = -vector_part[2];
    qdot(3,0) = vector_part[2];

    qdot.set_block(1, 1, 3, 3, -rcross(vector_part));

    qdot(1,1) += scalar_part;
    qdot(2,2) += scalar_part;
    qdot(3,3) += scalar_part;

  }

  matrix quaternion::left_mult_mat() const {

    matrix qdot(4,4);
    this->quaternion::left_mult_mat(qdot);
    return qdot;

  }

  void quaternion::right_mult_mat(matrix& dotq) const {

    dotq(0,0) = scalar_part;

    dotq(0,1) = -vector_part[0];
    dotq(1,0) = vector_part[0];

    dotq(0,2) = -vector_part[1];
    dotq(2,0) = vector_part[1];

    dotq(0,3) = -vector_part[2];
    dotq(3,0) = vector_part[2];

    dotq.set_block(1, 1, 3, 3, rcross(vector_part));

    dotq(1,1) += scalar_part;
    dotq(2,2) += scalar_part;
    dotq(3,3) += scalar_part;

  }

  matrix quaternion::right_mult_mat() const {

    matrix dotq(4,4);
    this->quaternion::right_mult_mat(dotq);
    return dotq;

  }

  //
  // BINARY OPERATOR OVERLOADING
  //

  quaternion operator /(quaternion q, const double s){

    q /= s;
    return q;

  }

  quaternion operator *(quaternion q, const double s){

    q *= s;
    return q;

  }

  quaternion operator *(const double s, quaternion q){

    q *= s;
    return q;

  }

  quaternion operator *(quaternion p, const quaternion& q){

    p *= q;
    return p;

  }

  quaternion operator +(quaternion p, const quaternion& q){

    p += q;
    return p;

  }

  quaternion operator -(quaternion p, const quaternion& q){

    p -= q;
    return p;

  }

  std::ostream& operator <<(std::ostream& stream, const quaternion& q){

    stream << q(0) << " " << q(1) << " " << q(2) << " " << q(3) << " " << endl;

    return stream;

  }

  //
  // OTHER FUNCTIONS ASSOCIATED WITH QUATERNIONS
  //

  void lie_exp(quaternion& q, const double *const u){

    const double theta = sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);

    if (theta < 1e-10){ // No trig involved in the norm(u)-based denominator so we don't have to use SMALL_ANGLE.

      q.scalar_part = 1.0;
      q.vector_part[0] = 0.0;
      q.vector_part[1] = 0.0;
      q.vector_part[2] = 0.0;

    } else {

      const double cs = cos(0.5*theta);
      const double sn = sin(0.5*theta)/theta;

      q.scalar_part = cs;
      q.vector_part[0] = sn*u[0];
      q.vector_part[1] = sn*u[1];
      q.vector_part[2] = sn*u[2];

    }

  }

  quaternion lie_exp(const double *const u){

    quaternion q;
    lie_exp(q, u);
    return q;

  }

  void midpoint_quaternion(quaternion& qmid, const quaternion& q1, const quaternion& q2){

    // qmid = (q2 * conj(q1))^0.5 * q1
    qmid = q1;
    qmid.conj_in_place();
    qmid = q2*qmid;
    qmid.sqrt_in_place();
    qmid *= q1;

  }

  quaternion midpoint_quaternion(const quaternion& q1, const quaternion& q2){

    quaternion qmid;
    midpoint_quaternion(qmid, q1, q2);
    return qmid;

  }

  void dexp(double *const out, const double *const u, const double *const v){

    const double theta2 = u[0]*u[0] + u[1]*u[1] + u[2]*u[2];
    const double theta = sqrt(theta2);

    double alpha, beta;

    if (theta < SMALL_ANGLE){

      // Rather than use a Taylor series approx., this is constructed to be the inverse of the
      // Taylor series approx. we use for dexpinv at small angles; i.e. to ensure that dexp(u, dexpinv(u,v))
      // = dexpinv(u, dexp(u,v)) = v holds for any choice of u and v.
      const double theta4 = theta2*theta2;

      beta = (1.0/3.0 + theta2/72.0)/(2.0 + theta2/6.0 + theta4/72.0);
      alpha = 1.0/6.0 + 2.0*beta*(1.0 - theta2/12.0);

    } else {

      const double theta3 = theta*theta2;

      alpha = sin(0.5*theta);
      alpha *= 2.0*alpha/theta2;

      beta = (theta - sin(theta))/theta3;

    }

    const double u_cross_v[3] = {u[1]*v[2] - u[2]*v[1], u[2]*v[0] - u[0]*v[2], u[0]*v[1] - u[1]*v[0]};

    out[0] = v[0] + alpha*u_cross_v[0] + beta*(u[1]*u_cross_v[2] - u[2]*u_cross_v[1]);
    out[1] = v[1] + alpha*u_cross_v[1] + beta*(u[2]*u_cross_v[0] - u[0]*u_cross_v[2]);
    out[2] = v[2] + alpha*u_cross_v[2] + beta*(u[0]*u_cross_v[1] - u[1]*u_cross_v[0]);

  }

  void dexpinv(double *const out, const double *const u, const double *const v){

    const double theta2 = u[0]*u[0] + u[1]*u[1] + u[2]*u[2];
    const double theta = sqrt(theta2);

    double fac;

    if (theta < SMALL_ANGLE){

      fac = -1.0/12.0;

    } else {

      fac = (0.5*theta*cos(0.5*theta)/sin(0.5*theta) - 1.0)/theta2;

    }

    const double u_cross_v[3] = {u[1]*v[2] - u[2]*v[1], u[2]*v[0] - u[0]*v[2], u[0]*v[1] - u[1]*v[0]};

    out[0] = v[0] - 0.5*u_cross_v[0] - fac*(u[1]*u_cross_v[2] - u[2]*u_cross_v[1]);
    out[1] = v[1] - 0.5*u_cross_v[1] - fac*(u[2]*u_cross_v[0] - u[0]*u_cross_v[2]);
    out[2] = v[2] - 0.5*u_cross_v[2] - fac*(u[0]*u_cross_v[1] - u[1]*u_cross_v[0]);

  }

  void dexpinv_transpose(double *const out, const double *const u, const double *const v){

    const double theta2 = u[0]*u[0] + u[1]*u[1] + u[2]*u[2];
    const double theta = sqrt(theta2);

    double fac;

    if (theta < SMALL_ANGLE){

      fac = -1.0/12.0;

    } else {

      fac = (0.5*theta*cos(0.5*theta)/sin(0.5*theta) - 1.0)/theta2;

    }

    const double u_cross_v[3] = {u[1]*v[2] - u[2]*v[1], u[2]*v[0] - u[0]*v[2], u[0]*v[1] - u[1]*v[0]};

    out[0] = v[0] + 0.5*u_cross_v[0] - fac*(u[1]*u_cross_v[2] - u[2]*u_cross_v[1]);
    out[1] = v[1] + 0.5*u_cross_v[1] - fac*(u[2]*u_cross_v[0] - u[0]*u_cross_v[2]);
    out[2] = v[2] + 0.5*u_cross_v[2] - fac*(u[0]*u_cross_v[1] - u[1]*u_cross_v[0]);

  }

  void bch(double *const out, const double *const u, const double *const v){

    // Implements the closed-form expression for the BCH formula from Engo (2001).
    // As with everything Lie algebra related, this expression is only valid locally
    // and it's going to hit the fan really quickly if either of these rotations are
    // through more than pi/2.

    const double theta = sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);
    const double phi = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);

    if (theta < SMALL_ANGLE){

      out[0] = v[0];
      out[1] = v[1];
      out[2] = v[2];

    } else if (phi < SMALL_ANGLE){

      out[0] = u[0];
      out[1] = u[1];
      out[2] = u[2];

    } else {

      const double dir_dot = (u[0]*v[0] + u[1]*v[1] + u[2]*v[2])/(theta*phi);
      const double ang = acos(dir_dot);

      const double st = sin(theta);
      const double sp = sin(phi);
      const double sa = sin(ang);

      double csq_half_phi = cos(0.5*phi);
      csq_half_phi *= csq_half_phi;
      const double ssq_half_phi = 1.0 - csq_half_phi;
      double csq_half_theta = cos(0.5*theta);
      csq_half_theta *= csq_half_theta;
      const double ssq_half_theta = 1.0 - csq_half_theta;

      const double a = st*csq_half_phi - sp*ssq_half_theta*dir_dot;
      const double b = sp*csq_half_theta - st*ssq_half_phi*dir_dot;
      const double c = 0.5*st*sp - 2.0*ssq_half_theta*ssq_half_phi*dir_dot;
      const double d = sqrt(a*a + b*b + 2.0*a*b*dir_dot + c*c*sa*sa);

      double alpha = asin(d)/d;
      double beta = alpha;
      double gamma = alpha;
      alpha *= a/theta;
      beta *= b/phi;
      gamma *= c/(theta*phi);

      out[0] = alpha*u[0] + beta*v[0] + gamma*(u[1]*v[2] - u[2]*v[1]);
      out[1] = alpha*u[1] + beta*v[1] + gamma*(u[2]*v[0] - u[0]*v[2]);
      out[2] = alpha*u[2] + beta*v[2] + gamma*(u[0]*v[1] - u[1]*v[0]);

    }

  }
