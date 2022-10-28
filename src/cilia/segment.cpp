// segment.cpp

#include "segment.hpp"
#include "matrix.hpp"
#include "config.hpp"

segment::~segment(){}

segment::segment(){}

void segment::initial_setup(const double *const x_in, const quaternion& q_in, const double *const u_in){

  x[0] = x_in[0];
  x[1] = x_in[1];
  x[2] = x_in[2];

  xm1[0] = x_in[0];
  xm1[1] = x_in[1];
  xm1[2] = x_in[2];

  xm2[0] = x_in[0];
  xm2[1] = x_in[1];
  xm2[2] = x_in[2];

  q = q_in;
  qm1 = q_in;

  u[0] = u_in[0];
  u[1] = u_in[1];
  u[2] = u_in[2];

  um1[0] = u_in[0];
  um1[1] = u_in[1];
  um1[2] = u_in[2];

}

void segment::initial_guess(const int nt){

  // This code only considers tethered filaments, which combined with the robot-arm constraint means that
  // segment positions are never true degrees of freedom in the system. As such, we don't make initial guesses
  // for them. However, this function doesn't just make an initial guess at the solution during this time-step
  // but also handles the time-index shifting, which should still happen for positions.

  xm2[0] = xm1[0];
  xm2[1] = xm1[1];
  xm2[2] = xm1[2];

  xm1[0] = x[0];
  xm1[1] = x[1];
  xm1[2] = x[2];

  if (nt <= NUM_EULER_STEPS){

    um1[0] = u[0];
    um1[1] = u[1];
    um1[2] = u[2];

  } else {

    double temp = 2.0*u[0] - um1[0];
    um1[0] = u[0];
    u[0] = temp;

    temp = 2.0*u[1] - um1[1];
    um1[1] = u[1];
    u[1] = temp;

    temp = 2.0*u[2] - um1[2];
    um1[2] = u[2];
    u[2] = temp;

  }

  qm1 = q;
  lie_exp(q, u);
  q *= qm1;

}

void segment::update(const double *const u_update){

  u[0] += u_update[0];
  u[1] += u_update[1];
  u[2] += u_update[2];

  lie_exp(q, u);
  q *= qm1;

}

void segment::tangent(double *const t) const {

  q.tangent(t);

}

void segment::normal(double *const n) const {

  q.normal(n);

}

void segment::binormal(double *const b) const {

  q.binormal(b);

}

void segment::tangent(matrix& t) const {

  q.tangent(t);

}

void segment::normal(matrix& n) const {

  q.normal(n);

}

void segment::binormal(matrix& b) const {

  q.binormal(b);

}
