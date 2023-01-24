// flow_field_evaluator.hpp

#ifndef MY_FLOW_FIELD_EVALUATOR_HEADER_INCLUDED
#define MY_FLOW_FIELD_EVALUATOR_HEADER_INCLUDED

#include "matrix.hpp"

class flow_field_evaluator {

public:

  // Values to be read in from saved data
  int NBOD;
  int NFIL;
  int NSEG;
  int NBLOB;
  int NFLOW;
  bool QUATERNION_BASED_DATA;
  bool REMOVE_NET_FORCE_STOKESLET;
  bool ADVANCE_TRACER_POSITIONS;
  Real MU;
  Real RSEG;
  Real DT; // Simulation time between data saves, not necessarily the DT of the original simulation.
  matrix fil_refs;
  matrix blob_refs;

  // GPU info
  int num_gpus;

  // CUDA variables
  Real *v_flow_host;
  Real **v_flow_device;

  Real *x_flow_host;
  Real **x_flow_device;

  Real *x_segs_host;
  Real **x_segs_device;

  Real *x_bods_host;
  Real **x_bods_device;

  Real *x_blobs_host;
  Real **x_blobs_device;

  Real *f_segs_host;
  Real **f_segs_device;

  Real *f_blobs_host;
  Real **f_blobs_device;

  int *num_flow_points;

  ~flow_field_evaluator();
  flow_field_evaluator();

  void setup(FILE *fil_ref_file, FILE *blob_ref_file);
  void read_positions_and_forces(FILE *body_state_file, FILE *seg_state_file, FILE *seg_force_file, FILE *blob_force_file);
  void velocities();

};

#endif // MY_FLOW_FIELD_EVALUATOR_HEADER_INCLUDED
