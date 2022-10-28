#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdio>
#include <string>
#include "flow_field_evaluator.hpp"

int main(int argc, char** argv){

  //////////////////////////////////////////
  // Read in the parameters
  //////////////////////////////////////////
  FILE *par_file = std::fopen("flow_field_simulation_parameters.input", "r"); // Read-only permission
  char buffer[100];

  std::cout << std::endl << ">>> Reading in parameters..." << std::endl << std::endl;

  char SIMULATION_NAME_CHAR_ARRAY[100];
  std::fscanf(par_file, "%s %s", &buffer, &SIMULATION_NAME_CHAR_ARRAY);
  const std::string SIMULATION_NAME = SIMULATION_NAME_CHAR_ARRAY;
  std::cout << ">>> Input file says SIMULATION_NAME = " << SIMULATION_NAME << std::endl << std::endl;

  // Find the simulation's parameters so we know how big arrays should be etc.
  FILE *sim_par_file = std::fopen((SIMULATION_NAME+".par").c_str(), "r");
  if (sim_par_file == NULL){
    std::cout << ">>> ERROR: Cannot open " << SIMULATION_NAME+".par" << std::endl;
    return -1;
  }
  std::cout << ">>> Found " << SIMULATION_NAME+".par" << std::endl;

  flow_field_evaluator flow;
  int TOTAL_TIME_STEPS, PLOT_FREQUENCY_IN_STEPS;
  double DT;

  std::fscanf(sim_par_file, "%i %s %s", &flow.NFIL, &buffer, &buffer);
  std::cout << "NFIL = " << flow.NFIL << std::endl;
  std::fscanf(sim_par_file, "%i %s %s", &flow.NSEG, &buffer, &buffer);
  std::cout << "NSEG = " << flow.NSEG << std::endl;
  std::fscanf(sim_par_file, "%i %s %s", &flow.NBOD, &buffer, &buffer);
  std::cout << "NBOD = " << flow.NBOD << std::endl;
  std::fscanf(sim_par_file, "%i %s %s", &flow.NBLOB, &buffer, &buffer);
  std::cout << "NBLOB = " << flow.NBLOB;
  if (flow.NBLOB == 0){
    std::cout << " (i.e. a half-space simulation)";
  }
  std::cout << std::endl;
  std::fscanf(sim_par_file, "%lf %s %s", &flow.MU, &buffer, &buffer);
  std::cout << "MU = " << flow.MU << std::endl;
  std::fscanf(sim_par_file, "%lf %s %s", &flow.RSEG, &buffer, &buffer);
  std::cout << "RSEG = " << flow.RSEG << std::endl;
  for (int n = 0; n < 6; n++){ // Throw away the next 6 values, all of which are doubles
    double temp;
    std::fscanf(sim_par_file, "%lf %s %s", &temp, &buffer, &buffer);
  }
  std::fscanf(sim_par_file, "%i %s %s", &TOTAL_TIME_STEPS, &buffer, &buffer);
  std::cout << "TOTAL_TIME_STEPS = " << TOTAL_TIME_STEPS << std::endl;
  std::fscanf(sim_par_file, "%lf %s %s", &DT, &buffer, &buffer);
  std::cout << "DT = " << DT << std::endl;
  std::fscanf(sim_par_file, "%i %s %s", &PLOT_FREQUENCY_IN_STEPS, &buffer, &buffer);
  std::cout << "PLOT_FREQUENCY_IN_STEPS = " << PLOT_FREQUENCY_IN_STEPS << std::endl;
  flow.DT = DT*PLOT_FREQUENCY_IN_STEPS; // Time between saves.
  const int max_saves = 1 + (TOTAL_TIME_STEPS-1)/PLOT_FREQUENCY_IN_STEPS;
  std::cout << std::endl << ">>> Data contains at most " << max_saves << " simulation states." << std::endl;
  std::fclose(sim_par_file);

  // Open the data files we need
  FILE *body_state_file = std::fopen((SIMULATION_NAME+"_body_states.dat").c_str(), "r");
  if (body_state_file == NULL){
    std::cout << std::endl << ">>> ERROR: Cannot open " << SIMULATION_NAME+"_body_states.dat" << std::endl;
    return -1;
  }
  std::cout << std::endl << ">>> Found " << SIMULATION_NAME+"_body_states.dat" << std::endl;

  FILE *seg_state_file = std::fopen((SIMULATION_NAME+"_seg_states.dat").c_str(), "r");
  if (seg_state_file == NULL){
    std::cout << std::endl << ">>> ERROR: Cannot open " << SIMULATION_NAME+"_seg_states.dat" << std::endl;
    return -1;
  }
  std::cout << std::endl << ">>> Found " << SIMULATION_NAME+"_seg_states.dat" << std::endl;

  FILE *seg_force_file = std::fopen((SIMULATION_NAME+"_seg_forces.dat").c_str(), "r");
  if (seg_force_file == NULL){
    std::cout << std::endl << ">>> ERROR: Cannot open " << SIMULATION_NAME+"_seg_forces.dat" << std::endl;
    return -1;
  }
  std::cout << std::endl << ">>> Found " << SIMULATION_NAME+"_seg_forces.dat" << std::endl;

  FILE *fil_ref_file = std::fopen((SIMULATION_NAME+"_fil_references.dat").c_str(), "r");
  if (fil_ref_file == NULL){
    std::cout << std::endl << ">>> ERROR: Cannot open " << SIMULATION_NAME+"_fil_references.dat" << std::endl;
    return -1;
  }
  std::cout << std::endl << ">>> Found " << SIMULATION_NAME+"_fil_references.dat" << std::endl;

  FILE *blob_ref_file, *blob_force_file;
  int blob_force_line_length;
  if (flow.NBLOB != 0){
    blob_force_file = std::fopen((SIMULATION_NAME+"_blob_forces.dat").c_str(), "r");
    if (blob_force_file == NULL){
      std::cout << std::endl << ">>> ERROR: Cannot open " << SIMULATION_NAME+"_blob_forces.dat" << std::endl;
      return -1;
    }
    std::cout << std::endl << ">>> Found " << SIMULATION_NAME+"_blob_forces.dat" << std::endl;

    blob_ref_file = std::fopen((SIMULATION_NAME+"_blob_references.dat").c_str(), "r");
    if (blob_ref_file == NULL){
      std::cout << std::endl << ">>> ERROR: Cannot open " << SIMULATION_NAME+"_blob_references.dat" << std::endl;
      return -1;
    }
    std::cout << std::endl << ">>> Found " << SIMULATION_NAME+"_blob_references.dat" << std::endl;
  }

  // Read in the number of points we should evaluate the flow field at
  std::fscanf(par_file, "%s %i", &buffer, &flow.NFLOW);
  std::cout << std::endl << ">>> Evaluating the flow field at NFLOW = " << flow.NFLOW << " points." << std::endl;

  // Check whether the filament data is quaternion-based or is already given as positions
  int temp;
  std::fscanf(par_file, "%s %i", &buffer, &temp);
  flow.QUATERNION_BASED_DATA = bool(temp);
  if (flow.QUATERNION_BASED_DATA){
    std::cout << std::endl << ">>> Using quaternion-based data." << std::endl;
  } else {
    std::cout << std::endl << ">>> Using direct position data." << std::endl;
  }

  // Check whether we're treating the data as periodic or not
  std::fscanf(par_file, "%s %i", &buffer, &temp);
  const bool TREAT_DATA_AS_PERIODIC = bool(temp);
  int NUMBER_OF_CYCLES_FOR_PERIODIC_DATA;
  std::fscanf(par_file, "%s %i", &buffer, &NUMBER_OF_CYCLES_FOR_PERIODIC_DATA);

  if (!TREAT_DATA_AS_PERIODIC){
    NUMBER_OF_CYCLES_FOR_PERIODIC_DATA = 1;
  }

  // Check if we want to remove the Stokeslet contribution of any net force
  std::fscanf(par_file, "%s %i", &buffer, &temp);
  flow.REMOVE_NET_FORCE_STOKESLET = bool(temp);

  // Check if we should advance the tracer positions.
  // turning this off is good for computing time-averaged flows on a known grid.
  std::fscanf(par_file, "%s %i", &buffer, &temp);
  flow.ADVANCE_TRACER_POSITIONS = bool(temp);

  std::fclose(par_file);

  std::cout << std::endl << ">>> Finished reading in parameters!" << std::endl << std::endl;

  //////////////////////////////////////////
  // Main loop over the saved data
  //////////////////////////////////////////
  flow.setup(fil_ref_file, blob_ref_file);

  for (int cycle = 0; cycle < NUMBER_OF_CYCLES_FOR_PERIODIC_DATA; cycle++){

    for (int n = 0; n < max_saves; n++){

      // Read the timestep at which the simulation saved out of each file.
      // This way, the flow field evaluator class only sees actual state data.
      int timestep_of_save;
      std::fscanf(body_state_file, "%i", &timestep_of_save);
      std::fscanf(seg_state_file, "%i", &timestep_of_save);
      std::fscanf(seg_force_file, "%i", &timestep_of_save);
      if (flow.NBLOB != 0){
        std::fscanf(blob_force_file, "%i", &timestep_of_save);
      }

      flow.read_positions_and_forces(body_state_file, seg_state_file, seg_force_file, blob_force_file);
      flow.velocities();

      // Write the flow data to file and we're done
      std::ofstream flow_location_file(SIMULATION_NAME+"_flow_locations.dat", std::ios::app);
      flow_location_file << timestep_of_save << " ";
      flow_location_file << std::scientific << std::setprecision(10);
      std::ofstream flow_velocities_file(SIMULATION_NAME+"_flow_velocities.dat", std::ios::app);
      flow_velocities_file << timestep_of_save << " ";
      flow_velocities_file << std::scientific << std::setprecision(10);

      for (int m = 0; m < 3*flow.NFLOW; m++){

        flow_location_file << flow.x_flow_host[m] << " ";
        flow_velocities_file << flow.v_flow_host[m] << " ";

      }

      flow_location_file << std::endl;
      flow_location_file.close();
      flow_velocities_file << std::endl;
      flow_velocities_file.close();

      std::cout << "Completed " << n+1 + cycle*max_saves << "/" << max_saves*NUMBER_OF_CYCLES_FOR_PERIODIC_DATA << std::endl;

    }

    // Close and re-open the data files so we loop back through
    std::fclose(body_state_file);
    std::fclose(seg_state_file);
    std::fclose(seg_force_file);

    body_state_file = std::fopen((SIMULATION_NAME+"_body_states.dat").c_str(), "r");
    seg_state_file = std::fopen((SIMULATION_NAME+"_seg_states.dat").c_str(), "r");
    seg_force_file = std::fopen((SIMULATION_NAME+"_seg_forces.dat").c_str(), "r");

    if (flow.NBLOB != 0){

      std::fclose(blob_force_file);

      blob_force_file = std::fopen((SIMULATION_NAME+"_blob_forces.dat").c_str(), "r");

    }

  }

  std::fclose(body_state_file);
  std::fclose(seg_state_file);
  std::fclose(seg_force_file);

  if (flow.NBLOB != 0){

    std::fclose(blob_force_file);

  }

  return 0;

}
