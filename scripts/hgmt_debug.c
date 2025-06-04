
// including standard files
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// including custom files
#include "../src/compton_chain_ordering.h"
#include "../src/helper_functions.h"
#include "../src/hgmt_structs.h"
#include "../src/llist.h"
#include "../src/vector_ops.h"
FILE *debug_out;
void print_data(double data, FILE *output) {
  fwrite(&data, sizeof(double), 1, output);
}
double get_phi(double cosX, double cosY, vec3d location) {

  double alpha = cosX;
  double beta = cosY;
  int sign_of_gamma = -1;
  double rad_to_deg = 180.0 / M_PI;

  // dir is the momentum unit vector of the gamma; has to check that the square
  // root is not imaginary (happens with floats)
  double operand = 1.0 - (alpha * alpha) + (beta * beta);
  if (operand < 0.0) {
    operand = 0.0;
  }
  vec3d dir = three_vec(alpha, beta, -sqrt(operand) * (float)sign_of_gamma);

  // constructing the normalized normal vector to the detector
  vec3d normal = location;
  normal.z = 0;
  normal = vec_norm(normal);

  // uses vector ops to find the angle between various vectors
  double phi =
      vec_angle(three_vec(0, 0, 1), vec_rejection(dir, normal)) * rad_to_deg;
  return phi;
}
int write_incidence_angle(FILE *source) {
  float junk;
  float cosX;
  float cosY;
  float energy;
  int worked = 0;
  int particle_id;
  float x;
  float y;
  float z;
  worked += fread(&x, sizeof(float), 1, source);
  worked += fread(&y, sizeof(float), 1, source);
  worked += fread(&z, sizeof(float), 1, source);
  worked += fread(&cosX, sizeof(float), 1, source);
  worked += fread(&cosY, sizeof(float), 1, source);
  worked += fread(&energy, sizeof(float), 1, source);
  worked += fread(&junk, sizeof(float), 1, source);
  worked += fread(&particle_id, sizeof(int), 1, source);
  fread(&junk, 2, 1, source);
  if (worked != 8) {
    return 0;
  } else {
    if (particle_id == 22 && energy > 0.510 && energy < 0.511) {
      print_data(get_phi(cosX, cosY, three_vec(x, y, z)), debug_out);
    }
    return 1;
  }
}
event *read_event(FILE *source) {

  uint event_id;
  double energy_deposit;
  float x;
  float y;
  float z;
  // float mx;
  // float my;
  // float mz;
  float tof;
  int particle_type;
  int track_id;
  int worked = 0;

  worked += fread(&event_id, sizeof(uint), 1, source);
  worked += fread(&energy_deposit, sizeof(double), 1, source);
  worked += fread(&x, sizeof(float), 1, source);
  worked += fread(&y, sizeof(float), 1, source);
  worked += fread(&z, sizeof(float), 1, source);
  // worked += fread(&mx, sizeof(float), 1, source);
  // worked += fread(&my, sizeof(float), 1, source);
  // worked += fread(&mz, sizeof(float), 1, source);
  worked += fread(&tof, sizeof(float), 1, source);
  worked += fread(&particle_type, sizeof(int), 1, source);
  worked += fread(&track_id, sizeof(int), 1, source);

  if (worked != 8) {
    return NULL;
  }

  // make a new event to be passed out
  event *new_event = (event *)malloc(sizeof(event));
  new_event->event_id = event_id;
  new_event->energy_deposit = energy_deposit;
  new_event->location = three_vec((double)x, (double)y, (double)z);
  // new_event->momentum = three_vec((double)mx, (double)my, (double)mz);
  new_event->tof = (double)tof;
  new_event->track_id = track_id;
  return new_event;
}
void read_angles(FILE *source) {
  while (write_incidence_angle(source) == 1) {
    continue;
  }
}
void hist_debug(FILE *input, float max_value, int num_bins) {
  histogram *hist = new_histogram(0.0, max_value, num_bins);

  double num;
  bool worked = fread(&num, sizeof(double), 1, input);
  double tot = 0;
  while (worked == 1) {
    tot += num;
    add_to_histogram(num, hist);
    worked = fread(&num, sizeof(double), 1, input);
  }
  printf("number of data points: %i\n", hist->count);
  printf("average: %lf\n", (double)tot / hist->count);
  print_histogram(hist);
}
event *read_history(int event_id, FILE *source) {
  event *new_event = read_event(source);
  while (new_event != NULL && new_event->event_id == event_id) {
    free(new_event);
    new_event = read_event(source);
  }
  return new_event;
}
void phsp_diagnostics(FILE *source) {
  event *new_event = read_event(source);
  int num_events = 0;
  while (new_event != NULL) {
    int id = new_event->event_id;
    free(new_event);
    new_event = read_history(id, source);
    num_events++;
  }
  printf("number of events: %i\n", num_events);
}
int main(int argc, char **argv) {
  char **flags = get_flags(argc, argv);
  char **args = get_args(argc, argv);
  // defines the help function and how to call it (by using -h or --help)
  for (int i = 0; i < num_flags(argc, argv); i++) {
    if (strcmp(flags[i], "-h") == 0) {
      printf("Usage: ./hgmt_debug -task debug_file_loc extra_args\n");
      printf("-h: print this help\n");
      printf("-hi: run with histogram, requires two extra args [max_value] and "
             "[num_bins]\n");
      printf("-p: run diagnostics on phsp file\n");
      printf("-i: run diagnosts on DetectorIn.phsp file, requires arg "
             "[DetectorIn.phsp]\n");
      exit(0);
    }
    if (strcmp(flags[i], "-hi") == 0) {
      printf("running histogram\n");
      FILE *input = fopen(args[0], "rb");
      hist_debug(input, atoi(args[1]), atoi(args[2]));
    }
    if (strcmp(flags[i], "-p") == 0) {
      printf("diagnostics on phsp file\n");
      FILE *input = fopen(args[0], "rb");
      phsp_diagnostics(input);
    }
    if (strcmp(flags[i], "-i") == 0) {
      printf("information about detector in file\n");
      FILE *input = fopen(args[0], "rb");
      debug_out = fopen("debug.data", "wb");
      printf("writing to debug.data\n");
      read_angles(input);
      printf("done\n");
    }
  }
}
