#ifndef hgmt_lor_creator_h
#define hgmt_lor_creator_h

#include "linear_algebra.h"
#include <stdbool.h>
#include <stdio.h>

#define COLS 105
#define SPD_LGHT 29.9792458 // cm/ns
#define REST_ENERGY 511.0   // KeV
#define LONG_UNC 0.1        // cm
#define CIRC_UNC 0.1        // cm
#define RAD_UNC 0.1         // cm
#define TIME_UNC 0.1        // ns
#define DETECTOR_THICKNESS 2.54
#define DETECTOR_SEGMENTATION 0
#define E_MAX 520.0
#define E_MIN 0.0
typedef unsigned int uint;

#define LONG_VAR LONG_UNC *LONG_UNC
#define CIRC_VAR CIRC_UNC *CIRC_UNC
#define RAD_VAR RAD_UNC *RAD_UNC
#define TIME_VAR TIME_UNC *TIME_UNC
static const double detector_positions[12] = {
    45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100}; // MUST BE SORTED
typedef struct event_ {
  double tof;
  double energy;
  vec3d position;
  vec3d direction;
  int detector_id;
  uint primary; // created by a primary scattering (1=gamma_1 and 2=gamma_2)
  uint number;  // 0=first scatter, 1 = second, etc, only valid if above nonzero
  bool detected;
} event;
typedef struct hit_ {
  vec3d position;
  double tof;
  event *source;
} hit;

typedef struct prim_lor_ {
  hit hit1;
  hit hit2;
} prim_lor;

typedef struct _lor {
  vec3d center;
  sym_matrix covariance;
} lor;
typedef struct _photon_path {
  int num_events;
  event **events;
  int num_hits;
  hit **hits;
} photon_path;

typedef struct _annihilation {
  vec3d origin;
  vec3d center;
  double time;
  uint num_events;
  event *events;
  uint num_hits;
  hit *hits;
  photon_path photon1_path;
  photon_path photon2_path;
} annihilation;
hit *event_to_hit(event *single_event, double eff_by_energy[COLS]);
bool is_detected(event *single_event, double eff_by_energy[COLS]);
int get_detector(vec3d position);
void free_annihilation(annihilation *annihilation_pointer);
#endif
