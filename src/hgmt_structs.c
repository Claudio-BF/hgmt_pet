#include <math.h>
#include <stdio.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// including custom files
#include "compton_chain_ordering.h"
#include "helper_functions.h"
#include "hgmt_structs.h"
#include "read_write.h"
hit *event_to_hit(event *single_event, double eff_by_energy[COLS]) {
  vec_mag(three_vec(single_event->position.x, single_event->position.y, 0.0));
  vec3d z_hat = three_vec(0.0, 0.0, 1.0);
  vec3d circ_hat = vec_norm(vec_cross(z_hat, single_event->position));
  vec3d offset = vec_add(vec_scale(z_hat, gaussian(LONG_UNC, 30)),
                         vec_scale(circ_hat, gaussian(CIRC_UNC, 30)));
  hit *new_hit = (hit *)malloc(sizeof(hit));
  new_hit->source = single_event;
  new_hit->position = single_event->position;
  new_hit->tof = single_event->tof + gaussian(TIME_UNC, 30);
  double rad_dist = radial_dist(new_hit->position);
  if (DETECTOR_SEGMENTATION) {
    // we move the radial component to the midpoint of the detector which it hit
    new_hit->position = radial_scale(
        new_hit->position, (detector_positions[single_event->detector_id] +
                            DETECTOR_THICKNESS / 2) /
                               rad_dist);
  } else {
    vec3d r_hat =
        vec_scale(three_vec(new_hit->position.x, new_hit->position.y, 0),
                  1.0 / radial_dist(new_hit->position));
    offset = vec_add(offset, vec_scale(r_hat, gaussian(RAD_UNC, 30)));
  }
  new_hit->position = vec_add(new_hit->position, offset);
  return new_hit;
}
bool is_detected(event *single_event, double eff_by_energy[COLS]) {
  if (single_event->detector_id != -1 &&
      drand48() < linear_interpolation(eff_by_energy, E_MIN, E_MAX,
                                       single_event->energy)) {
    return true;
  }
  return false;
}
int get_detector(vec3d position) {
  double rad_dist = radial_dist(position);
  for (int i = 0; i < sizeof(detector_positions) / sizeof(double); i++) {
    if (rad_dist > detector_positions[i] &&
        rad_dist < detector_positions[i] + DETECTOR_THICKNESS) {
      return i;
    }
  }
  return -1;
}
void free_annihilation(annihilation *annihilation_pointer) {
  free(annihilation_pointer->events);
  free(annihilation_pointer->hits);
  free(annihilation_pointer->photon1_path.events);
  free(annihilation_pointer->photon1_path.hits);
  free(annihilation_pointer->photon2_path.events);
  free(annihilation_pointer->photon2_path.hits);
}
