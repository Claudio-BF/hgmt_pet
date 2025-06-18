// including standard files
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>

// including custom files
#include "../neural_network/model_wrapper.hh"
#include "compton_chain_ordering.h"
#include "helper_functions.h"
#include "hgmt_structs.h"
#include "linear_algebra.h"
#define TIME_SD_TOLERANCE 3
#define OPTIMIZE 1
// This is just PCA but instead of using the mean i'm just using 0
vec3d principal_direction(hit *hits, uint num_hits) {
  sym_matrix mat = sym_zero();
  for (int i = 0; i < num_hits; i++) {
    vec3d position = hits[i].position;
    mat.xx += position.x * position.x;
    mat.xy += position.x * position.y;
    mat.xz += position.x * position.z;
    mat.yy += position.y * position.y;
    mat.yz += position.y * position.z;
    mat.zz += position.z * position.z;
  }
  double max_eigenvalue = sym_max_eigenvalue(mat);
  return sym_eigenvector(&mat, max_eigenvalue);
}
hit_split create_hit_split(hit *hits, uint num_hits) {
  hit_split split;
  vec3d dir = principal_direction(hits, num_hits);
  uint num_hits1 = 0;
  uint num_hits2 = 0;
  for (int i = 0; i < num_hits; i++) {
    if (vec_dot(hits[i].position, dir) > 0)
      num_hits1++;
    else
      num_hits2++;
  }
  split.hits1 = (hit **)malloc(sizeof(hit *) * num_hits1);
  split.hits2 = (hit **)malloc(sizeof(hit *) * num_hits2);
  num_hits1 = 0;
  num_hits2 = 0;
  for (int i = 0; i < num_hits; i++) {
    if (vec_dot(hits[i].position, dir) > 0) {
      split.hits1[num_hits1] = &hits[i];
      num_hits1++;
    } else {
      split.hits2[num_hits2] = &hits[i];
      num_hits2++;
    }
  }
  split.num_hits1 = num_hits1;
  split.num_hits2 = num_hits2;
  return split;
}
void free_hit_split(hit_split *split) {
  free(split->hits1);
  free(split->hits2);
}
hit *initial_by_neural_network(hit **hits, uint num_hits) {
  float *input = malloc(sizeof(float) * num_hits * 4);
  double least_tof = hits[0]->tof;
  vec3d average_pos = hits[0]->position;
  for (int i = 1; i < num_hits; i++) {
    least_tof = MIN(least_tof, hits[i]->tof);
    average_pos = vec_add(average_pos, hits[i]->position);
  }
  average_pos = vec_norm(average_pos);
  double sin_theta = average_pos.y;
  double cos_theta = average_pos.x;
  for (int i = 0; i < num_hits; i++) {
    input[4 * i] = hits[i]->tof - least_tof;
    // rotate by -theta about z axis
    vec3d pos = hits[i]->position;
    double new_x = cos_theta * pos.x + sin_theta * pos.y;
    double new_y = -sin_theta * pos.x + cos_theta * pos.y;
    input[4 * i + 1] = new_x;
    input[4 * i + 2] = new_y;
    input[4 * i + 3] = pos.z;
  }
  int index = predict(input, num_hits, 4);
  free(input);
  return hits[index];
}
hit *initial_by_least_radial(hit **hits, uint num_hits) {
  double best_rad = radial_dist(hits[0]->position);
  hit *initial = hits[0];
  for (int i = 1; i < num_hits; i++) {
    double new_rad = radial_dist(hits[i]->position);
    if (new_rad < best_rad) {
      best_rad = new_rad;
      initial = hits[i];
    }
  }
  return initial;
}
hit *initial_by_least_time(hit **hits, uint num_hits) {
  double best_time = hits[0]->tof;
  hit *initial = hits[0];
  for (int i = 1; i < num_hits; i++) {
    double new_time = hits[i]->tof;
    if (new_time < best_time) {
      best_time = new_time;
      initial = hits[i];
    }
  }
  return initial;
}
hit *initial_by_truth(hit **hits, uint num_hits) {
  uint best_number = hits[0]->source->number;
  hit *initial = hits[0];
  for (int i = 1; i < num_hits; i++) {
    uint new_number = hits[i]->source->number;
    if (new_number < best_number) {
      best_number = new_number;
      initial = hits[i];
    }
  }
  return initial;
}
