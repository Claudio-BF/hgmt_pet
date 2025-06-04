
// including standard files
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>

// including custom files
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
double variance_dist(vec3d loc1, vec3d loc2, double d) {
  double delta_x = loc1.x - loc2.x;
  double delta_y = loc1.y - loc2.y;
  double delta_z = loc1.z - loc2.z;
  double rad1 = radial_dist(loc1);
  double rad2 = radial_dist(loc2);
  double var_tan = CIRC_UNC * CIRC_UNC;
  double cos_theta1 = loc1.x / rad1;
  double sin_theta1 = loc1.y / rad1;
  double cos_theta2 = loc2.x / rad2;
  double sin_theta2 = loc1.y / rad2;
  double var_delta_x =
      cos_theta1 * cos_theta1 * RAD_VAR + sin_theta1 * sin_theta1 * var_tan +
      cos_theta2 * cos_theta2 * RAD_VAR + sin_theta2 * sin_theta2 * var_tan;
  double var_delta_y =
      sin_theta1 * sin_theta1 * RAD_VAR + cos_theta1 * cos_theta1 * var_tan +
      sin_theta2 * sin_theta2 * RAD_VAR + cos_theta2 * cos_theta2 * var_tan;
  double var_delta_z = 2 * LONG_UNC * LONG_UNC;
  double var_d = delta_x * delta_x * var_delta_x / (d * d) +
                 delta_y * delta_y * var_delta_y / (d * d) +
                 delta_z * delta_z * var_delta_z / (d * d);
  return var_d;
}
double time_FOM_cum(hit **hits, int *order, uint num_hits) {
  double FOM = 0;
  double var_t = TIME_UNC * TIME_UNC;
  double var_dt = 0;
  double dt = hits[order[0]]->tof;
  for (int i = 0; i < num_hits - 1; i++) {
    hit *hit1 = hits[order[i]];
    hit *hit2 = hits[order[i + 1]];
    // calculating sigma values
    double d = vec_dist(hit1->position, hit2->position);
    dt += d / SPD_LGHT;
    double var_d = variance_dist(hit1->position, hit2->position, d);
    var_dt += var_d / (SPD_LGHT * SPD_LGHT);
    // calculating FOM
    FOM += fabs(hit2->tof - dt) / sqrt(var_t + var_dt);
  }
  return FOM;
}

double chi_square_FOM(hit **hits, int *order, uint num_hits,
                      double eff_by_energy[COLS]) {
  double FOM = 0;
  double E_i = REST_ENERGY; // KeV
  for (int i = 0; i < num_hits - 1; i++) {
    vec3d hit0 = three_vec(0, 0, 0);
    if (i >= 1)
      hit0 = hits[order[i - 1]]->position;
    vec3d hit1 = hits[order[i]]->position;
    vec3d hit2 = hits[order[i + 1]]->position;
    vec3d a = vec_sub(hit1, hit0);
    vec3d b = vec_sub(hit2, hit1);
    double cos_theta = vec_dot(a, b) / (vec_mag(a) * vec_mag(b));
    double E_f = E_i / (1 + (E_i / REST_ENERGY) * (1.0 - cos_theta));
    double E_dep = E_f - E_i;
    FOM += log(linear_interpolation(eff_by_energy, E_MIN, E_MAX, E_dep));
    // printf("%lf ",log(linear_interpolation(eff_by_energy, E_MIN, E_MAX,
    // E_dep)));
    E_i = E_f;
  }
  // printf("\n");
  return -FOM;
}
double time_FOM(hit **hits, int *order, uint num_hits) {
  double FOM = 0;
  for (int i = 0; i < num_hits - 1; i++) {
    hit *hit1 = hits[order[i]];
    hit *hit2 = hits[order[i + 1]];
    // calculating sigma values
    double d = vec_dist(hit1->position, hit2->position);
    double var_d = variance_dist(hit1->position, hit2->position, d);
    double var_delta_t = TIME_UNC * TIME_UNC * 2;
    double var_dt = var_d / (SPD_LGHT * SPD_LGHT);
    // calculating FOM
    double delta_t = hit2->tof - hit1->tof;
    FOM += fabs(delta_t - d / SPD_LGHT) / sqrt(var_delta_t + var_dt);
  }
  return FOM;
}
hit *initial_by_best_order(hit **hits, uint num_hits,
                           double eff_by_energy[COLS]) {
  if (OPTIMIZE) {
    for (int i = 0; i < num_hits; i++) {
      if (hits[i]->tof - hits[0]->tof > TIME_SD_TOLERANCE * TIME_UNC) {
        num_hits = i + 1;
        break;
      }
    }
  }
  int factor = factorial(num_hits);
  perm *permutation = first_perm(num_hits);
  hit *the_best_hit = hits[0];
  double bestFOM =
      chi_square_FOM(hits, permutation->perm, num_hits, eff_by_energy);
  for (int i = 0; i < factor - 1; i++) {
    increment_perm(permutation);
    double figure_of_merit =
        chi_square_FOM(hits, permutation->perm, num_hits, eff_by_energy);
    if (figure_of_merit < bestFOM) {
      bestFOM = figure_of_merit;
      the_best_hit = hits[permutation->perm[0]];
    } else if (figure_of_merit == bestFOM &&
               the_best_hit->tof > hits[permutation->perm[0]]->tof) {
      the_best_hit = hits[permutation->perm[0]];
    }
  }
  free_perm(permutation);
  return the_best_hit;
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
