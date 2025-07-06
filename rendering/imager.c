#include <math.h>
#include <pthread.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../src/read_write.h"
#include "imager.h"

#define THREADS 20
#define EXP_LOOKUP_RES 1000

double min_exponent = -20;
FILE *output;
FILE *lor_file;
double exp_lookup[EXP_LOOKUP_RES];
graph sensitivity;
uint iterations;

void add_grid(grid *cells, uint i, uint j, uint k) {
  cells->counts[i * Y_RES * Z_RES + j * Z_RES + k]++;
}
void add_graph(graph *cells, uint i, uint j, uint k, double val) {
  cells->values[i * Y_RES * Z_RES + j * Z_RES + k] += val;
}
void graph_add(graph *cells, graph *add, double val) {
  for (int i = 0; i < RES; i++)
    cells->values[i] += add->values[i] * val;
}
void graph_mult(graph *cells, graph *mult) {
  for (int i = 0; i < RES; i++)
    cells->values[i] *= mult->values[i];
}
bool in_bounds(vec3d pos) {
  return (fabs(pos.x) < X_LENGTH / 2 && fabs(pos.y) < Y_LENGTH / 2 &&
          fabs(pos.z) < Z_LENGTH / 2);
}
bool add_point(grid *cells, vec3d pos) {
  if (in_bounds(pos)) {
    uint i = floor((pos.x / X_LENGTH + 0.5) * X_RES);
    uint j = floor((pos.y / Y_LENGTH + 0.5) * Y_RES);
    uint k = floor((pos.z / Z_LENGTH + 0.5) * Z_RES);
    add_grid(cells, i, j, k);
    return true;
  }
  return false;
}
void generate_image(grid *cells, graph *voxels) {
  uint total = 0;
  for (int i = 0; i < RES; i++)
    total += cells->counts[i];
  if (total == 0)
    return;
  for (int i = 0; i < RES; i++) {
    voxels->values[i] = cells->counts[i] * (double)RES / total;
    voxels->values[i] /= sensitivity.values[i];
  }
  voxels->background = 1e-6;
}
int bound_below(int num, int min) { return (num < min) ? min : num; }
int bound_above(int num, int max) { return (num > max) ? max : num; }
comp_lor get_comp_lor(lor *new_lor) {
  comp_lor new_comp;
  double det = sym_det(&new_lor->covariance);
  new_comp.new_dot = sym_scale(sym_adjucate(&new_lor->covariance), 1 / det);
  new_comp.constant = pow(2 * M_PI, -1.5) * pow(det, -0.5);
  new_comp.center = new_lor->center;
  vec3d col1 = {new_lor->covariance.xx, new_lor->covariance.xy,
                new_lor->covariance.xz};
  vec3d col2 = {new_lor->covariance.xy, new_lor->covariance.yy,
                new_lor->covariance.yz};
  vec3d col3 = {new_lor->covariance.xz, new_lor->covariance.yz,
                new_lor->covariance.zz};
  double x_range = fabs(
      col1.x * sqrt(-2 * min_exponent /
                    vec_dot(col1, sym_transform(&new_comp.new_dot, col1))));
  double y_range = fabs(
      col2.y * sqrt(-2 * min_exponent /
                    vec_dot(col2, sym_transform(&new_comp.new_dot, col2))));
  double z_range = fabs(
      col3.z * sqrt(-2 * min_exponent /
                    vec_dot(col3, sym_transform(&new_comp.new_dot, col3))));
  int i_left = floor(((new_comp.center.x - x_range) / X_LENGTH + 0.5) * X_RES);
  int i_right = floor(((new_comp.center.x + x_range) / X_LENGTH + 0.5) * X_RES);
  int j_left = floor(((new_comp.center.y - y_range) / Y_LENGTH + 0.5) * Y_RES);
  int j_right = floor(((new_comp.center.y + y_range) / Y_LENGTH + 0.5) * Y_RES);
  int k_left = floor(((new_comp.center.z - z_range) / Z_LENGTH + 0.5) * Z_RES);
  int k_right = floor(((new_comp.center.z + z_range) / Z_LENGTH + 0.5) * Z_RES);
  new_comp.i_left = bound_below(i_left, 0);
  new_comp.i_right = bound_above(i_right, X_RES - 1);
  new_comp.j_left = bound_below(j_left, 0);
  new_comp.j_right = bound_above(j_right, Y_RES - 1);
  new_comp.k_left = bound_below(k_left, 0);
  new_comp.k_right = bound_above(k_right, Z_RES - 1);
  return new_comp;
}
double fastexp(double x) {
  if (x < min_exponent)
    return 0;
  double pos = -(x - min_exponent) * EXP_LOOKUP_RES / min_exponent;
  uint top = ceil(pos);
  uint bottom = floor(pos);
  double weight = pos - (int)pos;
  return weight * exp_lookup[top] + (1 - weight) * exp_lookup[bottom];
}
void normalize(graph *cells) {
  double total = 0;
  for (int i = 0; i < RES; i++)
    total += cells->values[i];
  if (total == 0)
    return;
  double scale = ((double)RES) / total;
  for (int i = 0; i < RES; i++)
    cells->values[i] *= scale;
}
void *worker(void *arg) {
  thread_data *data = (thread_data *)arg;
  graph *result = (graph *)calloc(1, sizeof(graph));
  double one_over_likelihood_sum = 0;

  double x_const = X_LENGTH * (0.5 / X_RES - 0.5);
  double y_const = Y_LENGTH * (0.5 / Y_RES - 0.5);
  double z_const = Z_LENGTH * (0.5 / Z_RES - 0.5);
  double x_increment = ((double)X_LENGTH / X_RES);
  double y_increment = ((double)Y_LENGTH / Y_RES);
  double z_increment = ((double)Z_LENGTH / Z_RES);

  for (int r = 0; r < data->num_lors; r++) {
    comp_lor *current_lor = &data->data[r];
    graph posterior = {0};
    double likelihood = 0;

    for (int i = current_lor->i_left; i <= current_lor->i_right; i++) {
      for (int j = current_lor->j_left; j <= current_lor->j_right; j++) {
        for (int k = current_lor->k_left; k <= current_lor->k_right; k++) {
          int h = i * Y_RES * Z_RES + j * Z_RES + k;
          double x = x_const + i * x_increment;
          double y = y_const + j * y_increment;
          double z = z_const + k * z_increment;
          vec3d displacement = vec_sub(three_vec(x, y, z), current_lor->center);
          double exponent =
              -0.5 * vec_dot(displacement, sym_transform(&current_lor->new_dot,
                                                         displacement));
          double derivative = current_lor->constant * fastexp(exponent);
          likelihood += derivative * data->current_image->values[h];
          posterior.values[h] = derivative;
        }
      }
    }

    likelihood += data->current_image->background;

    if (likelihood > 1e-9) {
      double one_over_likelihood = 1.0 / likelihood;
      graph_add(result, &posterior, one_over_likelihood);
      one_over_likelihood_sum += one_over_likelihood;
    }
  }
  result->background = one_over_likelihood_sum;

  return (void *)result;
}

int main(int argc, char **argv) {
  char **flags = get_flags(argc, argv);
  char **args = get_args(argc, argv);
  for (int i = 0; i < num_flags(argc, argv); i++) {
    if (strcmp(flags[i], "-h") == 0) {
      printf("Usage: ./imager [lor_file_location.lor] [output_directory] "
             "[iterations]\n");
      exit(0);
    }
  }
  if (num_args(argc, argv) != 3) {
    printf("Incorrect number of arguments, three arguments required.\n");
    printf("Use the -h command to get options.\n");
    exit(1);
  }
  iterations = strtoul(args[2], NULL, 10);

  lor_file = fopen(args[0], "rb");
  char *filename;
  asprintf(&filename, "%simage.voxels", args[1]);
  output = fopen(filename, "wb");
  free(filename);
  grid cells = {0};
  graph voxels = {0};
  for (int i = 0; i < EXP_LOOKUP_RES; i++) {
    double val = min_exponent - i * min_exponent / (EXP_LOOKUP_RES - 1);
    exp_lookup[i] = exp(val);
  }
  for (int i = 0; i < RES; i++)
    sensitivity.values[i] = 1;

  print_double(X_LENGTH, output);
  print_double(Y_LENGTH, output);
  print_double(Z_LENGTH, output);
  print_int(iterations + 1, output);
  print_int(X_RES, output);
  print_int(Y_RES, output);
  print_int(Z_RES, output);

  lor new_lor;
  uint num_lors = 0;
  uint num_comps = 0;
  printf("First pass, creating heatmap as first guess...\n");
  while (true) {
    if (!read_lor(&new_lor, lor_file))
      break;
    if (add_point(&cells, new_lor.center))
      num_comps++;
    num_lors++;
  }
  generate_image(&cells, &voxels);

  printf("Initializing threads...\n");
  thread_data data[THREADS];

  for (int i = 0; i < THREADS; i++) {
    data[i].num_lors = num_comps / THREADS + (i < num_comps % THREADS);
    data[i].data = malloc(data[i].num_lors * sizeof(comp_lor));
    data[i].current_image = &voxels;
  }
  rewind(lor_file);
  uint idx = 0;
  for (uint i = 0; i < num_lors; i++) {
    read_lor(&new_lor, lor_file);
    if (in_bounds(new_lor.center)) {
      int thread = idx % THREADS;
      int index_in_thread = idx / THREADS;
      data[thread].data[index_in_thread] = get_comp_lor(&new_lor);
      idx++;
    }
  }

  int iteration = 0;
  double change = INFINITY;
  graph prev_voxels = {0};

  printf("Starting convergence loop...\n\n");
  fwrite(voxels.values, sizeof(double), RES, output);
  for (int i = 0; i < iterations; i++) {
    // --- E-Step ---
    memcpy(prev_voxels.values, voxels.values, RES * sizeof(double));

    pthread_t threads[THREADS];
    for (int j = 0; j < THREADS; j++)
      pthread_create(&threads[j], NULL, worker, &data[j]);

    graph *thread_results[THREADS] = {0};
    graph complete_expectation = {0};
    double total_background_update_sum = 0;

    for (int j = 0; j < THREADS; j++) {
      pthread_join(threads[j], (void *)&thread_results[j]);
      if (thread_results[j] != NULL) {
        graph_add(&complete_expectation, thread_results[j], 1.0);
        total_background_update_sum += thread_results[j]->background;
        free(thread_results[j]);
      }
    }

    // --- M-Step ---
    graph_mult(&voxels, &complete_expectation);
    graph_mult(&voxels, &sensitivity);

    voxels.background =
        voxels.background * total_background_update_sum / num_comps;

    normalize(&voxels);

    double total_diff = 0;
    for (int r = 0; r < RES; r++)
      total_diff += fabs(voxels.values[r] - prev_voxels.values[r]);
    change = total_diff / RES;
    iteration++;
    printf("Iteration %i/%i | Change (MAE): %e | Background: %e\n", iteration,
           iterations, change, voxels.background);

    fwrite(voxels.values, sizeof(double), RES, output);
  }

  printf("\n--- Run Finished ---\n");
  printf("Final image data saved.\n");

  for (int i = 0; i < THREADS; i++)
    free(data[i].data);
  fclose(lor_file);
  fclose(output);
  return 0;
}
