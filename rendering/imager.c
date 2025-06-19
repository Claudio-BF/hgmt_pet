#include <math.h> // Added for fabs and INFINITY
#include <pthread.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// EM algorithm, if you write out the math you find that you need to minimize
// the product of the dot products with the PDFs of the measurements.
#include "../src/read_write.h"
#include "imager.h"

#define THREADS 128
#define CONVERGENCE_THRESHOLD                                                  \
  1e-7 // Stop when the image change is below this value
#define EXP_LOOKUP_RES 1000

double min_exponent = -20;
FILE *output;
FILE *lor_file;
double exp_lookup[EXP_LOOKUP_RES];
graph sensitivity;
uint max_iterations;
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
  for (int i = 0; i < RES; i++) {
    voxels->values[i] = cells->counts[i] * (double)RES / total;
    voxels->values[i] /= sensitivity.values[i];
  }
}
comp_lor get_comp_lor(lor *new_lor) {
  comp_lor new_comp;
  double det = sym_det(&new_lor->covariance);
  // the matrix defining the new dot product is the inverse of the covariance
  new_comp.new_dot = sym_scale(sym_adjucate(&new_lor->covariance), 1 / det);
  new_comp.constant = pow(2 * M_PI, -1.5) * pow(det, -0.5);
  new_comp.center = new_lor->center;
  // now we need to compute the bounding box
  vec3d col1 = three_vec(new_lor->covariance.xx, new_lor->covariance.xy,
                         new_lor->covariance.xz);
  vec3d col2 = three_vec(new_lor->covariance.xy, new_lor->covariance.yy,
                         new_lor->covariance.yz);
  vec3d col3 = three_vec(new_lor->covariance.xz, new_lor->covariance.yz,
                         new_lor->covariance.zz);
  // directions of minimal ratio proj to (x,y,z)/variance
  new_comp.x_range = fabs(
      col1.x * sqrt(-2 * min_exponent /
                    vec_dot(col1, sym_transform(&new_comp.new_dot, col1))));
  new_comp.y_range = fabs(
      col2.y * sqrt(-2 * min_exponent /
                    vec_dot(col2, sym_transform(&new_comp.new_dot, col2))));
  new_comp.z_range = fabs(
      col3.z * sqrt(-2 * min_exponent /
                    vec_dot(col3, sym_transform(&new_comp.new_dot, col3))));
  return new_comp;
}
int bound_below(int num, int min) { return (num < min) ? min : num; }
int bound_above(int num, int max) { return (num > max) ? max : num; }
double fastexp(double x) {
  if (x > 0 || x < min_exponent)
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
  double scale = ((double)RES) / total;
  for (int i = 0; i < RES; i++)
    cells->values[i] *= scale;
}
void *worker(void *arg) {
  thread_data *data = (thread_data *)arg;
  graph *expectation = (graph *)calloc(1, sizeof(graph));
  double x_const = X_LENGTH * (0.5 / X_RES - 0.5);
  double y_const = Y_LENGTH * (0.5 / Y_RES - 0.5);
  double z_const = Z_LENGTH * (0.5 / Z_RES - 0.5);
  double x_increment = ((double)X_LENGTH / X_RES);
  double y_increment = ((double)Y_LENGTH / Y_RES);
  double z_increment = ((double)Z_LENGTH / Z_RES);
  for (int r = 0; r < data->num_lors; r++) {
    comp_lor *current_lor = &data->data[r];
    graph posterior = {0};
    int i_left = floor(
        ((current_lor->center.x - current_lor->x_range) / X_LENGTH + 0.5) *
        X_RES);
    int i_right = floor(
        ((current_lor->center.x + current_lor->x_range) / X_LENGTH + 0.5) *
        X_RES);
    int j_left = floor(
        ((current_lor->center.y - current_lor->y_range) / Y_LENGTH + 0.5) *
        Y_RES);
    int j_right = floor(
        ((current_lor->center.y + current_lor->y_range) / Y_LENGTH + 0.5) *
        Y_RES);
    int k_left = floor(
        ((current_lor->center.z - current_lor->z_range) / Z_LENGTH + 0.5) *
        Z_RES);
    int k_right = floor(
        ((current_lor->center.z + current_lor->z_range) / Z_LENGTH + 0.5) *
        Z_RES);
    i_left = bound_below(i_left, 0);
    i_right = bound_above(i_right, X_RES - 1);
    j_left = bound_below(j_left, 0);
    j_right = bound_above(j_right, Y_RES - 1);
    k_left = bound_below(k_left, 0);
    k_right = bound_above(k_right, Z_RES - 1);
    // printf("%i %i %i %i %i %i\n", i_left, i_right, j_left, j_right, k_left,
    //        k_right);
    double likelihood = 0;
    for (int i = i_left; i <= i_right; i++)
      for (int j = j_left; j <= j_right; j++)
        for (int k = k_left; k <= k_right; k++) {
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
    if (likelihood != 0)
      graph_add(expectation, &posterior, 1 / likelihood);
  }
  return (void *)expectation;
}
int main(int argc, char **argv) {
  char **flags = get_flags(argc, argv);
  char **args = get_args(argc, argv);
  // defines the help function and how to call it (by using -h or --help)
  for (int i = 0; i < num_flags(argc, argv); i++) {
    if (strcmp(flags[i], "-h") == 0) {
      printf("Usage: ./imager [lor_file_location.lor] [output_directory] [max "
             "iterations (0 for none)]\n");
      printf("-h: print this help\n");
      exit(0);
    }
  }
  if (num_args(argc, argv) != 3) {
    printf("Incorrect number of arguments, one argument required.\n");
    printf("Use the -h command to get options.\n\n");
    exit(1);
  }
  max_iterations = strtoul(args[2], NULL, 10);

  lor_file = fopen(args[0], "rb");
  char *filename;
  asprintf(&filename, "%simage.voxels", args[1]);
  output = fopen(filename, "wb");
  free(filename);
  grid cells = {0};
  graph voxels = {0};
  // initializing exp lookup table
  for (int i = 0; i < EXP_LOOKUP_RES; i++) {
    double val = min_exponent - i * min_exponent / (EXP_LOOKUP_RES - 1);
    exp_lookup[i] = exp(val);
  }
  // creating sensitivity index (prob that an LOR generated at that index is
  // detected). This is how one would do attenuation correction, but that's not
  // needed currently
  for (int i = 0; i < RES; i++)
    sensitivity.values[i] = 1;
  print_double(X_LENGTH, output);
  print_double(Y_LENGTH, output);
  print_double(Z_LENGTH, output);
  print_int(X_RES, output);
  print_int(Y_RES, output);
  print_int(Z_RES, output);

  lor new_lor;
  uint num_lors = 0;
  uint num_comps = 0;
  printf("First pass, creating heatmap as first guess...\n\n");
  while (true) {
    if (!read_lor(&new_lor, lor_file))
      break;
    if (add_point(&cells, new_lor.center))
      num_comps++;
    num_lors++;
  }
  generate_image(&cells, &voxels);

  printf("Initializing threads...\n\n");
  thread_data data[THREADS];
  graph *expectation[THREADS];

  for (int i = 0; i < THREADS; i++) {
    data[i].num_lors = num_comps / THREADS + (i < num_comps % THREADS);
    data[i].data = malloc(data[i].num_lors * sizeof(comp_lor));
  }
  rewind(lor_file);
  uint idx = 0;
  for (int i = 0; i < num_lors; i++) {
    int index = idx / THREADS;
    int thread = idx % THREADS;
    read_lor(&new_lor, lor_file);
    if (in_bounds(new_lor.center)) {
      idx++;
      data[thread].data[index] = get_comp_lor(&new_lor);
    }
  }
  for (int i = 0; i < THREADS; i++)
    data[i].current_image = &voxels;

  // --- CONVERGENCE LOOP ---
  int iteration = 0;
  double change = INFINITY;
  graph prev_voxels = {0};

  printf("Starting convergence loop...\n\n");
  while (change > CONVERGENCE_THRESHOLD && iteration < max_iterations) {
    memcpy(prev_voxels.values, voxels.values, RES * sizeof(double));
    pthread_t threads[THREADS];
    for (int j = 0; j < THREADS; j++)
      pthread_create(&threads[j], NULL, worker, &data[j]);
    for (int j = 0; j < THREADS; j++)
      pthread_join(threads[j], (void *)&expectation[j]);
    graph complete_expectation = {0};
    for (int j = 0; j < THREADS; j++)
      graph_add(&complete_expectation, expectation[j], 1.0 / num_lors);
    graph_mult(&voxels, &complete_expectation);
    graph_mult(&voxels, &sensitivity);
    normalize(&voxels);
    double total_diff = 0;
    for (int r = 0; r < RES; r++)
      total_diff += fabs(voxels.values[r] - prev_voxels.values[r]);
    change = total_diff / RES;
    iteration++;
    printf("Iteration %i/%i | Change (MAE): %e\n", iteration, max_iterations,
           change);
    // Free allocated memory for thread data
    for (int i = 0; i < THREADS; i++)
      free(expectation[i]); // expectation are calloc'd in worker
  }

  printf("\n--- Run Finished ---\n");
  if (change <= CONVERGENCE_THRESHOLD) {
    printf("Convergence reached after %d iterations.\n", iteration);
  } else {
    printf("Reached max iterations (%d) without converging.\n", max_iterations);
  }
  fwrite(voxels.values, sizeof(double), RES, output);
  printf("Final image data saved.\n");

  for (int i = 0; i < THREADS; i++)
    free(data[i].data);
  fclose(lor_file);
  fclose(output);
}
