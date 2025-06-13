
// including standard files
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// currently very simple algorithm, divides the image into boxes and adds up
// hits in each box. This works surprisingly well (better than the previous
// solution), however an EM algorithm utilizing the covariances would be ideal.
#include "../src/helper_functions.h"
#include "../src/read_write.h"
#include "imager.h"
void add_grid(grid *cells, uint i, uint j, uint k) {
  cells->counts[i * Y_RES * Z_RES + j * Z_RES + k]++;
}
bool add_point(grid *cells, vec3d pos) {
  if (fabs(pos.x) < X_LENGTH / 2 && fabs(pos.y) < Y_LENGTH / 2 &&
      fabs(pos.z) < Z_LENGTH / 2) {
    uint i = floor((pos.x / X_LENGTH + 0.5) * X_RES);
    uint j = floor((pos.y / Y_LENGTH + 0.5) * Y_RES);
    uint k = floor((pos.z / Z_LENGTH + 0.5) * Z_RES);
    add_grid(cells, i, j, k);
    return true;
  }
  return false;
}
void generate_image(grid *cells, image *pixels) {
  uint min = cells->counts[0];
  uint total = cells->counts[0];
  for (int i = 1; i < RES; i++) {
    min = MIN(min, cells->counts[i]);
    total += cells->counts[i];
    // printf("%u\n", cells.counts[i]);
  }
  double scale = (double)RES / (total - RES * min);
  for (int i = 0; i < RES; i++) {
    pixels->values[i] = (cells->counts[i] - min) * scale;
    // printf("%lf\n", pixels.values[i]);
  }
}

int main(int argc, char **argv) {
  char **flags = get_flags(argc, argv);
  char **args = get_args(argc, argv);
  // defines the help function and how to call it (by using -h or --help)
  for (int i = 0; i < num_flags(argc, argv); i++) {
    if (strcmp(flags[i], "-h") == 0) {
      printf("Usage: ./imager [lor_file_location.lor] [output_directory]\n");
      printf("-h: print this help\n");
      exit(0);
    }
  }
  if (num_args(argc, argv) != 2) {
    printf("Incorrect number of arguments, one argument required.\n");
    printf("Use the -h command to get options.\n\n");
    exit(1);
  }

  FILE *lor_file = fopen(args[0], "rb");
  char *filename;
  asprintf(&filename, "%simage.pixels", args[1]);
  FILE *output = fopen(filename, "wb");
  free(filename);
  grid cells = {0};
  image pixels = {0};
  print_double(X_LENGTH, output);
  print_double(Y_LENGTH, output);
  print_double(Z_LENGTH, output);
  print_int(X_RES, output);
  print_int(Y_RES, output);
  print_int(Z_RES, output);

  lor new_lor;
  bool worked = read_lor(&new_lor, lor_file);
  while (worked) {
    add_point(&cells, new_lor.center);
    worked = read_lor(&new_lor, lor_file);
  }
  generate_image(&cells, &pixels);
  fwrite(pixels.values, sizeof(double), RES, output);
  fclose(lor_file);
  fclose(output);
}
