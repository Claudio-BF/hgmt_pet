
// This file reads in all the annihilations, manipulates the data a little bit,
// and organizes it into batches that the neural network can train off of
// efficiently.
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// including custom files
#include "../src/compton_chain_ordering.h"
#include "../src/helper_functions.h"
#include "../src/hgmt_structs.h"
#include "../src/read_write.h"
#define MAX_HITS 16
#define BATCH_SIZE 128
double eff_by_energy[COLS];
uint num_annihilations = 1;
uint capacity[MAX_HITS] = {0};
FILE *output;
void flush_cache(uint num_hits, hit *cache[MAX_HITS][BATCH_SIZE]) {
  uint index = num_hits - 1;
  print_int(num_hits, output);
  for (int i = 0; i < BATCH_SIZE; i++)
    for (int j = 0; j < num_hits; j++) {
      print_double(cache[index][i][j].tof, output);
      print_double(cache[index][i][j].position.x, output);
      print_double(cache[index][i][j].position.y, output);
      print_double(cache[index][i][j].position.z, output);
    }
}
// apply some transformations to make it easier for the network to learn
// additionally mask information the network shouldn't have
void modify_hits(hit *hits, uint num_hits) {
  double least_tof = hits[0].tof;
  vec3d average_pos = hits[0].position;
  for (int i = 1; i < num_hits; i++) {
    least_tof = MIN(least_tof, hits[i].tof);
    average_pos = vec_add(average_pos, hits[i].position);
  }
  average_pos = vec_norm(average_pos);
  double sin_theta = average_pos.y;
  double cos_theta = average_pos.x;
  for (int i = 0; i < num_hits; i++) {
    hits[i].tof -= least_tof;
    // rotate by -theta about z axis
    vec3d pos = hits[i].position;
    double new_x = cos_theta * pos.x + sin_theta * pos.y;
    double new_y = -sin_theta * pos.x + cos_theta * pos.y;
    hits[i].position.x = new_x;
    hits[i].position.y = new_y;
  }
}
void add_cache(hit **hits, uint num_hits, hit *cache[MAX_HITS][BATCH_SIZE]) {
  uint index = num_hits - 1;
  for (int i = 0; i < num_hits; i++)
    cache[index][capacity[index]][i] = *hits[i];
  modify_hits(cache[index][capacity[index]], num_hits);
  capacity[index]++;
  if (capacity[index] == BATCH_SIZE) {
    capacity[index] = 0;
    flush_cache(num_hits, cache);
  }
}
int triangular(int n) { return n * (n + 1) / 2; }
int main(int argc, char **argv) {
  char **flags = get_flags(argc, argv);
  char **args = get_args(argc, argv);
  // handling all flags and arguments
  for (int i = 0; i < num_flags(argc, argv); i++)
    if (strcmp(flags[i], "-h") == 0) {
      printf("Usage: ./training_data [TOPAS_file_position (not the phsp)] "
             "[efficiency_table_position.csv]\n");
      exit(0);
    }

  // checks to make sure you have correct number of args
  if (num_args(argc, argv) != 2) {
    printf("Incorrect number of arguments, two arguments required.\n");
    printf("Use the -h command to get options.\n\n");
    exit(1);
  }
  // reads in efficiency table into 2D array called eff_by_ang
  printf(
      "Neural Network Exam Maker\n\nLoading in '%s' as efficiencies table...\n",
      args[1]);
  FILE *eff_table_file = fopen(args[1], "r");
  read_eff(eff_table_file, eff_by_energy);
  fclose(eff_table_file);

  // opens up a .lor file to output each LOR into
  output = fopen("training.data", "wb");
  FILE *phsp_file = fopen(args[0], "rb");
  printf("Loading in '%s' as the annihilations...\n", args[0]);

  printf("Constructing the data...\n");
  // I want a jagged 3d stack allocated array to store the hits
  // this achieves that effect sort of
  hit all[BATCH_SIZE * MAX_HITS * (MAX_HITS + 1) / 2] = {0};
  hit *cache[MAX_HITS][BATCH_SIZE];
  for (int i = 0; i < MAX_HITS; i++)
    for (int j = 0; j < BATCH_SIZE; j++)
      cache[i][j] = &all[BATCH_SIZE * triangular(i) + j * (i + 1)];
  annihilation new_annihilation;
  bool worked = read_annihilation(&new_annihilation, phsp_file, eff_by_energy);
  while (worked) {
    if (new_annihilation.photon1_path.num_hits >= 1 &&
        new_annihilation.photon1_path.hits[0]->source->number == 0)
      add_cache(new_annihilation.photon1_path.hits,
                new_annihilation.photon1_path.num_hits, cache);
    if (new_annihilation.photon2_path.num_hits >= 1 &&
        new_annihilation.photon2_path.hits[0]->source->number == 0)
      add_cache(new_annihilation.photon2_path.hits,
                new_annihilation.photon2_path.num_hits, cache);
    free_annihilation(&new_annihilation);
    printm(num_annihilations, 1000000);
    num_annihilations++;
    worked = read_annihilation(&new_annihilation, phsp_file, eff_by_energy);
  }
  printf("done!\n");
  return 0;
}
