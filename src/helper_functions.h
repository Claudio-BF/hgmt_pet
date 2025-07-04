#ifndef helper_functions_h
#define helper_functions_h

#include "compton_chain_ordering.h"
#include "hgmt_structs.h"
#include <stdbool.h>
#include <stdio.h>
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))
// gives a random number following a guassian distribution
typedef struct perm_ {
  uint *perm;
  uint *index;
  uint *places;
  bool *parity;
  int length;
} perm;
double gaussian(double sd);
char **get_args(int argc, char **argv);
char **get_flags(int argc, char **argv);
int num_flags(int argc, char **argv);
int num_args(int argc, char **argv);
void printm(int num, uint mod);
typedef struct histogram_ {
  uint *counts;
  uint num_bars;
  double min;
  double max;
  uint count;
} histogram;
histogram *new_histogram(double min, double max, int num_bars);
void add_to_histogram(double value, histogram *hist);
double linear_interpolation(double nums[COLS], double min, double max,
                            double value);
void print_histogram(histogram *hist);
#endif
