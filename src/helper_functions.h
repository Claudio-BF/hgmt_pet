#ifndef helper_functions_h
#define helper_functions_h
typedef unsigned int uint;

#include <stdio.h>

#include "hgmt_structs.h"
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))
// gives a random number following a guassian distribution
typedef struct perm_ {
  int *perm;
  int *index;
  int *places;
  int *parity;
  int length;
} perm;
double gaussian(double sd, int num_additions);
char **get_args(int argc, char **argv);
char **get_flags(int argc, char **argv);
int num_flags(int argc, char **argv);
int num_args(int argc, char **argv);
int factorial(int num);
void printm(int num, uint mod);
perm *first_perm(int n);
void free_perm(perm *permutation);
void increment_perm(perm *permutation);
void print_perm(perm *permutation);
typedef struct histogram_ {
  uint *counts;
  uint num_bars;
  double min;
  double max;
  uint count;
} histogram;
histogram *new_histogram(double min, double max, int num_bars);
void add_to_histogram(double value, histogram *hist);
void print_histogram(histogram *hist);
double linear_interpolation(double nums[COLS], double min, double max,
                            double value);
void print_double(double numb, FILE *output);
void print_int(int numb, FILE *output);
#endif
