#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "helper_functions.h"
int factlookup[10] = {1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800};
// gives a random number following a guassian distribution
double gaussian(double sd, int num_additions) {
  double rand = 0.0;
  for (int i = 0; i < num_additions; i++) {
    rand += drand48();
  }
  rand -= (double)num_additions / 2.0;
  rand *= sqrt(12.0) * sd;
  rand /= sqrt((double)num_additions);
  return rand;
}
int num_args(int argc, char **argv) {
  int numargs = 0;
  for (int i = 1; i < argc; i++) {
    if ('-' != argv[i][0]) {
      numargs++;
    }
  }
  return numargs;
}
int num_flags(int argc, char **argv) {
  int numflags = 0;
  for (int i = 1; i < argc; i++) {
    if ('-' == argv[i][0]) {
      numflags++;
    }
  }
  return numflags;
}
char **get_args(int argc, char **argv) {
  char **args = malloc(sizeof(char *) * num_args(argc, argv));
  for (int i = 1, j = 0; i < argc; i++) {
    if ('-' != argv[i][0]) {
      int len = strlen(argv[i]) + 1;
      args[j] = malloc(sizeof(char) * len);
      strcpy(args[j], argv[i]);
      j++;
    }
  }
  return args;
}

char **get_flags(int argc, char **argv) {
  char **flags = malloc(sizeof(char *) * num_flags(argc, argv));
  for (int i = 1, j = 0; i < argc; i++) {
    if ('-' == argv[i][0]) {
      int len = strlen(argv[i]) + 1;
      flags[j] = malloc(sizeof(char) * len);
      strcpy(flags[j], argv[i]);
      j++;
    }
  }
  return flags;
}
int factorial(int num) {
  if (num <= 10) {
    return factlookup[num - 1];
  }
  return num * factorial(num - 1);
}
// both these functions could be sped up massively with lookup tables someday,
// they just have terrible time complexity
perm *first_perm(int n) {
  perm *new_perm = (perm *)malloc(sizeof(perm));
  new_perm->perm = (int *)malloc(sizeof(int) * n);
  new_perm->places = (int *)malloc(sizeof(int) * n);
  new_perm->index = (int *)malloc(sizeof(int) * n);
  new_perm->parity = (int *)calloc(sizeof(int), n);
  new_perm->length = n - 1;
  for (int i = 0; i < n; i++) {
    new_perm->perm[i] = i;
    new_perm->places[i] = i;
    new_perm->index[i] = 0;
  }
  return new_perm;
}
void free_perm(perm *permutation) {
  free(permutation->perm);
  free(permutation->index);
  free(permutation->places);
  free(permutation->parity);
  free(permutation);
}
int increment_perm_index(perm *permutation, int i) {
  if (permutation->index[i] == permutation->length - i) {
    permutation->index[i] = 0;
    permutation->parity[i] = !permutation->parity[i];
    return increment_perm_index(permutation, i + 1);
  }
  permutation->index[i]++;
  return i;
}
void swap(int *a, int *b) {
  int tmp = *a;
  *a = *b;
  *b = tmp;
}
void increment_perm(perm *permutation) {
  int digit = increment_perm_index(permutation, 0);
  int place = permutation->places[digit];
  int parity = permutation->parity[digit];
  if (parity) {
    permutation->places[digit]--;
    permutation->places[permutation->perm[place - 1]]++;
    swap(&permutation->perm[place], &permutation->perm[place - 1]);
  } else {
    permutation->places[digit]++;
    permutation->places[permutation->perm[place + 1]]--;
    swap(&permutation->perm[place], &permutation->perm[place + 1]);
  }
}
void print_perm(perm *permutation) {
  for (int i = 0; i < permutation->length + 1; i++) {
    printf("%i ", permutation->perm[i]);
  }
  printf("\n");
}
histogram *new_histogram(double min, double max, int num_bars) {
  histogram *hist = (histogram *)malloc(sizeof(histogram));
  hist->num_bars = num_bars;
  hist->counts = (uint *)calloc(num_bars, sizeof(uint));
  hist->min = min;
  hist->max = max;
  hist->count = 0;
  return hist;
}
void add_to_histogram(double value, histogram *hist) {
  if (value >= hist->min && value < hist->max)
    hist->counts[(int)(hist->num_bars * (value - hist->min) /
                       (hist->max - hist->min))]++;
  hist->count++;
}
void print_histogram(histogram *hist) {
  double increment = (hist->max - hist->min) / hist->num_bars;
  for (int i = 0; i < hist->num_bars; i++) {
    printf("%lf-%lf: %lf\n", hist->min + i * increment,
           hist->min + (i + 1) * increment,
           ((double)hist->counts[i]) / hist->count);
  }
}
double linear_interpolation(double nums[COLS], double min, double max,
                            double value) {
  double i = (COLS - 1) * (value - min) / (max - min);
  int i_l = (int)i;
  int i_r = i_l + 1;
  double i_space = i - i_l;
  return nums[i_l] * i_space + nums[i_r] * (1.0 - i_space);
}
void printm(int num, uint mod) {
  if (num % mod == 0)
    printf("%i\n", num);
}
void print_double(double numb, FILE *output) {
  fwrite(&numb, sizeof(double), 1, output);
}
void print_int(int numb, FILE *output) {
  fwrite(&numb, sizeof(int), 1, output);
}
