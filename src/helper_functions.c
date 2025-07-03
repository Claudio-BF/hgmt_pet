#include <math.h>
#include <pthread.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "helper_functions.h"
// gives a random number following a guassian distribution, uses Box-Muller
// Transform
double gaussian(double sd) {
  static __thread struct drand48_data rand_buf;
  static __thread int seeded = 0;
  static __thread int cached = 0;
  static __thread double cache;

  if (!seeded) {
    // Use a unique seed per thread: combine time and thread id
    srand48_r(time(NULL) ^ (intptr_t)pthread_self(), &rand_buf);
    seeded = 1;
  }

  if (cached) {
    cached = 0;
    return sd * cache;
  }

  double u1, u2;
  drand48_r(&rand_buf, &u1);
  drand48_r(&rand_buf, &u2);

  double r = sqrt(-2.0 * log(u1));
  double theta = 2.0 * M_PI * u2;

  cache = r * sin(theta);
  cached = 1;
  return sd * r * cos(theta);
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
// both these functions could be sped up massively with lookup tables someday,
// they just have terrible time complexity
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
