#ifndef imager_h
#define imager_h

#include <stdbool.h>
// including custom files
#include "../src/compton_chain_ordering.h"
#include "../src/helper_functions.h"
#include "../src/hgmt_lor_creator.h"
#include "../src/linear_algebra.h"
#define X_LENGTH 25.0
#define Y_LENGTH 25.0
#define Z_LENGTH 25.0
#define X_RES 200
#define Y_RES 200
#define Z_RES 1
#define RES X_RES *Y_RES *Z_RES
typedef struct grid {
  uint counts[X_RES * Y_RES * Z_RES];
} grid;
typedef struct image {
  double values[X_RES * Y_RES * Z_RES];
} image;
#endif
