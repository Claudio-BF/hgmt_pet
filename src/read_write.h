#ifndef read_write_h
#define read_write_h
#include "compton_chain_ordering.h"
#include "helper_functions.h"
#include "hgmt_structs.h"
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void read_event(event *new_event, FILE *source, double eff_by_energy[COLS]);
bool read_annihilation(annihilation *new_annihilation, FILE *source,
                       double eff_by_energy[COLS]);
void print_lor(lor *new_lor, FILE *output);
void read_eff(FILE *source, double eff_by_energy[COLS]);
void print_sym_matrix(sym_matrix *mat, FILE *output);
bool read_vec(vec3d *vec, FILE *source);
bool read_sym(sym_matrix *mat, FILE *source);
// prints a vector as the three values
void print_vec(vec3d a, FILE *source);
void print_double(double numb, FILE *output);
void print_int(int numb, FILE *output);
bool read_lor(lor *new_lor, FILE *input);
typedef struct debug_context {
  annihilation *annihil;
  hit_split *split;
  primitive_lor *prim_lor;
  lor *lor;
} debug_context;
void print_history(debug_context context, FILE *output);
#endif
