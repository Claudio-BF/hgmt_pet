#ifndef compton_chain_ordering_h
#define compton_chain_ordering_h

#include "hgmt_structs.h"
#include <stdio.h>
typedef struct _hit_split {
  uint num_hits1;
  uint num_hits2;
  hit **hits1;
  hit **hits2;
} hit_split;
hit_split create_hit_split(hit *hits, uint num_hits);
void free_hit_split(hit_split *split);
hit *initial_by_best_order(hit **hits, uint num_hits,
                           double eff_by_energy[COLS]);
hit *initial_by_least_radial(hit **hits, uint num_hits);
hit *initial_by_least_time(hit **hits, uint num_hits);
#endif
