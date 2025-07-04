// including standard files
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// including custom files
#include "compton_chain_ordering.h"
#include "helper_functions.h"
#include "hgmt_structs.h"
#include "linear_algebra.h"
#include "read_write.h"

#define THREADS 20
double eff_by_energy[COLS];
// params
bool writing_to_lor = true;
uint vis_events = 0;
uint counter = 0;
// global variables
#define NUM_CUTS 6
#define NUM_DEBUG_OPTIONS 5
// cuts are {occured, interacted with something, wasn't inpatient, detected,
// first or second detected}
char *cut_descriptions[] = {
    "occured",
    "interected with something",
    "something detected",
    "wasn't inpatient",
    "first scatter detected",
    "first scatter identified (requires lor to be made)"};
uint cuts[NUM_CUTS] = {0};
// dual cuts are the same but require both to happen in an annihilation
uint num_annihilations = 0;
uint max_annihilations = 0;
uint dual_cuts[NUM_CUTS] = {0};
uint num_scatters = 0;
uint num_hits = 0;
event *first_event;
bool debug_options[NUM_DEBUG_OPTIONS];
FILE *debug[NUM_DEBUG_OPTIONS];
FILE *visualization;
FILE *lor_output;
FILE *eff_table_file;
FILE *input_file;
pthread_mutex_t read_lock = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t write_lock = PTHREAD_MUTEX_INITIALIZER;
primitive_lor create_prim_lor(hit_split split) {
  primitive_lor prim_lor;
  hit *hit1 = initial_by_neural_network(split.hits1, split.num_hits1);
  hit *hit2 = initial_by_neural_network(split.hits2, split.num_hits2);
  // hit *hit1 = initial_by_truth(split.hits1, split.num_hits1);
  // hit *hit2 = initial_by_truth(split.hits2, split.num_hits2);
  // hit *hit1 = initial_by_least_time(split.hits1, split.num_hits1);
  // hit *hit2 = initial_by_least_time(split.hits2, split.num_hits2);
  prim_lor.hit1 = *hit1;
  prim_lor.hit2 = *hit2;
  return prim_lor;
}

double impact_parameter(vec3d loc1, vec3d loc2, double tof1, double tof2,
                        vec3d true_center) {
  vec3d c = vec_sub(loc1, loc2);
  vec3d geometric_center = vec_add(loc2, vec_scale(c, 0.5));
  vec3d c_hat = vec_norm(c);
  double delta_t = tof2 - tof1;
  vec3d displacement_from_center = vec_scale(c_hat, SPD_LGHT * delta_t * 0.5);
  vec3d estimated_loc = vec_add(geometric_center, displacement_from_center);
  return vec_mag(vec_rejection(vec_sub(estimated_loc, true_center), c));
}

// gets the detector an event happened in. return -1 if it didn't happen in a
// detector
sym_matrix hit_covariance(vec3d position) {
  sym_matrix covariance;
  double rad = vec_mag(position);
  double cos = position.x / rad;
  double sin = position.y / rad;
  covariance.xx = RAD_VAR * cos * cos + CIRC_VAR * sin * sin;
  covariance.xy = (RAD_VAR - CIRC_VAR) * sin * cos;
  covariance.xz = 0;
  covariance.yy = CIRC_VAR * cos * cos + RAD_VAR * sin * sin;
  covariance.yz = 0;
  covariance.zz = LONG_VAR;
  return covariance;
}
lor create_lor(primitive_lor *prim_lor) {

  vec3d a = prim_lor->hit1.position;
  vec3d b = prim_lor->hit2.position;
  vec3d c = vec_sub(a, b);
  vec3d geometric_center = vec_scale(vec_add(a, b), 0.5);
  vec3d c_hat = vec_norm(c);
  double delta_t = -(prim_lor->hit1.tof - prim_lor->hit2.tof);
  vec3d displacement_from_center = vec_scale(c_hat, SPD_LGHT * delta_t * 0.5);
  vec3d annihilation_loc = vec_add(geometric_center, displacement_from_center);

  lor new_lor;
  new_lor.center = annihilation_loc;
  new_lor.covariance = sym_add(hit_covariance(prim_lor->hit1.position),
                               hit_covariance(prim_lor->hit2.position));
  new_lor.covariance =
      sym_add(new_lor.covariance,
              sym_scale(sym_proj(c_hat), SPD_LGHT * TIME_VAR * 0.5));
  new_lor.covariance = sym_add(new_lor.covariance, sym_id(DIFFUSION_VARIANCE));
  return new_lor;
}

// provide debug statistics
int debug_path(photon_path *path, debug_context context) {
  if (path->num_events == 0)
    return 0;
  // getting all the important statistics
  if (debug_options[1])
    print_double(path->events[0]->detector_id, debug[1]);
  // figure out which cut the photon got to, format is: if (not cut n) cut=n-1
  int cut;
  if (path->num_hits == 0)
    cut = 1;
  else if (path->events[0]->detector_id == -1)
    cut = 2;
  else if (!path->events[0]->detected)
    cut = 3;
  else if (context.prim_lor == NULL ||
           (path->events[0] != context.prim_lor->hit1.source &&
            path->events[0] != context.prim_lor->hit2.source))
    cut = 4;
  else
    cut = 5;
  return cut;
}
int debug_annihilation(debug_context context) {
  annihilation *annihil = context.annihil;
  num_scatters += annihil->num_events;
  num_hits += annihil->num_hits;
  if (debug_options[0])
    for (int j = 0; j < annihil->num_events; j++)
      print_double(annihil->events[j].detector_id, debug[0]);

  if (vis_events > 0) {
    print_annihilation(annihil, visualization);
    vis_events--;
  }
  int cut1 = debug_path(&annihil->photon1, context);
  int cut2 = debug_path(&annihil->photon2, context);
  int cut = MIN(cut1, cut2);
  cuts[cut1]++;
  cuts[cut2]++;
  dual_cuts[cut]++;

  if (debug_options[4] && cut >= 1) {
    for (int i = 0; i < MIN(annihil->photon1.num_events, 4); i++)
      for (int j = 0; j < MIN(annihil->photon2.num_events, 4); j++) {
        vec3d true_center = annihil->center;
        vec3d loc1 = annihil->photon1.events[i]->position;
        vec3d loc2 = annihil->photon2.events[j]->position;
        double tof1 = annihil->photon1.events[i]->tof;
        double tof2 = annihil->photon2.events[j]->tof;
        print_int(i + 1, debug[4]);
        print_int(j + 1, debug[4]);
        print_double(impact_parameter(loc1, loc2, tof1, tof2, true_center),
                     debug[4]);
      }
  }
  return cut;
}
void debug_lor(debug_context context) {
  if (debug_options[2] || debug_options[3]) {
    lor *new_lor = context.lor;
    vec3d truecenter = context.annihil->center;
    vec3d dir = vec_norm(vec_sub(context.prim_lor->hit1.position,
                                 context.prim_lor->hit2.position));
    if (debug_options[2])
      print_double(
          vec_mag(vec_rejection(vec_sub(new_lor->center, truecenter), dir)),
          debug[2]);
    if (debug_options[3])
      print_double(
          vec_mag(vec_projection(vec_sub(new_lor->center, truecenter), dir)),
          debug[3]);
  }
}
void debug_prim_lor(debug_context context) {}
void debug_all(debug_context context) {
  debug_annihilation(context);
  num_annihilations++;
  printm(num_annihilations, 1000000);
  if (context.prim_lor != NULL) {
    debug_prim_lor(context);
    debug_lor(context);
  }
}
void *worker(void *arg) {
  annihilation annihil;
  while (true) {
    pthread_mutex_lock(&read_lock);
    bool worked =
        read_annihilation(&annihil, input_file, eff_by_energy) &&
        (max_annihilations == 0 || num_annihilations < max_annihilations);
    pthread_mutex_unlock(&read_lock);
    if (!worked)
      return NULL;
    hit_split split = create_hit_split(annihil.hits, annihil.num_hits);
    debug_context context = {0};
    context.annihil = &annihil;
    context.split = &split;
    primitive_lor prim_lor;
    lor new_lor;
    if (split.num_hits1 >= 1 && split.num_hits2 >= 1) {
      prim_lor = create_prim_lor(split);
      new_lor = create_lor(&prim_lor);
      context.prim_lor = &prim_lor;
      context.lor = &new_lor;
    }
    pthread_mutex_lock(&write_lock);
    debug_all(context);
    if (context.prim_lor != NULL) {
      if (writing_to_lor)
        print_lor(&new_lor, lor_output);
    }
    pthread_mutex_unlock(&write_lock);
    free_annihilation(&annihil);
    free_hit_split(&split);
  }
}
int main(int argc, char **argv) {
  char **flags = get_flags(argc, argv);
  char **args = get_args(argc, argv);
  // handling all flags and arguments
  for (int i = 0; i < num_flags(argc, argv); i++) {
    if (strcmp(flags[i], "-h") == 0) {
      printf("Usage: ./hgmt_lor_creator [TOPAS_file_position (not the phsp)] "
             "[efficiency_table_position.csv] [output_directory]\n");
      printf("-d: run in debug mode, do not write to lor file\n");
      printf("-v#: visualize # events\n");
      printf("-m#: only read a maximum of # annihilations\n");
      printf("-e#: run with debug option #\n");
      printf("\t0: histogram of detector vs number of scatters\n");
      printf("\t1: histogram of detector vs number of first scatters\n");
      printf("\t2: lor reconstruction error to real center (transverse)\n");
      printf("\t3: lor reconstruction error to real center (longitudinal)\n");
      printf("\t4: scatter number vs reconstruction error truth study\n");
      exit(0);
    } else if (strcmp(flags[i], "-d") == 0) {
      printf("running in debug mode, won't write to a lor file\n");
      writing_to_lor = false;
    } else if (strncmp(flags[i], "-e", 2) == 0) {
      uint debug_option;
      sscanf(flags[i], "-e%u", &debug_option);
      debug_options[debug_option] = true;
    } else if (strncmp(flags[i], "-m", 2) == 0) {
      sscanf(flags[i], "-m%u", &max_annihilations);
    } else if (strncmp(flags[i], "-v", 2) == 0) {
      sscanf(flags[i], "-v%u", &vis_events);
      printf("outputting data to visualize %u events\n", vis_events);
    }
  }

  // checks to make sure you have correct number of args
  if (num_args(argc, argv) != 3) {
    printf("Incorrect number of arguments, three arguments required.\n");
    printf("Use the -h command to get options.\n");
    exit(1);
  }
  // opens files for debug output
  printf("running with debug options:"); // Output: 42
  for (int i = 0; i < NUM_DEBUG_OPTIONS; i++)
    if (debug_options[i]) {
      printf(" %i", i);
      char *filename;
      asprintf(&filename, "%sdebug%d.data", args[2], i);
      debug[i] = fopen(filename, "wb");
      free(filename);
    }
  printf("\n");
  // reads in efficiency table into 2D array called eff_by_ang
  printf("HGMT LOR Creator\n\nLoading in '%s' as efficiencies table...\n",
         args[1]);
  eff_table_file = fopen(args[1], "r");
  read_eff(eff_table_file, eff_by_energy);
  fclose(eff_table_file);

  // opens up a .lor file to output each LOR into
  if (writing_to_lor) {
    char *lor_file_loc;
    asprintf(&lor_file_loc, "%sHGMTDerenzo.lor", args[2]);
    lor_output = fopen(lor_file_loc, "wb");
    free(lor_file_loc);
  }
  if (vis_events) {
    char *filename;
    asprintf(&filename, "%svisualization.data", args[2]);
    visualization = fopen(filename, "w");
    free(filename);
  }
  input_file = fopen(args[0], "rb");
  printf("Loading in '%s' as the annihilations...\n", args[0]);

  printf("Constructing the hits...\n");
  // everythings happens here, process and read the annihilations
  pthread_t threads[THREADS];
  for (int i = 0; i < THREADS; i++)
    pthread_create(&threads[i], NULL, worker, NULL);
  for (int i = 0; i < THREADS; i++)
    pthread_join(threads[i], NULL);
  printf("\n");
  // fixing cuts formating to be cumulative
  for (int i = NUM_CUTS - 2; i >= 0; i--) {
    cuts[i] += cuts[i + 1];
    dual_cuts[i] += dual_cuts[i + 1];
  }
  printf("total annihilations: %u\n", num_annihilations);
  printf("total scatters: %u\n", num_scatters);
  printf("total hits: %u\n\n", num_hits);
  printf(
      "(DUAL)CUT 'N': 'num' 'percent passing' 'cumulative percent passing'\n");
  for (int i = 1; i < NUM_CUTS; i++)
    printf("%u: %s\n", i, cut_descriptions[i]);
  printf("\n");
  for (int i = 1; i < NUM_CUTS; i++)
    printf("CUT %u: %u %lf %lf\n", i, cuts[i], (double)cuts[i] / cuts[i - 1],
           (double)cuts[i] / cuts[0]);
  printf("\n");
  for (int i = 1; i < NUM_CUTS; i++)
    printf("DUALCUT %u: %u %lf %lf\n", i, dual_cuts[i],
           (double)dual_cuts[i] / dual_cuts[i - 1],
           (double)dual_cuts[i] / dual_cuts[0]);
  // closing stuff out
  pthread_mutex_destroy(&read_lock);
  pthread_mutex_destroy(&write_lock);
  for (int i = 0; i < NUM_DEBUG_OPTIONS; i++)
    if (debug_options[i])
      fclose(debug[i]);
  if (visualization != NULL)
    fclose(visualization);
  if (lor_output != NULL)
    fclose(lor_output);
  return 0;
}
