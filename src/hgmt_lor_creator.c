// including standard files
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// including custom files
#include "compton_chain_ordering.h"
#include "helper_functions.h"
#include "hgmt_structs.h"
#include "linear_algebra.h"

// params
bool writing_to_lor = true;
uint vis_events = 0;
uint counter = 0;
double detector_positions[12] = {45, 50, 55, 60, 65, 70,
                                 75, 80, 85, 90, 95, 100}; // MUST BE SORTED

// global variables
#define NUM_CUTS 6
#define NUM_DEBUG_OPTIONS 5
// cuts are {occured, interacted with something, wasn't inpatient, detected,
// first or second detected}
char *cut_descriptions[] = {"occured",
                            "interected with something",
                            "detected",
                            "wasn't inpatient",
                            "first scatter detected",
                            "first scatter identified"};
uint cuts[NUM_CUTS] = {0};
// dual cuts are the same but require both to happen in an annihilation
uint num_annihilations = 0;
uint max_annihilations = 0;
uint dual_cuts[NUM_CUTS] = {0};
uint num_scatters = 0;
uint num_hits = 0;
event *first_event;
double eff_by_energy[COLS];
bool debug_options[NUM_DEBUG_OPTIONS];
FILE *debug[NUM_DEBUG_OPTIONS];
FILE *visualization;

void print_lor(lor *new_lor, FILE *output) {
  print_vec(new_lor->center, output);
  print_sym_matrix(&new_lor->covariance, output);
}
prim_lor create_prim_lor(hit_split split) {
  prim_lor primitive_lor;
  hit *hit1 = split.hits1[0];
  hit *hit2 = split.hits2[0];
  // hit *hit1 =
  //     initial_by_best_order(split.hits1, split.num_hits1, eff_by_energy);
  // hit *hit2 =
  //     initial_by_best_order(split.hits2, split.num_hits2, eff_by_energy);

  primitive_lor.hit1 = *hit1;
  primitive_lor.hit2 = *hit2;
  return primitive_lor;
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

void read_eff(FILE *source) {
  // loops through all the entries in a row
  for (int i = 0; i < COLS; i++) {
    fscanf(source, "%lf,", &eff_by_energy[i]);
  }
}
// gets the detector an event happened in. return -1 if it didn't happen in a
// detector
int get_detector(vec3d position) {
  double rad_dist = radial_dist(position);
  for (int i = 0; i < sizeof(detector_positions) / sizeof(double); i++) {
    if (rad_dist > detector_positions[i] &&
        rad_dist < detector_positions[i] + DETECTOR_THICKNESS) {
      return i;
    }
  }
  return -1;
}
bool is_detected(event *single_event) {
  if (single_event->detector_id != -1 &&
      drand48() < linear_interpolation(eff_by_energy, E_MIN, E_MAX,
                                       single_event->energy)) {
    return true;
  }
  return false;
}
hit *event_to_hit(event *single_event) {
  vec_mag(three_vec(single_event->position.x, single_event->position.y, 0.0));
  vec3d z_hat = three_vec(0.0, 0.0, 1.0);
  vec3d circ_hat = vec_norm(vec_cross(z_hat, single_event->position));
  vec3d offset = vec_add(vec_scale(z_hat, gaussian(LONG_UNC, 30)),
                         vec_scale(circ_hat, gaussian(CIRC_UNC, 30)));
  hit *new_hit = (hit *)malloc(sizeof(hit));
  new_hit->source = single_event;
  new_hit->position = single_event->position;
  new_hit->tof = single_event->tof + gaussian(TIME_UNC, 30);
  double rad_dist = radial_dist(new_hit->position);
  if (DETECTOR_SEGMENTATION) {
    // we move the radial component to the midpoint of the detector which it hit
    new_hit->position = radial_scale(
        new_hit->position, (detector_positions[single_event->detector_id] +
                            DETECTOR_THICKNESS / 2) /
                               rad_dist);
  } else {
    vec3d r_hat =
        vec_scale(three_vec(new_hit->position.x, new_hit->position.y, 0),
                  1.0 / radial_dist(new_hit->position));
    offset = vec_add(offset, vec_scale(r_hat, gaussian(RAD_UNC, 30)));
  }
  new_hit->position = vec_add(new_hit->position, offset);
  return new_hit;
}
sym_matrix hit_covariance(vec3d position) {
  sym_matrix covariance;
  double rad = vec_mag(position);
  double cos = position.x / rad;
  double sin = position.y / rad;
  covariance.xx = RAD_VAR + (CIRC_VAR - RAD_VAR) * sin * sin;
  covariance.xy = (RAD_VAR - CIRC_VAR) * sin * cos;
  covariance.xz = 0;
  covariance.yy = CIRC_VAR + (RAD_VAR - CIRC_VAR) * sin * sin;
  covariance.yz = 0;
  covariance.zz = LONG_VAR;
  return covariance;
}
lor create_lor(prim_lor *primitive_lor) {

  vec3d a = primitive_lor->hit1.position;
  vec3d b = primitive_lor->hit2.position;
  vec3d c = vec_sub(a, b);
  vec3d geometric_center = vec_scale(vec_add(a, b), 0.5);
  vec3d c_hat = vec_norm(c);
  double delta_t = -(primitive_lor->hit1.tof - primitive_lor->hit2.tof);
  vec3d displacement_from_center = vec_scale(c_hat, SPD_LGHT * delta_t * 0.5);
  vec3d annihilation_loc = vec_add(geometric_center, displacement_from_center);

  lor new_lor;
  new_lor.center = annihilation_loc;
  new_lor.covariance = sym_add(hit_covariance(primitive_lor->hit1.position),
                               hit_covariance(primitive_lor->hit2.position));
  new_lor.covariance =
      sym_add(new_lor.covariance,
              sym_scale(sym_proj(c_hat), SPD_LGHT * TIME_VAR * 0.5));
  return new_lor;
}
void read_event(event *new_event, FILE *source) {
  fread(&new_event->tof, sizeof(double), 1, source);
  fread(&new_event->energy, sizeof(double), 1, source);
  read_vec(&new_event->position, source);
  read_vec(&new_event->direction, source);
  fread(&new_event->primary, sizeof(uint), 1, source);
  new_event->detector_id = get_detector(new_event->position);
  new_event->detected = is_detected(new_event);
}
int compare_hits(const void *hit1, const void *hit2) {
  return ((hit *)hit1)->tof > ((hit *)hit2)->tof;
}
bool read_annihilation(annihilation *new_annihilation, FILE *source) {
  if (!fread(&new_annihilation->time, sizeof(double), 1, source))
    return false;
  read_vec(&new_annihilation->origin, source);
  read_vec(&new_annihilation->center, source);
  fread(&new_annihilation->num_events, sizeof(uint), 1, source);
  uint num_primary1 = 0;
  uint num_primary2 = 0;
  new_annihilation->events =
      (event *)malloc(sizeof(event) * new_annihilation->num_events);
  uint num_hits = 0;
  uint num_primary1_hits = 0;
  uint num_primary2_hits = 0;
  for (int i = 0; i < new_annihilation->num_events; i++) {
    event *next_event = &new_annihilation->events[i];
    read_event(next_event, source);
    if (next_event->primary == 1) {
      next_event->number = num_primary1;
      num_primary1++;
    } else if (next_event->primary == 2) {
      next_event->number = num_primary2;
      num_primary2++;
    }
    if (next_event->detected) {
      num_hits++;
      if (next_event->primary == 1)
        num_primary1_hits++;
      else if (next_event->primary == 2)
        num_primary2_hits++;
    }
  }
  new_annihilation->hits = (hit *)malloc(sizeof(hit) * num_hits);
  new_annihilation->photon1_path.events =
      (event **)malloc(sizeof(event *) * num_primary1);
  new_annihilation->photon2_path.events =
      (event **)malloc(sizeof(event *) * num_primary2);
  new_annihilation->photon1_path.hits =
      (hit **)malloc(sizeof(hit *) * num_primary1_hits);
  new_annihilation->photon2_path.hits =
      (hit **)malloc(sizeof(hit *) * num_primary2_hits);
  num_primary1 = 0;
  num_primary2 = 0;
  num_hits = 0;
  num_primary1_hits = 0;
  num_primary2_hits = 0;
  for (int i = 0; i < new_annihilation->num_events; i++) {
    event *current_event = &new_annihilation->events[i];
    if (current_event->primary == 1) {
      new_annihilation->photon1_path.events[num_primary1] = current_event;
      num_primary1++;
    } else if (current_event->primary == 2) {
      new_annihilation->photon2_path.events[num_primary2] = current_event;
      num_primary2++;
    }
    if (current_event->detected) {
      hit *detector_hit = event_to_hit(current_event);
      new_annihilation->hits[num_hits] = *detector_hit;
      free(detector_hit);
      num_hits++;
    }
  }
  qsort(new_annihilation->hits, num_hits, sizeof(hit), compare_hits);
  for (int i = 0; i < num_hits; i++) {
    hit *current_hit = &new_annihilation->hits[i];
    if (current_hit->source->primary == 1) {
      new_annihilation->photon1_path.hits[num_primary1_hits] = current_hit;
      num_primary1_hits++;
    } else if (current_hit->source->primary == 2) {
      new_annihilation->photon2_path.hits[num_primary2_hits] = current_hit;
      num_primary2_hits++;
    }
  }
  new_annihilation->num_hits = num_hits;
  new_annihilation->photon1_path.num_events = num_primary1;
  new_annihilation->photon2_path.num_events = num_primary2;
  new_annihilation->photon1_path.num_hits = num_primary1_hits;
  new_annihilation->photon2_path.num_hits = num_primary2_hits;
  return true;
}

void free_annihilation(annihilation *annihilation_pointer) {
  free(annihilation_pointer->events);
  free(annihilation_pointer->hits);
  free(annihilation_pointer->photon1_path.events);
  free(annihilation_pointer->photon1_path.hits);
  free(annihilation_pointer->photon2_path.events);
  free(annihilation_pointer->photon2_path.hits);
}
void print_path(photon_path *path) {
  for (int i = 0; i < path->num_events; i++)
    // format: x,y,z, energy deposit, detected
    fprintf(visualization, "%lf %lf %lf %lf %d \n", path->events[i]->position.x,
            path->events[i]->position.y, path->events[i]->position.z,
            path->events[i]->energy, path->events[i]->detected ? 1 : 0);
}
void print_annihilation(annihilation *new_annihilation) {
  fprintf(visualization, "%lf %lf %lf \n\n", new_annihilation->center.x,
          new_annihilation->center.y, new_annihilation->center.z);
  print_path(&new_annihilation->photon1_path);
  fprintf(visualization, "\n");
  print_path(&new_annihilation->photon2_path);
  fprintf(visualization, "\n\n");
}
// provide debug statistics
int debug_path(photon_path *path) {
  if (path->num_events == 0) {
    cuts[0]++;
    return 0;
  }
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
  else
    cut = 4;
  cuts[cut]++;
  return cut;
}
int debug_annihilation(annihilation *new_annihilation) {
  num_scatters += new_annihilation->num_events;
  num_hits += new_annihilation->num_hits;
  if (debug_options[0])
    for (int j = 0; j < new_annihilation->num_events; j++)
      print_double(new_annihilation->events[j].detector_id, debug[0]);

  if (vis_events > 0) {
    print_annihilation(new_annihilation);
    vis_events--;
  }
  int cut1 = debug_path(&new_annihilation->photon1_path);
  int cut2 = debug_path(&new_annihilation->photon2_path);
  int cut = MIN(cut1, cut2);
  dual_cuts[cut]++;

  if (debug_options[4] && cut >= 1) {
    for (int i = 0; i < MIN(new_annihilation->photon1_path.num_events, 4); i++)
      for (int j = 0; j < MIN(new_annihilation->photon2_path.num_events, 4);
           j++) {
        vec3d true_center = new_annihilation->center;
        vec3d loc1 = new_annihilation->photon1_path.events[i]->position;
        vec3d loc2 = new_annihilation->photon2_path.events[j]->position;
        double tof1 = new_annihilation->photon1_path.events[i]->tof;
        double tof2 = new_annihilation->photon2_path.events[j]->tof;
        print_int(i + 1, debug[4]);
        print_int(j + 1, debug[4]);
        print_double(impact_parameter(loc1, loc2, tof1, tof2, true_center),
                     debug[4]);
      }
  }
  return cut;
}
void debug_lor(lor *new_lor, vec3d truecenter) {
  if (debug_options[2]) {
    vec3d dir = sym_eigenvector(&new_lor->covariance,
                                sym_max_eigenvalue(new_lor->covariance));
    print_double(
        vec_mag(vec_rejection(vec_sub(new_lor->center, truecenter), dir)),
        debug[2]);
  }
  if (debug_options[3]) {
    vec3d dir = sym_eigenvector(&new_lor->covariance,
                                sym_max_eigenvalue(new_lor->covariance));
    print_double(
        vec_mag(vec_projection(vec_sub(new_lor->center, truecenter), dir)),
        debug[3]);
  }
}
void debug_prim_lor(prim_lor *new_lor) {
  int a = 0;
  if (new_lor->hit1.source->primary && new_lor->hit1.source->number == 0) {
    a++;
    cuts[4]--;
    cuts[5]++;
  }
  if (new_lor->hit2.source->primary && new_lor->hit2.source->number == 0) {
    a++;
    cuts[4]--;
    cuts[5]++;
  }
  if (a == 2) {
    dual_cuts[4]--;
    dual_cuts[5]++;
  }
}
int main(int argc, char **argv) {
  char **flags = get_flags(argc, argv);
  char **args = get_args(argc, argv);
  // handling all flags and arguments
  for (int i = 0; i < num_flags(argc, argv); i++) {
    if (strcmp(flags[i], "-h") == 0) {
      printf("Usage: ./hgmt_lor_creator [TOPAS_file_position.phsp] "
             "[efficiency_table_position.csv] [output_directory]\n");
      printf("-h: print this help\n");
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
    printf("Use the -h command to get options.\n\n");
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
  FILE *eff_table_file = fopen(args[1], "r");
  read_eff(eff_table_file);
  fclose(eff_table_file);

  // opens up a .lor file to output each LOR into
  FILE *lor_output = NULL;
  if (writing_to_lor) {
    printf("Unable to open output file for writing\n");
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
  FILE *phsp_file = fopen(args[0], "rb");
  printf("Loading in '%s' as the phsp file...\n", args[0]);

  printf("Constructing the hits...\n");
  annihilation new_annihilation;
  bool worked = read_annihilation(&new_annihilation, phsp_file);
  while (worked) {
    hit_split split =
        create_hit_split(new_annihilation.hits, new_annihilation.num_hits);
    if (split.num_hits1 && split.num_hits2) {
      debug_annihilation(&new_annihilation);
      prim_lor primitive_lor = create_prim_lor(split);
      debug_prim_lor(&primitive_lor);
      lor new_lor = create_lor(&primitive_lor);
      if (writing_to_lor)
        print_lor(&new_lor, lor_output);
      debug_lor(&new_lor, new_annihilation.center);
    }
    free_annihilation(&new_annihilation);
    free_hit_split(&split);
    num_annihilations++;
    printm(num_annihilations, 1000000);
    if (max_annihilations != 0 && num_annihilations >= max_annihilations)
      break;
    worked = read_annihilation(&new_annihilation, phsp_file);
  }
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

  for (int i = 0; i < NUM_DEBUG_OPTIONS; i++)
    if (debug_options[i])
      fclose(debug[i]);
  if (visualization != NULL)
    fclose(visualization);
  if (lor_output != NULL)
    fclose(lor_output);
  return 0;
}
