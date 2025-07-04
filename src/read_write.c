#include "read_write.h"
#include "hgmt_structs.h"

void read_event(event *new_event, FILE *source, double eff_by_energy[COLS]) {
  fread(&new_event->tof, sizeof(double), 1, source);
  fread(&new_event->energy, sizeof(double), 1, source);
  read_vec(&new_event->position, source);
  read_vec(&new_event->direction, source);
  fread(&new_event->primary, sizeof(uint), 1, source);
  new_event->detector_id = get_detector(new_event->position);
  new_event->detected = is_detected(new_event, eff_by_energy);
}
bool read_annihilation(annihilation *annihil, FILE *source,
                       double eff_by_energy[COLS]) {
  if (!fread(&annihil->time, sizeof(double), 1, source))
    return false;
  read_vec(&annihil->origin, source);
  read_vec(&annihil->center, source);
  fread(&annihil->num_events, sizeof(uint), 1, source);
  uint num_primary1 = 0;
  uint num_primary2 = 0;
  annihil->events = (event *)malloc(sizeof(event) * annihil->num_events);
  uint num_hits = 0;
  uint num_primary1_hits = 0;
  uint num_primary2_hits = 0;
  for (int i = 0; i < annihil->num_events; i++) {
    event *next_event = &annihil->events[i];
    read_event(next_event, source, eff_by_energy);
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
  annihil->hits = (hit *)malloc(sizeof(hit) * num_hits);
  annihil->photon1.events = (event **)malloc(sizeof(event *) * num_primary1);
  annihil->photon2.events = (event **)malloc(sizeof(event *) * num_primary2);
  annihil->photon1.hits = (hit **)malloc(sizeof(hit *) * num_primary1_hits);
  annihil->photon2.hits = (hit **)malloc(sizeof(hit *) * num_primary2_hits);
  annihil->photon1.num_events = 0;
  annihil->photon2.num_events = 0;
  annihil->num_hits = 0;
  annihil->photon1.num_hits = 0;
  annihil->photon2.num_hits = 0;
  for (int i = 0; i < annihil->num_events; i++) {
    event *current_event = &annihil->events[i];
    if (current_event->primary == 1) {
      annihil->photon1.events[annihil->photon1.num_events] = current_event;
      annihil->photon1.num_events++;
    } else if (current_event->primary == 2) {
      annihil->photon2.events[annihil->photon2.num_events] = current_event;
      annihil->photon2.num_events++;
    }
    if (current_event->detected) {
      annihil->hits[annihil->num_hits] = event_to_hit(current_event);
      hit *current_hit = &annihil->hits[annihil->num_hits];
      annihil->num_hits++;
      if (current_event->primary == 1) {
        annihil->photon1.hits[annihil->photon1.num_hits] = current_hit;
        annihil->photon1.num_hits++;
      } else if (current_event->primary == 2) {
        annihil->photon2.hits[annihil->photon2.num_hits] = current_hit;
        annihil->photon2.num_hits++;
      }
    }
  }
  return true;
}
void print_lor(lor *new_lor, FILE *output) {
  print_vec(new_lor->center, output);
  print_sym_matrix(&new_lor->covariance, output);
}
void read_eff(FILE *source, double eff_by_energy[COLS]) {
  // loops through all the entries in a row
  for (int i = 0; i < COLS; i++) {
    fscanf(source, "%lf,", &eff_by_energy[i]);
  }
}
void print_double(double numb, FILE *output) {
  fwrite(&numb, sizeof(double), 1, output);
}
void print_int(int numb, FILE *output) {
  fwrite(&numb, sizeof(int), 1, output);
}
void print_vec(vec3d vec, FILE *output) {
  fwrite(&vec.x, sizeof(double), 1, output);
  fwrite(&vec.y, sizeof(double), 1, output);
  fwrite(&vec.z, sizeof(double), 1, output);
}
void print_sym_matrix(sym_matrix *mat, FILE *output) {
  fwrite(&mat->xx, sizeof(double), 1, output);
  fwrite(&mat->xy, sizeof(double), 1, output);
  fwrite(&mat->xz, sizeof(double), 1, output);
  fwrite(&mat->yy, sizeof(double), 1, output);
  fwrite(&mat->yz, sizeof(double), 1, output);
  fwrite(&mat->zz, sizeof(double), 1, output);
}
bool read_vec(vec3d *vec, FILE *source) {
  bool worked = fread(&vec->x, sizeof(double), 1, source);
  fread(&vec->y, sizeof(double), 1, source);
  fread(&vec->z, sizeof(double), 1, source);
  return worked;
}
bool read_sym(sym_matrix *mat, FILE *source) {
  bool worked = fread(&mat->xx, sizeof(double), 1, source);
  fread(&mat->xy, sizeof(double), 1, source);
  fread(&mat->xz, sizeof(double), 1, source);
  fread(&mat->yy, sizeof(double), 1, source);
  fread(&mat->yz, sizeof(double), 1, source);
  fread(&mat->zz, sizeof(double), 1, source);
  return worked;
}
bool read_lor(lor *new_lor, FILE *input) {
  read_vec(&new_lor->center, input);
  return read_sym(&new_lor->covariance, input);
}

void print_path(photon_path *path, FILE *output) {
  fwrite(&path->num_events, sizeof(uint), 1, output);
  for (int i = 0; i < path->num_events; i++) {
    print_vec(path->events[i]->position, output);
    fwrite(&path->events[i]->energy, sizeof(double), 1, output);
    fwrite(&path->events[i]->detected, sizeof(bool), 1, output);
  }
}
void print_annihilation(annihilation *annihil, FILE *output) {
  print_vec(annihil->center, output);
  print_path(&annihil->photon1, output);
  print_path(&annihil->photon2, output);
}
