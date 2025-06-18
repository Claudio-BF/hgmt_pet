#include "read_write.h"

void read_event(event *new_event, FILE *source, double eff_by_energy[COLS]) {
  fread(&new_event->tof, sizeof(double), 1, source);
  fread(&new_event->energy, sizeof(double), 1, source);
  read_vec(&new_event->position, source);
  read_vec(&new_event->direction, source);
  fread(&new_event->primary, sizeof(uint), 1, source);
  new_event->detector_id = get_detector(new_event->position);
  new_event->detected = is_detected(new_event, eff_by_energy);
}
int compare_hits(const void *hit1, const void *hit2) {
  return ((hit *)hit1)->tof > ((hit *)hit2)->tof;
}
bool read_annihilation(annihilation *new_annihilation, FILE *source,
                       double eff_by_energy[COLS]) {
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
      hit detector_hit = event_to_hit(current_event);
      new_annihilation->hits[num_hits] = detector_hit;
      num_hits++;
    }
  }
  // qsort(new_annihilation->hits, num_hits, sizeof(hit), compare_hits);
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
