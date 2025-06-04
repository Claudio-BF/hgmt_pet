#include "linear_algebra.h"
#include "helper_functions.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define MAX(a, b) ((a) > (b) ? (a) : (b))
// defines a new 3-vector
vec3d three_vec(double x, double y, double z) {
  vec3d new;
  new.x = x;
  new.y = y;
  new.z = z;
  return new;
}

// determines the distance between two vectors
double vec_dist(vec3d a, vec3d b) {
  // subtract b from a
  vec3d new = vec_sub(a, b);
  // get the distance of the subtracted vector
  double distance = vec_mag(new);
  return distance;
}

// returns the magnitude of a 3-vector
double vec_mag(vec3d vector) {
  double x2 = vector.x * vector.x;
  double y2 = vector.y * vector.y;
  double z2 = vector.z * vector.z;
  return sqrt(x2 + y2 + z2);
}

double vec_mag2(vec3d vector) {
  double x2 = vector.x * vector.x;
  double y2 = vector.y * vector.y;
  double z2 = vector.z * vector.z;
  return x2 + y2 + z2;
}
// returns the dot product of two 3-vectors
double vec_dot(vec3d a, vec3d b) {
  double first = a.x * b.x;
  double second = a.y * b.y;
  double third = a.z * b.z;
  return first + second + third;
}

// returns the addition of two 3-vectors as a new 3-vector
vec3d vec_add(vec3d a, vec3d b) {
  vec3d sum;
  sum.x = a.x + b.x;
  sum.y = a.y + b.y;
  sum.z = a.z + b.z;
  return sum;
}

// returns the subtraction of the second vector from the first
vec3d vec_sub(vec3d a, vec3d b) {
  vec3d sum;
  sum.x = a.x - b.x;
  sum.y = a.y - b.y;
  sum.z = a.z - b.z;
  return sum;
}

// returns the cross product of the two vectors a X b
vec3d vec_cross(vec3d a, vec3d b) {
  vec3d cross;
  cross.x = (a.y * b.z) - (a.z * b.y);
  cross.y = (a.z * b.x) - (a.x * b.z);
  cross.z = (a.x * b.y) - (a.y * b.x);
  return cross;
}

// calculates the angle between two 3 vectors in radians
double vec_angle(vec3d a, vec3d b) {
  // arccos((a * b) / (||a|| ||b||)) = theta
  double a_dot_b = vec_dot(a, b);
  double norms = vec_mag(a) * (vec_mag(b));
  double operand = a_dot_b / norms;

  // checks to make sure we're not going out-of-bounds in the domain of arccos
  // sometimes happens when working with floats
  if (operand > 1.0) {
    return acos(1);
  } else if (operand < -1.0) {
    return acos(-1);
  } else {
    return acos(a_dot_b / norms);
  }
}

void print_vec(vec3d vec, FILE *output) {
  fwrite(&vec.x, sizeof(double), 1, output);
  fwrite(&vec.y, sizeof(double), 1, output);
  fwrite(&vec.z, sizeof(double), 1, output);
}

// normalizes the given vector, returns as a new vector structure
vec3d vec_norm(vec3d a) {
  double mag = vec_mag(a);
  return three_vec(a.x / mag, a.y / mag, a.z / mag);
}

// multiplies a vector by a scalar
vec3d vec_scale(vec3d a, double b) {
  return three_vec(b * a.x, b * a.y, b * a.z);
}

// projects a onto b
vec3d vec_projection(vec3d a, vec3d b) {
  return vec_scale(b, (vec_dot(a, b) / vec_dot(b, b)));
}

// rejects b from a (i.e. projects a onto the plane perpendicular to b)
vec3d vec_rejection(vec3d a, vec3d b) {
  // find the projection of a onto b
  vec3d projection = vec_projection(a, b);
  return vec_sub(a, projection);
}
// returns the radial distance from axis of the detecter
double radial_dist(vec3d a) { return vec_mag(three_vec(a.x, a.y, 0.0)); }
// scales only the radial component
vec3d radial_scale(vec3d a, double b) {
  return three_vec(a.x * b, a.y * b, a.z);
}
void vec_print(vec3d vec) { printf("%f %f %f\n", vec.x, vec.y, vec.z); }
sym_matrix sym_zero() {
  sym_matrix mat;
  mat.xx = 0;
  mat.xy = 0;
  mat.xz = 0;
  mat.yy = 0;
  mat.yz = 0;
  mat.zz = 0;
  return mat;
}
sym_matrix sym_proj(vec3d vec) {
  sym_matrix mat;
  mat.xx = vec.x * vec.x;
  mat.xy = vec.x * vec.y;
  mat.xz = vec.x * vec.z;
  mat.yy = vec.y * vec.y;
  mat.yz = vec.y * vec.z;
  mat.zz = vec.z * vec.z;
  return mat;
}
sym_matrix sym_add(sym_matrix a, sym_matrix b) {
  sym_matrix c;
  c.xx = a.xx + b.xx;
  c.xy = a.xy + b.xy;
  c.xz = a.xz + b.xz;
  c.yy = a.yy + b.yy;
  c.yz = a.yz + b.yz;
  c.zz = a.zz + b.zz;
  return c;
}
sym_matrix sym_scale(sym_matrix mat, double factor) {
  mat.xx *= factor;
  mat.xy *= factor;
  mat.xz *= factor;
  mat.yy *= factor;
  mat.yz *= factor;
  mat.zz *= factor;
  return mat;
}
double sym_det(sym_matrix *mat) {
  return mat->xx * mat->yy * mat->zz + 2 * mat->xy * mat->xz * mat->yz -
         mat->xx * mat->yz * mat->yz - mat->yy * mat->xz * mat->xz -
         mat->zz * mat->xy * mat->xy;
}
sym_matrix square_coefficients(sym_matrix *mat) {
  sym_matrix new_mat;
  new_mat.xx = mat->xx * mat->xx;
  new_mat.xy = mat->xy * mat->xy;
  new_mat.xz = mat->xz * mat->xz;
  new_mat.yy = mat->yy * mat->yy;
  new_mat.yz = mat->yz * mat->yz;
  new_mat.zz = mat->zz * mat->zz;
  return new_mat;
}
double sym_max_eigenvalue(sym_matrix mat) {
  // geometrictools.com/Documentation/RobustEigenSymmetric3x3.pdf
  double q = (mat.xx + mat.yy + mat.zz) / 3;
  mat.xx -= q;
  mat.yy -= q;
  mat.zz -= q;
  sym_matrix math = square_coefficients(&mat);
  double p = sqrt(
      (math.xx + math.yy + math.zz + 2 * (math.xy + math.xz + math.yz)) / 6);
  double det = sym_det(&mat) / (p * p * p);
  double theta = fabs(det) >= 2 ? (det > 0 ? 0 : M_PI) : acos(det / 2) / 3;
  double root1 = 2 * cos(theta + 2 * M_PI / 3);
  double root2 = 2 * cos(theta + 4 * M_PI / 3);
  double root3 = 2 * cos(theta + 6 * M_PI / 3);
  return p * MAX(root1, MAX(root2, root3)) + q;
}
vec3d sym_row(sym_matrix *mat, uint row) {
  switch (row) {
  case 0:
    return three_vec(mat->xx, mat->xy, mat->xz);
  case 1:
    return three_vec(mat->xy, mat->yy, mat->yz);
  case 2:
    return three_vec(mat->xz, mat->yz, mat->zz);
  }
}
vec3d sym_eigenvector(sym_matrix *mat, double eigenvalue) {
  vec3d row1 = sym_row(mat, 0);
  vec3d row2 = sym_row(mat, 1);
  vec3d row3 = sym_row(mat, 2);
  row1.x -= eigenvalue;
  row2.y -= eigenvalue;
  row3.z -= eigenvalue;
  double dot1 = fabs(vec_dot(row1, row2));
  double dot2 = fabs(vec_dot(row1, row3));
  double dot3 = fabs(vec_dot(row2, row3));
  // taking the cross product so that the dot product has the minumum value
  // (avoids some nullity problems)
  if (dot1 < dot2) {
    if (dot1 < dot3)
      return vec_cross(row1, row2);
    else
      return vec_cross(row2, row3);
  } else {
    if (dot2 < dot3)
      return vec_cross(row1, row3);
    else
      return vec_cross(row2, row3);
  }
}
lower_matrix chopesky(sym_matrix *mat) {
  lower_matrix lower;
  lower.l11 = sqrt(mat->xx);
  lower.l21 = mat->xy / lower.l11;
  lower.l31 = mat->xz / lower.l11;
  lower.l22 = sqrt(mat->yy - lower.l21 * lower.l21);
  lower.l32 = (mat->yz - lower.l21 * lower.l31) / lower.l22;
  return lower;
}
vec3d lower_transform(lower_matrix *lower, vec3d vec) {
  return three_vec(vec.x * lower->l11, vec.x * lower->l21 + vec.y * lower->l22,
                   vec.x * lower->l31 + vec.y * lower->l32 +
                       vec.z * lower->l33);
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
