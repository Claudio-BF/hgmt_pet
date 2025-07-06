#ifndef linear_algebra_h
#define linear_algebra_h
#include <stdbool.h>
#include <stdio.h>

typedef struct vec3d {
  double x, y, z;
} vec3d;
typedef struct sym_matrix {
  double xx, xy, xz, yy, yz, zz;
} sym_matrix;
typedef struct lower_matrix {
  double l11, l21, l31, l22, l32, l33;
} lower_matrix;
// defines a new 3-vector
vec3d three_vec(double x, double y, double z);

// returns the magnitude of a 3-vector
double vec_mag(vec3d vector);
// same but sqaured cause I don't like square roots
double vec_mag2(vec3d vector);

// returns the dot product of two 3-vectors
double vec_dot(vec3d a, vec3d b);

// returns the addition of two 3-vectors as a new 3-vector
vec3d vec_add(vec3d a, vec3d b);

// returns the subtraction of the second vector from the first
// if the first vector is NULL acts as negating
vec3d vec_sub(vec3d a, vec3d b);

// returns the cross product of the two vectors a X b as a new vector
vec3d vec_cross(vec3d a, vec3d b);

// calculates the angle between two 3 vectors in radians
double vec_angle(vec3d a, vec3d b);

// determines the distance between two vectors
double vec_dist(vec3d a, vec3d b);

// // makes a new copy of a vector
// vec3d vec_copy(vec3d a);

// normalizes the given vector, returns as a new vector structure
vec3d vec_norm(vec3d a);

// multiplies a vector by a scalar
vec3d vec_scale(vec3d a, double b);

// projects a onto b
vec3d vec_projection(vec3d a, vec3d b);

// rejects a from b (i.e. projects a onto the plane perpendicular to b)
vec3d vec_rejection(vec3d a, vec3d b);
double radial_dist(vec3d a);
vec3d radial_scale(vec3d a, double b);
void vec_print(vec3d vec);
sym_matrix sym_add(sym_matrix a, sym_matrix b);
double sym_max_eigenvalue(sym_matrix mat);
vec3d sym_eigenvector(sym_matrix *mat, double eigenvalue);
sym_matrix sym_zero();
double sym_det(sym_matrix *mat);
sym_matrix sym_scale(sym_matrix mat, double factor);
sym_matrix sym_proj(vec3d vec);         // gives the outer square of a vector
lower_matrix cholesky(sym_matrix *mat); // find L so that LL^T = S
sym_matrix sym_adjucate(sym_matrix *sym);
vec3d sym_transform(sym_matrix *sym, vec3d vec);
sym_matrix sym_id(double diagonal);
#endif
