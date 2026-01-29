#pragma once

#include<cmath>
#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include <climits> 
#include <cfloat>

class Geometry
{
public:
    Geometry();
    ~Geometry();
    double get_point_angle(double x, double y);

    double distance(size_t num_dim, double* x, double* y);

    double cos_angle(size_t num_dim, double* x, double* y, double* z);

	bool get_power_vertex(size_t num_dim, size_t num_points, double** centers, double* radii, double* pv);

	int get_normal_component(size_t num_dim, size_t num_basis, double** basis, double* vect, double &norm);

	double dot_product(size_t num_dim, double* v1, double* v2);

	bool normalize_vector(size_t num_dim, double* vec);

	bool overlapping_spheres(size_t num_dim, double* sphere_i, double* sphere_j);

	double get_3d_triangle_area(double** corners);

	bool get_3d_triangle_normal(double** corners, double* normal);

	int project_to_3d_triangle(double* p, double** corners, double* q, double &proj_dist);

	int project_to_3d_line_segment(double* p, double* xo, double* xn, double* q, double& proj_dist);

	int project_to_3d_line(double* p, double* xo, double* edir, double* q, double &proj_dist);
};