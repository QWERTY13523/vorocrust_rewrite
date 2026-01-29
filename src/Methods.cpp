#include "Methods.h"

int Methods::get_overlapping_spheres(double* sphere, MeshingTree* spheres,
	                                        size_t &num_overlapping_spheres, size_t* &overlapping_spheres)
{
	#pragma region Get Overlapping Spheres:
	if (spheres->get_num_tree_points() == 0)
	{
		num_overlapping_spheres = 0;
		overlapping_spheres = 0;
		return 0;
	}

	size_t iclosest; double hclosest(DBL_MAX);
	spheres->get_closest_tree_point(sphere, iclosest, hclosest);
	double rclosest = spheres->get_tree_point_attrib(iclosest, 0);

	double r = rclosest + 3 * hclosest;
	if (r < sphere[3]) r = sphere[3];

	double R = r * 5;
	R += r;

	size_t num(0), cap(100);
	size_t* neighbor_spheres = new size_t[cap];
	spheres->get_tree_points_in_sphere(sphere, R, num, neighbor_spheres, cap);

	num_overlapping_spheres = 0;
	for (size_t i = 0; i < num; i++)
	{
		double* neighbor_sphere = new double[4];
		spheres->get_tree_point(neighbor_spheres[i], 4, neighbor_sphere);
		if (_geom.overlapping_spheres(3, sphere, neighbor_sphere)) num_overlapping_spheres++;
		delete[] neighbor_sphere;
	}

	overlapping_spheres = new size_t[num_overlapping_spheres];
	num_overlapping_spheres = 0;
	for (size_t i = 0; i < num; i++)
	{
		double* neighbor_sphere = new double[4];
		spheres->get_tree_point(neighbor_spheres[i], 4, neighbor_sphere);
		if (_geom.overlapping_spheres(3, sphere, neighbor_sphere))
		{
			overlapping_spheres[num_overlapping_spheres] = neighbor_spheres[i];
			num_overlapping_spheres++;
		}
		delete[] neighbor_sphere;
	}
	delete[] neighbor_spheres;
	return 0;
	#pragma endregion
}