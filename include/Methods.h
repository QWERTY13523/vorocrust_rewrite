#include <cmath>
#include "MeshingTree.h"
#include "Geometry.h"

class Methods
{
public:
    int get_overlapping_spheres(double* sphere, MeshingTree* spheres,
	                                    size_t &num_overlapping_spheres, size_t* &overlapping_spheres);
    bool point_covered(double* point, MeshingTree* spheres, double alpha_coverage, size_t si, size_t sj, size_t sk,
	                         size_t &num_covering_spheres, size_t& cap_covering_spheres, size_t* &covering_spheres);
private:
    Geometry _geom;
};