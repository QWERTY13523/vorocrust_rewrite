#include <cmath>
#include "MeshingTree.h"
#include "Geometry.h"

class Methods
{
public:
    int get_overlapping_spheres(double* sphere, MeshingTree* spheres,
	                                    size_t &num_overlapping_spheres, size_t* &overlapping_spheres);
private:
    Geometry _geom;
};