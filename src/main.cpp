#include "Generator.h"
#include "MeshingTree.h"
#include <vector>

extern void optimizer(MeshingTree* seeds, MeshingTree* spheres, std::vector<int> face_flat);

double **points;
size_t **faces;
double *seeds_sizing, *seedes;
size_t *seeds_region_id;
std::vector<int> face_flat;
int main()
{
    Generator generator;
    MeshingTree *spheres = new MeshingTree();
    MeshingTree *upper_seeds = new MeshingTree();
    MeshingTree *lower_seeds = new MeshingTree();
    MeshingTree *seeds = new MeshingTree();
    size_t num_points, num_faces, num_faces1;
    generator.read_input_obj_file("../data/obj/mobius1.obj",num_points,points,num_faces,faces);
    if (num_points == 0) {
        std::cerr << "[Error] Failed to read input points or empty file." << std::endl;
        return 1;
    }

    generator.generate_spheres("../data/spheres/Sphere_6000_55.csv", spheres);
    if (spheres->get_num_tree_points() == 0) {
        std::cerr << "[Error] Failed to load spheres or empty file." << std::endl;
        return 1;
    }

    generator.read_obj_faces("../data/obj/Ours_6000_mobius1_Remesh.obj",face_flat, num_faces1);
    if (face_flat.empty()) {
        std::cerr << "[Error] Failed to read face indices or empty file." << std::endl;
        return 1;
    }

    // generator.generate_surface_seeds(num_points, points, num_faces, faces,
    //      spheres, upper_seeds, lower_seeds);
    generator.generate_surface_seeds(num_points, points, num_faces, faces,
        spheres, upper_seeds, lower_seeds, num_faces1, face_flat);
    generator.color_surface_seeds(num_faces1, spheres,upper_seeds, lower_seeds, seeds,
        face_flat, points, seedes, seeds_region_id, seeds_sizing);
    optimizer(seeds, spheres, face_flat);
    generator.generate_surface_mesh(seeds, spheres, "surface_mesh1.obj");
}