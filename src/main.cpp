#include "Generator.h"
#include "MeshingTree.h"
#include <vector>
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
    generator.read_input_obj_file("./data/obj/block.obj",num_points,points,num_faces,faces);
    generator.generate_spheres("./data/spheres/Sphere_6000_55_block.csv", spheres);
    generator.read_obj_faces("./data/obj/Ours_6000_block_Remesh.obj",face_flat, num_faces1);
    generator.generate_surface_seeds(num_points, points, num_faces, faces,
         spheres, upper_seeds, lower_seeds);
    generator.color_surface_seeds(num_faces1, spheres,upper_seeds, lower_seeds, seeds,
        face_flat, points, seedes, seeds_region_id, seeds_sizing);
    //generator.generate_seed_csv("seeds.csv",3,seeds->get_num_tree_points(),seedes, seeds_sizing, seeds_region_id);
    generator.generate_surface_mesh(seeds, "surface_mesh1.obj");
}