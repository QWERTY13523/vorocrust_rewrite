#pragma once
#include "Geometry.h"
#include "MeshingTree.h"
#include "Methods.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <iterator>
#include <climits>
#include <vector>
#include <iomanip>
#include <cstring>
#include <tuple>

class Generator
{
public:
    Generator();
    ~Generator();
    int read_input_obj_file(std::string filename, size_t &num_points, double** &points, size_t &num_faces, size_t** &faces);
    void generate_spheres(const char* filename, MeshingTree *&spheres);
    void read_obj_faces(const char* filename, std::vector<int>& faces_flat, size_t& num_faces);
    //void generate_connections(const char* filename, std::vector<std::vector<std::tuple<size_t,size_t,size_t>>> &connections);
    void generate_surface_seeds(size_t num_points, double **points, size_t num_faces, size_t **faces,
         MeshingTree *surface_spheres,MeshingTree *upper_seeds, MeshingTree *lower_seeds);
    void color_surface_seeds(
        int num_faces, 
        MeshingTree *surface_spheres, 
        MeshingTree *upper_seeds, 
        MeshingTree *lower_seeds, 
        MeshingTree *seeds,
        std::vector<int> face, 
        double** points,
        double* &seedes, 
        size_t* &seeds_region_id, 
        double* &seeds_sizing
    );
    void generate_surface_mesh(MeshingTree* seeds, const char* output_filename);
    void generate_surface_mesh1(MeshingTree* seeds, const char* output_filename);
    void generate_seed_csv(
        const char* filename, 
        int num_dim, 
        size_t num_seeds, 
        double* spheres, 
        double* spheres_sizing, 
        size_t* spheres_region_id
    );
private:
    Methods _methods;
    Geometry _geom;
};