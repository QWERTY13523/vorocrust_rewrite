#pragma once

#include <cstddef>
#include <string>

class Inout
{
public:
    static bool analyze_seed_csv(
        const std::string& filename,
        const char* cloud_obj = "seed_cloud.obj",
        const char* query_obj = "error_query_points.obj",
        const char* combined_obj = "error_query_with_nearest.obj",
        const char* nearest_obj = "error_nearest_points.obj",
        const char* reference_mesh_obj = "../data/obj/block.obj",
        std::size_t num_tests = 10000
    );
};
