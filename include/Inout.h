#pragma once

#include <cstddef>
#include <string>

class Inout
{
public:
    static bool benchmark_seed_nearest_csv(
        const std::string& filename,
        std::size_t num_tests = 1000000,
        const std::string& csv_output = ""
    );

    static bool analyze_seed_csv(
        const std::string& filename,
        const char* cloud_obj = "seed_cloud.obj",
        const char* query_obj = "error_query_points.obj",
        const char* combined_obj = "error_query_with_nearest.obj",
        const char* nearest_obj = "error_nearest_points.obj",
        const char* reference_mesh_obj = "../data/obj/man.obj",
        std::size_t num_tests = 10000
    );

    static bool benchmark_mesh_inside_outside(
        const std::string& reconstructed_mesh_obj,
        const std::string& reference_mesh_obj,
        const std::string& seed_pair_csv,
        std::size_t num_tests = 100000,
        const std::string& csv_output = ""
    );
};
