#include "Inout.h"
#include "MeshingTree.h"

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#ifdef USE_CGAL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/helpers.h>
#endif

namespace {

struct SeedPoint {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    double r = 0.0;
    size_t region = 0;
};

#ifdef USE_CGAL
using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point3 = Kernel::Point_3;
using SurfaceMesh = CGAL::Surface_mesh<Point3>;
using SideTester = CGAL::Side_of_triangle_mesh<SurfaceMesh, Kernel>;

static bool load_reference_mesh_with_fallback(
    const char* reference_mesh_obj,
    SurfaceMesh& mesh,
    std::string& used_path
)
{
    if (!reference_mesh_obj || reference_mesh_obj[0] == '\0') {
        return false;
    }

    std::vector<std::string> candidates;
    const auto add_candidate = [&](const std::filesystem::path& p) {
        const std::string s = p.string();
        if (s.empty()) return;
        if (std::find(candidates.begin(), candidates.end(), s) == candidates.end()) {
            candidates.push_back(s);
        }
    };

    const std::filesystem::path ref_path(reference_mesh_obj);
    add_candidate(ref_path);
    if (!ref_path.is_absolute()) {
        add_candidate(std::filesystem::path("build") / ref_path);
        add_candidate(std::filesystem::path("../build") / ref_path);
    }

    for (const auto& candidate : candidates) {
        SurfaceMesh candidate_mesh;
        if (!CGAL::Polygon_mesh_processing::IO::read_polygon_mesh(candidate, candidate_mesh)) {
            continue;
        }
        if (!CGAL::is_triangle_mesh(candidate_mesh)) {
            continue;
        }
        if (candidate_mesh.number_of_faces() == 0) {
            continue;
        }

        mesh = std::move(candidate_mesh);
        used_path = candidate;
        return true;
    }

    return false;
}
#endif

static inline void trim_in_place(std::string& s)
{
    size_t start = 0;
    while (start < s.size() && std::isspace(static_cast<unsigned char>(s[start]))) {
        start++;
    }
    size_t end = s.size();
    while (end > start && std::isspace(static_cast<unsigned char>(s[end - 1]))) {
        end--;
    }
    if (start == 0 && end == s.size()) return;
    s = s.substr(start, end - start);
}

static bool parse_csv_line(const std::string& line, SeedPoint& out)
{
    std::string s = line;
    trim_in_place(s);
    if (s.empty()) return false;
    if (s[0] == '#') return false;

    std::vector<std::string> tokens;
    std::stringstream ss(s);
    std::string token;
    while (std::getline(ss, token, ',')) {
        trim_in_place(token);
        if (!token.empty()) tokens.push_back(token);
    }

    if (tokens.size() < 4) return false;

    if (!tokens[0].empty() && std::isalpha(static_cast<unsigned char>(tokens[0][0]))) {
        return false;
    }

    try {
        out.x = std::stod(tokens[0]);
        out.y = std::stod(tokens[1]);
        out.z = std::stod(tokens[2]);
        out.r = std::stod(tokens[3]);
        if (tokens.size() >= 5) {
            out.region = static_cast<size_t>(std::stoull(tokens[4]));
        } else {
            out.region = 1;
        }
    } catch (...) {
        return false;
    }

    return true;
}

} // namespace

bool Inout::analyze_seed_csv(
    const std::string& filename,
    const char* cloud_obj,
    const char* query_obj,
    const char* combined_obj,
    const char* nearest_obj,
    const char* reference_mesh_obj,
    std::size_t num_tests
)
{
    if (num_tests == 0) {
        std::cerr << "Error: num_tests must be greater than 0." << std::endl;
        return false;
    }

    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: cannot open seed csv: " << filename << std::endl;
        return false;
    }

    std::vector<SeedPoint> points;
    std::string line;
    while (std::getline(file, line)) {
        SeedPoint p;
        if (parse_csv_line(line, p)) {
            points.push_back(p);
        }
    }

    if (points.empty()) {
        std::cerr << "Error: no valid points parsed from " << filename << std::endl;
        return false;
    }

    if (!reference_mesh_obj || reference_mesh_obj[0] == '\0') {
        std::cerr << "Error: reference_mesh_obj is empty." << std::endl;
        return false;
    }

    double min_x = std::numeric_limits<double>::infinity();
    double min_y = std::numeric_limits<double>::infinity();
    double min_z = std::numeric_limits<double>::infinity();
    double max_x = -std::numeric_limits<double>::infinity();
    double max_y = -std::numeric_limits<double>::infinity();
    double max_z = -std::numeric_limits<double>::infinity();

    for (const auto& p : points) {
        min_x = std::min(min_x, p.x);
        min_y = std::min(min_y, p.y);
        min_z = std::min(min_z, p.z);
        max_x = std::max(max_x, p.x);
        max_y = std::max(max_y, p.y);
        max_z = std::max(max_z, p.z);
    }

    const double span_x = max_x - min_x;
    const double span_y = max_y - min_y;
    const double span_z = max_z - min_z;
    const double max_span = std::max({span_x, span_y, span_z, 1e-9});
    const double padding = 0.02 * max_span;
    min_x -= padding;
    min_y -= padding;
    min_z -= padding;
    max_x += padding;
    max_y += padding;
    max_z += padding;

    MeshingTree seed_tree;
    for (const auto& p : points) {
        double x[3] = {p.x, p.y, p.z};
        seed_tree.add_tree_point(3, x, nullptr, nullptr);
    }

    SurfaceMesh mesh;
    std::string loaded_mesh_path;
    if (!load_reference_mesh_with_fallback(reference_mesh_obj, mesh, loaded_mesh_path)) {
        std::cerr << "Error: cannot read a valid triangular reference mesh obj from: "
                  << reference_mesh_obj
                  << " (also tried build/ and ../build/ fallbacks)"
                  << std::endl;
        return false;
    }
    std::cout << "[Validation] reference mesh: " << loaded_mesh_path
              << " (faces=" << mesh.number_of_faces() << ")" << std::endl;

    SideTester side_tester(mesh);

    std::mt19937_64 gen(std::random_device{}());
    std::uniform_real_distribution<double> dist_x(min_x, max_x);
    std::uniform_real_distribution<double> dist_y(min_y, max_y);
    std::uniform_real_distribution<double> dist_z(min_z, max_z);

    struct ErrorPair {
        double qx;
        double qy;
        double qz;
        double nx;
        double ny;
        double nz;
    };

    std::vector<ErrorPair> errors;
    std::size_t mismatch_count = 0;
    std::size_t boundary_count = 0;
    std::size_t valid_count = 0;
    std::chrono::nanoseconds cgal_time_ns(0);
    std::chrono::nanoseconds method_time_ns(0);

    for (std::size_t i = 0; i < num_tests; i++) {
        const double qx = dist_x(gen);
        const double qy = dist_y(gen);
        const double qz = dist_z(gen);

        const auto cgal_start = std::chrono::steady_clock::now();
        const CGAL::Bounded_side side = side_tester(Point3(qx, qy, qz));
        const auto cgal_end = std::chrono::steady_clock::now();
        cgal_time_ns += std::chrono::duration_cast<std::chrono::nanoseconds>(cgal_end - cgal_start);

        if (side == CGAL::ON_BOUNDARY) {
            boundary_count++;
            continue;
        }

        valid_count++;

        const auto method_start = std::chrono::steady_clock::now();
        double query[3] = {qx, qy, qz};
        size_t nearest_idx = points.size();
        double nearest_dist = std::numeric_limits<double>::max();
        if (seed_tree.get_closest_tree_point(query, nearest_idx, nearest_dist) != 0 || nearest_idx >= points.size()) {
            std::cerr << "Error: nearest-neighbor query failed at test " << i << std::endl;
            return false;
        }

        const bool predicted_inside = points[nearest_idx].region != 0;
        const bool reference_inside = (side == CGAL::ON_BOUNDED_SIDE);
        const auto method_end = std::chrono::steady_clock::now();
        method_time_ns += std::chrono::duration_cast<std::chrono::nanoseconds>(method_end - method_start);

        if (predicted_inside != reference_inside) {
            mismatch_count++;
            const SeedPoint& nearest = points[nearest_idx];
            errors.push_back({qx, qy, qz, nearest.x, nearest.y, nearest.z});
        }

        if ((i + 1) % 1000 == 0) {
            std::cout << "[Validation] processed " << (i + 1)
                      << " / " << num_tests
                      << ", mismatches=" << mismatch_count
                      << ", boundary_skipped=" << boundary_count
                      << std::endl;
        }
    }

    const double accuracy = (valid_count > 0)
        ? (1.0 - static_cast<double>(mismatch_count) / static_cast<double>(valid_count))
        : 1.0;

    const double cgal_total_ms = std::chrono::duration<double, std::milli>(cgal_time_ns).count();
    const double method_total_ms = std::chrono::duration<double, std::milli>(method_time_ns).count();
    const double cgal_avg_us = (num_tests > 0)
        ? (cgal_total_ms * 1000.0 / static_cast<double>(num_tests))
        : 0.0;
    const double method_avg_us = (valid_count > 0)
        ? (method_total_ms * 1000.0 / static_cast<double>(valid_count))
        : 0.0;

    std::cout << std::fixed << std::setprecision(6)
              << "[Validation Summary] total_tests=" << num_tests
              << ", effective_tests=" << valid_count
              << ", boundary_skipped=" << boundary_count
              << ", mismatches=" << mismatch_count
              << ", accuracy=" << (accuracy * 100.0) << "%"
              << ", cgal_total_ms=" << cgal_total_ms
              << ", cgal_avg_us=" << cgal_avg_us
              << ", method_total_ms=" << method_total_ms
              << ", method_avg_us=" << method_avg_us
              << std::endl;

    {
        if (cloud_obj && cloud_obj[0] != '\0') {
            std::ofstream out(cloud_obj);
            if (out.is_open()) {
                out << std::fixed << std::setprecision(16);
                for (const auto& p : points) {
                    out << "v " << p.x << " " << p.y << " " << p.z << "\n";
                }
            } else {
                std::cerr << "Warning: cannot write " << cloud_obj << std::endl;
            }
        }
    }

    {
        if (query_obj && query_obj[0] != '\0') {
            std::ofstream out(query_obj);
            if (out.is_open()) {
                out << std::fixed << std::setprecision(16);
                for (const auto& e : errors) {
                    out << "v " << e.qx << " " << e.qy << " " << e.qz << "\n";
                }
            } else {
                std::cerr << "Warning: cannot write " << query_obj << std::endl;
            }
        }
    }

    {
        if (nearest_obj && nearest_obj[0] != '\0') {
            std::ofstream out(nearest_obj);
            if (out.is_open()) {
                out << std::fixed << std::setprecision(16);
                for (const auto& e : errors) {
                    out << "v " << e.nx << " " << e.ny << " " << e.nz << "\n";
                }
            } else {
                std::cerr << "Warning: cannot write " << nearest_obj << std::endl;
            }
        }
    }

    {
        if (combined_obj && combined_obj[0] != '\0') {
            std::ofstream out(combined_obj);
            if (out.is_open()) {
                const std::size_t num_errors = errors.size();
                out << std::fixed << std::setprecision(16);
                for (const auto& e : errors) {
                    out << "v " << e.qx << " " << e.qy << " " << e.qz << " 1 0 0\n";
                }
                for (const auto& e : errors) {
                    out << "v " << e.nx << " " << e.ny << " " << e.nz << " 0 1 0\n";
                }
                for (std::size_t i = 0; i < num_errors; i++) {
                    out << "l " << (i + 1) << " " << (num_errors + i + 1) << "\n";
                }
            } else {
                std::cerr << "Warning: cannot write " << combined_obj << std::endl;
            }
        }
    }

    if (errors.empty()) {
        std::cout << "[Validation] No mismatches found. Error OBJ files are empty." << std::endl;
    } else {
        std::cout << "[Validation] Wrote " << errors.size() << " mismatch pairs to OBJ files." << std::endl;
    }

    return true;
}
