#include "Generator.h"
#include <algorithm>
#include <cmath>
#include <map>
#include <fstream>
#include <random>
#include <iomanip>
#include <limits>
#include <set>
#include <string>
#include <array>
#include <tuple>
#include <vector>
#include <sstream>
#include <unordered_map>
#include <unordered_set>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>

#include <CGAL/Delaunay_triangulation_cell_base_3.h>
#define M_PI 3.14159265358979323846
typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_with_info_3<size_t, K> Vb;
typedef CGAL::Delaunay_triangulation_cell_base_3<K> Cb;
typedef CGAL::Triangulation_data_structure_3<Vb, Cb> Tds;
typedef CGAL::Delaunay_triangulation_3<K, Tds> Delaunay;
typedef Delaunay::Point Point_3;
typedef Delaunay::Vertex_handle Vertex_handle;
typedef Delaunay::Cell_handle Cell_handle;
typedef Delaunay::Edge Edge;
typedef K::Point_3 Point;

// 用于 std::tuple<size_t, size_t, size_t> 的哈希
struct TupleHash {
    inline size_t operator()(const std::tuple<size_t, size_t, size_t>& t) const {
        std::hash<size_t> hasher;
        size_t seed = 0;
        seed ^= hasher(std::get<0>(t)) + 0x9e3779b9 + (seed<<6) + (seed>>2);
        seed ^= hasher(std::get<1>(t)) + 0x9e3779b9 + (seed<<6) + (seed>>2);
        seed ^= hasher(std::get<2>(t)) + 0x9e3779b9 + (seed<<6) + (seed>>2);
        return seed;
    }
};

void optimizer(MeshingTree* seeds, MeshingTree* spheres, std::vector<int> face_flat)
{
    struct InsideFourthEvent {
        int iter;
        size_t seed_index;
        size_t s0;
        size_t s1;
        size_t s2;
        size_t fourth_sid;
        double px;
        double py;
        double pz;
    };

    Geometry geom;
    size_t num_spheres = spheres->get_num_tree_points();
    size_t num_faces_total = face_flat.size() / 3;

    std::cout << "=== Optimizer: Seed-to-Face Distance Driven Shrink (One Sphere Per Face) ===" << std::endl;

    auto make_face_key = [](size_t a, size_t b, size_t c) -> std::tuple<size_t, size_t, size_t> {
        size_t arr[3] = {a, b, c};
        if (arr[0] > arr[1]) std::swap(arr[0], arr[1]);
        if (arr[1] > arr[2]) std::swap(arr[1], arr[2]);
        if (arr[0] > arr[1]) std::swap(arr[0], arr[1]);
        return std::make_tuple(arr[0], arr[1], arr[2]);
    };

    // 辅助函数：计算种子点到某个面的距离
    auto dist_to_face = [&](double* seed_pos, size_t s0, size_t s1, size_t s2) -> double {
        double* corners[3];
        double sp0[4], sp1[4], sp2[4];
        spheres->get_tree_point(s0, 4, sp0);
        spheres->get_tree_point(s1, 4, sp1);
        spheres->get_tree_point(s2, 4, sp2);
        corners[0] = sp0; corners[1] = sp1; corners[2] = sp2;
        
        double q[3];
        double dist = 0.0;
        geom.project_to_3d_triangle(seed_pos, corners, q, dist);
        return dist;
    };

    // 1. 构建球到面的连接关系
    std::vector<std::vector<size_t>> sphere_to_faces(num_spheres);
    for (size_t fi = 0; fi < num_faces_total; fi++) {
        for (int k = 0; k < 3; k++) {
            int sid = face_flat[fi * 3 + k];
            if (sid >= 0 && (size_t)sid < num_spheres)
                sphere_to_faces[sid].push_back(fi);
        }
    }

    // Build 1-ring sphere adjacency from face connectivity.
    std::vector<std::vector<size_t>> sphere_face_neighbors(num_spheres);
    auto add_sphere_neighbor = [&](size_t a, size_t b) {
        if (a >= num_spheres || b >= num_spheres || a == b) return;
        sphere_face_neighbors[a].push_back(b);
        sphere_face_neighbors[b].push_back(a);
    };
    for (size_t fi = 0; fi < num_faces_total; ++fi) {
        int a = face_flat[fi * 3 + 0];
        int b = face_flat[fi * 3 + 1];
        int c = face_flat[fi * 3 + 2];
        if (a < 0 || b < 0 || c < 0) continue;
        add_sphere_neighbor(static_cast<size_t>(a), static_cast<size_t>(b));
        add_sphere_neighbor(static_cast<size_t>(b), static_cast<size_t>(c));
        add_sphere_neighbor(static_cast<size_t>(c), static_cast<size_t>(a));
    }
    for (size_t sid = 0; sid < num_spheres; ++sid) {
        auto& nbrs = sphere_face_neighbors[sid];
        std::sort(nbrs.begin(), nbrs.end());
        nbrs.erase(std::unique(nbrs.begin(), nbrs.end()), nbrs.end());
    }

    // ==========================================================
    // [极致优化]：预计算所有球的曲率，避免在循环中重复计算！
    // ==========================================================
    std::vector<double> sphere_curvatures(num_spheres, DBL_MAX);
    for (size_t sid = 0; sid < num_spheres; sid++) {
        const auto& connected_faces = sphere_to_faces[sid];
        if (connected_faces.empty()) continue;
        
        std::vector<size_t> n_sids;
        for (size_t fi : connected_faces) {
            for (int k = 0; k < 3; k++) {
                size_t other = (size_t)face_flat[fi * 3 + k];
                if (other != sid && other < num_spheres) {
                    n_sids.push_back(other);
                }
            }
        }
        std::sort(n_sids.begin(), n_sids.end());
        n_sids.erase(std::unique(n_sids.begin(), n_sids.end()), n_sids.end());

        double n0[3];
        double* n0_ptr = spheres->get_tree_point_normal(sid);
        n0[0] = n0_ptr[0]; n0[1] = n0_ptr[1]; n0[2] = n0_ptr[2];

        std::vector<double> n_normals;
        n_normals.reserve(n_sids.size() * 3);
        for (size_t other : n_sids) {
            double* ni_ptr = spheres->get_tree_point_normal(other);
            n_normals.push_back(ni_ptr[0]);
            n_normals.push_back(ni_ptr[1]);
            n_normals.push_back(ni_ptr[2]);
        }

        sphere_curvatures[sid] = geom.estimate_curvature_normal_variation(
            n0, n_normals.data(), n_sids.size(), 3);
    }

    auto compute_face_seeds = [&](size_t si, size_t sj, size_t sk,
                                   double* seed_a, double* seed_b,
                                   double* normal_out) -> bool {
        double sp_i[4], sp_j[4], sp_k[4];
        spheres->get_tree_point(si, 4, sp_i);
        spheres->get_tree_point(sj, 4, sp_j);
        spheres->get_tree_point(sk, 4, sp_k);

        if (geom.distance(3, sp_i, sp_j) > sp_i[3] + sp_j[3] + 1e-10) return false;
        if (geom.distance(3, sp_i, sp_k) > sp_i[3] + sp_k[3] + 1e-10) return false;
        if (geom.distance(3, sp_j, sp_k) > sp_j[3] + sp_k[3] + 1e-10) return false;

        double* centers[3] = {sp_i, sp_j, sp_k};
        double radii[3]    = {sp_i[3], sp_j[3], sp_k[3]};
        double c_ijk[3];
        if (!geom.get_power_vertex(3, 3, centers, radii, c_ijk)) return false;

        double hi = geom.distance(3, c_ijk, sp_i);
        if (hi > sp_i[3] - 1e-10) return false;

        double vi = sqrt(fmax(0.0, sp_i[3] * sp_i[3] - hi * hi));

        if (geom.get_3d_triangle_area(centers) < 1e-12) return false;
        geom.get_3d_triangle_normal(centers, normal_out);

        double min_r = fmin(sp_i[3], fmin(sp_j[3], sp_k[3]));
        for (int d = 0; d < 3; d++) {
            seed_a[d] = c_ijk[d] + vi * normal_out[d];
            seed_b[d] = c_ijk[d] - vi * normal_out[d];
        }
        seed_a[3] = min_r;
        seed_b[3] = min_r;
        return true;
    };

    auto compute_seed_dot = [&](double* seed_pos, size_t si, size_t sj, size_t sk) -> double {
        double dot_sum = 0.0;
        int count = 0;
        size_t sids[3] = {si, sj, sk};
        for (int i = 0; i < 3; i++) {
            if (sids[i] >= num_spheres) continue;
            double center[4];
            spheres->get_tree_point(sids[i], 4, center);
            double* normal = spheres->get_tree_point_normal(sids[i]);
            if (std::abs(normal[0]) < 1e-9 && std::abs(normal[1]) < 1e-9 && std::abs(normal[2]) < 1e-9) continue;
            
            dot_sum += (seed_pos[0] - center[0]) * normal[0]
                     + (seed_pos[1] - center[1]) * normal[1]
                     + (seed_pos[2] - center[2]) * normal[2];
            count++;
        }
        return (count > 0) ? (dot_sum / count) : 0.0;
    };

    const double kCandidateSphereCenterDistanceFactor = 4.0;

    auto collect_candidate_spheres = [&](size_t p0, size_t p1, size_t p2, double* seed_pos) {
        std::vector<size_t> candidate_spheres = {p0, p1, p2};

        double sp0[4], sp1[4], sp2[4];
        double local_search_radius = 0.0;
        if (p0 < num_spheres) {
            spheres->get_tree_point(p0, 4, sp0);
            local_search_radius = std::max(local_search_radius, sp0[3]);
        }
        if (p1 < num_spheres) {
            spheres->get_tree_point(p1, 4, sp1);
            local_search_radius = std::max(local_search_radius, sp1[3]);
        }
        if (p2 < num_spheres) {
            spheres->get_tree_point(p2, 4, sp2);
            local_search_radius = std::max(local_search_radius, sp2[3]);
        }
        local_search_radius = std::max(1e-6, kCandidateSphereCenterDistanceFactor * local_search_radius);

        size_t num_points_in_sphere = 0;
        size_t capacity = 0;
        size_t* points_in_sphere = nullptr;
        spheres->get_tree_points_in_sphere(seed_pos, local_search_radius,
                                           num_points_in_sphere, points_in_sphere, capacity);
        for (size_t i = 0; i < num_points_in_sphere; ++i) {
            candidate_spheres.push_back(points_in_sphere[i]);
        }
        delete[] points_in_sphere;

        std::sort(candidate_spheres.begin(), candidate_spheres.end());
        candidate_spheres.erase(std::unique(candidate_spheres.begin(), candidate_spheres.end()), candidate_spheres.end());
        return candidate_spheres;
    };

    const int max_iterations = 10;
    const double kCurvatureShrinkThresholdRadians = 0.5; 
    const int binary_search_iterations = 30;

    std::unordered_map<std::tuple<size_t,size_t,size_t>, std::vector<size_t>, TupleHash> face_key_to_seeds;
    std::vector<std::string> shrink_blocked_logs;
    shrink_blocked_logs.reserve(1024);
    std::vector<InsideFourthEvent> inside_fourth_events_iter0;
    inside_fourth_events_iter0.reserve(1024);
    std::unordered_set<std::string> inside_fourth_event_keys_iter0;
    inside_fourth_event_keys_iter0.reserve(2048);
    bool inside_fourth_iter0_dumped = false;

    auto append_uv_sphere_mesh = [&](std::ofstream& out, size_t& vertex_offset,
                                     double cx, double cy, double cz, double r) {
        if (r <= 0.0) return;
        const int slices = 20;
        const int stacks = 10;
        const int rings = stacks - 1;

        const size_t north = vertex_offset++;
        out << "v " << cx << " " << cy << " " << (cz + r) << "\n";

        const size_t ring_start = vertex_offset;
        for (int i = 1; i < stacks; ++i) {
            double phi = M_PI * static_cast<double>(i) / static_cast<double>(stacks);
            double sp = std::sin(phi);
            double cp = std::cos(phi);
            for (int j = 0; j < slices; ++j) {
                double theta = 2.0 * M_PI * static_cast<double>(j) / static_cast<double>(slices);
                double ct = std::cos(theta);
                double st = std::sin(theta);
                out << "v "
                    << (cx + r * sp * ct) << " "
                    << (cy + r * sp * st) << " "
                    << (cz + r * cp) << "\n";
                vertex_offset++;
            }
        }

        const size_t south = vertex_offset++;
        out << "v " << cx << " " << cy << " " << (cz - r) << "\n";

        if (rings <= 0) return;

        for (int j = 0; j < slices; ++j) {
            const size_t v1 = ring_start + static_cast<size_t>(j);
            const size_t v2 = ring_start + static_cast<size_t>((j + 1) % slices);
            out << "f " << north << " " << v1 << " " << v2 << "\n";
        }

        for (int i = 0; i < rings - 1; ++i) {
            const size_t cur = ring_start + static_cast<size_t>(i * slices);
            const size_t nxt = cur + static_cast<size_t>(slices);
            for (int j = 0; j < slices; ++j) {
                const size_t jn = static_cast<size_t>((j + 1) % slices);
                const size_t a = cur + static_cast<size_t>(j);
                const size_t b = nxt + static_cast<size_t>(j);
                const size_t c = nxt + jn;
                const size_t d = cur + jn;
                out << "f " << a << " " << b << " " << c << "\n";
                out << "f " << a << " " << c << " " << d << "\n";
            }
        }

        const size_t last_ring_start = ring_start + static_cast<size_t>((rings - 1) * slices);
        for (int j = 0; j < slices; ++j) {
            const size_t v1 = last_ring_start + static_cast<size_t>(j);
            const size_t v2 = last_ring_start + static_cast<size_t>((j + 1) % slices);
            out << "f " << v1 << " " << south << " " << v2 << "\n";
        }
    };

    auto dump_inside_fourth_iter0 = [&](const std::vector<InsideFourthEvent>& events) {
        const char* points_filename = "inside_fourth_points_iter0.obj";
        const char* spheres_filename = "inside_fourth_spheres_iter0.obj";

        // 1) Points file
        {
            std::ofstream out_points(points_filename);
            if (!out_points.is_open()) {
                std::cerr << "Failed to write " << points_filename << std::endl;
            } else {
                out_points << std::fixed << std::setprecision(16);
                out_points << "# iter0 inside-fourth points\n";
                out_points << "# event_count=" << events.size() << "\n";
                size_t vertex_offset = 1;
                for (size_t ei = 0; ei < events.size(); ++ei) {
                    const auto& ev = events[ei];
                    out_points << "# event " << (ei + 1)
                               << " iter=" << ev.iter
                               << " seed=" << ev.seed_index
                               << " face=(" << ev.s0 << "," << ev.s1 << "," << ev.s2 << ")"
                               << " inside_sphere=" << ev.fourth_sid << "\n";
                    out_points << "v " << ev.px << " " << ev.py << " " << ev.pz << "\n";
                    out_points << "p " << vertex_offset << "\n";
                    vertex_offset++;
                }
                out_points.close();
            }
        }

        // 2) Spheres file
        {
            std::set<size_t> involved_sphere_ids;
            for (const auto& ev : events) {
                involved_sphere_ids.insert(ev.s0);
                involved_sphere_ids.insert(ev.s1);
                involved_sphere_ids.insert(ev.s2);
                involved_sphere_ids.insert(ev.fourth_sid);
            }

            std::ofstream out_spheres(spheres_filename);
            if (!out_spheres.is_open()) {
                std::cerr << "Failed to write " << spheres_filename << std::endl;
            } else {
                out_spheres << std::fixed << std::setprecision(16);
                out_spheres << "# iter0 inside-fourth involved spheres\n";
                out_spheres << "# sphere_count=" << involved_sphere_ids.size() << "\n";
                size_t vertex_offset = 1;
                for (size_t sid : involved_sphere_ids) {
                    if (sid >= num_spheres) continue;
                    double sp[4];
                    spheres->get_tree_point(sid, 4, sp);
                    out_spheres << "g sphere_sid_" << sid << "\n";
                    append_uv_sphere_mesh(out_spheres, vertex_offset, sp[0], sp[1], sp[2], sp[3]);
                }
                out_spheres.close();
            }

            std::cout << "Iter0 inside-fourth debug exported: "
                      << points_filename << " (points=" << events.size() << "), "
                      << spheres_filename << " (spheres=" << involved_sphere_ids.size() << ")"
                      << std::endl;
        }
    };

    // Deadlock reporting disabled per request; keep logic below commented out.
    // std::vector<std::array<double, 3>> deadlock_positions;

    for (int iter = 0; iter < max_iterations; iter++) {

        std::unordered_set<size_t> processed_this_iter;
        size_t num_seeds = seeds->get_num_tree_points();
        size_t active_seed_count = 0;

        // Pre-Delaunay illegal intrusion detection:
        // A seed is illegal if it lies inside a non-own candidate sphere.
        // Intrusion ratio = (r - h) / r, where h = distance(seed, sphere center).
        {
            size_t illegal_seed_count = 0;
            size_t illegal_hit_count = 0;
            std::vector<double> intrusion_ratios;
            intrusion_ratios.reserve(num_seeds / 4 + 1);

            for (size_t i = 0; i < num_seeds; ++i) {
                if (!seeds->tree_point_is_active(i)) continue;
                active_seed_count++;

                size_t* attrib = seeds->get_tree_point_attrib(i);
                size_t p0 = attrib[2], p1 = attrib[3], p2 = attrib[4];
                if (p0 >= num_spheres || p1 >= num_spheres || p2 >= num_spheres) continue;

                double* seed_pos = seeds->get_tree_point(i);
                std::vector<size_t> candidate_spheres = collect_candidate_spheres(p0, p1, p2, seed_pos);

                bool this_seed_illegal = false;
                for (size_t c_sid : candidate_spheres) {
                    if (c_sid >= num_spheres) continue;
                    if (c_sid == p0 || c_sid == p1 || c_sid == p2) continue;

                    double c_sp[4];
                    spheres->get_tree_point(c_sid, 4, c_sp);
                    double h = geom.distance(3, seed_pos, c_sp);
                    double r = c_sp[3];
                    if (r <= 1e-14) continue;
                    if (h < r - 1e-10) {
                        this_seed_illegal = true;
                        illegal_hit_count++;
                        intrusion_ratios.push_back((r - h) / r);
                    }
                }

                if (this_seed_illegal) illegal_seed_count++;
            }

            if (!intrusion_ratios.empty()) {
                double sum = 0.0;
                for (double v : intrusion_ratios) sum += v;
                double mean = sum / static_cast<double>(intrusion_ratios.size());

                std::vector<double> tmp = intrusion_ratios;
                const size_t n = tmp.size();
                const size_t mid = n / 2;
                std::nth_element(tmp.begin(), tmp.begin() + mid, tmp.end());
                double median = tmp[mid];
                if (n % 2 == 0) {
                    double lower = *std::max_element(tmp.begin(), tmp.begin() + mid);
                    median = 0.5 * (lower + median);
                }

                const size_t p95_idx = static_cast<size_t>(0.95 * static_cast<double>(n - 1));
                std::nth_element(tmp.begin(), tmp.begin() + p95_idx, tmp.end());
                double p95 = tmp[p95_idx];
                double max_ratio = *std::max_element(intrusion_ratios.begin(), intrusion_ratios.end());

                std::cout << "  [Pre-Delaunay illegal check] iter=" << iter
                          << ", active_seeds=" << active_seed_count
                          << ", illegal_seeds=" << illegal_seed_count
                          << ", illegal_hits=" << illegal_hit_count
                          << ", ratio(mean/median/p95/max)="
                          << mean << "/" << median << "/" << p95 << "/" << max_ratio
                          << std::endl;
            } else {
                std::cout << "  [Pre-Delaunay illegal check] iter=" << iter
                          << ", active_seeds=" << active_seed_count
                          << ", illegal_seeds=0, illegal_hits=0"
                          << std::endl;
            }
        }

        face_key_to_seeds.clear();
        face_key_to_seeds.reserve(num_seeds / 2);
        for (size_t i = 0; i < num_seeds; i++) {
            if (!seeds->tree_point_is_active(i)) continue;
            size_t* attrib = seeds->get_tree_point_attrib(i);
            auto key = make_face_key(attrib[2], attrib[3], attrib[4]);
            face_key_to_seeds[key].push_back(i);
        }

        // ==========================================================
        // 枚举每个种子点，寻找距离越界的“侵略面”
        // ==========================================================
        std::set<std::tuple<size_t, size_t, size_t>> target_faces;
        std::unordered_set<std::tuple<size_t, size_t, size_t>, TupleHash> inside_fourth_target_faces;
        std::unordered_map<std::tuple<size_t, size_t, size_t>, std::vector<size_t>, TupleHash> target_face_shrink_pool;
        auto add_unique_sphere_to_pool = [&](const std::tuple<size_t, size_t, size_t>& key, size_t sid) {
            if (sid >= num_spheres) return;
            std::vector<size_t>& pool = target_face_shrink_pool[key];
            if (std::find(pool.begin(), pool.end(), sid) == pool.end()) {
                pool.push_back(sid);
            }
        };
        auto add_face_spheres_to_pool = [&](const std::tuple<size_t, size_t, size_t>& key, size_t s0, size_t s1, size_t s2) {
            add_unique_sphere_to_pool(key, s0);
            add_unique_sphere_to_pool(key, s1);
            add_unique_sphere_to_pool(key, s2);
        };

        for (size_t i = 0; i < num_seeds; i++) {
            if (!seeds->tree_point_is_active(i)) continue;
            size_t* attrib = seeds->get_tree_point_attrib(i);
            size_t p0 = attrib[2], p1 = attrib[3], p2 = attrib[4];
            if (p0 >= num_spheres || p1 >= num_spheres || p2 >= num_spheres) continue;

            double* seed_pos = seeds->get_tree_point(i);
            std::vector<size_t> candidate_spheres = collect_candidate_spheres(p0, p1, p2, seed_pos);

            auto key_own = make_face_key(p0, p1, p2);
            bool has_inside_fourth = false;
            for (size_t c_sid : candidate_spheres) {
                if (c_sid >= num_spheres) continue;
                if (c_sid == p0 || c_sid == p1 || c_sid == p2) continue;

                double c_sp[4];
                spheres->get_tree_point(c_sid, 4, c_sp);
                double h = geom.distance(3, seed_pos, c_sp);
                if (h >= c_sp[3] - 1e-10) continue;

                bool found_other_face = false;
                double best_other_dist = DBL_MAX;
                size_t best_f0 = 0, best_f1 = 0, best_f2 = 0;
                for (size_t fi : sphere_to_faces[c_sid]) {
                    size_t f0 = face_flat[fi*3];
                    size_t f1 = face_flat[fi*3+1];
                    size_t f2 = face_flat[fi*3+2];
                    
                    auto key_other = make_face_key(f0, f1, f2);
                    if (key_own == key_other) continue;
                    double dist_other = dist_to_face(seed_pos, f0, f1, f2);

                    if (!found_other_face || dist_other < best_other_dist) {
                        found_other_face = true;
                        best_other_dist = dist_other;
                        best_f0 = f0;
                        best_f1 = f1;
                        best_f2 = f2;
                    }
                }

                if (!found_other_face) continue;
                has_inside_fourth = true;
                add_face_spheres_to_pool(key_own, p0, p1, p2);
                add_face_spheres_to_pool(key_own, best_f0, best_f1, best_f2);

                // Export request: collect first-iteration detections only.
                if (iter == 0) {
                    const long long qx = llround(seed_pos[0] * 1e9);
                    const long long qy = llround(seed_pos[1] * 1e9);
                    const long long qz = llround(seed_pos[2] * 1e9);
                    std::string event_key = std::to_string(c_sid) + "|" +
                                            std::to_string(qx) + "|" +
                                            std::to_string(qy) + "|" +
                                            std::to_string(qz);
                    if (inside_fourth_event_keys_iter0.insert(event_key).second) {
                        inside_fourth_events_iter0.push_back(
                            InsideFourthEvent{
                                iter, i, p0, p1, p2, c_sid,
                                seed_pos[0], seed_pos[1], seed_pos[2]
                            }
                        );
                    }
                }
            }

            if (has_inside_fourth) {
                target_faces.insert(key_own);
                inside_fourth_target_faces.insert(key_own);
            }
        }

        if (iter == 0 && !inside_fourth_iter0_dumped) {
            dump_inside_fourth_iter0(inside_fourth_events_iter0);
            inside_fourth_iter0_dumped = true;
        }

        std::cout << "  Iter " << iter
                  << ": active seeds = " << active_seed_count
                  << ", faces targeted for optimization = " << target_faces.size() << std::endl;

        if (target_faces.empty()) break;

        size_t fixes_this_iter = 0;

        for (const auto& target_face : target_faces) {
            size_t face_sids[3] = {
                std::get<0>(target_face),
                std::get<1>(target_face),
                std::get<2>(target_face)
            };
            
            if (face_sids[0] >= num_spheres || face_sids[1] >= num_spheres || face_sids[2] >= num_spheres) continue;
            bool is_inside_fourth_target = (inside_fourth_target_faces.find(target_face) != inside_fourth_target_faces.end());

            bool has_sharp_vertex = false;
            for (int k = 0; k < 3; k++) {
                size_t sid = face_sids[k];
                double curvature = sphere_curvatures[sid];
                if (curvature != DBL_MAX && curvature >= kCurvatureShrinkThresholdRadians) {
                    has_sharp_vertex = true;
                    break;
                }
            }
            if (!has_sharp_vertex && !is_inside_fourth_target) continue;

            std::vector<size_t> shrink_pool = {face_sids[0], face_sids[1], face_sids[2]};
            auto pool_it = target_face_shrink_pool.find(target_face);
            if (pool_it != target_face_shrink_pool.end() && !pool_it->second.empty()) {
                shrink_pool = pool_it->second;
            }
            std::sort(shrink_pool.begin(), shrink_pool.end());
            shrink_pool.erase(std::unique(shrink_pool.begin(), shrink_pool.end()), shrink_pool.end());

            size_t best_sid = SIZE_MAX;
            double best_valid_r = 0.0;
            double best_orig_r = 0.0;
            double best_shrink_ratio = 0.0;
            double best_min_seed_pair_dist = -DBL_MAX;
            std::vector<size_t> best_seeds_to_delete;
            auto log_shrink_blocked = [&](size_t sid, const char* reason) {
                std::ostringstream oss;
                oss << "iter=" << iter
                    << ", target_face=(" << face_sids[0] << "," << face_sids[1] << "," << face_sids[2] << ")"
                    << ", sphere=" << sid
                    << ", reason=" << reason;
                shrink_blocked_logs.push_back(oss.str());
            };

            for (size_t sid : shrink_pool) {
                if (sid >= num_spheres) continue;
                // 【保护机制】：如果该球在本轮大迭代中已经被缩小过一次了，不能再动它
                if (processed_this_iter.count(sid)) {
                    log_shrink_blocked(sid, "sphere already shrunk once in this iteration");
                    continue;
                }

                double* sph_data = spheres->get_tree_point(sid);
                double  orig_r   = sph_data[3];

                const auto& connected_faces = sphere_to_faces[sid];
                if (connected_faces.empty()) {
                    log_shrink_blocked(sid, "sphere has no connected faces");
                    continue;
                }
                double r_min = 0.0;

                for (size_t fi : connected_faces) {
                    for (int k = 0; k < 3; k++) {
                        size_t other = (size_t)face_flat[fi * 3 + k];
                        if (other == sid) continue;
                        double other_sp[4];
                        spheres->get_tree_point(other, 4, other_sp);
                        double d   = geom.distance(3, sph_data, other_sp);
                        double req = d - other_sp[3] + 1e-6;
                        if (req > r_min) r_min = req;
                    }
                }
                if (r_min <= 0.0) r_min = 1e-6;
                // 拦截：如果再缩会导致断裂，放弃缩小
                if (r_min >= orig_r - 1e-10) {
                    log_shrink_blocked(sid, "minimum required radius reached; further shrink would break overlap");
                    continue;
                }

                std::vector<size_t> seeds_to_del;
                for (size_t fi : connected_faces) {
                    auto key = make_face_key(
                        (size_t)face_flat[fi*3],
                        (size_t)face_flat[fi*3+1],
                        (size_t)face_flat[fi*3+2]);
                    
                    auto it = face_key_to_seeds.find(key);
                    if (it != face_key_to_seeds.end()) {
                        for (size_t s : it->second) {
                            if (seeds->tree_point_is_active(s))
                                seeds_to_del.push_back(s);
                        }
                    }
                }
                if (!seeds_to_del.empty()) {
                    std::sort(seeds_to_del.begin(), seeds_to_del.end());
                    seeds_to_del.erase(std::unique(seeds_to_del.begin(), seeds_to_del.end()), seeds_to_del.end());
                }

                double lo   = 0;
                double hi_r = orig_r;
                double valid_r = orig_r;
                bool   found   = false;

                for (int bs = 0; bs < binary_search_iterations; bs++) {
                    double mid = (lo + hi_r) / 2.0;
                    sph_data[3] = mid;

                    bool all_ok = true;

                    for (size_t fi : connected_faces) {
                        double sa[4], sb[4], nm[3];
                        if (!compute_face_seeds(
                                (size_t)face_flat[fi*3],
                                (size_t)face_flat[fi*3+1],
                                (size_t)face_flat[fi*3+2],
                                sa, sb, nm)) {
                            all_ok = false;
                            break;
                        }
                    }

                    if (all_ok) {
                        valid_r = mid;
                        found   = true;
                        hi_r    = mid;
                    } else {
                        lo = mid;
                    }
                }

                sph_data[3] = orig_r;  // 先还原，供下一个球尝试

                if (!found) {
                    log_shrink_blocked(sid, "no valid shrink radius found by binary search");
                    continue;
                }

                double shrink_ratio = (orig_r - valid_r) / orig_r;

                // Evaluate the sphere by the minimum seed-pair distance across all faces touching this sphere.
                sph_data[3] = valid_r;
                double min_seed_pair_dist = DBL_MAX;
                bool metric_ok = false;
                for (size_t fi : connected_faces) {
                    double sa[4], sb[4], nm[3];
                    if (!compute_face_seeds(
                            (size_t)face_flat[fi*3],
                            (size_t)face_flat[fi*3+1],
                            (size_t)face_flat[fi*3+2],
                            sa, sb, nm)) {
                        metric_ok = false;
                        break;
                    }

                    double pair_dist = geom.distance(3, sa, sb);
                    if (pair_dist < min_seed_pair_dist) {
                        min_seed_pair_dist = pair_dist;
                    }
                    metric_ok = true;
                }
                sph_data[3] = orig_r;

                if (!metric_ok) {
                    log_shrink_blocked(sid, "seed-pair metric invalid after trial shrink");
                    continue;
                }

                // Previous logic (commented by request): choose sphere with largest shrink ratio.
                /*
                if (shrink_ratio > best_shrink_ratio) {
                    best_sid = sid;
                    best_valid_r = valid_r;
                    best_orig_r = orig_r;
                    best_shrink_ratio = shrink_ratio;
                    best_seeds_to_delete = std::move(seeds_to_del);
                }
                */

                // Default path: maximize minimum seed-pair distance.
                // Tie-break by larger shrink ratio.
                if (best_sid == SIZE_MAX ||
                    min_seed_pair_dist > best_min_seed_pair_dist ||
                    (std::abs(min_seed_pair_dist - best_min_seed_pair_dist) <= 1e-12 &&
                     shrink_ratio > best_shrink_ratio)) {
                    best_sid = sid;
                    best_valid_r = valid_r;
                    best_orig_r = orig_r;
                    best_shrink_ratio = shrink_ratio;
                    best_min_seed_pair_dist = min_seed_pair_dist;
                    best_seeds_to_delete = std::move(seeds_to_del);
                }
            }

            // ==========================================================
            // 如果找到了最佳球，仅对该球执行替换和拓扑刷新
            // ==========================================================
            if (best_sid != SIZE_MAX) {
                double* sph_data = spheres->get_tree_point(best_sid);
                sph_data[3] = best_valid_r;
                
                // 将被选中的最佳球加入黑名单，本轮不再碰它
                processed_this_iter.insert(best_sid);

                // std::cout << "    Targeted Face (" << face_sids[0] << "," << face_sids[1] << "," << face_sids[2] 
                //           << ") -> ONLY shrunk best sphere " << best_sid
                //           << ": radius " << best_orig_r << " -> " << best_valid_r
                //           << " (shrink " << (best_shrink_ratio * 100.0) << "%)"
                //           << ", min-seed-pair-dist " << best_min_seed_pair_dist << std::endl;

                for (size_t s : best_seeds_to_delete) {
                    seeds->lazy_delete_tree_point(s);
                }

                const auto& connected_faces = sphere_to_faces[best_sid];
                for (size_t fi : connected_faces) {
                    size_t s0 = (size_t)face_flat[fi*3];
                    size_t s1 = (size_t)face_flat[fi*3+1];
                    size_t s2 = (size_t)face_flat[fi*3+2];

                    double seed_a[4], seed_b[4], normal[3];
                    if (!compute_face_seeds(s0, s1, s2, seed_a, seed_b, normal))
                        continue;

                    double dot_a = compute_seed_dot(seed_a, s0, s1, s2);
                    double dot_b = compute_seed_dot(seed_b, s0, s1, s2);

                    size_t region_a = (dot_a >= dot_b) ? 0 : 1;
                    size_t region_b = (dot_a >= dot_b) ? 1 : 0;

                    size_t idx_a = seeds->get_num_tree_points();
                    size_t idx_b = idx_a + 1;

                    size_t att_a[6] = {6, idx_b, s0, s1, s2, region_a};
                    size_t att_b[6] = {6, idx_a, s0, s1, s2, region_b};

                    double normal_a[4] = { normal[0],  normal[1],  normal[2], 0.0};
                    double normal_b[4] = {-normal[0], -normal[1], -normal[2], 0.0};

                    seeds->add_tree_point(4, seed_a, normal_a, att_a);
                    seeds->add_tree_point(4, seed_b, normal_b, att_b);

                    auto key = make_face_key(s0, s1, s2);
                    face_key_to_seeds[key].push_back(idx_a);
                    face_key_to_seeds[key].push_back(idx_b);
                }

                fixes_this_iter++;
            }
        }

        std::cout << "  Applied " << fixes_this_iter << " sphere shrinks" << std::endl;
        
        // ==========================================================
        // 死锁反馈输出逻辑：如果卡死，则计算所有高曲率死锁面几何中心，并写入 xyz 文件
        // ==========================================================
        if (fixes_this_iter == 0) {
            if (!target_faces.empty()) {
                size_t high_curvature_deadlocks_count = 0;
                
                for (const auto& target_face : target_faces) {
                    size_t s0 = std::get<0>(target_face);
                    size_t s1 = std::get<1>(target_face);
                    size_t s2 = std::get<2>(target_face);

                    if (s0 >= num_spheres || s1 >= num_spheres || s2 >= num_spheres) continue;

                    bool has_sharp_vertex = false;
                    size_t face_sids[3] = {s0, s1, s2};
                    for (int k = 0; k < 3; k++) {
                        size_t sid = face_sids[k];
                        double curvature = sphere_curvatures[sid];
                        if (curvature != DBL_MAX && curvature >= kCurvatureShrinkThresholdRadians) {
                            has_sharp_vertex = true;
                            break;
                        }
                    }
                    
                    if (!has_sharp_vertex) continue;

                    double p0[4], p1[4], p2[4];
                    spheres->get_tree_point(s0, 4, p0);
                    spheres->get_tree_point(s1, 4, p1);
                    spheres->get_tree_point(s2, 4, p2);

                    double cx = (p0[0] + p1[0] + p2[0]) / 3.0;
                    double cy = (p0[1] + p1[1] + p2[1]) / 3.0;
                    double cz = (p0[2] + p1[2] + p2[2]) / 3.0;

                    //deadlock_positions.push_back({cx, cy, cz});
                    high_curvature_deadlocks_count++;
                }
                
                if (high_curvature_deadlocks_count > 0) {
                    std::cout << "\n=== [DEADLOCK DETECTED] ===" << std::endl;
                    std::cout << "Optimizer stalled: Found " << high_curvature_deadlocks_count 
                              << " high-curvature faces that cannot be shrunk further." << std::endl;
                    
                    // std::ofstream out_file("../deadlocks.xyz");
                    // if (out_file.is_open()) {
                    //     for (const auto& pos : deadlock_positions) {
                    //         out_file << pos[0] << " " << pos[1] << " " << pos[2] << "\n";
                    //     }
                    //     out_file.close();
                    //     std::cout << "Successfully exported " << high_curvature_deadlocks_count 
                    //               << " deadlock points to '../deadlocks.xyz'." << std::endl;
                    // } else {
                    //     std::cerr << "Error: Failed to open ../deadlocks.xyz for writing." << std::endl;
                    // }
                    std::cout << "===========================\n" << std::endl;
                }
            }
            break; // 优化器卡死，结束迭代
        }
    }

    std::ofstream blocked_log_file("shrink_blocked_report.txt");
    if (blocked_log_file.is_open()) {
        blocked_log_file << "# sphere shrink blocked report\n";
        blocked_log_file << "# format: iter, target_face, sphere, reason\n";
        for (const auto& row : shrink_blocked_logs) {
            blocked_log_file << row << "\n";
        }
        blocked_log_file.close();
        std::cout << "Shrink blocked report written to shrink_blocked_report.txt ("
                  << shrink_blocked_logs.size() << " records)" << std::endl;
    } else {
        std::cerr << "Failed to write shrink_blocked_report.txt" << std::endl;
    }

    // Debug export: sphere-id-indexed mesh for MeshLab picking.
    // In this OBJ:
    // - Vertex index (1-based) == sphere sid + 1
    // - Face indices are the sphere triplets from face_flat
    {
        std::ofstream sphere_obj("debug_sphere_faces.obj");
        if (sphere_obj.is_open()) {
            sphere_obj << std::fixed << std::setprecision(16);
            sphere_obj << "# Debug sphere-face mesh\n";
            sphere_obj << "# Vertex index (1-based) maps to sphere sid: sid = vid - 1\n";
            sphere_obj << "# Face indices are sphere sids + 1 from face_flat\n";

            for (size_t sid = 0; sid < num_spheres; sid++) {
                double sp[4];
                spheres->get_tree_point(sid, 4, sp);
                sphere_obj << "v " << sp[0] << " " << sp[1] << " " << sp[2] << "\n";
            }

            for (size_t fi = 0; fi < num_faces_total; fi++) {
                int a = face_flat[fi * 3 + 0];
                int b = face_flat[fi * 3 + 1];
                int c = face_flat[fi * 3 + 2];
                if (a < 0 || b < 0 || c < 0) continue;
                if ((size_t)a >= num_spheres || (size_t)b >= num_spheres || (size_t)c >= num_spheres) continue;
                sphere_obj << "f " << (a + 1) << " " << (b + 1) << " " << (c + 1) << "\n";
            }

            sphere_obj.close();
            std::cout << "Debug sphere-face mesh written to debug_sphere_faces.obj" << std::endl;
        } else {
            std::cerr << "Failed to write debug_sphere_faces.obj" << std::endl;
        }

        std::ofstream sphere_map("debug_spheres.txt");
        if (sphere_map.is_open()) {
            sphere_map << "# sid, x, y, z, r\n";
            sphere_map << std::fixed << std::setprecision(16);
            for (size_t sid = 0; sid < num_spheres; sid++) {
                double sp[4];
                spheres->get_tree_point(sid, 4, sp);
                sphere_map << sid << ", " << sp[0] << ", " << sp[1] << ", " << sp[2] << ", " << sp[3] << "\n";
            }
            sphere_map.close();
            std::cout << "Debug sphere map written to debug_spheres.txt" << std::endl;
        } else {
            std::cerr << "Failed to write debug_spheres.txt" << std::endl;
        }
    }

    std::cout << "=== Optimizer finished ===" << std::endl;
}

void Generator::generate_surface_mesh(MeshingTree* seeds, MeshingTree* spheres, const char* output_filename)
{
    std::cout << "Generating surface mesh using CGAL Voronoi..." << std::endl;

    
    size_t num_seeds = seeds->get_num_tree_points();
    if (num_seeds == 0) {
        std::cerr << "Error: No seeds provided." << std::endl;
        return;
    }

    // 1. Build Delaunay triangulation with seed index info
    Delaunay dt;
    std::vector<Vertex_handle> vertex_handles(num_seeds);
    std::vector<std::pair<size_t, size_t>> all_seed_pairs;
    std::set<std::pair<size_t, size_t>> strict_seed_pair_keys;
    std::vector<char> reconstruction_seed_indices(num_seeds, false);

    auto pair_key = [](size_t a, size_t b) {
        return (a < b) ? std::make_pair(a, b) : std::make_pair(b, a);
    };

    {
        std::set<std::pair<size_t, size_t>> seen_pairs;
        for (size_t i = 0; i < num_seeds; i++) {
            if (!seeds->tree_point_is_active(i)) continue;

            size_t* attrib = seeds->get_tree_point_attrib(i);
            size_t pair_idx = attrib[1];
            if (pair_idx >= num_seeds || !seeds->tree_point_is_active(pair_idx)) continue;

            size_t* pair_attrib = seeds->get_tree_point_attrib(pair_idx);
            if (attrib[5] == pair_attrib[5]) continue;

            auto key = pair_key(i, pair_idx);
            if (seen_pairs.insert(key).second) {
                all_seed_pairs.push_back(key);
            }
        }

        for (const auto& seed_pair : all_seed_pairs) {
            strict_seed_pair_keys.insert(seed_pair);
            reconstruction_seed_indices[seed_pair.first] = true;
            reconstruction_seed_indices[seed_pair.second] = true;
        }

        std::cout << "  * Selected all " << all_seed_pairs.size()
                  << " seed pairs for Voronoi reconstruction" << std::endl;

        if (all_seed_pairs.empty()) {
            std::cerr << "Error: No seed pairs available for reconstruction." << std::endl;
            return;
        }
    }

    auto build_delaunay = [&]() {
        dt.clear();
        std::fill(vertex_handles.begin(), vertex_handles.end(), Vertex_handle());
        for (size_t i = 0; i < num_seeds; i++) {
            if (!seeds->tree_point_is_active(i)) continue;
            if (!reconstruction_seed_indices[i]) continue;
            double* pt = seeds->get_tree_point(i);
            Vertex_handle vh = dt.insert(Point_3(pt[0], pt[1], pt[2]));
            vh->info() = i;
            vertex_handles[i] = vh;
        }
    };

    // Rebuild once to ensure dt/vertex_handles match final adjusted seed positions.
    build_delaunay();

    std::cout << "  * Delaunay triangulation built with " << dt.number_of_vertices() << " vertices" << std::endl;

    // // Export global Voronoi diagram
    // {
    //     std::cout << "  * Exporting global Voronoi diagram to global_voronoi.obj..." << std::endl;
        
    //     std::ofstream voronoi_file("global_voronoi.obj");
    //     if (!voronoi_file.is_open()) {
    //         std::cerr << "Error: Cannot create global_voronoi.obj" << std::endl;
    //     } else {
    //         voronoi_file << std::fixed << std::setprecision(16);
    //         voronoi_file << "# Global Voronoi Diagram\n";
    //         voronoi_file << "# Generated from Delaunay triangulation with " << dt.number_of_vertices() << " vertices\n\n";
            
    //         size_t vertex_offset = 1;
    //         size_t facet_count = 0;
            
    //         // Export all finite Voronoi facets
    //         for (auto eit = dt.finite_edges_begin(); eit != dt.finite_edges_end(); ++eit) {
    //             Cell_handle c = eit->first;
    //             int i1 = eit->second;
    //             int i2 = eit->third;
                
    //             // Get the dual Voronoi facet for this Delaunay edge
    //             std::vector<Point_3> facet_vertices;
                
    //             Delaunay::Cell_circulator cc = dt.incident_cells(*eit);
    //             Delaunay::Cell_circulator done = cc;
                
    //             if (cc != nullptr) {
    //                 do {
    //                     if (!dt.is_infinite(cc)) {
    //                         Point_3 center = dt.dual(cc);
    //                         facet_vertices.push_back(center);
    //                     }
    //                     ++cc;
    //                 } while (cc != done);
    //             }
                
    //             // Only export facets with at least 3 vertices and all coords within [-100, 100]
    //             if (facet_vertices.size() >= 3) {
    //                 bool all_in_range = true;
    //                 for (const auto& pt : facet_vertices) {
    //                     double px = CGAL::to_double(pt.x());
    //                     double py = CGAL::to_double(pt.y());
    //                     double pz = CGAL::to_double(pt.z());
    //                     if (fabs(px) > 100.0 || fabs(py) > 100.0 || fabs(pz) > 100.0) {
    //                         all_in_range = false;
    //                         break;
    //                     }
    //                 }
    //                 if (!all_in_range) continue;

    //                 // Write vertices
    //                 for (const auto& pt : facet_vertices) {
    //                     voronoi_file << "v " << CGAL::to_double(pt.x()) << " " 
    //                                << CGAL::to_double(pt.y()) << " " 
    //                                << CGAL::to_double(pt.z()) << "\n";
    //                 }
                    
    //                 // Write face as polygon
    //                 voronoi_file << "f";
    //                 for (size_t i = 0; i < facet_vertices.size(); ++i) {
    //                     voronoi_file << " " << (vertex_offset + i);
    //                 }
    //                 voronoi_file << "\n";
                    
    //                 vertex_offset += facet_vertices.size();
    //                 facet_count++;
    //             }
    //         }
            
    //         voronoi_file.close();
    //         std::cout << "  * Global Voronoi diagram exported: " << facet_count << " facets (global_voronoi.obj)" << std::endl;
    //     }
    // }




    {
        std::cout << "  * Exporting inside/outside seed points to OBJ files..." << std::endl;

        std::ofstream inside_seeds("inside_seeds.obj");
        std::ofstream outside_seeds("outside_seeds.obj");

        if (!inside_seeds.is_open() || !outside_seeds.is_open()) {
            std::cerr << "Error: Cannot create seed point output files" << std::endl;
        }
        else {
            inside_seeds << std::fixed << std::setprecision(16);
            outside_seeds << std::fixed << std::setprecision(16);

            inside_seeds << "# Inside Seed Points (seed + sphere centers + links)\n";
            outside_seeds << "# Outside Seed Points (seed + sphere centers + links)\n";

            size_t inside_count = 0;
            size_t outside_count = 0;
            size_t inside_vertex_offset = 1;
            size_t outside_vertex_offset = 1;
            size_t num_spheres = spheres ? static_cast<size_t>(spheres->get_num_tree_points()) : 0;

            auto write_seed_with_spheres = [&](std::ofstream& out,
                                               size_t& vertex_offset,
                                               double* seed_pt,
                                               size_t* attrib) {
                size_t seed_vertex_index = vertex_offset;
                out << "v " << seed_pt[0] << " " << seed_pt[1] << " " << seed_pt[2] << "\n";
                vertex_offset++;

                if (!spheres) return;

                size_t sphere_ids[3] = { attrib[2], attrib[3], attrib[4] };
                for (size_t si : sphere_ids) {
                    if (si >= num_spheres) continue;
                    double sphere_pt[4];
                    spheres->get_tree_point(si, 4, sphere_pt);
                    size_t sphere_vertex_index = vertex_offset;
                    out << "v " << sphere_pt[0] << " " << sphere_pt[1] << " " << sphere_pt[2] << "\n";
                    out << "l " << seed_vertex_index << " " << sphere_vertex_index << "\n";
                    vertex_offset++;
                }
            };

            for (size_t i = 0; i < num_seeds; i++) {
                if (!seeds->tree_point_is_active(i)) continue;

                double* pt = seeds->get_tree_point(i);
                size_t* attrib = seeds->get_tree_point_attrib(i);

                // attrib[5] is the region label: 0 = outside, 1 = inside.
                if (attrib[5] == 1) {
                    write_seed_with_spheres(inside_seeds, inside_vertex_offset, pt, attrib);
                    inside_count++;
                }
                else {
                    write_seed_with_spheres(outside_seeds, outside_vertex_offset, pt, attrib);
                    outside_count++;
                }
            }

            inside_seeds.close();
            outside_seeds.close();

            std::cout << "  * Inside seeds: " << inside_count << " (saved to inside_seeds.obj)" << std::endl;
            std::cout << "  * Outside seeds: " << outside_count << " (saved to outside_seeds.obj)" << std::endl;
        }
    }

    {
        std::cout << "  * Exporting seed pair connections to seed_pairs.obj ..." << std::endl;
        
        std::ofstream pairs_out("seed_pairs.obj");
        if (pairs_out.is_open()) {
            pairs_out << std::fixed << std::setprecision(16);
            pairs_out << "# Seed Pair Connections\n";
            
            // 棣栧厛杈撳嚭鎵€鏈夌瀛愮偣浣滀负椤剁偣
            for (size_t i = 0; i < num_seeds; i++) {
                if (!seeds->tree_point_is_active(i)) continue;
                double* pt = seeds->get_tree_point(i);
                pairs_out << "v " << pt[0] << " " << pt[1] << " " << pt[2] << "\n";
            }
            
            // 鐒跺悗杈撳嚭閰嶅杩炴帴绾?
            size_t pair_count = 0;
            std::set<std::pair<size_t, size_t>> processed_pairs;
            std::vector<double> pair_dists;
            pair_dists.reserve(num_seeds / 2);
            
            for (size_t i = 0; i < num_seeds; i++) {
                if (!seeds->tree_point_is_active(i)) continue;
                
                size_t* attrib = seeds->get_tree_point_attrib(i);
                size_t pair_idx = attrib[1]; // attrib[1] 鏄厤瀵圭瀛愮储寮?
                
                if (pair_idx < num_seeds && seeds->tree_point_is_active(pair_idx)) {
                    // 閬垮厤閲嶅杈撳嚭鍚屼竴瀵?
                    auto pair_key = std::minmax(i, pair_idx);
                    if (processed_pairs.find(pair_key) == processed_pairs.end()) {
                        // OBJ浣跨敤1-based绱㈠紩
                        pairs_out << "l " << (i + 1) << " " << (pair_idx + 1) << "\n";
                        processed_pairs.insert(pair_key);
                        pair_count++;

                        double* p0 = seeds->get_tree_point(i);
                        double* p1 = seeds->get_tree_point(pair_idx);
                        double dx = p0[0] - p1[0];
                        double dy = p0[1] - p1[1];
                        double dz = p0[2] - p1[2];
                        pair_dists.push_back(std::sqrt(dx * dx + dy * dy + dz * dz));
                    }
                }
            }
            
            pairs_out.close();
            std::cout << "  * Seed pairs: " << pair_count << " connections (saved to seed_pairs.obj)" << std::endl;

            if (!pair_dists.empty()) {
                double sum = 0.0;
                for (double d : pair_dists) sum += d;
                double mean = sum / static_cast<double>(pair_dists.size());

                const size_t n = pair_dists.size();
                const size_t mid = n / 2;
                std::nth_element(pair_dists.begin(), pair_dists.begin() + mid, pair_dists.end());
                double median = pair_dists[mid];
                if (n % 2 == 0) {
                    double lower = *std::max_element(pair_dists.begin(), pair_dists.begin() + mid);
                    median = 0.5 * (lower + median);
                }

                std::cout << "  * Seed pair distance: mean = " << mean
                          << ", median = " << median << std::endl;
            }
        } else {
            std::cerr << "Error: Cannot write to seed_pairs.obj" << std::endl;
        }
    }
    // 2. Surface extraction is now handled by the topological Delaunay-Voronoi dual
    //    in exportBoundarySurfaceImproved() below.  We intentionally do not build
    //    a separate point-only facet list here, because that was the source of
    //    per-face duplicated OBJ vertices.

    // 3. Write Surface Mesh OBJ file
    /*
    std::ofstream obj_file(output_filename);
    if (!obj_file.is_open()) {
        std::cerr << "Error: Cannot open output file " << output_filename << std::endl;
        return;
    }
    
    obj_file << "# Voronoi surface mesh generated by VoroCrust" << std::endl;
    obj_file << "# Number of facets: " << voronoi_facets.size() << std::endl;
    
    // Use a map to avoid duplicate vertices
    // Note: Using Exact kernel points as map keys might be slow or tricky with rounding.
    // Your rounding strategy is good for merging close points.
    std::map<std::tuple<double, double, double>, size_t> vertex_map;
    std::vector<Point_3> unique_vertices;
    std::vector<std::vector<size_t>> face_indices;
    
    auto get_vertex_index = [&](const Point_3& p) -> size_t {
        // Round to avoid floating point precision issues
        double x = std::round(CGAL::to_double(p.x()) * 1e10) / 1e10;
        double y = std::round(CGAL::to_double(p.y()) * 1e10) / 1e10;
        double z = std::round(CGAL::to_double(p.z()) * 1e10) / 1e10;
        auto key = std::make_tuple(x, y, z);
        
        auto it = vertex_map.find(key);
        if (it != vertex_map.end()) {
            return it->second;
        }
        size_t idx = unique_vertices.size();
        vertex_map[key] = idx;
        unique_vertices.push_back(p);
        return idx;
    };
    
    // Process facets and triangulate polygons
    for (const auto& facet : voronoi_facets) {
        if (facet.size() < 3) continue;
        
        std::vector<size_t> indices;
        for (const auto& pt : facet) {
            indices.push_back(get_vertex_index(pt));
        }
        
        // Fan triangulation for polygon with more than 3 vertices
        // Improve triangulation: use centroid fan to create better triangles
        Point_3 centroid = CGAL::ORIGIN;
        for(const auto& pt : facet) centroid = centroid + (pt - CGAL::ORIGIN);
        double s = static_cast<double>(facet.size());
        centroid = Point_3(centroid.x() / s, centroid.y() / s, centroid.z() / s);
        size_t centroid_idx = get_vertex_index(centroid);

        for (size_t i = 0; i < indices.size(); i++) {
            size_t idx0 = indices[i];
            size_t idx1 = indices[(i + 1) % indices.size()];
            
            // 璺宠繃閫€鍖栦笁瑙掑舰
            if (idx0 == idx1 || idx0 == centroid_idx || idx1 == centroid_idx) continue;

            face_indices.push_back({centroid_idx, idx0, idx1});
        }
    }
    
    // Write vertices
    for (const auto& v : unique_vertices) {
        obj_file << "v " << std::setprecision(16) 
                 << CGAL::to_double(v.x()) << " " 
                 << CGAL::to_double(v.y()) << " " 
                 << CGAL::to_double(v.z()) << std::endl;
    }
    
    // Write faces (OBJ uses 1-based indexing)
    for (const auto& face : face_indices) {
        obj_file << "f " << (face[0] + 1) << " " << (face[1] + 1) << " " << (face[2] + 1) << std::endl;
    }
    
    obj_file.close();
    
    std::cout << "  * Surface mesh saved to " << output_filename << std::endl;
    std::cout << "  * Total vertices: " << unique_vertices.size() << ", triangles: " << face_indices.size() << std::endl;
    */

    {
        std::map<Vertex_handle, int> vertex_labels;
        for (size_t i = 0; i < num_seeds; i++) {
            if (!seeds->tree_point_is_active(i)) continue;
            Vertex_handle vh = vertex_handles[i];
            if (vh == Vertex_handle()) continue;
            size_t* attrib = seeds->get_tree_point_attrib(i);
            vertex_labels[vh] = static_cast<int>(attrib[5]);
        }

        auto write_voronoi_facets_original_obj = [&](
            const std::string& filename,
            const std::vector<std::vector<Point_3>>& voronoi_facets,
            const std::vector<std::array<size_t, 8>>* facet_sources = nullptr,
            const std::string& sources_filename = "voronoi_original_map.csv"
        ) {
            std::ofstream out(filename);
            if (!out) {
                std::cerr << "Cannot open file " << filename << std::endl;
                return;
            }

            out << "# Original strict-pair Voronoi polygon facets OBJ\n";
            out << "# Same printed coordinates share one OBJ vertex index.\n";
            out << "# No epsilon merge, no coordinate-neighborhood deduplication, no triangulation.\n";
            out << "# Face duplicate indices are cleaned only by removing consecutive and closing duplicates; ring order is preserved.\n";

            // 注意：这里不是 epsilon 合并。
            // 只把最终 OBJ 中打印出来完全相同的坐标复用为同一个 vertex index。
            // 这可以解决每个面重复创建同一坐标顶点的问题。
            std::map<std::string, int> vertex_key_to_index;
            std::vector<Point_3> vertices;
            std::vector<std::vector<int>> faces;
            std::vector<size_t> face_src_facet;

            vertices.reserve(voronoi_facets.size() * 4 + 16);
            faces.reserve(voronoi_facets.size());
            face_src_facet.reserve(voronoi_facets.size());

            auto make_vertex_key = [&](const Point_3& p) -> std::string {
                std::ostringstream oss;
                oss << std::fixed << std::setprecision(16)
                    << CGAL::to_double(p.x()) << " "
                    << CGAL::to_double(p.y()) << " "
                    << CGAL::to_double(p.z());
                return oss.str();
            };

            auto get_vertex_index = [&](const Point_3& p) -> int {
                std::string key = make_vertex_key(p);
                auto it = vertex_key_to_index.find(key);
                if (it != vertex_key_to_index.end()) {
                    return it->second;
                }

                int idx = static_cast<int>(vertices.size()) + 1; // OBJ is 1-based
                vertex_key_to_index[key] = idx;
                vertices.push_back(p);
                return idx;
            };

            size_t input_facets = 0;
            size_t exported_faces = 0;
            size_t dropped_lt3 = 0;
            size_t duplicate_refs_removed = 0;

            for (size_t facet_idx = 0; facet_idx < voronoi_facets.size(); ++facet_idx) {
                const auto& facet = voronoi_facets[facet_idx];
                input_facets++;

                if (facet.size() < 3) {
                    dropped_lt3++;
                    continue;
                }

                std::vector<int> face;
                face.reserve(facet.size());
                for (const auto& p : facet) {
                    face.push_back(get_vertex_index(p));
                }

                const size_t before_cleanup = face.size();

                // 只删除连续重复 index 和首尾重复 index，保留 incident_cells 给出的环序。
                // 不做 sort，不删除非连续重复 index，也不因为 face 内存在重复 index 就跳过整张面。
                std::vector<int> cleaned_face;
                cleaned_face.reserve(face.size());
                for (int idx : face) {
                    if (idx <= 0) continue;
                    if (!cleaned_face.empty() && cleaned_face.back() == idx) {
                        duplicate_refs_removed++;
                        continue;
                    }
                    cleaned_face.push_back(idx);
                }
                if (cleaned_face.size() > 1 && cleaned_face.front() == cleaned_face.back()) {
                    cleaned_face.pop_back();
                    duplicate_refs_removed++;
                }

                if (cleaned_face.size() < 3) {
                    dropped_lt3++;
                    continue;
                }

                faces.push_back(std::move(cleaned_face));
                face_src_facet.push_back(facet_idx);
                exported_faces++;
            }

            out << std::fixed << std::setprecision(16);

            for (const auto& p : vertices) {
                out << "v "
                    << CGAL::to_double(p.x()) << " "
                    << CGAL::to_double(p.y()) << " "
                    << CGAL::to_double(p.z()) << "\n";
            }

            for (const auto& f : faces) {
                out << "f";
                for (int idx : f) {
                    out << " " << idx;
                }
                out << "\n";
            }

            std::cout << "  * Original strict-pair Voronoi OBJ exported: " << filename << std::endl;
            std::cout << "    input_facets = " << input_facets << std::endl;
            std::cout << "    exported_faces = " << exported_faces << std::endl;
            std::cout << "    vertices = " << vertices.size() << std::endl;
            std::cout << "    dropped_lt3 = " << dropped_lt3 << std::endl;
            std::cout << "    duplicate_refs_removed = " << duplicate_refs_removed << std::endl;

            if (facet_sources && facet_sources->size() == voronoi_facets.size()) {
                std::ofstream meta_out(sources_filename);
                if (meta_out.is_open()) {
                    meta_out << "face_index_0based,face_index_1based,src_facet_index,seed_a,seed_b,"
                             << "seed_a_s0,seed_a_s1,seed_a_s2,seed_b_s0,seed_b_s1,seed_b_s2\n";
                    for (size_t face_idx = 0; face_idx < face_src_facet.size(); ++face_idx) {
                        size_t src = face_src_facet[face_idx];
                        const auto& m = (*facet_sources)[src];
                        meta_out << face_idx << "," << (face_idx + 1) << "," << src << ","
                                 << m[0] << "," << m[1] << ","
                                 << m[2] << "," << m[3] << "," << m[4] << ","
                                 << m[5] << "," << m[6] << "," << m[7] << "\n";
                    }
                }
            }
        };

        auto write_voronoi_facets_to_obj = [&](
            const std::string& filename,
            const std::vector<std::vector<Point_3>>& voronoi_facets
        ) {
            std::ofstream out(filename);
            if (!out) {
                std::cerr << "Cannot open file " << filename << std::endl;
                return;
            }

            out << "# Voronoi facets OBJ\n";
            size_t vertex_offset = 1;

            for (const auto& facet : voronoi_facets) {
                if (facet.size() < 3) continue;

                for (const auto& p : facet) {
                    out << "v " << p.x() << " " << p.y() << " " << p.z() << "\n";
                }

                out << "f";
                for (size_t i = 0; i < facet.size(); ++i) {
                    out << " " << (vertex_offset + i);
                }
                out << "\n";

                vertex_offset += facet.size();
            }
        };

        auto write_voronoi_facets_triangulated_obj = [&](
            const std::string& filename,
            const std::vector<std::vector<Point_3>>& voronoi_facets,
            const std::vector<std::array<size_t, 8>>* facet_sources = nullptr,
            const std::string& sources_filename = "voronoi_triangulated_map.csv"
        ) {
            std::ofstream out(filename);
            if (!out) {
                std::cerr << "Cannot open file " << filename << std::endl;
                return;
            }

            out << "# Triangulated Voronoi facets\n";

            std::map<Point_3, int> vertex_index;
            std::vector<Point_3> vertices;
            std::vector<std::array<int, 3>> triangles;
            std::vector<size_t> triangle_src_facet;

            auto get_index = [&](const Point_3& p) {
                auto it = vertex_index.find(p);
                if (it != vertex_index.end()) return it->second;
                int idx = static_cast<int>(vertices.size()) + 1;
                vertex_index[p] = idx;
                vertices.push_back(p);
                return idx;
            };

            for (size_t facet_idx = 0; facet_idx < voronoi_facets.size(); facet_idx++) {
                const auto& facet = voronoi_facets[facet_idx];
                if (facet.size() < 3) continue;
                int v0 = get_index(facet[0]);
                for (size_t i = 1; i + 1 < facet.size(); ++i) {
                    int v1 = get_index(facet[i]);
                    int v2 = get_index(facet[i + 1]);
                    triangles.push_back({ v0, v1, v2 });
                    triangle_src_facet.push_back(facet_idx);
                }
            }

            for (const auto& p : vertices) {
                out << "v " << p.x() << " " << p.y() << " " << p.z() << "\n";
            }
            for (const auto& t : triangles) {
                out << "f " << t[0] << " " << t[1] << " " << t[2] << "\n";
            }

            if (facet_sources && facet_sources->size() == voronoi_facets.size()) {
                std::ofstream meta_out(sources_filename);
                if (meta_out.is_open()) {
                    meta_out << "tri_index_0based,tri_index_1based,src_facet_index,seed_a,seed_b,"
                             << "seed_a_s0,seed_a_s1,seed_a_s2,seed_b_s0,seed_b_s1,seed_b_s2\n";
                    for (size_t tri_idx = 0; tri_idx < triangle_src_facet.size(); tri_idx++) {
                        size_t src = triangle_src_facet[tri_idx];
                        const auto& m = (*facet_sources)[src];
                        meta_out << tri_idx << "," << (tri_idx + 1) << "," << src << ","
                                 << m[0] << "," << m[1] << ","
                                 << m[2] << "," << m[3] << "," << m[4] << ","
                                 << m[5] << "," << m[6] << "," << m[7] << "\n";
                    }
                }
            }
        };

        auto write_voronoi_facets_to_obj_dedup = [&](
            const std::string& filename,
            const std::vector<std::vector<Point_3>>& voronoi_facets,
            const std::vector<std::array<size_t, 8>>* facet_sources = nullptr,
            const std::string& sources_filename = "voronoi_dedup_map.csv",
            double epsilon = 1.05e-5, // 1e-5 的安全放大，适合单精度坐标
            double angle_threshold_deg = 170.0,
            bool debug = false
        ) {
            std::cout << "  * Starting optimized deduplication export to " << filename << "..." << std::endl;
            // 优化 IO：使用 buffer 稍微大一点，虽然 ofstream 默认有 buffer
            std::ofstream out(filename);
            if (!out) {
                std::cerr << "Cannot open file " << filename << std::endl;
                return;
            }
            out << "# Voronoi facets (epsilon=" << epsilon
                << ", angle_threshold=" << angle_threshold_deg << "deg)\n";

            std::vector<Point_3> vertices;
            // 预分配内存，避免频繁 realloc。假设平均每个面由 4 个顶点组成
            vertices.reserve(voronoi_facets.size() * 2); 
            
            std::vector<std::vector<int>> faces;
            faces.reserve(voronoi_facets.size());
            std::vector<size_t> face_src_facet;
            face_src_facet.reserve(voronoi_facets.size());
            
            // --- 优化 1：使用 unordered_map 替代 map (O(logN) -> O(1)) ---
            // 定义 Key 结构体以避免 std::tuple 的开销
            struct VertexKey {
                long long x, y, z;
                
                // 重载 == 操作符
                bool operator==(const VertexKey& other) const {
                    return x == other.x && y == other.y && z == other.z;
                }
            };

            // 自定义 Hash 函数
            struct VertexKeyHash {
                std::size_t operator()(const VertexKey& k) const {
                    // 使用简单的位运算混合 hash，比 boost::hash_combine 快且无依赖
                    size_t h1 = std::hash<long long>{}(k.x);
                    size_t h2 = std::hash<long long>{}(k.y);
                    size_t h3 = std::hash<long long>{}(k.z);
                    return h1 ^ (h2 << 1) ^ (h3 << 2); // 简单混合
                }
            };

            // 使用 unordered_map
            std::unordered_map<VertexKey, int, VertexKeyHash> vertex_map;
            // 预分配 map 空间以减少 rehash
            vertex_map.reserve(voronoi_facets.size() * 3);

            double scale_factor = 1.0 / (epsilon > 0 ? epsilon : 1e-6);

            int total_input_facets = 0;
            int total_output_facets = 0;

            auto find_or_add_vertex = [&](const Point_3& p) -> int {
                double x = CGAL::to_double(p.x());
                double y = CGAL::to_double(p.y());
                double z = CGAL::to_double(p.z());
                
                // 坐标离散化逻辑保持不变
                long long ix = std::llround(x * scale_factor);
                long long iy = std::llround(y * scale_factor);
                long long iz = std::llround(z * scale_factor);
                
                VertexKey key{ix, iy, iz};
                
                auto it = vertex_map.find(key);
                if (it != vertex_map.end()) {
                    return it->second;
                }
                
                vertices.push_back(p);
                // vertices.size() 返回的是 size_t，转换为 int
                int idx = static_cast<int>(vertices.size()); 
                vertex_map[key] = idx; // OBJ 索引通常从 1 开始，这里存的是内部索引，写入时可能要注意
                return idx;
            };

            // 角度计算逻辑保持不变，但为了性能，内联并减少对象创建
            auto compute_angle = [&](const Point_3& p1, const Point_3& p2, const Point_3& p3) -> double {
                // p1 -> prev, p2 -> curr, p3 -> next
                // Vector_3 构造可能涉及精确运算，这里直接转 double 计算以加速 (逻辑保持一致)
                double v1x = CGAL::to_double(p1.x() - p2.x());
                double v1y = CGAL::to_double(p1.y() - p2.y());
                double v1z = CGAL::to_double(p1.z() - p2.z());

                double v2x = CGAL::to_double(p3.x() - p2.x());
                double v2y = CGAL::to_double(p3.y() - p2.y());
                double v2z = CGAL::to_double(p3.z() - p2.z());

                double len1_sq = v1x*v1x + v1y*v1y + v1z*v1z;
                double len2_sq = v2x*v2x + v2y*v2y + v2z*v2z;

                if (len1_sq < 1e-20 || len2_sq < 1e-20) return 0.0;

                double dot = v1x*v2x + v1y*v2y + v1z*v2z;
                double cos_angle = dot / std::sqrt(len1_sq * len2_sq);
                
                // Clamp
                if (cos_angle > 1.0) cos_angle = 1.0;
                else if (cos_angle < -1.0) cos_angle = -1.0;

                return std::acos(cos_angle) * 180.0 / M_PI;
            };

            // --- 优化 2：减少循环内的内存分配 ---
            // 将临时容器提到循环外
            std::vector<Point_3> corner_points;
            corner_points.reserve(16); // 预估多边形顶点数
            std::vector<int> face_indices;
            face_indices.reserve(16);
            std::vector<int> cleaned_face;
            cleaned_face.reserve(16);

            // 处理所有面片
            for (size_t facet_idx = 0; facet_idx < voronoi_facets.size(); facet_idx++) {
                const auto& facet = voronoi_facets[facet_idx];
                total_input_facets++;
                if (facet.size() < 3) continue;

                // --- 提取角点逻辑 (内联以复用 corner_points) ---
                corner_points.clear();
                int n = static_cast<int>(facet.size());
                
                // 快速路径：如果是三角形，不需要计算角度（通常都是角点，或者逻辑上保留原样）
                // 原逻辑是 <=3 直接返回。
                if (n <= 3) {
                    corner_points = facet; 
                } else {
                    bool has_corner = false;
                    for (int i = 0; i < n; ++i) {
                        const Point_3& prev = facet[(i - 1 + n) % n];
                        const Point_3& curr = facet[i];
                        const Point_3& next = facet[(i + 1) % n];
                        
                        if (compute_angle(prev, curr, next) < angle_threshold_deg) {
                            corner_points.push_back(curr);
                            has_corner = true;
                        }
                    }
                    // 如果没有检测到角点（例如是个圆），保留原始多边形以防丢失
                    if (corner_points.size() < 3) {
                        corner_points = facet; 
                    }
                }
                // -------------------------------------------

                if (corner_points.size() < 3) continue;

                face_indices.clear();
                for (const auto& p : corner_points) {
                    face_indices.push_back(find_or_add_vertex(p));
                }

                // 清理连续重复索引
                cleaned_face.clear();
                if (!face_indices.empty()) {
                    cleaned_face.push_back(face_indices[0]);
                    for (size_t i = 1; i < face_indices.size(); ++i) {
                        if (face_indices[i] != face_indices[i - 1]) {
                            cleaned_face.push_back(face_indices[i]);
                        }
                    }
                    // 检查首尾是否闭合
                    if (cleaned_face.size() > 1 && cleaned_face.front() == cleaned_face.back()) {
                        cleaned_face.pop_back();
                    }
                }

                if (cleaned_face.size() >= 3) {
                    // 检查非连续重复顶点（epsilon合并可能导致）
                    bool has_duplicate = false;
                    for (size_t i = 0; i < cleaned_face.size() && !has_duplicate; ++i) {
                        for (size_t j = i + 1; j < cleaned_face.size() && !has_duplicate; ++j) {
                            if (cleaned_face[i] == cleaned_face[j]) has_duplicate = true;
                        }
                    }
                    if (!has_duplicate) {
                        faces.push_back(cleaned_face);
                        face_src_facet.push_back(facet_idx);
                        total_output_facets++;
                    }
                }
            }

            // 写入文件
            // 优化 3：使用 '\n' 替代 std::endl 以避免频繁 flush
            out << std::fixed << std::setprecision(16);
            for (const auto& p : vertices) {
                out << "v " << CGAL::to_double(p.x()) << " " 
                            << CGAL::to_double(p.y()) << " " 
                            << CGAL::to_double(p.z()) << "\n";
            }
            // OBJ 索引从 1 开始
            for (const auto& f : faces) {
                out << "f";
                for (int idx : f) out << " " << (idx); // 你的 find_or_add 返回的是 size，OBJ 是 1-based，如果你之前的逻辑是直接用 size 做索引，那这里需要注意。通常 OBJ 索引需要 +1，但如果你之前的代码生成的 OBJ 能用，那就保持原样。
                // *注意*：你原来的代码 idx 是 vertices.size()，如果它是 1, 2, 3... 那就是 1-based。
                // 这里的实现 vertices.size() 在 push 之后返回的是 1, 2, 3... 所以是 1-based。无需修改。
                out << "\n";
            }

            std::cout << "  * Exported " << total_output_facets << " facets (" 
                      << vertices.size() << " vertices) to " << filename << std::endl;

            if (facet_sources && facet_sources->size() == voronoi_facets.size()) {
                std::ofstream meta_out(sources_filename);
                if (meta_out.is_open()) {
                    meta_out << "face_index_0based,face_index_1based,src_facet_index,seed_a,seed_b,"
                             << "seed_a_s0,seed_a_s1,seed_a_s2,seed_b_s0,seed_b_s1,seed_b_s2\n";
                    for (size_t face_idx = 0; face_idx < face_src_facet.size(); face_idx++) {
                        size_t src = face_src_facet[face_idx];
                        const auto& m = (*facet_sources)[src];
                        meta_out << face_idx << "," << (face_idx + 1) << "," << src << ","
                                 << m[0] << "," << m[1] << ","
                                 << m[2] << "," << m[3] << "," << m[4] << ","
                                 << m[5] << "," << m[6] << "," << m[7] << "\n";
                    }
                }
            }
        };

        auto export_single_voronoi_polygon = [&](Vertex_handle v1, Vertex_handle v2, const std::string& filename) {
            std::cout << "=== Exporting single Voronoi polygon ===" << std::endl;
            std::cout << "  * Edge: seed " << v1->info() << " <-> seed " << v2->info() << std::endl;

            std::ofstream obj_file(filename);
            if (!obj_file.is_open()) {
                std::cerr << "Error: Cannot open " << filename << std::endl;
                return;
            }

            obj_file << std::fixed << std::setprecision(16);
            obj_file << "# Voronoi polygon for Delaunay edge\n";
            obj_file << "# Seed 1 (index " << v1->info() << ")\n";
            obj_file << "# Seed 2 (index " << v2->info() << ")\n\n";

            std::vector<Point_3> polygon_vertices;
            std::vector<Cell_handle> cells_v1;
            dt.incident_cells(v1, std::back_inserter(cells_v1));

            for (auto cell : cells_v1) {
                if (dt.is_infinite(cell)) continue;

                bool contains_v2 = false;
                for (int i = 0; i < 4; i++) {
                    if (cell->vertex(i) == v2) {
                        contains_v2 = true;
                        break;
                    }
                }

                if (contains_v2) {
                    Point_3 center = dt.dual(cell);
                    polygon_vertices.push_back(center);
                }
            }

            std::cout << "  * Polygon has " << polygon_vertices.size() << " vertices" << std::endl;

            if (polygon_vertices.size() < 3) {
                std::cout << "  * Warning: Not enough vertices to form a polygon" << std::endl;
                obj_file.close();
                return;
            }

            for (size_t i = 0; i < polygon_vertices.size(); i++) {
                const auto& v = polygon_vertices[i];
                obj_file << "v " << CGAL::to_double(v.x()) << " "
                    << CGAL::to_double(v.y()) << " "
                    << CGAL::to_double(v.z()) << "\n";
            }

            obj_file << "\n# Polygon edges\n";
            for (size_t i = 0; i < polygon_vertices.size(); i++) {
                size_t next_i = (i + 1) % polygon_vertices.size();
                obj_file << "l " << (i + 1) << " " << (next_i + 1) << "\n";
            }

            obj_file << "\n# Polygon face\n";
            obj_file << "f";
            for (size_t i = 0; i < polygon_vertices.size(); i++) {
                obj_file << " " << (i + 1);
            }
            obj_file << "\n";

            obj_file << "\n# Corresponding seed points\n";
            obj_file << "# v " << CGAL::to_double(v1->point().x()) << " "
                << CGAL::to_double(v1->point().y()) << " "
                << CGAL::to_double(v1->point().z()) << " # Seed 1\n";
            obj_file << "# v " << CGAL::to_double(v2->point().x()) << " "
                << CGAL::to_double(v2->point().y()) << " "
                << CGAL::to_double(v2->point().z()) << " # Seed 2\n";

            obj_file.close();
            std::cout << "  * Saved to " << filename << std::endl;
        };

        // Directly compute the common Voronoi-cell face for each strict seed pair.
        // This implements:
        //   V_i ∩ V_j = Π_ij ∩ ⋂_{k≠i,j}{ x | ||x-s_i|| <= ||x-s_k|| }
        // where Π_ij is the perpendicular bisector plane of the seed pair.
        // It does NOT use Delaunay tetrahedron circumcenters to form the polygon vertices.
        auto exportBoundarySurfaceImproved = [&]() {
            std::cout << "\n=== Extracting Boundary Surface (Direct Voronoi Cell Intersection, Strict Pair) ===" << std::endl;

            struct Vec3 {
                double x, y, z;
            };

            auto make_vec = [](double x, double y, double z) -> Vec3 {
                return Vec3{x, y, z};
            };
            auto add = [](const Vec3& a, const Vec3& b) -> Vec3 {
                return Vec3{a.x + b.x, a.y + b.y, a.z + b.z};
            };
            auto sub = [](const Vec3& a, const Vec3& b) -> Vec3 {
                return Vec3{a.x - b.x, a.y - b.y, a.z - b.z};
            };
            auto mul = [](const Vec3& a, double s) -> Vec3 {
                return Vec3{a.x * s, a.y * s, a.z * s};
            };
            auto dot = [](const Vec3& a, const Vec3& b) -> double {
                return a.x * b.x + a.y * b.y + a.z * b.z;
            };
            auto cross = [](const Vec3& a, const Vec3& b) -> Vec3 {
                return Vec3{
                    a.y * b.z - a.z * b.y,
                    a.z * b.x - a.x * b.z,
                    a.x * b.y - a.y * b.x
                };
            };
            auto norm2 = [&](const Vec3& a) -> double {
                return dot(a, a);
            };
            auto norm = [&](const Vec3& a) -> double {
                return std::sqrt(norm2(a));
            };
            auto normalize = [&](const Vec3& a) -> Vec3 {
                double l = norm(a);
                if (l <= 1e-300) return Vec3{0.0, 0.0, 0.0};
                return mul(a, 1.0 / l);
            };
            auto lerp = [&](const Vec3& a, const Vec3& b, double t) -> Vec3 {
                return Vec3{
                    a.x + (b.x - a.x) * t,
                    a.y + (b.y - a.y) * t,
                    a.z + (b.z - a.z) * t
                };
            };
            auto sqdist = [&](const Vec3& a, const Vec3& b) -> double {
                return norm2(sub(a, b));
            };

            // ------------------------------------------------------------------
            // Step 1. Cache all reconstruction seed positions.
            // The Voronoi diagram here is defined only by reconstruction seeds,
            // i.e. active seeds participating in strict inside/outside pairs.
            // ------------------------------------------------------------------
            std::vector<Vec3> seed_pos(num_seeds, Vec3{0.0, 0.0, 0.0});
            std::vector<char> seed_pos_valid(num_seeds, false);
            std::vector<size_t> reconstruction_seed_ids;
            reconstruction_seed_ids.reserve(static_cast<size_t>(dt.number_of_vertices()));

            Vec3 bb_min{ DBL_MAX, DBL_MAX, DBL_MAX };
            Vec3 bb_max{-DBL_MAX,-DBL_MAX,-DBL_MAX };

            for (size_t sid = 0; sid < num_seeds; ++sid) {
                if (!seeds->tree_point_is_active(sid)) continue;
                if (!reconstruction_seed_indices[sid]) continue;

                double* p = seeds->get_tree_point(sid);
                Vec3 q{p[0], p[1], p[2]};
                seed_pos[sid] = q;
                seed_pos_valid[sid] = true;
                reconstruction_seed_ids.push_back(sid);

                bb_min.x = std::min(bb_min.x, q.x);
                bb_min.y = std::min(bb_min.y, q.y);
                bb_min.z = std::min(bb_min.z, q.z);
                bb_max.x = std::max(bb_max.x, q.x);
                bb_max.y = std::max(bb_max.y, q.y);
                bb_max.z = std::max(bb_max.z, q.z);
            }

            if (reconstruction_seed_ids.empty()) {
                std::cerr << "Error: no reconstruction seed positions available." << std::endl;
                return;
            }

            const double bbox_diag = std::max(1e-12, norm(sub(bb_max, bb_min)));
            const double initial_half_size = 8.0 * bbox_diag;
            const double clip_tol = 1e-12 * std::max(1.0, bbox_diag * bbox_diag);
            const double cleanup_eps = 1e-12 * bbox_diag;
            const double cleanup_eps2 = cleanup_eps * cleanup_eps;

            // ------------------------------------------------------------------
            // Clip a convex polygon by one Voronoi halfspace:
            //   ||x-a||^2 <= ||x-k||^2
            // equivalent to:
            //   2(k-a)·x <= |k|^2 - |a|^2
            // We represent the inside value as:
            //   value(x) = rhs - 2(k-a)·x >= 0
            // ------------------------------------------------------------------
            auto clip_polygon_by_seed_halfspace = [&](const std::vector<Vec3>& poly,
                                                       const Vec3& a,
                                                       const Vec3& k) -> std::vector<Vec3> {
                if (poly.empty()) return {};

                const Vec3 ka = sub(k, a);
                const double rhs = norm2(k) - norm2(a);

                auto value = [&](const Vec3& x) -> double {
                    return rhs - 2.0 * dot(ka, x);
                };

                std::vector<Vec3> out;
                out.reserve(poly.size() + 2);

                Vec3 prev = poly.back();
                double f_prev = value(prev);
                bool prev_inside = (f_prev >= -clip_tol);

                for (const Vec3& curr : poly) {
                    double f_curr = value(curr);
                    bool curr_inside = (f_curr >= -clip_tol);

                    if (prev_inside && curr_inside) {
                        out.push_back(curr);
                    } else if (prev_inside && !curr_inside) {
                        double denom = f_prev - f_curr;
                        if (std::abs(denom) > 1e-300) {
                            double t = f_prev / denom;
                            if (t < 0.0) t = 0.0;
                            if (t > 1.0) t = 1.0;
                            out.push_back(lerp(prev, curr, t));
                        }
                    } else if (!prev_inside && curr_inside) {
                        double denom = f_prev - f_curr;
                        if (std::abs(denom) > 1e-300) {
                            double t = f_prev / denom;
                            if (t < 0.0) t = 0.0;
                            if (t > 1.0) t = 1.0;
                            out.push_back(lerp(prev, curr, t));
                        }
                        out.push_back(curr);
                    }

                    prev = curr;
                    f_prev = f_curr;
                    prev_inside = curr_inside;
                }

                return out;
            };

            auto cleanup_polygon_points = [&](const std::vector<Vec3>& poly) -> std::vector<Vec3> {
                std::vector<Vec3> out;
                out.reserve(poly.size());
                for (const Vec3& p : poly) {
                    if (!out.empty() && sqdist(out.back(), p) <= cleanup_eps2) {
                        continue;
                    }
                    out.push_back(p);
                }
                if (out.size() > 1 && sqdist(out.front(), out.back()) <= cleanup_eps2) {
                    out.pop_back();
                }
                return out;
            };

            auto compute_pair_voronoi_face = [&](size_t seed_i,
                                                  size_t seed_j,
                                                  std::vector<Point_3>& out_face) -> bool {
                out_face.clear();
                if (seed_i >= num_seeds || seed_j >= num_seeds) return false;
                if (!seed_pos_valid[seed_i] || !seed_pos_valid[seed_j]) return false;

                const Vec3 a = seed_pos[seed_i];
                const Vec3 b = seed_pos[seed_j];
                Vec3 ab = sub(b, a);
                const double ab_len = norm(ab);
                if (ab_len <= 1e-14 * bbox_diag) return false;

                Vec3 n = mul(ab, 1.0 / ab_len);
                Vec3 mid = mul(add(a, b), 0.5);

                // Build an orthonormal basis on the bisector plane.
                Vec3 tmp = (std::abs(n.x) < 0.85) ? Vec3{1.0, 0.0, 0.0} : Vec3{0.0, 1.0, 0.0};
                Vec3 u = normalize(cross(n, tmp));
                if (norm2(u) <= 1e-24) {
                    tmp = Vec3{0.0, 0.0, 1.0};
                    u = normalize(cross(n, tmp));
                }
                Vec3 v = normalize(cross(n, u));

                // Large initial square on Π_ij.  If the true face is unbounded,
                // this produces a bounded diagnostic clipping; for the intended
                // reconstruction pairs the face should normally be bounded by other seeds.
                const double R = initial_half_size;
                std::vector<Vec3> poly;
                poly.reserve(16);
                poly.push_back(add(mid, add(mul(u, -R), mul(v, -R))));
                poly.push_back(add(mid, add(mul(u,  R), mul(v, -R))));
                poly.push_back(add(mid, add(mul(u,  R), mul(v,  R))));
                poly.push_back(add(mid, add(mul(u, -R), mul(v,  R))));

                // Clip by all other reconstruction seeds.
                for (size_t sid_k : reconstruction_seed_ids) {
                    if (sid_k == seed_i || sid_k == seed_j) continue;
                    poly = clip_polygon_by_seed_halfspace(poly, a, seed_pos[sid_k]);
                    if (poly.size() < 3) return false;
                }

                poly = cleanup_polygon_points(poly);
                if (poly.size() < 3) return false;

                out_face.reserve(poly.size());
                for (const Vec3& p : poly) {
                    out_face.emplace_back(p.x, p.y, p.z);
                }
                return true;
            };

            std::vector<std::vector<Point_3>> direct_facets;
            std::vector<std::array<size_t, 8>> facet_sources;
            direct_facets.reserve(all_seed_pairs.size());
            facet_sources.reserve(all_seed_pairs.size());

            size_t strict_pair_total = 0;
            size_t reciprocal_pair_total = 0;
            size_t direct_faces_exported = 0;
            size_t direct_faces_empty = 0;
            size_t invalid_pair_count = 0;

            for (const auto& seed_pair : all_seed_pairs) {
                size_t seed_idx1 = seed_pair.first;
                size_t seed_idx2 = seed_pair.second;
                strict_pair_total++;

                if (seed_idx1 >= num_seeds || seed_idx2 >= num_seeds) {
                    invalid_pair_count++;
                    continue;
                }
                if (!seeds->tree_point_is_active(seed_idx1) || !seeds->tree_point_is_active(seed_idx2)) {
                    invalid_pair_count++;
                    continue;
                }

                size_t* attrib1 = seeds->get_tree_point_attrib(seed_idx1);
                size_t* attrib2 = seeds->get_tree_point_attrib(seed_idx2);

                if (attrib1[5] == attrib2[5]) {
                    invalid_pair_count++;
                    continue;
                }

                bool is_reciprocal_pair =
                    attrib1[1] == seed_idx2 && attrib2[1] == seed_idx1;
                if (!is_reciprocal_pair) {
                    invalid_pair_count++;
                    continue;
                }
                reciprocal_pair_total++;

                std::vector<Point_3> face;
                if (!compute_pair_voronoi_face(seed_idx1, seed_idx2, face)) {
                    direct_faces_empty++;
                    continue;
                }

                direct_facets.push_back(std::move(face));
                facet_sources.push_back({
                    seed_idx1, seed_idx2,
                    attrib1[2], attrib1[3], attrib1[4],
                    attrib2[2], attrib2[3], attrib2[4]
                });
                direct_faces_exported++;
            }

            write_voronoi_facets_original_obj(
                output_filename,
                direct_facets,
                &facet_sources,
                "voronoi_direct_pair_map.csv"
            );

            // Keep the old dedup pipeline available on the directly computed faces.
            write_voronoi_facets_to_obj_dedup(
                "voronoi_dedup.obj",
                direct_facets,
                &facet_sources,
                "voronoi_dedup_map.csv",
                1.05e-5,
                170.0,
                false
            );

            {
                std::ofstream stats_out("vorocrust_surface_stats.csv");
                if (stats_out.is_open()) {
                    stats_out << "mode,total_seed_pairs,reciprocal_strict_pairs,reconstruction_seed_count,"
                              << "direct_faces_exported,direct_faces_empty,invalid_pairs,"
                              << "bbox_diag,initial_half_size,clip_tol,cleanup_eps\n";
                    stats_out << "direct_pair_cell_intersection," 
                              << strict_pair_total << ','
                              << reciprocal_pair_total << ','
                              << reconstruction_seed_ids.size() << ','
                              << direct_faces_exported << ','
                              << direct_faces_empty << ','
                              << invalid_pair_count << ','
                              << bbox_diag << ','
                              << initial_half_size << ','
                              << clip_tol << ','
                              << cleanup_eps << '\n';
                }
            }

            std::cout << "  * Direct strict-pair Voronoi cell intersection exported" << std::endl;
            std::cout << "    reconstruction seeds      = " << reconstruction_seed_ids.size() << std::endl;
            std::cout << "    strict pairs total        = " << strict_pair_total << std::endl;
            std::cout << "    reciprocal strict pairs   = " << reciprocal_pair_total << std::endl;
            std::cout << "    exported polygon faces    = " << direct_faces_exported << std::endl;
            std::cout << "    empty / clipped-out faces = " << direct_faces_empty << std::endl;
            std::cout << "    invalid pairs             = " << invalid_pair_count << std::endl;
            std::cout << "    bbox_diag                 = " << bbox_diag << std::endl;
        };

        auto exportDelaunayBoundarySurface = [&](const std::string& filename) {
            std::ofstream out(filename);
            if (!out) {
                std::cerr << "Error: Cannot create file " << filename << std::endl;
                return;
            }

            out << "# Delaunay-based boundary surface" << std::endl;
            std::cout << "\nExtracting Delaunay boundary surface..." << std::endl;

            std::set<Point> unique_points;
            std::vector<std::array<Point, 3>> triangles;

            for (auto fit = dt.finite_facets_begin(); fit != dt.finite_facets_end(); ++fit) {
                Cell_handle c1 = fit->first;
                int idx = fit->second;
                Cell_handle c2 = c1->neighbor(idx);

                if (dt.is_infinite(c1) || dt.is_infinite(c2)) continue;

                std::set<int> labels;
                for (int i = 0; i < 4; ++i) {
                    Vertex_handle v = c1->vertex(i);
                    auto it = vertex_labels.find(v);
                    if (it != vertex_labels.end()) labels.insert(it->second);
                }
                for (int i = 0; i < 4; ++i) {
                    Vertex_handle v = c2->vertex(i);
                    auto it = vertex_labels.find(v);
                    if (it != vertex_labels.end()) labels.insert(it->second);
                }

                if (labels.size() >= 2) {
                    std::array<Point, 3> triangle;
                    int tri_idx = 0;
                    for (int i = 0; i < 4; ++i) {
                        if (i != idx) triangle[tri_idx++] = c1->vertex(i)->point();
                    }
                    triangles.push_back(triangle);
                    for (const auto& p : triangle) unique_points.insert(p);
                }
            }

            std::map<Point, int> point_to_index;
            int idx = 1;
            for (const auto& p : unique_points) point_to_index[p] = idx++;

            for (const auto& p : unique_points) {
                out << "v " << p.x() << " " << p.y() << " " << p.z() << std::endl;
            }
            for (const auto& tri : triangles) {
                out << "f " << point_to_index[tri[0]]
                    << " " << point_to_index[tri[1]]
                    << " " << point_to_index[tri[2]] << std::endl;
            }
            std::cout << "Delaunay surface saved to: " << filename << std::endl;
        };

        // exportBoundarySurfaceImproved("voronoi_surface.obj");
        exportBoundarySurfaceImproved();
        // exportDelaunayBoundarySurface("delaunay_surface.obj");

    }
}
