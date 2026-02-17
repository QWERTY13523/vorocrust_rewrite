#include "Generator.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <fstream>
#include <random>
#include <iomanip>
#include <limits>
#include <set>
#include <string>
#include <vector>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Delaunay_triangulation_cell_base_3.h> // 蹇呴』鍖呭惈
#define M_PI 3.14159265358979323846
typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_with_info_3<size_t, K> Vb;
typedef CGAL::Delaunay_triangulation_cell_base_3<K> Cb; // 瀹氫箟 Cell Base
typedef CGAL::Triangulation_data_structure_3<Vb, Cb> Tds; // 浼犲叆 Vb 鍜?Cb
typedef CGAL::Delaunay_triangulation_3<K, Tds> Delaunay;
typedef Delaunay::Point Point_3;
typedef Delaunay::Vertex_handle Vertex_handle;
typedef Delaunay::Cell_handle Cell_handle;
typedef Delaunay::Edge Edge;
typedef K::Point_3 Point;

void optimizer(MeshingTree* seeds, MeshingTree* spheres, std::vector<int> face_flat)
{
    Geometry geom;
    size_t num_spheres = spheres->get_num_tree_points();
    size_t num_faces_total = face_flat.size() / 3;

    std::cout << "=== Optimizer: fixing non-adjacent seed pairs ===" << std::endl;

    auto pair_key = [](size_t a, size_t b) -> std::pair<size_t, size_t> {
        return (a < b) ? std::make_pair(a, b) : std::make_pair(b, a);
    };

    auto make_face_key = [](size_t a, size_t b, size_t c) -> std::tuple<size_t, size_t, size_t> {
        size_t arr[3] = {a, b, c};
        if (arr[0] > arr[1]) std::swap(arr[0], arr[1]);
        if (arr[1] > arr[2]) std::swap(arr[1], arr[2]);
        if (arr[0] > arr[1]) std::swap(arr[0], arr[1]);
        return std::make_tuple(arr[0], arr[1], arr[2]);
    };

    // Build sphere-to-faces map: for each sphere, which face indices contain it
    std::vector<std::vector<size_t>> sphere_to_faces(num_spheres);
    for (size_t fi = 0; fi < num_faces_total; fi++) {
        for (int k = 0; k < 3; k++) {
            int sid = face_flat[fi * 3 + k];
            if (sid >= 0 && (size_t)sid < num_spheres)
                sphere_to_faces[sid].push_back(fi);
        }
    }

    // Helper: compute seed pair for a face given 3 sphere indices
    // Returns true if valid; fills seed_a[4], seed_b[4], normal_out[3]
    auto compute_face_seeds = [&](size_t si, size_t sj, size_t sk,
                                   double* seed_a, double* seed_b,
                                   double* normal_out) -> bool {
        double sp_i[4], sp_j[4], sp_k[4];
        spheres->get_tree_point(si, 4, sp_i);
        spheres->get_tree_point(sj, 4, sp_j);
        spheres->get_tree_point(sk, 4, sp_k);

        // Pairwise overlap check (two spheres intersect iff dist < r1+r2)
        if (geom.distance(3, sp_i, sp_j) > sp_i[3] + sp_j[3] + 1e-10) return false;
        if (geom.distance(3, sp_i, sp_k) > sp_i[3] + sp_k[3] + 1e-10) return false;
        if (geom.distance(3, sp_j, sp_k) > sp_j[3] + sp_k[3] + 1e-10) return false;

        // Power vertex (radical center)
        double* centers[3] = {sp_i, sp_j, sp_k};
        double radii[3]    = {sp_i[3], sp_j[3], sp_k[3]};
        double c_ijk[3];
        if (!geom.get_power_vertex(3, 3, centers, radii, c_ijk)) return false;

        double hi = geom.distance(3, c_ijk, sp_i);
        if (hi > sp_i[3] - 1e-10) return false;

        double vi = sqrt(fmax(0.0, sp_i[3] * sp_i[3] - hi * hi));
        if (vi < 1e-5) vi = 1e-5;

        if (geom.get_3d_triangle_area(centers) < 1e-10) return false;
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

    // Helper: compute dot product with sphere normals for region determination
    auto compute_seed_dot = [&](double* seed_pos, size_t si, size_t sj, size_t sk) -> double {
        double dot_sum = 0.0;
        int count = 0;
        size_t sids[3] = {si, sj, sk};
        for (int i = 0; i < 3; i++) {
            if (sids[i] >= num_spheres) continue;
            double center[4];
            spheres->get_tree_point(sids[i], 4, center);
            double* normal = spheres->get_tree_point_normal(sids[i]);
            double nlen_sq = normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2];
            if (nlen_sq < 1e-20) continue;
            dot_sum += (seed_pos[0] - center[0]) * normal[0]
                     + (seed_pos[1] - center[1]) * normal[1]
                     + (seed_pos[2] - center[2]) * normal[2];
            count++;
        }
        return (count > 0) ? (dot_sum / count) : 0.0;
    };

    const int max_iterations = 10;

    for (int iter = 0; iter < max_iterations; iter++) {

        size_t num_seeds = seeds->get_num_tree_points();

        // ---- Build Delaunay triangulation ----
        Delaunay dt;
        std::vector<Vertex_handle> vertex_handles(num_seeds);
        for (size_t i = 0; i < num_seeds; i++) {
            if (!seeds->tree_point_is_active(i)) continue;
            double* pt = seeds->get_tree_point(i);
            Vertex_handle vh = dt.insert(Point_3(pt[0], pt[1], pt[2]));
            vh->info() = i;
            vertex_handles[i] = vh;
        }

        // ---- Collect Delaunay edges ----
        std::set<std::pair<size_t, size_t>> delaunay_edges;
        for (auto eit = dt.finite_edges_begin(); eit != dt.finite_edges_end(); ++eit) {
            Cell_handle c = eit->first;
            Vertex_handle v1 = c->vertex(eit->second);
            Vertex_handle v2 = c->vertex(eit->third);
            if (dt.is_infinite(v1) || dt.is_infinite(v2)) continue;
            size_t s1 = v1->info();
            size_t s2 = v2->info();
            if (s1 >= num_seeds || s2 >= num_seeds || s1 == s2) continue;
            delaunay_edges.insert(pair_key(s1, s2));
        }

        // ---- Find seed pairs whose Voronoi cells are NOT adjacent ----
        std::vector<std::pair<size_t, size_t>> non_adjacent_pairs;
        for (size_t i = 0; i < num_seeds; i++) {
            if (!seeds->tree_point_is_active(i)) continue;
            size_t* attrib = seeds->get_tree_point_attrib(i);
            size_t j = attrib[1];
            if (j <= i || j >= num_seeds) continue;
            if (!seeds->tree_point_is_active(j)) continue;
            size_t* pair_attrib = seeds->get_tree_point_attrib(j);
            if (pair_attrib[1] != i) continue;
            if (delaunay_edges.find(pair_key(i, j)) == delaunay_edges.end()) {
                non_adjacent_pairs.emplace_back(i, j);
            }
        }

        std::cout << "  Iter " << iter << ": " << non_adjacent_pairs.size()
                  << " non-adjacent pairs" << std::endl;
        if (non_adjacent_pairs.empty()) break;

        // ---- Build face-key -> seeds map ----
        std::map<std::tuple<size_t,size_t,size_t>, std::vector<size_t>> face_key_to_seeds;
        for (size_t i = 0; i < num_seeds; i++) {
            if (!seeds->tree_point_is_active(i)) continue;
            size_t* attrib = seeds->get_tree_point_attrib(i);
            auto key = make_face_key(attrib[2], attrib[3], attrib[4]);
            face_key_to_seeds[key].push_back(i);
        }

        std::set<size_t> processed_spheres;
        size_t fixes_this_iter = 0;

        for (size_t pi = 0; pi < non_adjacent_pairs.size(); pi++) {
            size_t seed_i = non_adjacent_pairs[pi].first;
            size_t seed_j = non_adjacent_pairs[pi].second;
            if (!seeds->tree_point_is_active(seed_i) ||
                !seeds->tree_point_is_active(seed_j)) continue;

            size_t* attrib_i = seeds->get_tree_point_attrib(seed_i);
            size_t sph_ids[3] = {attrib_i[2], attrib_i[3], attrib_i[4]};

            // ---- Choose sphere with fewest faces (least impact) ----
            int    best_k   = -1;
            size_t best_sid = SIZE_MAX;
            size_t min_fc   = SIZE_MAX;
            for (int k = 0; k < 3; k++) {
                size_t sid = sph_ids[k];
                if (sid >= num_spheres) continue;
                if (processed_spheres.count(sid)) continue;
                size_t fc = sphere_to_faces[sid].size();
                if (fc < min_fc) { min_fc = fc; best_k = k; best_sid = sid; }
            }
            if (best_k < 0) continue;
            processed_spheres.insert(best_sid);

            double* sph_data = spheres->get_tree_point(best_sid);
            double  orig_r   = sph_data[3];

            // ---- Compute minimum viable radius ----
            // For every face containing this sphere, the sphere must still
            // overlap with the other two spheres (dist < r_new + r_other).
            // => r_new > dist - r_other   (plus small margin)
            double r_min = 0.0;
            for (size_t fi : sphere_to_faces[best_sid]) {
                for (int k = 0; k < 3; k++) {
                    size_t other = (size_t)face_flat[fi * 3 + k];
                    if (other == best_sid) continue;
                    double other_sp[4];
                    spheres->get_tree_point(other, 4, other_sp);
                    double d   = geom.distance(3, sph_data, other_sp);
                    double req = d - other_sp[3] + 1e-6;
                    if (req > r_min) r_min = req;
                }
            }
            if (r_min <= 0.0) r_min = 1e-6;
            if (r_min >= orig_r - 1e-10) continue;  // cannot shrink

            // ---- Collect old seeds from all affected faces ----
            std::set<std::tuple<size_t,size_t,size_t>> affected_keys;
            for (size_t fi : sphere_to_faces[best_sid]) {
                affected_keys.insert(make_face_key(
                    (size_t)face_flat[fi*3],
                    (size_t)face_flat[fi*3+1],
                    (size_t)face_flat[fi*3+2]));
            }

            std::vector<size_t> seeds_to_delete;
            for (auto& fk : affected_keys) {
                auto it = face_key_to_seeds.find(fk);
                if (it == face_key_to_seeds.end()) continue;
                for (size_t s : it->second) {
                    if (seeds->tree_point_is_active(s))
                        seeds_to_delete.push_back(s);
                }
            }

            // ---- Binary search for smallest valid radius ----
            // Constraints checked:
            //   C1. Every affected face still produces valid seeds
            //       (pairwise overlap, power vertex inside all 3 spheres).
            //   C2. New seed positions > 1e-3 from all non-deleted active seeds.
            double lo   = r_min;
            double hi_r = orig_r;
            double valid_r = orig_r;
            bool   found   = false;

            for (int bs = 0; bs < 40; bs++) {
                double mid = (lo + hi_r) / 2.0;
                sph_data[3] = mid;

                bool all_ok = true;

                // -- C1: face validity --
                for (size_t fi : sphere_to_faces[best_sid]) {
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

                // -- C2: distance check --
                if (all_ok) {
                    for (size_t fi : sphere_to_faces[best_sid]) {
                        double sa[4], sb[4], nm[3];
                        compute_face_seeds(
                            (size_t)face_flat[fi*3],
                            (size_t)face_flat[fi*3+1],
                            (size_t)face_flat[fi*3+2],
                            sa, sb, nm);

                        // Constraint: distance between new seed pair > 1e-4
                        if (geom.distance(3, sa, sb) <= 1e-4) { all_ok = false; break; }

                        size_t ic; double hc = DBL_MAX;
                        seeds->get_closest_tree_point(
                            sa, seeds_to_delete.size(),
                            seeds_to_delete.data(), ic, hc);
                        if (hc < 1e-3) { all_ok = false; break; }

                        hc = DBL_MAX;
                        seeds->get_closest_tree_point(
                            sb, seeds_to_delete.size(),
                            seeds_to_delete.data(), ic, hc);
                        if (hc < 1e-3) { all_ok = false; break; }
                    }
                }

                if (all_ok) {
                    valid_r = mid;
                    found   = true;
                    hi_r    = mid;   // try even smaller
                } else {
                    lo = mid;        // too small, increase
                }
            }

            sph_data[3] = orig_r;  // restore temporarily
            if (!found) continue;

            // ---- Apply the shrink ----
            sph_data[3] = valid_r;
            std::cout << "    Sphere " << best_sid
                      << ": radius " << orig_r << " -> " << valid_r << std::endl;

            // ---- Delete old seeds from affected faces ----
            for (size_t s : seeds_to_delete) {
                seeds->lazy_delete_tree_point(s);
            }

            // ---- Recompute and add new seeds for each affected face ----
            for (size_t fi : sphere_to_faces[best_sid]) {
                size_t s0 = (size_t)face_flat[fi*3];
                size_t s1 = (size_t)face_flat[fi*3+1];
                size_t s2 = (size_t)face_flat[fi*3+2];

                double seed_a[4], seed_b[4], normal[3];
                if (!compute_face_seeds(s0, s1, s2, seed_a, seed_b, normal))
                    continue;

                // Determine inside/outside region for each seed using dot product comparison
                // Larger dot product -> Outside (0)
                // Smaller dot product -> Inside (1)
                double dot_a = compute_seed_dot(seed_a, s0, s1, s2);
                double dot_b = compute_seed_dot(seed_b, s0, s1, s2);

                size_t region_a = (dot_a >= dot_b) ? 0 : 1;
                size_t region_b = (dot_a >= dot_b) ? 1 : 0;

                size_t idx_a = seeds->get_num_tree_points();
                size_t idx_b = idx_a + 1;

                // attrib layout: [0]=count, [1]=pair, [2..4]=spheres, [5]=region
                size_t att_a[6] = {6, idx_b, s0, s1, s2, region_a};
                size_t att_b[6] = {6, idx_a, s0, s1, s2, region_b};

                // normal stored with num_dim=4 elements
                double normal_a[4] = { normal[0],  normal[1],  normal[2], 0.0};
                double normal_b[4] = {-normal[0], -normal[1], -normal[2], 0.0};

                seeds->add_tree_point(4, seed_a, normal_a, att_a);
                seeds->add_tree_point(4, seed_b, normal_b, att_b);
            }

            fixes_this_iter++;
        }

        std::cout << "  Applied " << fixes_this_iter << " sphere shrinks" << std::endl;
        if (fixes_this_iter == 0) break;
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

    auto pair_key = [](size_t a, size_t b) {
        return (a < b) ? std::make_pair(a, b) : std::make_pair(b, a);
    };

    auto build_delaunay = [&]() {
        dt.clear();
        std::fill(vertex_handles.begin(), vertex_handles.end(), Vertex_handle());
        for (size_t i = 0; i < num_seeds; i++) {
            if (!seeds->tree_point_is_active(i)) continue;
            double* pt = seeds->get_tree_point(i);
            Vertex_handle vh = dt.insert(Point_3(pt[0], pt[1], pt[2]));
            vh->info() = i;
            vertex_handles[i] = vh;
        }
    };

    // auto calc_bbox_diag = [&]() {
    //     double min_x = std::numeric_limits<double>::max();
    //     double min_y = std::numeric_limits<double>::max();
    //     double min_z = std::numeric_limits<double>::max();
    //     double max_x = -std::numeric_limits<double>::max();
    //     double max_y = -std::numeric_limits<double>::max();
    //     double max_z = -std::numeric_limits<double>::max();
    //     bool has_active_seed = false;

    //     for (size_t i = 0; i < num_seeds; i++) {
    //         if (!seeds->tree_point_is_active(i)) continue;
    //         double* pt = seeds->get_tree_point(i);
    //         min_x = std::min(min_x, pt[0]);
    //         min_y = std::min(min_y, pt[1]);
    //         min_z = std::min(min_z, pt[2]);
    //         max_x = std::max(max_x, pt[0]);
    //         max_y = std::max(max_y, pt[1]);
    //         max_z = std::max(max_z, pt[2]);
    //         has_active_seed = true;
    //     }

    //     if (!has_active_seed) return 1.0;
    //     double dx = max_x - min_x;
    //     double dy = max_y - min_y;
    //     double dz = max_z - min_z;
    //     return std::sqrt(dx * dx + dy * dy + dz * dz);
    // };

    // // Enforce adjacency for paired seeds by shrinking non-adjacent pair distance iteratively.
    // const int max_pair_adjust_iters = 20;
    // const double bbox_diag = calc_bbox_diag();
    // const double move_step_start = std::max(6e-4, bbox_diag * 6e-4);
    // const double move_step_end = std::max(2e-4, bbox_diag * 2e-4);
    // size_t last_non_adjacent_count = 0;

    // for (int iter = 0; iter < max_pair_adjust_iters; ++iter) {
    //     build_delaunay();

    //     // Learning-rate style schedule: larger step early, smaller step later.
    //     double t = (max_pair_adjust_iters > 1)
    //         ? static_cast<double>(iter) / static_cast<double>(max_pair_adjust_iters - 1)
    //         : 1.0;
    //     double move_step = move_step_start + (move_step_end - move_step_start) * t;

    //     std::set<std::pair<size_t, size_t>> delaunay_edges;
    //     for (auto eit = dt.finite_edges_begin(); eit != dt.finite_edges_end(); ++eit) {
    //         Cell_handle c = eit->first;
    //         int i1 = eit->second;
    //         int i2 = eit->third;
    //         Vertex_handle v1 = c->vertex(i1);
    //         Vertex_handle v2 = c->vertex(i2);
    //         if (dt.is_infinite(v1) || dt.is_infinite(v2)) continue;

    //         size_t s1 = v1->info();
    //         size_t s2 = v2->info();
    //         if (s1 >= num_seeds || s2 >= num_seeds || s1 == s2) continue;
    //         delaunay_edges.insert(pair_key(s1, s2));
    //     }

    //     std::vector<std::pair<size_t, size_t>> non_adjacent_pairs;
    //     for (size_t i = 0; i < num_seeds; i++) {
    //         if (!seeds->tree_point_is_active(i)) continue;
    //         size_t* attrib = seeds->get_tree_point_attrib(i);
    //         size_t j = attrib[1];
    //         if (j <= i || j >= num_seeds) continue;
    //         if (!seeds->tree_point_is_active(j)) continue;
    //         size_t* pair_attrib = seeds->get_tree_point_attrib(j);
    //         if (pair_attrib[1] != i) continue;

    //         if (delaunay_edges.find(pair_key(i, j)) == delaunay_edges.end()) {
    //             non_adjacent_pairs.emplace_back(i, j);
    //         }
    //     }

    //     std::cout << "iter : " << iter << " " << last_non_adjacent_count << "\n";
    //     last_non_adjacent_count = non_adjacent_pairs.size();
    //     if (last_non_adjacent_count == 0) {
    //         if (iter > 0) {
    //             std::cout << "  * Seed-pair adjacency converged after " << iter
    //                       << " iterations" << std::endl;
    //         }
    //         break;
    //     }

    //     for (const auto& pr : non_adjacent_pairs) {
    //         double* p1 = seeds->get_tree_point(pr.first);
    //         double* p2 = seeds->get_tree_point(pr.second);

    //         double dx = p2[0] - p1[0];
    //         double dy = p2[1] - p1[1];
    //         double dz = p2[2] - p1[2];
    //         double dist = std::sqrt(dx * dx + dy * dy + dz * dz);
    //         if (dist < 1e-12) continue;

    //         double move_len = std::min(move_step, 0.45 * dist);
    //         if (move_len <= 0.0) continue;

    //         double inv_dist = 1.0 / dist;
    //         p1[0] += dx * inv_dist * move_len;
    //         p1[1] += dy * inv_dist * move_len;
    //         p1[2] += dz * inv_dist * move_len;
    //         p2[0] -= dx * inv_dist * move_len;
    //         p2[1] -= dy * inv_dist * move_len;
    //         p2[2] -= dz * inv_dist * move_len;
    //     }
    // }

    // Rebuild once to ensure dt/vertex_handles match final adjusted seed positions.
    build_delaunay();
    // if (last_non_adjacent_count > 0) {
    //     std::cout << "  * Warning: still " << last_non_adjacent_count
    //               << " non-adjacent seed pairs after max iterations" << std::endl;
    // }

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

                // 鏍规嵁鍖哄煙ID鍒ゆ柇鍐呭
                // attrib[5] 鏄尯鍩烮D: 0 = 澶栭儴, 闈? = 鍐呴儴
                if (attrib[5] == 1) {
                    // 澶栭儴绉嶅瓙鐐?
                    write_seed_with_spheres(outside_seeds, outside_vertex_offset, pt, attrib);
                    outside_count++;
                }
                else {
                    // 鍐呴儴绉嶅瓙鐐?
                    write_seed_with_spheres(inside_seeds, inside_vertex_offset, pt, attrib);
                    inside_count++;
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
                }
            }
        }
        
        pairs_out.close();
        std::cout << "  * Seed pairs: " << pair_count << " connections (saved to seed_pairs.obj)" << std::endl;
    } else {
        std::cerr << "Error: Cannot write to seed_pairs.obj" << std::endl;
    }
}
    // 2. Collect Voronoi facets for inside/outside seed pairs
    std::vector<std::vector<Point_3>> voronoi_facets;
    
    for (auto eit = dt.finite_edges_begin(); eit != dt.finite_edges_end(); ++eit) {
        Cell_handle c = eit->first;
        int i1 = eit->second;
        int i2 = eit->third;
        
        Vertex_handle v1 = c->vertex(i1);
        Vertex_handle v2 = c->vertex(i2);
        
        size_t seed_idx1 = v1->info();
        size_t seed_idx2 = v2->info();
        
        // Check if this edge connects a seed pair (inside/outside)
        size_t* attrib1 = seeds->get_tree_point_attrib(seed_idx1);
        size_t* attrib2 = seeds->get_tree_point_attrib(seed_idx2);
        
        // attrib[1] is the pair seed index
        // 杩欓噷浣跨敤涔嬪墠鐢熸垚鐨勯厤瀵逛俊鎭紝鎴栬€呭彲浠ユ敼鐢ㄥ尯鍩?ID (attrib[5]) 鍒ゆ柇: 
        //bool is_interface = (attrib1[5] != attrib2[5]);
        bool is_interface = (attrib1[5] != attrib2[5]) && attrib1[1] == seed_idx2 && attrib2[1] == seed_idx1;
        if (!is_interface) continue;
      //  if (!is_pair) continue;
        int label1 = seed_idx1;
        int label2 = seed_idx2;
        // This edge connects a seed pair - extract the dual Voronoi facet
        // The Voronoi facet is the convex polygon formed by circumcenters of 
        // cells incident to this edge
        std::vector<Point_3> facet_vertices;
        
        Delaunay::Cell_circulator cc = dt.incident_cells(*eit);
        Delaunay::Cell_circulator done = cc;
        
        if (cc == nullptr) continue;
        
        do {
            if (!dt.is_infinite(cc)) {
                // Compute circumcenter manually or use CGAL dual
                // Using dual() is safer and usually cached/optimized in some kernels, 
                // but manual construction as you did is also fine for Exact kernel.
                Point_3 center = dt.dual(cc); 
                facet_vertices.push_back(center);
            }
            ++cc;
        } while (cc != done);
        
        if (facet_vertices.size() >= 3) {
            voronoi_facets.push_back(facet_vertices);
        }
    }
    
    std::cout << "  * Found " << voronoi_facets.size() << " Voronoi facets for seed pairs" << std::endl;

    // 3. Write Surface Mesh OBJ file
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

    {
        std::map<Vertex_handle, int> vertex_labels;
        for (size_t i = 0; i < num_seeds; i++) {
            if (!seeds->tree_point_is_active(i)) continue;
            Vertex_handle vh = vertex_handles[i];
            if (vh == Vertex_handle()) continue;
            size_t* attrib = seeds->get_tree_point_attrib(i);
            vertex_labels[vh] = static_cast<int>(attrib[5]);
        }

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
            const std::vector<std::vector<Point_3>>& voronoi_facets
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

            auto get_index = [&](const Point_3& p) {
                auto it = vertex_index.find(p);
                if (it != vertex_index.end()) return it->second;
                int idx = static_cast<int>(vertices.size()) + 1;
                vertex_index[p] = idx;
                vertices.push_back(p);
                return idx;
            };

            for (const auto& facet : voronoi_facets) {
                if (facet.size() < 3) continue;
                int v0 = get_index(facet[0]);
                for (size_t i = 1; i + 1 < facet.size(); ++i) {
                    int v1 = get_index(facet[i]);
                    int v2 = get_index(facet[i + 1]);
                    triangles.push_back({ v0, v1, v2 });
                }
            }

            for (const auto& p : vertices) {
                out << "v " << p.x() << " " << p.y() << " " << p.z() << "\n";
            }
            for (const auto& t : triangles) {
                out << "f " << t[0] << " " << t[1] << " " << t[2] << "\n";
            }
        };

        auto write_voronoi_facets_to_obj_dedup = [&](
            const std::string& filename,
            const std::vector<std::vector<Point_3>>& voronoi_facets,
            double epsilon = 5e-2,
            double angle_threshold_deg = 172.0,
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
            for (const auto& facet : voronoi_facets) {
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

        auto exportBoundarySurfaceImproved = [&](const std::string& filename) {
            std::ofstream out(filename);
            if (!out) {
                std::cerr << "Error: Cannot create file " << filename << std::endl;
                return;
            }

            out << "# Voronoi boundary surface (improved with validation)" << std::endl;
            std::cout << "\n=== Extracting Boundary Surface (Improved) ===" << std::endl;

            int total_edges = 0;
            int boundary_edges_count = 0;
            int inner_edges = 0;
            int outer_edges = 0;
            int no_label_edges = 0;

            std::vector<std::vector<Point_3>> facets;
            bool did_export_single = false;

            for (auto eit = dt.finite_edges_begin(); eit != dt.finite_edges_end(); ++eit) {
                total_edges++;

                Cell_handle c = eit->first;
                int i = eit->second;
                int j = eit->third;

                Vertex_handle v1 = c->vertex(i);
                Vertex_handle v2 = c->vertex(j);

                auto it1 = vertex_labels.find(v1);
                auto it2 = vertex_labels.find(v2);
                if (it1 == vertex_labels.end() || it2 == vertex_labels.end()) {
                    no_label_edges++;
                    continue;
                }

                int label1 = it1->second;
                int label2 = it2->second;

                if (label1 != label2) {
                    std::vector<Point_3> facet_vertices;
                    Delaunay::Cell_circulator cc = dt.incident_cells(*eit);
                    Delaunay::Cell_circulator done = cc;
                    if (cc == nullptr) continue;

                    do {
                        if (!dt.is_infinite(cc)) {
                            Point_3 center = dt.dual(cc);
                            facet_vertices.push_back(center);
                        }
                        ++cc;
                    } while (cc != done);

                    if (facet_vertices.size() >= 3) {
                        facets.push_back(facet_vertices);
                        // if (!did_export_single) {
                        //     export_single_voronoi_polygon(v1, v2, "single_voronoi_polygon.obj");
                        //     did_export_single = true;
                        // }
                    }

                    boundary_edges_count++;
                }
                else {
                    if (label1 == 1) inner_edges++;
                    else if (label1 == 2) outer_edges++;
                }
            }

            write_voronoi_facets_to_obj("voronoi.obj", facets);
            write_voronoi_facets_to_obj_dedup("voronoi_dedup.obj", facets);
            write_voronoi_facets_triangulated_obj("voronoi_dedup_triangulated.obj", facets);
            std::cout << "  * Found " << facets.size() << " Voronoi facets for seed pairs" << std::endl;

            {
                std::map<std::tuple<double, double, double>, size_t> vertex_map;
                std::vector<Point_3> unique_vertices;
                std::vector<std::vector<size_t>> face_indices;

                auto get_vertex_index = [&](const Point_3& p) -> size_t {
                    double x = std::round(CGAL::to_double(p.x()) * 1e10) / 1e10;
                    double y = std::round(CGAL::to_double(p.y()) * 1e10) / 1e10;
                    double z = std::round(CGAL::to_double(p.z()) * 1e10) / 1e10;
                    auto key = std::make_tuple(x, y, z);

                    auto it = vertex_map.find(key);
                    if (it != vertex_map.end()) return it->second;
                    size_t idx = unique_vertices.size();
                    vertex_map[key] = idx;
                    unique_vertices.push_back(p);
                    return idx;
                };

                for (const auto& facet : facets) {
                    if (facet.size() < 3) continue;

                    std::vector<size_t> indices;
                    for (const auto& pt : facet) indices.push_back(get_vertex_index(pt));

                    Point_3 centroid = CGAL::ORIGIN;
                    for (const auto& pt : facet) centroid = centroid + (pt - CGAL::ORIGIN);
                    double s = static_cast<double>(facet.size());
                    centroid = Point_3(centroid.x() / s, centroid.y() / s, centroid.z() / s);
                    size_t centroid_idx = get_vertex_index(centroid);

                    for (size_t i = 0; i < indices.size(); i++) {
                        size_t idx0 = indices[i];
                        size_t idx1 = indices[(i + 1) % indices.size()];
                        if (idx0 == idx1 || idx0 == centroid_idx || idx1 == centroid_idx) continue;
                        face_indices.push_back({ centroid_idx, idx0, idx1 });
                    }
                }

                for (const auto& v : unique_vertices) {
                    out << "v " << std::setprecision(16)
                        << CGAL::to_double(v.x()) << " "
                        << CGAL::to_double(v.y()) << " "
                        << CGAL::to_double(v.z()) << std::endl;
                }
                for (const auto& face : face_indices) {
                    out << "f " << (face[0] + 1) << " " << (face[1] + 1) << " " << (face[2] + 1) << std::endl;
                }

                std::cout << "  * Total vertices: " << unique_vertices.size() << ", triangles: " << face_indices.size() << std::endl;
            }

            std::cout << "\n--- Edge Classification ---" << std::endl;
            std::cout << "Total Delaunay edges: " << total_edges << std::endl;
            std::cout << "|- Boundary edges (1-2): " << boundary_edges_count
                << " (" << (100.0 * boundary_edges_count / std::max(1, total_edges)) << "%)" << std::endl;
            std::cout << "|- Inner edges (1-1): " << inner_edges
                << " (" << (100.0 * inner_edges / std::max(1, total_edges)) << "%)" << std::endl;
            std::cout << "|- Outer edges (2-2): " << outer_edges
                << " (" << (100.0 * outer_edges / std::max(1, total_edges)) << "%)" << std::endl;
            std::cout << "'- No label edges: " << no_label_edges << std::endl;
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

        exportBoundarySurfaceImproved("voronoi_surface.obj");
        exportDelaunayBoundarySurface("delaunay_surface.obj");

        // // ====== Detect protrusions and export nearby Voronoi diagram ======
        // {
        //     std::cout << "\n=== Detecting Protrusions and Exporting Nearby Voronoi ===" << std::endl;

        //     std::set<size_t> protrusion_seeds;
        //     double protrusion_threshold = 3.0; // flag if Voronoi vertex > N * pair_distance from midpoint

        //     for (auto eit = dt.finite_edges_begin(); eit != dt.finite_edges_end(); ++eit) {
        //         Cell_handle c = eit->first;
        //         Vertex_handle v1 = c->vertex(eit->second);
        //         Vertex_handle v2 = c->vertex(eit->third);

        //         size_t idx1 = v1->info();
        //         size_t idx2 = v2->info();

        //         if (!seeds->tree_point_is_active(idx1) || !seeds->tree_point_is_active(idx2)) continue;

        //         size_t* attrib1 = seeds->get_tree_point_attrib(idx1);
        //         size_t* attrib2 = seeds->get_tree_point_attrib(idx2);
        //         bool is_interface = (attrib1[5] != attrib2[5]) && attrib1[1] == idx2 && attrib2[1] == idx1;
        //         if (!is_interface) continue;

        //         double* pt1 = seeds->get_tree_point(idx1);
        //         double* pt2 = seeds->get_tree_point(idx2);
        //         double mid[3] = {(pt1[0]+pt2[0])*0.5, (pt1[1]+pt2[1])*0.5, (pt1[2]+pt2[2])*0.5};
        //         double pair_dist = std::sqrt((pt1[0]-pt2[0])*(pt1[0]-pt2[0]) +
        //                                      (pt1[1]-pt2[1])*(pt1[1]-pt2[1]) +
        //                                      (pt1[2]-pt2[2])*(pt1[2]-pt2[2]));
        //         if (pair_dist < 1e-15) continue;

        //         Delaunay::Cell_circulator cc = dt.incident_cells(*eit);
        //         Delaunay::Cell_circulator done = cc;
        //         if (cc == nullptr) continue;

        //         bool is_protrusion = false;
        //         do {
        //             if (!dt.is_infinite(cc)) {
        //                 Point_3 center = dt.dual(cc);
        //                 double cx = CGAL::to_double(center.x());
        //                 double cy = CGAL::to_double(center.y());
        //                 double cz = CGAL::to_double(center.z());
        //                 double dist = std::sqrt((cx-mid[0])*(cx-mid[0]) +
        //                                         (cy-mid[1])*(cy-mid[1]) +
        //                                         (cz-mid[2])*(cz-mid[2]));
        //                 if (dist > protrusion_threshold * pair_dist) {
        //                     is_protrusion = true;
        //                     break;
        //                 }
        //             }
        //             ++cc;
        //         } while (cc != done);

        //         if (is_protrusion) {
        //             protrusion_seeds.insert(idx1);
        //             protrusion_seeds.insert(idx2);
        //         }
        //     }

        //     if (protrusion_seeds.empty()) {
        //         std::cout << "  * No protrusions detected." << std::endl;
        //     } else {
        //         std::cout << "  * Found " << protrusion_seeds.size()
        //                   << " seeds near protrusions." << std::endl;

        //         // Expand to Delaunay neighbors for a fuller picture
        //         std::set<size_t> expanded_seeds = protrusion_seeds;
        //         for (auto vit = dt.finite_vertices_begin(); vit != dt.finite_vertices_end(); ++vit) {
        //             if (protrusion_seeds.count(vit->info())) {
        //                 std::vector<Vertex_handle> neighbors;
        //                 dt.adjacent_vertices(vit, std::back_inserter(neighbors));
        //                 for (auto& nv : neighbors) {
        //                     if (!dt.is_infinite(nv)) {
        //                         expanded_seeds.insert(nv->info());
        //                     }
        //                 }
        //             }
        //         }
        //         std::cout << "  * Expanded to " << expanded_seeds.size()
        //                   << " seeds (with Delaunay neighbors)." << std::endl;

        //         // Collect ALL Voronoi facets where at least one endpoint is in expanded_seeds
        //         std::vector<std::vector<Point_3>> all_nearby_facets;
        //         std::vector<std::vector<Point_3>> interface_nearby_facets;

        //         for (auto eit = dt.finite_edges_begin(); eit != dt.finite_edges_end(); ++eit) {
        //             Cell_handle c = eit->first;
        //             Vertex_handle v1 = c->vertex(eit->second);
        //             Vertex_handle v2 = c->vertex(eit->third);
        //             size_t idx1 = v1->info();
        //             size_t idx2 = v2->info();

        //             if (!expanded_seeds.count(idx1) && !expanded_seeds.count(idx2)) continue;

        //             std::vector<Point_3> facet_verts;
        //             Delaunay::Cell_circulator cc = dt.incident_cells(*eit);
        //             Delaunay::Cell_circulator done = cc;
        //             if (cc == nullptr) continue;

        //             do {
        //                 if (!dt.is_infinite(cc)) {
        //                     facet_verts.push_back(dt.dual(cc));
        //                 }
        //                 ++cc;
        //             } while (cc != done);

        //             if (facet_verts.size() >= 3) {
        //                 all_nearby_facets.push_back(facet_verts);

        //                 if (seeds->tree_point_is_active(idx1) && seeds->tree_point_is_active(idx2)) {
        //                     size_t* a1 = seeds->get_tree_point_attrib(idx1);
        //                     size_t* a2 = seeds->get_tree_point_attrib(idx2);
        //                     bool iface = (a1[5] != a2[5]) && a1[1] == idx2 && a2[1] == idx1;
        //                     if (iface) interface_nearby_facets.push_back(facet_verts);
        //                 }
        //             }
        //         }

        //         write_voronoi_facets_to_obj("protrusion_voronoi_all.obj", all_nearby_facets);
        //         write_voronoi_facets_to_obj("protrusion_voronoi_interface.obj", interface_nearby_facets);

        //         // Export seed points near protrusions (colored by region)
        //         {
        //             std::ofstream psf("protrusion_seeds.obj");
        //             psf << std::fixed << std::setprecision(16);
        //             psf << "# Seeds near protrusions (red=region1, blue=region0)\n";
        //             for (size_t sidx : expanded_seeds) {
        //                 if (!seeds->tree_point_is_active(sidx)) continue;
        //                 double* pt = seeds->get_tree_point(sidx);
        //                 size_t* attr = seeds->get_tree_point_attrib(sidx);
        //                 if (attr[5] == 1) {
        //                     psf << "v " << pt[0] << " " << pt[1] << " " << pt[2]
        //                         << " 1.0 0.0 0.0\n";
        //                 } else {
        //                     psf << "v " << pt[0] << " " << pt[1] << " " << pt[2]
        //                         << " 0.0 0.0 1.0\n";
        //                 }
        //             }
        //             psf.close();
        //         }

        //         std::cout << "  * Exported " << all_nearby_facets.size()
        //                   << " total Voronoi facets near protrusions -> protrusion_voronoi_all.obj" << std::endl;
        //         std::cout << "  * Exported " << interface_nearby_facets.size()
        //                   << " interface facets near protrusions -> protrusion_voronoi_interface.obj" << std::endl;
        //         std::cout << "  * Exported seed points -> protrusion_seeds.obj" << std::endl;
        //     }
        // }
    }
}


struct LabeledPoint {
    Point point;
    int label;

    LabeledPoint(double x, double y, double z, int l)
        : point(x, y, z), label(l) {
    }
};


// Get Voronoi face vertices for a Delaunay edge
std::vector<Point> getVoronoiFace(const Delaunay& dt,
    Vertex_handle v1,
    Vertex_handle v2) {
    std::vector<Point> face_vertices;

    // Get all cells incident to both vertices
    std::vector<Cell_handle> cells;
    dt.incident_cells(v1, std::back_inserter(cells));

    std::vector<Cell_handle> common_cells;
    for (const auto& cell : cells) {
        if (!dt.is_infinite(cell)) {
            for (int i = 0; i < 4; ++i) {
                if (cell->vertex(i) == v2) {
                    common_cells.push_back(cell);
                    break;
                }
            }
        }
    }

    // Get circumcenters of common cells (Voronoi vertices)
    for (const auto& cell : common_cells) {
        face_vertices.push_back(dt.dual(cell));
    }

    return face_vertices;
}

// Export boundary surface as triangulated mesh
void exportBoundarySurface(const Delaunay& dt,
    const std::map<Vertex_handle, int>& vertex_labels,
    const std::string& filename) {
    std::ofstream out(filename);
    if (!out) {
        std::cerr << "Error: Cannot create file " << filename << std::endl;
        return;
    }

    out << "# Voronoi boundary surface mesh" << std::endl;
    out << "# Triangulated surface between different labels" << std::endl;

    std::cout << "\nExtracting boundary surface..." << std::endl;

    // Find all edges between vertices with different labels
    std::set<std::pair<Vertex_handle, Vertex_handle>> boundary_edges;

    for (auto eit = dt.finite_edges_begin(); eit != dt.finite_edges_end(); ++eit) {
        Cell_handle c = eit->first;
        int i = eit->second;
        int j = eit->third;

        Vertex_handle v1 = c->vertex(i);
        Vertex_handle v2 = c->vertex(j);

        if (vertex_labels.find(v1) != vertex_labels.end() &&
            vertex_labels.find(v2) != vertex_labels.end()) {
            if (vertex_labels.at(v1) != vertex_labels.at(v2)) {
                // Normalize edge (smaller pointer first)
                if (v1 < v2) {
                    boundary_edges.insert({ v1, v2 });
                }
                else {
                    boundary_edges.insert({ v2, v1 });
                }
            }
        }
    }

    std::cout << "Found " << boundary_edges.size() << " boundary edges" << std::endl;

    // Collect all Voronoi vertices and faces
    std::map<Point, int> point_indices;
    int vertex_count = 0;
    std::vector<std::vector<int>> faces;

    for (const auto& edge_pair : boundary_edges) {
        Vertex_handle v1 = edge_pair.first;
        Vertex_handle v2 = edge_pair.second;

        std::vector<Point> face_points = getVoronoiFace(dt, v1, v2);

        if (face_points.size() >= 3) {
            std::vector<int> face_indices;
            for (const auto& p : face_points) {
                if (point_indices.find(p) == point_indices.end()) {
                    point_indices[p] = ++vertex_count;
                }
                face_indices.push_back(point_indices[p]);
            }
            faces.push_back(face_indices);
        }
    }

    // Write vertices
    std::vector<Point> ordered_points(vertex_count);
    for (const auto& pair : point_indices) {
        ordered_points[pair.second - 1] = pair.first;
    }

    for (const auto& p : ordered_points) {
        out << "v " << p.x() << " " << p.y() << " " << p.z() << std::endl;
    }

    // Triangulate and write faces
    int triangle_count = 0;
    for (const auto& face : faces) {
        if (face.size() >= 3) {
            // Simple fan triangulation from first vertex
            for (size_t i = 1; i < face.size() - 1; ++i) {
                out << "f " << face[0] << " " << face[i] << " " << face[i + 1] << std::endl;
                triangle_count++;
            }
        }
    }

    out.close();

    std::cout << "Surface mesh exported:" << std::endl;
    std::cout << "  Vertices: " << vertex_count << std::endl;
    std::cout << "  Voronoi faces: " << faces.size() << std::endl;
    std::cout << "  Triangles: " << triangle_count << std::endl;
    std::cout << "  Saved to: " << filename << std::endl;
}

// Export Delaunay facets as triangulated surface
void exportDelaunayBoundarySurface(const Delaunay& dt,
    const std::map<Vertex_handle, int>& vertex_labels,
    const std::string& filename) {
    std::ofstream out(filename);
    if (!out) {
        std::cerr << "Error: Cannot create file " << filename << std::endl;
        return;
    }

    out << "# Delaunay-based boundary surface" << std::endl;

    std::cout << "\nExtracting Delaunay boundary surface..." << std::endl;

    // Find boundary facets (Delaunay triangles with different labels on each side)
    std::set<Point> unique_points;
    std::vector<std::array<Point, 3>> triangles;

    for (auto fit = dt.finite_facets_begin(); fit != dt.finite_facets_end(); ++fit) {
        Cell_handle c1 = fit->first;
        int idx = fit->second;
        Cell_handle c2 = c1->neighbor(idx);

        if (dt.is_infinite(c1) || dt.is_infinite(c2)) continue;

        // Check if this facet separates different labels
        std::set<int> labels;
        for (int i = 0; i < 4; ++i) {
            Vertex_handle v = c1->vertex(i);
            if (vertex_labels.find(v) != vertex_labels.end()) {
                labels.insert(vertex_labels.at(v));
            }
        }
        for (int i = 0; i < 4; ++i) {
            Vertex_handle v = c2->vertex(i);
            if (vertex_labels.find(v) != vertex_labels.end()) {
                labels.insert(vertex_labels.at(v));
            }
        }

        if (labels.size() >= 2) {
            // This is a boundary facet - get its 3 vertices
            std::array<Point, 3> triangle;
            int tri_idx = 0;
            for (int i = 0; i < 4; ++i) {
                if (i != idx) {
                    triangle[tri_idx++] = c1->vertex(i)->point();
                }
            }

            triangles.push_back(triangle);
            for (const auto& p : triangle) {
                unique_points.insert(p);
            }
        }
    }

    std::cout << "Found " << triangles.size() << " boundary triangles" << std::endl;
    std::cout << "With " << unique_points.size() << " unique vertices" << std::endl;

    // Map points to indices
    std::map<Point, int> point_to_index;
    int idx = 1;
    for (const auto& p : unique_points) {
        point_to_index[p] = idx++;
    }

    // Write vertices
    for (const auto& p : unique_points) {
        out << "v " << p.x() << " " << p.y() << " " << p.z() << std::endl;
    }

    // Write faces
    for (const auto& tri : triangles) {
        out << "f " << point_to_index[tri[0]]
            << " " << point_to_index[tri[1]]
            << " " << point_to_index[tri[2]] << std::endl;
    }

    out.close();
    std::cout << "Delaunay surface saved to: " << filename << std::endl;
}
// 鏀硅繘鐗堬細甯﹁缁嗙粺璁″拰楠岃瘉鐨勮竟鐣岄潰鎻愬彇
// void export_single_voronoi_polygon(
//     const Delaunay& dt,
//     Vertex_handle v1,
//     Vertex_handle v2,
//     const std::string& filename)
// {
//     std::cout << "=== Exporting single Voronoi polygon ===" << std::endl;
//     std::cout << "  * Edge: seed " << v1->info() << " <-> seed " << v2->info() << std::endl;

//     std::ofstream obj_file(filename);
//     if (!obj_file.is_open()) {
//         std::cerr << "Error: Cannot open " << filename << std::endl;
//         return;
//     }

//     obj_file << std::fixed << std::setprecision(16);
//     obj_file << "# Voronoi polygon for Delaunay edge\n";
//     obj_file << "# Seed 1 (index " << v1->info() << ")\n";
//     obj_file << "# Seed 2 (index " << v2->info() << ")\n\n";

//     // 鏀堕泦鍥寸粫璇ヨ竟鐨勬墍鏈夊洓闈綋鐨勫鎺ョ悆蹇?
//     std::vector<Point_3> polygon_vertices;

//     // 闇€瑕佷粠杈规瀯閫爀dge
//     // 鍦–GAL涓紝edge鏄€氳繃(cell, i, j)涓夊厓缁勮〃绀虹殑
//     // 杩欓噷鎴戜滑鐢ㄥ彟涓€绉嶆柟娉曪細鎵惧埌鍖呭惈v1鍜寁2鐨勬墍鏈塩ell

//     std::vector<Cell_handle> cells_v1;
//     dt.incident_cells(v1, std::back_inserter(cells_v1));

//     for (auto cell : cells_v1) {
//         if (dt.is_infinite(cell)) continue;

//         // 妫€鏌ヨ繖涓猚ell鏄惁涔熷寘鍚玽2
//         bool contains_v2 = false;
//         for (int i = 0; i < 4; i++) {
//             if (cell->vertex(i) == v2) {
//                 contains_v2 = true;
//                 break;
//             }
//         }

//         if (contains_v2) {
//             // 杩欎釜cell鍚屾椂鍖呭惈v1鍜寁2锛屽叾澶栨帴鐞冨績鏄杈瑰舰鐨勪竴涓《鐐?
//             Point_3 center = dt.dual(cell);
//             polygon_vertices.push_back(center);
//         }
//     }

//     std::cout << "  * Polygon has " << polygon_vertices.size() << " vertices" << std::endl;

//     if (polygon_vertices.size() < 3) {
//         std::cout << "  * Warning: Not enough vertices to form a polygon" << std::endl;
//         obj_file.close();
//         return;
//     }

//     // 杈撳嚭澶氳竟褰㈤《鐐?
//     for (size_t i = 0; i < polygon_vertices.size(); i++) {
//         const auto& v = polygon_vertices[i];
//         obj_file << "v " << CGAL::to_double(v.x()) << " "
//             << CGAL::to_double(v.y()) << " "
//             << CGAL::to_double(v.z()) << "\n";
//     }

//     // 杈撳嚭澶氳竟褰㈢殑杈癸紙绾挎锛?
//     obj_file << "\n# Polygon edges\n";
//     for (size_t i = 0; i < polygon_vertices.size(); i++) {
//         size_t next_i = (i + 1) % polygon_vertices.size();
//         obj_file << "l " << (i + 1) << " " << (next_i + 1) << "\n";
//     }

//     // 鍙€夛細涔熻緭鍑哄杈瑰舰鐨勯潰锛堝鏋滈《鐐瑰叡闈級
//     obj_file << "\n# Polygon face\n";
//     obj_file << "f";
//     for (size_t i = 0; i < polygon_vertices.size(); i++) {
//         obj_file << " " << (i + 1);
//     }
//     obj_file << "\n";

//     // 杈撳嚭瀵瑰簲鐨勭瀛愮偣浣嶇疆锛堢敤浜庡弬鑰冿級
//     obj_file << "\n# Corresponding seed points\n";
//     obj_file << "# v " << CGAL::to_double(v1->point().x()) << " "
//         << CGAL::to_double(v1->point().y()) << " "
//         << CGAL::to_double(v1->point().z()) << " # Seed 1\n";
//     obj_file << "# v " << CGAL::to_double(v2->point().x()) << " "
//         << CGAL::to_double(v2->point().y()) << " "
//         << CGAL::to_double(v2->point().z()) << " # Seed 2\n";

//     obj_file.close();
//     std::cout << "  * Saved to " << filename << std::endl;
// }

void write_voronoi_facets_to_obj_dedup(
    const std::string& filename,
    const std::vector<std::vector<Point_3>>& voronoi_facets,
    double epsilon = 1e-1,
    double angle_threshold_deg = 170.0,  // 瑙掑害闃堝€硷紙搴︼級锛屽ぇ浜庢鍊艰涓哄叡绾?
    bool debug = false  // 鏄惁杈撳嚭璋冭瘯淇℃伅
) {
    std::ofstream out(filename);
    if (!out) {
        std::cerr << "Cannot open file " << filename << std::endl;
        return;
    }
    out << "# Voronoi facets (epsilon=" << epsilon
        << ", angle_threshold=" << angle_threshold_deg << "掳)\n";

    std::vector<Point_3> vertices;
    std::vector<std::vector<int>> faces;
    double epsilon_squared = epsilon * epsilon;
    double cos_threshold = std::cos(angle_threshold_deg * M_PI / 180.0);

    int total_input_facets = 0;
    int total_output_facets = 0;

    // ---------- 杈呭姪鍑芥暟锛氭煡鎵炬垨娣诲姞椤剁偣 ----------
    auto find_or_add_vertex = [&](const Point_3& p) -> int {
        for (size_t i = 0; i < vertices.size(); ++i) {
            double dist_sq = CGAL::to_double(CGAL::squared_distance(vertices[i], p));
            if (dist_sq < epsilon_squared) {
                return static_cast<int>(i) + 1;
            }
        }
        vertices.push_back(p);
        return static_cast<int>(vertices.size());
        };

    // ---------- 杈呭姪鍑芥暟锛氳绠椾笁鐐瑰舰鎴愮殑瑙掑害锛堝害鏁帮級----------
    auto compute_angle = [&](const Point_3& p1, const Point_3& p2, const Point_3& p3) -> double {
        typedef CGAL::Vector_3<K> Vector_3;
        Vector_3 v1 = p1 - p2;  // 浠巔2鎸囧悜p1
        Vector_3 v2 = p3 - p2;  // 浠巔2鎸囧悜p3

        double len1_sq = CGAL::to_double(v1.squared_length());
        double len2_sq = CGAL::to_double(v2.squared_length());

        if (len1_sq < 1e-20 || len2_sq < 1e-20) return 0.0;

        double dot = CGAL::to_double(v1 * v2);
        double cos_angle = dot / std::sqrt(len1_sq * len2_sq);

        // 闄愬埗cos鍊煎湪[-1, 1]鑼冨洿鍐咃紝閬垮厤娴偣璇樊
        cos_angle = std::max(-1.0, std::min(1.0, cos_angle));

        double angle_rad = std::acos(cos_angle);
        return angle_rad * 180.0 / M_PI;
        };

    // ---------- 鏀硅繘鐨勮鐐规彁鍙栧嚱鏁?----------
    auto extract_corner_vertices = [&](const std::vector<Point_3>& polygon) -> std::vector<Point_3> {
        if (polygon.size() <= 3) return polygon;

        std::vector<Point_3> corners;
        std::vector<bool> is_corner(polygon.size(), false);
        int n = polygon.size();

        // 绗竴閬嶏細鏍囪鎵€鏈夎鐐?
        for (int i = 0; i < n; ++i) {
            const Point_3& prev = polygon[(i - 1 + n) % n];
            const Point_3& curr = polygon[i];
            const Point_3& next = polygon[(i + 1) % n];

            double angle = compute_angle(prev, curr, next);

            // 濡傛灉瑙掑害灏忎簬闃堝€硷紙涓嶆帴杩?80搴︼級锛屽垯鏄鐐?
            if (angle < angle_threshold_deg) {
                is_corner[i] = true;
            }
        }

        // 鏀堕泦瑙掔偣
        for (int i = 0; i < n; ++i) {
            if (is_corner[i]) {
                corners.push_back(polygon[i]);
            }
        }

        // 濡傛灉妫€娴嬪埌鐨勮鐐瑰お灏戯紝鐩存帴杩斿洖鍘熷澶氳竟褰?
        if (corners.size() < 3) {
            if (debug) {
                std::cout << "Warning: Only found " << corners.size()
                    << " corners, using original polygon with "
                    << polygon.size() << " vertices" << std::endl;
            }
            return polygon;
        }

        return corners;
        };

    // ---------- 鏋勫缓椤剁偣鍒楄〃鍜岄潰绱㈠紩 ----------
    for (const auto& facet : voronoi_facets) {
        total_input_facets++;

        if (facet.size() < 3) {
            if (debug) {
                std::cout << "Skipping facet with < 3 vertices" << std::endl;
            }
            continue;
        }

        // 鎻愬彇瑙掔偣
        std::vector<Point_3> corner_points = extract_corner_vertices(facet);

        if (debug && corner_points.size() != facet.size()) {
            std::cout << "Facet " << total_input_facets << ": "
                << facet.size() << " vertices -> "
                << corner_points.size() << " corners" << std::endl;
        }

        if (corner_points.size() < 3) continue;

        std::vector<int> face;
        for (const auto& p : corner_points) {
            int idx = find_or_add_vertex(p);
            face.push_back(idx);
        }

        // 绉婚櫎閲嶅鐨勮繛缁储寮?
        std::vector<int> cleaned_face;
        for (size_t i = 0; i < face.size(); ++i) {
            if (i == 0 || face[i] != face[i - 1]) {
                cleaned_face.push_back(face[i]);
            }
        }

        // 妫€鏌ラ灏炬槸鍚︾浉鍚?
        if (cleaned_face.size() > 1 && cleaned_face.front() == cleaned_face.back()) {
            cleaned_face.pop_back();
        }

        if (cleaned_face.size() >= 3) {
            faces.push_back(cleaned_face);
            total_output_facets++;
        }
    }

    // ---------- 杈撳嚭椤剁偣 ----------
    for (const auto& p : vertices) {
        out << "v "
            << p.x() << " "
            << p.y() << " "
            << p.z() << "\n";
    }

    // ---------- 杈撳嚭闈?----------
    for (const auto& f : faces) {
        out << "f";
        for (int idx : f) {
            out << " " << idx;
        }
        out << "\n";
    }

    out.close();

    std::cout << "=== Summary ===" << std::endl;
    std::cout << "Input facets: " << total_input_facets << std::endl;
    std::cout << "Output facets: " << total_output_facets << std::endl;
    std::cout << "Unique vertices: " << vertices.size() << std::endl;
    std::cout << "Wrote to: " << filename << std::endl;
}

void write_voronoi_facets_triangulated_obj(
    const std::string& filename,
    const std::vector<std::vector<Point_3>>& voronoi_facets
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

    auto get_index = [&](const Point_3& p) {
        auto it = vertex_index.find(p);
        if (it != vertex_index.end())
            return it->second;

        int idx = static_cast<int>(vertices.size()) + 1; // OBJ 1-based
        vertex_index[p] = idx;
        vertices.push_back(p);
        return idx;
        };

    // ---------- triangulate each facet ----------
    for (const auto& facet : voronoi_facets) {
        if (facet.size() < 3) continue;

        int v0 = get_index(facet[0]);

        for (size_t i = 1; i + 1 < facet.size(); ++i) {
            int v1 = get_index(facet[i]);
            int v2 = get_index(facet[i + 1]);

            triangles.push_back({ v0, v1, v2 });
        }
    }

    // ---------- write vertices ----------
    for (const auto& p : vertices) {
        out << "v "
            << p.x() << " "
            << p.y() << " "
            << p.z() << "\n";
    }

    // ---------- write triangles ----------
    for (const auto& t : triangles) {
        out << "f "
            << t[0] << " "
            << t[1] << " "
            << t[2] << "\n";
    }

    out.close();
}
void write_voronoi_facets_to_obj(
    const std::string& filename,
    const std::vector<std::vector<Point_3>>& voronoi_facets
) {
    std::ofstream out(filename);
    if (!out) {
        std::cerr << "Cannot open file " << filename << std::endl;
        return;
    }

    out << "# Voronoi facets OBJ\n";

    size_t vertex_offset = 1; // OBJ index starts from 1

    for (const auto& facet : voronoi_facets) {
        if (facet.size() < 3) continue;

        // 1. write vertices
        for (const auto& p : facet) {
            out << "v "
                << p.x() << " "
                << p.y() << " "
                << p.z() << "\n";
        }

        // 2. write face
        out << "f";
        for (size_t i = 0; i < facet.size(); ++i) {
            out << " " << (vertex_offset + i);
        }
        out << "\n";

        vertex_offset += facet.size();
    }

    out.close();
}
void exportBoundarySurfaceImproved(const Delaunay& dt,
    const std::map<Vertex_handle, int>& vertex_labels,
    const std::string& filename) {
    std::ofstream out(filename);
    if (!out) {
        std::cerr << "Error: Cannot create file " << filename << std::endl;
        return;
    }

    out << "# Voronoi boundary surface (improved with validation)" << std::endl;
    std::cout << "\n=== Extracting Boundary Surface (Improved) ===" << std::endl;

    // ========================================
    // 缁熻锛氬府鍔╃悊瑙ｅ摢浜涢潰琚繚鐣?涓㈠純
    // ========================================
    int total_edges = 0;
    int boundary_edges_count = 0;  // label 1-2
    int inner_edges = 0;            // label 1-1
    int outer_edges = 0;            // label 2-2
    int no_label_edges = 0;         // 娌℃湁鏍囩鐨勮竟

    std::set<std::pair<Vertex_handle, Vertex_handle>> boundary_edges;

    // ========================================
    // 姝ラ1锛氶亶鍘嗘墍鏈塂elaunay杈?
    // 姣忔潯杈瑰搴斾竴涓猇oronoi闈?
    // ========================================
    std::vector<std::vector<Point_3>> voronoi_facets;
    int fla = 0;
    for (auto eit = dt.finite_edges_begin(); eit != dt.finite_edges_end(); ++eit) {
        total_edges++;

        Cell_handle c = eit->first;
        int i = eit->second;
        int j = eit->third;

        // 鑾峰彇杈圭殑涓や釜绔偣锛堢珯鐐癸級
        Vertex_handle v1 = c->vertex(i);
        Vertex_handle v2 = c->vertex(j);

        // 妫€鏌ユ槸鍚﹂兘鏈夋爣绛?
        if (vertex_labels.find(v1) == vertex_labels.end() ||
            vertex_labels.find(v2) == vertex_labels.end()) {
            no_label_edges++;
            continue;
        }

        int label1 = vertex_labels.at(v1);
        int label2 = vertex_labels.at(v2);

        // ========================================
        // 鏍稿績鍒ゆ柇锛氫繚鐣欐潯浠?
        // 
        // Voronoi闈㈢敱2涓珯鐐圭敓鎴?
        // 淇濈暀鏉′欢锛氳繖2涓珯鐐圭殑label涓嶅悓
        // ========================================
        if (label1 != label2) {
        
            std::vector<Point_3> facet_vertices;

            Delaunay::Cell_circulator cc = dt.incident_cells(*eit);
            Delaunay::Cell_circulator done = cc;

            if (cc == nullptr) continue;

            do {
                if (!dt.is_infinite(cc)) {
                    // Compute circumcenter manually or use CGAL dual
                    // Using dual() is safer and usually cached/optimized in some kernels, 
                    // but manual construction as you did is also fine for Exact kernel.
                    Point_3 center = dt.dual(cc);
                    facet_vertices.push_back(center);
                }
                ++cc;
            } while (cc != done);

            if (facet_vertices.size() >= 3) {
              /*  if (fla == 0)
                {
                    export_single_voronoi_polygon(dt, v1, v2, "single_voronoi_polygon.obj");
                    fla = 1;
                }*/
                voronoi_facets.push_back(facet_vertices);
            }
       
            // [OK] 杈圭晫闈細涓€涓唴鐐?+ 涓€涓鐐?
            boundary_edges_count++;

            // 瑙勮寖鍖栬竟锛堜繚璇佸敮涓€鎬э級
            if (v1 < v2) {
                boundary_edges.insert({ v1, v2 });
            }
            else {
                boundary_edges.insert({ v2, v1 });
            }
        }
        else {
            // [X] 闈炶竟鐣岄潰
            if (label1 == 1) {
                inner_edges++;  // 鍐?鍐?
            }
            else if (label1 == 2) {
                outer_edges++;  // 澶?澶?
            }
        }
    }

    write_voronoi_facets_to_obj("voronoi.obj", voronoi_facets);
    write_voronoi_facets_to_obj_dedup("voronoi_dedup.obj", voronoi_facets);
    write_voronoi_facets_triangulated_obj("voronoi_dedup_triangulated.obj", voronoi_facets);
std::cout << "  * Found " << voronoi_facets.size() << " Voronoi facets for seed pairs" << std::endl;

// 3. Write Surface Mesh OBJ file
std::ofstream obj_file("output_filename.obj");
if (!obj_file.is_open()) {
    //  std::cerr << "Error: Cannot open output file " << output_filename << std::endl;
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
    for (const auto& pt : facet) centroid = centroid + (pt - CGAL::ORIGIN);
    double s = static_cast<double>(facet.size());
    centroid = Point_3(centroid.x() / s, centroid.y() / s, centroid.z() / s);
    size_t centroid_idx = get_vertex_index(centroid);

    for (size_t i = 0; i < indices.size(); i++) {
        size_t idx0 = indices[i];
        size_t idx1 = indices[(i + 1) % indices.size()];

        // 璺宠繃閫€鍖栦笁瑙掑舰
        if (idx0 == idx1 || idx0 == centroid_idx || idx1 == centroid_idx) continue;

        face_indices.push_back({ centroid_idx, idx0, idx1 });
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

//  std::cout << "  * Surface mesh saved to " << output_filename << std::endl;
std::cout << "  * Total vertices: " << unique_vertices.size() << ", triangles: " << face_indices.size() << std::endl;
    // ========================================
    // 鎵撳嵃缁熻锛氶獙璇佺粨鏋滄槸鍚﹀悎鐞?
    // ========================================
    std::cout << "\n--- Edge Classification ---" << std::endl;
    std::cout << "Total Delaunay edges: " << total_edges << std::endl;
    std::cout << "|- Boundary edges (1-2): " << boundary_edges_count
        << " (" << (100.0 * boundary_edges_count / total_edges) << "%)" << std::endl;
    std::cout << "|- Inner edges (1-1): " << inner_edges
        << " (" << (100.0 * inner_edges / total_edges) << "%)" << std::endl;
    std::cout << "|- Outer edges (2-2): " << outer_edges
        << " (" << (100.0 * outer_edges / total_edges) << "%)" << std::endl;
    std::cout << "'- No label edges: " << no_label_edges << std::endl;





    // ========================================
    // 姝ラ2锛氫负姣忎釜杈圭晫杈规彁鍙朧oronoi闈?
    // ========================================
    std::map<Point, int> point_indices;
    int vertex_count = 0;
    std::vector<std::vector<int>> faces;

    int valid_faces = 0;
    int degenerate_faces = 0;  // 閫€鍖栫殑闈紙椤剁偣<3鎴栭噸澶嶏級

    std::cout << "\n--- Extracting Voronoi Faces ---" << std::endl;

    for (const auto& edge_pair : boundary_edges) {
        Vertex_handle v1 = edge_pair.first;
        Vertex_handle v2 = edge_pair.second;

        // 鑾峰彇杩欐潯杈瑰搴旂殑Voronoi闈㈢殑椤剁偣
        std::vector<Point> face_points = getVoronoiFace(dt, v1, v2);

        // 杩囨护閫€鍖栫殑闈?
        if (face_points.size() < 3) {
            degenerate_faces++;
            continue;
        }

        // 妫€鏌ユ槸鍚︽湁閲嶅椤剁偣
        std::set<Point> unique_points(face_points.begin(), face_points.end());
        if (unique_points.size() < 3) {
            degenerate_faces++;
            continue;
        }

        valid_faces++;

        // 鏀堕泦椤剁偣绱㈠紩
        std::vector<int> face_indices;
        for (const auto& p : face_points) {
            if (point_indices.find(p) == point_indices.end()) {
                point_indices[p] = ++vertex_count;
            }
            face_indices.push_back(point_indices[p]);
        }
        faces.push_back(face_indices);
    }

    std::cout << "Valid faces: " << valid_faces << std::endl;
    std::cout << "Degenerate faces (filtered out): " << degenerate_faces << std::endl;

    // ========================================
    // 姝ラ3锛氬啓鍏BJ鏂囦欢
    // ========================================

    // 鍐欏叆椤剁偣
    std::vector<Point> ordered_points(vertex_count);
    for (const auto& pair : point_indices) {
        ordered_points[pair.second - 1] = pair.first;
    }

    for (const auto& p : ordered_points) {
        out << "v " << p.x() << " " << p.y() << " " << p.z() << std::endl;
    }

    // 涓夎鍖栧苟鍐欏叆闈?
    int triangle_count = 0;
    std::map<int, int> face_vertex_distribution;

    for (const auto& face : faces) {
        int num_vertices = face.size();
        face_vertex_distribution[num_vertices]++;

        if (face.size() >= 3) {
            // 鎵囧舰涓夎鍖栦粠绗竴涓《鐐瑰紑濮?
            for (size_t i = 1; i < face.size() - 1; ++i) {
                out << "f " << face[0] << " "
                    << face[i] << " " << face[i + 1] << std::endl;
                triangle_count++;
            }
        }
    }

    out.close();

    // ========================================
    // 鏈€缁堟姤鍛?
    // ========================================
    std::cout << "\n--- Final Statistics ---" << std::endl;
    std::cout << "Output file: " << filename << std::endl;
    std::cout << "|- Vertices: " << vertex_count << std::endl;
    std::cout << "|- Voronoi faces: " << valid_faces << std::endl;
    std::cout << "'- Triangles: " << triangle_count << std::endl;

    std::cout << "\nFace complexity (vertices per face):" << std::endl;
    for (const auto& pair : face_vertex_distribution) {
        std::cout << "  " << pair.first << " vertices: "
            << pair.second << " faces" << std::endl;
    }

    // ========================================
    // 楠岃瘉缁撴灉
    // ========================================
    std::cout << "\n--- Validation ---" << std::endl;

    double boundary_percentage = 100.0 * boundary_edges_count / total_edges;

    std::cout << "Boundary edges percentage: " << boundary_percentage << "%" << std::endl;

    // 鏈熸湜鍊硷細濡傛灉鍐呭鐐规暟閲忔帴杩戯紝杈圭晫杈瑰簲璇ュ崰50-70%
    if (boundary_percentage < 10.0) {
        std::cout << "[WARN]  WARNING: Very few boundary edges!" << std::endl;
        std::cout << "   -> Check if labels are correct" << std::endl;
        std::cout << "   -> Both labels should be present in the data" << std::endl;
    }
    else if (boundary_percentage > 90.0) {
        std::cout << "[WARN]  WARNING: Almost all edges are boundary!" << std::endl;
        std::cout << "   -> This is unusual, check your labels" << std::endl;
    }
    else if (boundary_percentage >= 40.0 && boundary_percentage <= 70.0) {
        std::cout << "[OK] Boundary percentage looks reasonable" << std::endl;
        std::cout << "   -> Labels appear to be well-distributed" << std::endl;
    }

    // 妫€鏌ラ潰鐨勬湁鏁堟€?
    double valid_ratio = 100.0 * valid_faces / boundary_edges_count;
    std::cout << "\nValid face ratio: " << valid_ratio << "%" << std::endl;

    if (valid_ratio < 80.0) {
        std::cout << "[WARN]  WARNING: Many faces were filtered as degenerate" << std::endl;
        std::cout << "   -> This might indicate issues with the triangulation" << std::endl;
    }
    else {
        std::cout << "[OK] Most faces are valid (good!)" << std::endl;
    }

    std::cout << "\n=== Extraction Complete ===" << std::endl;
}
void Generator::generate_surface_mesh1(MeshingTree* seeds, const char* output_filename)
{
    std::cout << "Generating surface mesh using CGAL Voronoi1..." << std::endl;

    size_t num_seeds = seeds->get_num_tree_points();
    if (num_seeds == 0) {
        std::cerr << "Error: No seeds provided." << std::endl;
        return;
    }

    // 1. Build Delaunay triangulation with seed index info
  //  Delaunay dt;
    std::vector<Vertex_handle> vertex_handles(num_seeds);
    std::vector<LabeledPoint> labeled_points;
    for (size_t i = 0; i < num_seeds; i++) {
        if (!seeds->tree_point_is_active(i)) continue;
        double* pt = seeds->get_tree_point(i);
      /*  Vertex_handle vh = dt.insert(Point_3(pt[0], pt[1], pt[2]));
        vh->info() = i;
        vertex_handles[i] = vh;*/ 
        size_t* attrib = seeds->get_tree_point_attrib(i);

        // 鏍规嵁鍖哄煙ID鍒ゆ柇鍐呭
        // attrib[5] 鏄尯鍩烮D: 0 = 澶栭儴, 闈? = 鍐呴儴
       // if (attrib[5] == 1)
        labeled_points.push_back(LabeledPoint(
           pt[0],
            pt[1],
            pt[2],
            attrib[5]
        ));
    }

    if (labeled_points.empty()) {
        std::cerr << "Error: No points loaded!" << std::endl;
        return ;
    }

    std::cout << "Loaded " << labeled_points.size() << " points" << std::endl;

    // Count labels
    std::map<int, int> label_counts;
    for (const auto& lp : labeled_points) {
        label_counts[lp.label]++;
    }

    std::cout << "\nLabel distribution:" << std::endl;
    for (const auto& pair : label_counts) {
        std::cout << "  Label " << pair.first << ": " << pair.second << " points" << std::endl;
    }

    // Build Delaunay triangulation
    std::cout << "\nBuilding Delaunay triangulation..." << std::endl;
    Delaunay dt;
    std::map<Vertex_handle, int> vertex_labels;

    for (const auto& lp : labeled_points) {
        Vertex_handle vh = dt.insert(lp.point);
        vertex_labels[vh] = lp.label;
    }

    std::cout << "Number of vertices: " << dt.number_of_vertices() << std::endl;
    std::cout << "Number of cells: " << dt.number_of_cells() << std::endl;
    std::cout << "Number of facets: " << dt.number_of_facets() << std::endl;

    // Export results
    std::cout << "\n=== Exporting Surface Meshes ===" << std::endl;

    // Method 1: Voronoi-based surface (from Voronoi faces)
    exportBoundarySurfaceImproved(dt, vertex_labels, "voronoi_surface.obj");

    // Method 2: Delaunay-based surface (from Delaunay facets)
    exportDelaunayBoundarySurface(dt, vertex_labels, "delaunay_surface.obj");

    // Export input points for reference
    std::ofstream points_out("input_points.obj");
    std::ofstream points_out1("output_points.obj");
    if (points_out) {
        points_out << "# Input points" << std::endl;
        for (const auto& lp : labeled_points) {
            if (lp.label != 1)  // Skip label 0 if needed
                points_out << "v " << lp.point.x() << " "
                << lp.point.y() << " " << lp.point.z() << std::endl;

            else
                points_out1 << "v " << lp.point.x() << " "
                << lp.point.y() << " " << lp.point.z() << std::endl;
        }
        points_out.close();
        points_out1.close();
        std::cout << "\nInput points exported to input_points.obj" << std::endl;
    }

    std::cout << "\n=== Complete ===" << std::endl;
    std::cout << "\nGenerated surface files:" << std::endl;
    std::cout << "  - voronoi_surface.obj    : Voronoi-based surface mesh (triangulated)" << std::endl;
    std::cout << "  - delaunay_surface.obj   : Delaunay-based surface mesh (faster)" << std::endl;
    std::cout << "  - input_points.obj       : Original input points" << std::endl;
    std::cout << "\nRecommendation:" << std::endl;
    std::cout << "  - Use delaunay_surface.obj for faster results" << std::endl;
    std::cout << "  - Use voronoi_surface.obj for true Voronoi surface" << std::endl;
} 
