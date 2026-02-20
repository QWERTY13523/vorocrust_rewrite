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


struct PairHash {
    inline size_t operator()(const std::pair<size_t, size_t> &v) const {
        std::hash<size_t> hasher;
        size_t seed = 0;
        // Hash combine 算法
        seed ^= hasher(v.first) + 0x9e3779b9 + (seed<<6) + (seed>>2);
        seed ^= hasher(v.second) + 0x9e3779b9 + (seed<<6) + (seed>>2);
        return seed;
    }
};

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
    Geometry geom;
    size_t num_spheres = spheres->get_num_tree_points();
    size_t num_faces_total = face_flat.size() / 3;

    std::cout << "=== Optimizer: fixing non-adjacent seed pairs (Single-thread Optimized) ===" << std::endl;

    auto pair_key = [](size_t a, size_t b) -> std::pair<size_t, size_t> {
        return (a < b) ? std::make_pair(a, b) : std::make_pair(b, a);
    };

    auto make_face_key = [](size_t a, size_t b, size_t c) -> std::tuple<size_t, size_t, size_t> {
        size_t arr[3] = {a, b, c};
        // 简单的排序网络
        if (arr[0] > arr[1]) std::swap(arr[0], arr[1]);
        if (arr[1] > arr[2]) std::swap(arr[1], arr[2]);
        if (arr[0] > arr[1]) std::swap(arr[0], arr[1]);
        return std::make_tuple(arr[0], arr[1], arr[2]);
    };

    // [优化1] 将拓扑结构的构建移出循环。
    // 球和面的连接关系在优化过程中是不变的。
    std::vector<std::vector<size_t>> sphere_to_faces(num_spheres);
    for (size_t fi = 0; fi < num_faces_total; fi++) {
        for (int k = 0; k < 3; k++) {
            int sid = face_flat[fi * 3 + k];
            if (sid >= 0 && (size_t)sid < num_spheres)
                sphere_to_faces[sid].push_back(fi);
        }
    }

    // 辅助函数：计算种子点
    auto compute_face_seeds = [&](size_t si, size_t sj, size_t sk,
                                   double* seed_a, double* seed_b,
                                   double* normal_out) -> bool {
        double sp_i[4], sp_j[4], sp_k[4];
        spheres->get_tree_point(si, 4, sp_i);
        spheres->get_tree_point(sj, 4, sp_j);
        spheres->get_tree_point(sk, 4, sp_k);

        // 快速拒绝
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

    // 辅助函数：点积判断区域
    auto compute_seed_dot = [&](double* seed_pos, size_t si, size_t sj, size_t sk) -> double {
        double dot_sum = 0.0;
        int count = 0;
        size_t sids[3] = {si, sj, sk};
        for (int i = 0; i < 3; i++) {
            if (sids[i] >= num_spheres) continue;
            double center[4];
            spheres->get_tree_point(sids[i], 4, center);
            double* normal = spheres->get_tree_point_normal(sids[i]);
            // 简单的零向量检查
            if (std::abs(normal[0]) < 1e-9 && std::abs(normal[1]) < 1e-9 && std::abs(normal[2]) < 1e-9) continue;
            
            dot_sum += (seed_pos[0] - center[0]) * normal[0]
                     + (seed_pos[1] - center[1]) * normal[1]
                     + (seed_pos[2] - center[2]) * normal[2];
            count++;
        }
        return (count > 0) ? (dot_sum / count) : 0.0;
    };

    // 辅助函数：检查power vertex是否在三角形内部（使用重心坐标）
    auto is_point_inside_triangle = [&](double* p, double* a, double* b, double* c) -> bool {
        // 计算三角形法向量
        double v0[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]};
        double v1[3] = {c[0] - a[0], c[1] - a[1], c[2] - a[2]};
        double v2[3] = {p[0] - a[0], p[1] - a[1], p[2] - a[2]};
        
        // 计算点积
        double dot00 = v0[0]*v0[0] + v0[1]*v0[1] + v0[2]*v0[2];
        double dot01 = v0[0]*v1[0] + v0[1]*v1[1] + v0[2]*v1[2];
        double dot02 = v0[0]*v2[0] + v0[1]*v2[1] + v0[2]*v2[2];
        double dot11 = v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2];
        double dot12 = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
        
        // 计算重心坐标
        double inv_denom = 1.0 / (dot00 * dot11 - dot01 * dot01);
        double u = (dot11 * dot02 - dot01 * dot12) * inv_denom;
        double v = (dot00 * dot12 - dot01 * dot02) * inv_denom;
        
        // 检查是否在三角形内部
        return (u >= -1e-10) && (v >= -1e-10) && (u + v <= 1.0 + 1e-10);
    };

    const int max_iterations = 10;
    const int binary_search_iterations = 40; // 保持原有的迭代次数

    // [优化2] 内存复用，避免循环内反复分配
    std::vector<Vertex_handle> vertex_handles;
    std::vector<std::pair<size_t, size_t>> non_adjacent_pairs;
    std::vector<size_t> seeds_to_del;
    // 使用 unordered_map 替代 map
    std::unordered_map<std::tuple<size_t,size_t,size_t>, std::vector<size_t>, TupleHash> face_key_to_seeds;
    // 使用 unordered_set 替代 set
    std::unordered_set<std::pair<size_t, size_t>, PairHash> delaunay_edges;
    std::unordered_set<size_t> processed_spheres;

    for (int iter = 0; iter < max_iterations; iter++) {

        size_t num_seeds = seeds->get_num_tree_points();
        vertex_handles.assign(num_seeds, Vertex_handle()); 

        // ---- Build Delaunay triangulation ----
        Delaunay dt;
        for (size_t i = 0; i < num_seeds; i++) {
            if (!seeds->tree_point_is_active(i)) continue;
            double* pt = seeds->get_tree_point(i);
            Vertex_handle vh = dt.insert(Point_3(pt[0], pt[1], pt[2]));
            vh->info() = i;
            vertex_handles[i] = vh;
        }

        // ---- Collect Delaunay edges ----
        delaunay_edges.clear();
        // 预估大小，减少 rehash
        delaunay_edges.reserve(num_seeds * 7); 

        for (auto eit = dt.finite_edges_begin(); eit != dt.finite_edges_end(); ++eit) {
            Cell_handle c = eit->first;
            Vertex_handle v1 = c->vertex(eit->second);
            Vertex_handle v2 = c->vertex(eit->third);
            if (dt.is_infinite(v1) || dt.is_infinite(v2)) continue;
            size_t s1 = v1->info();
            size_t s2 = v2->info();
            if (s1 >= num_seeds || s2 >= num_seeds || s1 == s2) continue;
            delaunay_edges.insert(pair_key(s1, s2)); // O(1) 插入
        }

        // ---- Find non-adjacent seed pairs ----
        non_adjacent_pairs.clear();
        for (size_t i = 0; i < num_seeds; i++) {
            if (!seeds->tree_point_is_active(i)) continue;
            size_t* attrib = seeds->get_tree_point_attrib(i);
            size_t j = attrib[1];
            if (j <= i || j >= num_seeds) continue;
            if (!seeds->tree_point_is_active(j)) continue;
            size_t* pair_attrib = seeds->get_tree_point_attrib(j);
            if (pair_attrib[1] != i) continue;
            
            // O(1) 查找
            if (delaunay_edges.find(pair_key(i, j)) == delaunay_edges.end()) {
                non_adjacent_pairs.emplace_back(i, j);
            }
        }

        // ---- Find same-side close seed pairs ----
        std::vector<std::pair<size_t, size_t>> same_side_close_pairs;
        for (size_t i = 0; i < num_seeds; i++) {
            if (!seeds->tree_point_is_active(i)) continue;
            size_t* attrib_i = seeds->get_tree_point_attrib(i);
            size_t region_i = attrib_i[5];
            double* pos_i = seeds->get_tree_point(i);
            
            for (size_t j = i + 1; j < num_seeds; j++) {
                if (!seeds->tree_point_is_active(j)) continue;
                size_t* attrib_j = seeds->get_tree_point_attrib(j);
                size_t region_j = attrib_j[5];
                
                // 检查是否同侧（相同region ID）且不是配对种子
                if (region_i == region_j && attrib_i[1] != j) {
                    double* pos_j = seeds->get_tree_point(j);
                    double dist = geom.distance(3, pos_i, pos_j);
                    if (dist < 1e-2) {
                        same_side_close_pairs.emplace_back(i, j);
                    }
                }
            }
        }

        std::cout << "  Iter " << iter << ": " << non_adjacent_pairs.size()
                  << " non-adjacent pairs, " << same_side_close_pairs.size()
                  << " same-side close pairs" << std::endl;
        if (non_adjacent_pairs.empty() && same_side_close_pairs.empty()) break;

        // ---- Build face-key -> seeds map ----
        face_key_to_seeds.clear();
        // 预估 bucket 数量
        face_key_to_seeds.reserve(num_seeds / 2);
        for (size_t i = 0; i < num_seeds; i++) {
            if (!seeds->tree_point_is_active(i)) continue;
            size_t* attrib = seeds->get_tree_point_attrib(i);
            auto key = make_face_key(attrib[2], attrib[3], attrib[4]);
            face_key_to_seeds[key].push_back(i);
        }

        processed_spheres.clear();
        size_t fixes_this_iter = 0;

        // 合并两种问题种子对
        std::vector<std::tuple<size_t, size_t, bool>> problematic_pairs;
        for (const auto& pr : non_adjacent_pairs) {
            problematic_pairs.emplace_back(pr.first, pr.second, false); // false = non-adjacent
        }
        for (const auto& pr : same_side_close_pairs) {
            problematic_pairs.emplace_back(pr.first, pr.second, true); // true = same-side close
        }

        for (size_t pi = 0; pi < problematic_pairs.size(); pi++) {
            size_t seed_i = std::get<0>(problematic_pairs[pi]);
            size_t seed_j = std::get<1>(problematic_pairs[pi]);
            bool is_same_side = std::get<2>(problematic_pairs[pi]);
            
            if (!seeds->tree_point_is_active(seed_i) ||
                !seeds->tree_point_is_active(seed_j)) continue;

            size_t* attrib_i = seeds->get_tree_point_attrib(seed_i);
            size_t* attrib_j = seeds->get_tree_point_attrib(seed_j);
            
            size_t sph_ids[3];
            if (is_same_side) {
                // 对于同侧close pair，需要找到power vertex在三角形外的face
                // 检查两个face（seed_i的face和seed_j的face）
                size_t face_i[3] = {attrib_i[2], attrib_i[3], attrib_i[4]};
                size_t face_j[3] = {attrib_j[2], attrib_j[3], attrib_j[4]};
                
                bool use_face_i = false;
                
                // 检查face_i的power vertex是否在三角形外
                double sp_i[4], sp_j[4], sp_k[4];
                spheres->get_tree_point(face_i[0], 4, sp_i);
                spheres->get_tree_point(face_i[1], 4, sp_j);
                spheres->get_tree_point(face_i[2], 4, sp_k);
                
                double* centers[3] = {sp_i, sp_j, sp_k};
                double radii[3] = {sp_i[3], sp_j[3], sp_k[3]};
                double c_ijk[3];
                
                if (geom.get_power_vertex(3, 3, centers, radii, c_ijk)) {
                    if (!is_point_inside_triangle(c_ijk, sp_i, sp_j, sp_k)) {
                        use_face_i = true;
                    }
                }
                
                // 如果face_i的power vertex在内部，检查face_j
                if (!use_face_i) {
                    spheres->get_tree_point(face_j[0], 4, sp_i);
                    spheres->get_tree_point(face_j[1], 4, sp_j);
                    spheres->get_tree_point(face_j[2], 4, sp_k);
                    
                    centers[0] = sp_i; centers[1] = sp_j; centers[2] = sp_k;
                    radii[0] = sp_i[3]; radii[1] = sp_j[3]; radii[2] = sp_k[3];
                    
                    if (geom.get_power_vertex(3, 3, centers, radii, c_ijk)) {
                        if (!is_point_inside_triangle(c_ijk, sp_i, sp_j, sp_k)) {
                            sph_ids[0] = face_j[0];
                            sph_ids[1] = face_j[1];
                            sph_ids[2] = face_j[2];
                        } else {
                            continue; // 两个face的power vertex都在内部，跳过
                        }
                    } else {
                        continue; // power vertex计算失败
                    }
                } else {
                    sph_ids[0] = face_i[0];
                    sph_ids[1] = face_i[1];
                    sph_ids[2] = face_i[2];
                }
            } else {
                // 对于non-adjacent pair，使用原来的逻辑
                sph_ids[0] = attrib_i[2];
                sph_ids[1] = attrib_i[3];
                sph_ids[2] = attrib_i[4];
            }

            // ---- Try all 3 spheres ----
            size_t best_sid = SIZE_MAX;
            double best_valid_r = 0.0;
            double best_orig_r = 0.0;
            double best_shrink_ratio = 0.0; 
            std::vector<size_t> best_seeds_to_delete;

            for (int try_k = 0; try_k < 3; try_k++) {
                size_t sid = sph_ids[try_k];
                if (sid >= num_spheres) continue;
                if (processed_spheres.count(sid)) continue;

                double* sph_data = spheres->get_tree_point(sid);
                double  orig_r   = sph_data[3];

                // ---- Compute minimum viable radius ----
                double r_min = 0.0;
                // 直接使用外层计算好的 sphere_to_faces，无需搜索
                const auto& connected_faces = sphere_to_faces[sid];
                
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
                if (r_min >= orig_r - 1e-10) continue; 

                // ---- Collect old seeds ----
                seeds_to_del.clear();
                
                for (size_t fi : connected_faces) {
                    auto key = make_face_key(
                        (size_t)face_flat[fi*3],
                        (size_t)face_flat[fi*3+1],
                        (size_t)face_flat[fi*3+2]);
                    
                    // O(1) 查找
                    auto it = face_key_to_seeds.find(key);
                    if (it != face_key_to_seeds.end()) {
                        for (size_t s : it->second) {
                            if (seeds->tree_point_is_active(s))
                                seeds_to_del.push_back(s);
                        }
                    }
                }

                // ---- Binary search ----
                double lo   = r_min;
                double hi_r = orig_r;
                double valid_r = orig_r;
                bool   found   = false;

                for (int bs = 0; bs < binary_search_iterations; bs++) {
                    double mid = (lo + hi_r) / 2.0;
                    sph_data[3] = mid;

                    bool all_ok = true;

                    // C1: face validity check (较快)
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

                    // C2: distance check (较慢，只在 C1 通过后执行)
                    if (all_ok) {
                        for (size_t fi : connected_faces) {
                            double sa[4], sb[4], nm[3];
                            compute_face_seeds(
                                (size_t)face_flat[fi*3],
                                (size_t)face_flat[fi*3+1],
                                (size_t)face_flat[fi*3+2],
                                sa, sb, nm);

                            // [优化3] 几何剪枝
                            // 在调用昂贵的 get_closest_tree_point 之前，先检查生成的一对种子点是否太近
                            // 如果它们彼此都分不开，就更别提和别的种子分开了。
                            if (geom.distance(3, sa, sb) <= 1e-4) { all_ok = false; break; }

                            // 这里的搜索比较耗时，但它是必须的。
                            // 由于我们已经优化了外围结构，这部分是剩余的硬骨头。
                            size_t ic; double hc = DBL_MAX;
                            seeds->get_closest_tree_point(
                                sa, seeds_to_del.size(),
                                seeds_to_del.data(), ic, hc);
                            if (hc < 1e-2) { all_ok = false; break; }

                            hc = DBL_MAX;
                            seeds->get_closest_tree_point(
                                sb, seeds_to_del.size(),
                                seeds_to_del.data(), ic, hc);
                            if (hc < 1e-2) { all_ok = false; break; }
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

                sph_data[3] = orig_r;  // restore

                if (!found) continue;

                double shrink_ratio = (orig_r - valid_r) / orig_r;
                if (best_sid == SIZE_MAX || shrink_ratio > best_shrink_ratio) {
                    best_sid = sid;
                    best_valid_r = valid_r;
                    best_orig_r = orig_r;
                    best_shrink_ratio = shrink_ratio;
                    best_seeds_to_delete = std::move(seeds_to_del);
                }
            }

            if (best_sid == SIZE_MAX) continue; 
            processed_spheres.insert(best_sid);

            // ---- Apply the best shrink ----
            double* sph_data = spheres->get_tree_point(best_sid);
            sph_data[3] = best_valid_r;
            
            std::cout << "    Sphere " << best_sid
                      << ": radius " << best_orig_r << " -> " << best_valid_r
                      << " (shrink " << (best_shrink_ratio * 100.0) << "%)" << std::endl;

            for (size_t s : best_seeds_to_delete) {
                seeds->lazy_delete_tree_point(s);
            }

            for (size_t fi : sphere_to_faces[best_sid]) {
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

                // Update face_key_to_seeds map so subsequent iterations know about the new seeds
                auto key = make_face_key(s0, s1, s2);
                face_key_to_seeds[key].push_back(idx_a);
                face_key_to_seeds[key].push_back(idx_b);
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
        bool is_interface = (attrib1[5] != attrib2[5]) && (attrib1[1] == seed_idx2 || attrib2[1] == seed_idx1);
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
            double epsilon = 1e-2,
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
            int boundary_pair_edges = 0;
            int boundary_non_pair_edges = 0;
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
                    size_t seed_idx1 = v1->info();
                    size_t seed_idx2 = v2->info();
                    size_t* attrib1 = seeds->get_tree_point_attrib(seed_idx1);
                    size_t* attrib2 = seeds->get_tree_point_attrib(seed_idx2);
                    bool is_pair_connection = (attrib1[1] == seed_idx2) || (attrib2[1] == seed_idx1);
                    boundary_edges_count++;
                    if (!is_pair_connection) {
                        boundary_non_pair_edges++;
                        continue;
                    }

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
                        boundary_pair_edges++;
                        // if (!did_export_single) {
                        //     export_single_voronoi_polygon(v1, v2, "single_voronoi_polygon.obj");
                        //     did_export_single = true;
                        // }
                    }
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
        std::cout << "|- Pair boundary edges (1-2) preserved: " << boundary_pair_edges
            << " (" << (100.0 * boundary_pair_edges / std::max(1, total_edges)) << "%)" << std::endl;
        std::cout << "|- Boundary edges (1-2) skipped (no pair link): " << boundary_non_pair_edges
            << std::endl;
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

    }
}
