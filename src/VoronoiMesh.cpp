#include "Generator.h"
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <fstream>
#include <iomanip>
#include <array>
#include <queue>
#include <string>
#include<vector>
#include <algorithm>
#include <random>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Delaunay_triangulation_cell_base_3.h> // 必须包含
#include <CGAL/intersections.h>
#include <CGAL/Timer.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_with_info_3<size_t, K> Vb;
typedef CGAL::Delaunay_triangulation_cell_base_3<K> Cb; // 定义 Cell Base
typedef CGAL::Triangulation_data_structure_3<Vb, Cb> Tds; // 传入 Vb 和 Cb
typedef CGAL::Delaunay_triangulation_3<K, Tds> Delaunay;
typedef Delaunay::Point Point_3;
typedef Delaunay::Vertex_handle Vertex_handle;
typedef Delaunay::Cell_handle Cell_handle;
typedef Delaunay::Edge Edge;

namespace {
struct EdgeKey {
    size_t a;
    size_t b;
    bool operator==(const EdgeKey& o) const noexcept { return a == o.a && b == o.b; }
};

struct EdgeKeyHash {
    size_t operator()(const EdgeKey& k) const noexcept {
        size_t h = 0;
        h ^= std::hash<size_t>{}(k.a) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        h ^= std::hash<size_t>{}(k.b) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        return h;
    }
};

struct SeedPairKey {
    size_t a;
    size_t b;
    bool operator==(const SeedPairKey& o) const noexcept { return a == o.a && b == o.b; }
};

struct SeedPairKeyHash {
    size_t operator()(const SeedPairKey& k) const noexcept {
        size_t h = 0;
        h ^= std::hash<size_t>{}(k.a) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        h ^= std::hash<size_t>{}(k.b) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        return h;
    }
};

struct EdgeAdj {
    int count = 0;
    size_t tri0 = static_cast<size_t>(-1);
    size_t tri1 = static_cast<size_t>(-1);
};

struct VoronoiFacet {
    std::vector<Point_3> vertices;
    size_t seed1;
    size_t seed2;
};

bool is_dual_vertex_reliable(
        const Point_3& dual_pt, 
        const Point_3& seed_a, 
        const Point_3& seed_b, 
        double threshold_ratio = 50.0) 
    {
        // 计算种子对的间距平方
        double d_seeds_sq = CGAL::to_double(CGAL::squared_distance(seed_a, seed_b));
        if (d_seeds_sq < 1e-12) return false; // 种子重合，异常

        // 计算外心到种子对中点的距离平方
        Point_3 midpoint = CGAL::midpoint(seed_a, seed_b);
        double d_dual_sq = CGAL::to_double(CGAL::squared_distance(dual_pt, midpoint));

        // 如果外心距离 远大于 种子间距 (例如 50 倍)，则认为是 Sliver 造成的飞逸点
        // 阈值 50.0 可以根据实际数据的密度调整，通常 10-100 都是合理的
        return d_dual_sq <= (d_seeds_sq * threshold_ratio * threshold_ratio);
    }

}

void Generator::generate_surface_mesh(MeshingTree* seeds, const char* output_filename)
{
    std::cout << "Generating surface mesh using CGAL Voronoi..." << std::endl;
    
    size_t num_seeds = seeds->get_num_tree_points();
    if (num_seeds == 0) {
        std::cerr << "Error: No seeds provided." << std::endl;
        return;
    }

    // ==================================================================================
    // [新增部分] 0. 导出未删除的种子点 (Export Undeleted Seeds for Debugging)
    // ==================================================================================
    {
        const char* debug_seeds_filename = "debug_active_seeds_for_mesh.obj";
        std::cout << "  * Exporting active seeds to " << debug_seeds_filename << " ..." << std::endl;
        std::ofstream seeds_out(debug_seeds_filename);
        if (seeds_out.is_open()) {
            seeds_out << "# Active seeds used for meshing\n";
            seeds_out << "# Format: v x y z r g b\n";
            
            size_t count = 0;
            for (size_t i = 0; i < num_seeds; i++) {
                if (!seeds->tree_point_is_active(i)) continue;
                
                double* pt = seeds->get_tree_point(i);
                size_t* attrib = seeds->get_tree_point_attrib(i);
                size_t region_id = attrib[5]; // 获取区域ID用于着色

                // 简单着色：1=内(红), 2=外(蓝), 其他=绿/黄
                double r = 0.0, g = 0.0, b = 0.0;
                if (region_id == 1)      { r = 1.0; g = 0.0; b = 0.0; } // Inside: Red
                else if (region_id == 2) { r = 0.0; g = 0.0; b = 1.0; } // Outside: Blue
                else                     { r = 0.0; g = 1.0; b = 0.0; } // Others: Green

                seeds_out << "v " << pt[0] << " " << pt[1] << " " << pt[2] 
                          << " " << r << " " << g << " " << b << "\n";
                count++;
            }
            seeds_out.close();
            std::cout << "  * Exported " << count << " active seeds." << std::endl;
        } else {
            std::cerr << "[Warning] Could not write to " << debug_seeds_filename << std::endl;
        }
    }
    // ==================================================================================

    // 1. Build Delaunay triangulation with seed index info
    Delaunay dt;
    std::vector<Vertex_handle> vertex_handles(num_seeds);
    
    for (size_t i = 0; i < num_seeds; i++) {
        if (!seeds->tree_point_is_active(i)) continue;
        double* pt = seeds->get_tree_point(i);
        Vertex_handle vh = dt.insert(Point_3(pt[0], pt[1], pt[2]));
        vh->info() = i;
        vertex_handles[i] = vh;
    }
    
    std::cout << "  * Delaunay triangulation built with " << dt.number_of_vertices() << " vertices" << std::endl;

    // [保留原有] 输出全局维诺图线框 (Global Voronoi Diagram Wireframe)
    {
        std::cout << "  * Exporting global Voronoi diagram to global_voronoi.obj ..." << std::endl;
        std::ofstream vor_out("global_voronoi.obj");
        if (vor_out.is_open()) {
            vor_out << std::fixed << std::setprecision(16);
            vor_out << "# Global Voronoi Diagram (Finite Edges Only)\n";
            
            size_t edge_v_count = 1;
            for (auto fit = dt.finite_facets_begin(); fit != dt.finite_facets_end(); ++fit) {
                CGAL::Object o = dt.dual(*fit);
                if (const K::Segment_3* s = CGAL::object_cast<K::Segment_3>(&o)) {
                    if (s->squared_length() > 1e12) continue; 
                    vor_out << "v " << s->source().x() << " " << s->source().y() << " " << s->source().z() << "\n";
                    vor_out << "v " << s->target().x() << " " << s->target().y() << " " << s->target().z() << "\n";
                    vor_out << "l " << edge_v_count << " " << edge_v_count + 1 << "\n";
                    edge_v_count += 2;
                }
            }
            vor_out.close();
            std::cout << "  * Done." << std::endl;
        }
    }

    // 2. Collect Voronoi facets for inside/outside seed pairs
    std::vector<VoronoiFacet> voronoi_facets;
    
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
        
        //bool is_pair = (attrib1[1] == seed_idx2) || (attrib2[1] == seed_idx1);
        bool is_pair = attrib1[5] != attrib2[5];
        
        if (!is_pair) continue;
        
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
            voronoi_facets.push_back({facet_vertices, seed_idx1, seed_idx2});
        }
    }
    
    std::cout << "  * Found " << voronoi_facets.size() << " Voronoi facets for seed pairs" << std::endl;

    std::cout << "  * Generated " << voronoi_facets.size() << " facets." << std::endl;

    // 3. Write Surface Mesh OBJ file
    std::ofstream obj_file(output_filename);
    if (!obj_file.is_open()) {
        std::cerr << "Error: Cannot open output file " << output_filename << std::endl;
        return;
    }
    
    obj_file << "# Voronoi surface mesh generated by VoroCrust" << std::endl;
    obj_file << "# Number of facets: " << voronoi_facets.size() << std::endl;
    
    std::map<std::tuple<double, double, double>, size_t> vertex_map;
    std::vector<Point_3> unique_vertices;
    std::vector<std::vector<size_t>> face_indices;
    std::vector<SeedPairKey> face_seed_pairs;
    
    auto get_vertex_index = [&](const Point_3& p) -> size_t {
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
    
    for (const auto& facet : voronoi_facets) {
        if (facet.vertices.size() < 3) continue;
        
        std::vector<size_t> indices;
        for (const auto& pt : facet.vertices) {
            indices.push_back(get_vertex_index(pt));
        }
        
        Point_3 centroid = CGAL::ORIGIN;
        for(const auto& pt : facet.vertices) centroid = centroid + (pt - CGAL::ORIGIN);
        double s = static_cast<double>(facet.vertices.size());
        centroid = Point_3(centroid.x() / s, centroid.y() / s, centroid.z() / s);
        size_t centroid_idx = get_vertex_index(centroid);

        SeedPairKey sp{facet.seed1, facet.seed2};
        if (sp.a > sp.b) std::swap(sp.a, sp.b);

        for (size_t i = 0; i < indices.size(); i++) {
            size_t idx0 = indices[i];
            size_t idx1 = indices[(i + 1) % indices.size()];
            if (idx0 == idx1 || idx0 == centroid_idx || idx1 == centroid_idx) continue;
            face_indices.push_back({centroid_idx, idx0, idx1});
            face_seed_pairs.push_back(sp);
        }
    }

    if (face_indices.size() != face_seed_pairs.size()) {
        std::cerr << "[Warning] face_indices and face_seed_pairs size mismatch." << std::endl;
    }

    // std::unordered_map<EdgeKey, EdgeAdj, EdgeKeyHash> edge_adjacency;
    // edge_adjacency.reserve(face_indices.size() * 3);
    // std::unordered_map<size_t, std::vector<size_t>> boundary_vertex_graph;

    // auto add_edge = [&](size_t u, size_t v, size_t tri_id) {
    //     if (u == v) return;
    //     if (u > v) std::swap(u, v);
    //     EdgeKey k{u, v};
    //     auto& adj = edge_adjacency[k];
    //     if (adj.count == 0) adj.tri0 = tri_id;
    //     else if (adj.count == 1) adj.tri1 = tri_id;
    //     adj.count++;
    // };

    // for (size_t ti = 0; ti < face_indices.size(); ti++) {
    //     const auto& tri = face_indices[ti];
    //     if (tri.size() != 3) continue;
    //     add_edge(tri[0], tri[1], ti);
    //     add_edge(tri[1], tri[2], ti);
    //     add_edge(tri[2], tri[0], ti);
    // }

    // for (const auto& kv : edge_adjacency) {
    //     if (kv.second.count == 1) {
    //         const EdgeKey& e = kv.first;
    //         boundary_vertex_graph[e.a].push_back(e.b);
    //         boundary_vertex_graph[e.b].push_back(e.a);
    //     }
    // }

    // std::unordered_set<size_t> boundary_vertices_visited;
    // boundary_vertices_visited.reserve(boundary_vertex_graph.size() * 2);
    // std::vector<std::vector<size_t>> hole_boundary_components;

    // for (const auto& kv : boundary_vertex_graph) {
    //     size_t start_v = kv.first;
    //     if (boundary_vertices_visited.find(start_v) != boundary_vertices_visited.end()) continue;

    //     std::vector<size_t> component_vertices;
    //     std::queue<size_t> q;
    //     q.push(start_v);
    //     boundary_vertices_visited.insert(start_v);

    //     while (!q.empty()) {
    //         size_t v = q.front();
    //         q.pop();
    //         component_vertices.push_back(v);

    //         auto it = boundary_vertex_graph.find(v);
    //         if (it == boundary_vertex_graph.end()) continue;
    //         for (size_t nb : it->second) {
    //             if (boundary_vertices_visited.insert(nb).second) {
    //                 q.push(nb);
    //             }
    //         }
    //     }

    //     if (!component_vertices.empty()) {
    //         hole_boundary_components.push_back(std::move(component_vertices));
    //     }
    // }

    // auto collect_voronoi_facets_for_seed = [&](size_t seed_idx) {
    //     std::vector<std::vector<Point_3>> facets;
    //     if (seed_idx >= num_seeds) return facets;
    //     if (!seeds->tree_point_is_active(seed_idx)) return facets;

    //     for (auto eit = dt.finite_edges_begin(); eit != dt.finite_edges_end(); ++eit) {
    //         Cell_handle c = eit->first;
    //         int i1 = eit->second;
    //         int i2 = eit->third;

    //         Vertex_handle v1 = c->vertex(i1);
    //         Vertex_handle v2 = c->vertex(i2);

    //         size_t s1 = v1->info();
    //         size_t s2 = v2->info();

    //         if (s1 != seed_idx && s2 != seed_idx) continue;

    //         std::vector<Point_3> facet_vertices;
    //         Delaunay::Cell_circulator cc = dt.incident_cells(*eit);
    //         Delaunay::Cell_circulator done = cc;
    //         if (cc == nullptr) continue;

    //         do {
    //             if (!dt.is_infinite(cc)) {
    //                 facet_vertices.push_back(dt.dual(cc));
    //             }
    //             ++cc;
    //         } while (cc != done);

    //         if (facet_vertices.size() >= 3) {
    //             facets.push_back(std::move(facet_vertices));
    //         }
    //     }
    //     return facets;
    // };

    // auto export_voronoi_cells_obj = [&](
    //     const std::string& filename,
    //     const std::vector<size_t>& seeds_to_export
    // ) {
    //     std::ofstream out(filename);
    //     if (!out.is_open()) {
    //         std::cerr << "[Warning] Cannot open " << filename << " for writing." << std::endl;
    //         return;
    //     }

    //     std::map<std::tuple<double, double, double>, size_t> vmap;
    //     std::vector<Point_3> verts;
    //     std::vector<std::array<size_t, 3>> tris;

    //     auto get_vid = [&](const Point_3& p) -> size_t {
    //         double x = std::round(CGAL::to_double(p.x()) * 1e10) / 1e10;
    //         double y = std::round(CGAL::to_double(p.y()) * 1e10) / 1e10;
    //         double z = std::round(CGAL::to_double(p.z()) * 1e10) / 1e10;
    //         auto key = std::make_tuple(x, y, z);
    //         auto it = vmap.find(key);
    //         if (it != vmap.end()) return it->second;
    //         size_t idx = verts.size();
    //         vmap[key] = idx;
    //         verts.push_back(p);
    //         return idx;
    //     };

    //     for (size_t seed_idx : seeds_to_export) {
    //         if (seed_idx >= num_seeds) continue;
    //         if (!seeds->tree_point_is_active(seed_idx)) continue;

    //         auto cell_facets = collect_voronoi_facets_for_seed(seed_idx);
    //         for (const auto& facet : cell_facets) {
    //             if (facet.size() < 3) continue;

    //             std::vector<size_t> ids;
    //             ids.reserve(facet.size());
    //             for (const auto& p : facet) ids.push_back(get_vid(p));

    //             Point_3 centroid = CGAL::ORIGIN;
    //             for (const auto& p : facet) centroid = centroid + (p - CGAL::ORIGIN);
    //             double s = static_cast<double>(facet.size());
    //             centroid = Point_3(centroid.x() / s, centroid.y() / s, centroid.z() / s);
    //             size_t cid = get_vid(centroid);

    //             for (size_t i = 0; i < ids.size(); i++) {
    //                 size_t a = ids[i];
    //                 size_t b = ids[(i + 1) % ids.size()];
    //                 if (a == b || a == cid || b == cid) continue;
    //                 tris.push_back({cid, a, b});
    //             }
    //         }
    //     }

    //     out << "# Voronoi cells exported for hole-adjacent seed pairs\n";
    //     for (const auto& v : verts) {
    //         out << "v " << std::setprecision(16)
    //             << CGAL::to_double(v.x()) << " "
    //             << CGAL::to_double(v.y()) << " "
    //             << CGAL::to_double(v.z()) << "\n";
    //     }
    //     for (const auto& t : tris) {
    //         out << "f " << (t[0] + 1) << " " << (t[1] + 1) << " " << (t[2] + 1) << "\n";
    //     }
    //     out.close();
    // };

    // auto bbox_overlaps = [&](const CGAL::Bbox_3& a, const CGAL::Bbox_3& b) -> bool {
    //     if (a.xmax() < b.xmin() || b.xmax() < a.xmin()) return false;
    //     if (a.ymax() < b.ymin() || b.ymax() < a.ymin()) return false;
    //     if (a.zmax() < b.zmin() || b.zmax() < a.zmin()) return false;
    //     return true;
    // };

    // auto build_cell_tris_and_bbox = [&](
    //     size_t seed_idx,
    //     std::vector<K::Triangle_3>& tris,
    //     CGAL::Bbox_3& bbox
    // ) -> bool {
    //     tris.clear();
    //     bool have_bbox = false;
    //     CGAL::Bbox_3 local_bbox;

    //     if (seed_idx >= num_seeds) return false;
    //     if (!seeds->tree_point_is_active(seed_idx)) return false;

    //     auto cell_facets = collect_voronoi_facets_for_seed(seed_idx);
    //     for (const auto& facet : cell_facets) {
    //         if (facet.size() < 3) continue;

    //         std::vector<size_t> ids;
    //         ids.reserve(facet.size());
    //         for (const auto& p : facet) {
    //             if (!have_bbox) {
    //                 local_bbox = p.bbox();
    //                 have_bbox = true;
    //             } else {
    //                 local_bbox = local_bbox + p.bbox();
    //             }
    //         }

    //         Point_3 centroid = CGAL::ORIGIN;
    //         for (const auto& p : facet) centroid = centroid + (p - CGAL::ORIGIN);
    //         double s = static_cast<double>(facet.size());
    //         centroid = Point_3(centroid.x() / s, centroid.y() / s, centroid.z() / s);

    //         if (!have_bbox) {
    //             local_bbox = centroid.bbox();
    //             have_bbox = true;
    //         } else {
    //             local_bbox = local_bbox + centroid.bbox();
    //         }

    //         for (size_t i = 0; i < facet.size(); i++) {
    //             const Point_3& a = facet[i];
    //             const Point_3& b = facet[(i + 1) % facet.size()];
    //             if (a == b || a == centroid || b == centroid) continue;
    //             tris.push_back(K::Triangle_3(centroid, a, b));
    //         }
    //     }

    //     if (!have_bbox) return false;
    //     bbox = local_bbox;
    //     return true;
    // };

    // auto voronoi_cells_intersect = [&](size_t seed_a, size_t seed_b) -> bool {
    //     std::vector<K::Triangle_3> tris_a;
    //     std::vector<K::Triangle_3> tris_b;
    //     CGAL::Bbox_3 bbox_a;
    //     CGAL::Bbox_3 bbox_b;

    //     if (!build_cell_tris_and_bbox(seed_a, tris_a, bbox_a)) return false;
    //     if (!build_cell_tris_and_bbox(seed_b, tris_b, bbox_b)) return false;
    //     if (!bbox_overlaps(bbox_a, bbox_b)) return false;

    //     for (const auto& ta : tris_a) {
    //         if (!bbox_overlaps(ta.bbox(), bbox_b)) continue;
    //         for (const auto& tb : tris_b) {
    //             if (!bbox_overlaps(ta.bbox(), tb.bbox())) continue;
    //             if (CGAL::do_intersect(ta, tb)) return true;
    //         }
    //     }
    //     return false;
    // };

    // if (!hole_boundary_components.empty()) {
    //     std::string base(output_filename);
    //     size_t dot = base.find_last_of('.');
    //     if (dot != std::string::npos) base = base.substr(0, dot);

    //     for (size_t hi = 0; hi < hole_boundary_components.size(); hi++) {
    //         const auto& comp_vertices = hole_boundary_components[hi];
    //         std::unordered_set<EdgeKey, EdgeKeyHash> comp_edges;
    //         comp_edges.reserve(comp_vertices.size() * 4);
    //         std::unordered_set<size_t> boundary_tris;
    //         boundary_tris.reserve(comp_vertices.size() * 4);

    //         for (size_t v : comp_vertices) {
    //             auto it = boundary_vertex_graph.find(v);
    //             if (it == boundary_vertex_graph.end()) continue;
    //             for (size_t nb : it->second) {
    //                 size_t a = v, b = nb;
    //                 if (a > b) std::swap(a, b);
    //                 EdgeKey ek{a, b};
    //                 if (!comp_edges.insert(ek).second) continue;
    //                 auto adj_it = edge_adjacency.find(ek);
    //                 if (adj_it != edge_adjacency.end() && adj_it->second.count == 1 && adj_it->second.tri0 != static_cast<size_t>(-1)) {
    //                     boundary_tris.insert(adj_it->second.tri0);
    //                 }
    //             }
    //         }

    //         std::unordered_set<SeedPairKey, SeedPairKeyHash> seed_pairs;
    //         seed_pairs.reserve(boundary_tris.size() * 2);
    //         std::unordered_set<size_t> seeds_to_export_set;
    //         seeds_to_export_set.reserve(boundary_tris.size() * 4);

    //         for (size_t ti : boundary_tris) {
    //             if (ti >= face_seed_pairs.size()) continue;
    //             SeedPairKey sp = face_seed_pairs[ti];
    //             if (sp.a > sp.b) std::swap(sp.a, sp.b);
    //             if (seed_pairs.insert(sp).second) {
    //                 seeds_to_export_set.insert(sp.a);
    //                 seeds_to_export_set.insert(sp.b);
    //             }
    //         }

    //         std::vector<size_t> seeds_to_export;
    //         seeds_to_export.reserve(seeds_to_export_set.size());
    //         for (size_t s : seeds_to_export_set) seeds_to_export.push_back(s);
    //         std::sort(seeds_to_export.begin(), seeds_to_export.end());

    //         {
    //             size_t num_intersect = 0;
    //             size_t num_not_intersect = 0;
    //             for (const auto& sp : seed_pairs) {
    //                 bool is_intersect = voronoi_cells_intersect(sp.a, sp.b);
    //                 std::cout << "  * Hole " << hi << " Voronoi cells intersect? seeds (" << sp.a << ", " << sp.b << ") -> "
    //                           << (is_intersect ? "YES" : "NO") << std::endl;
    //                 if (is_intersect) num_intersect++;
    //                 else num_not_intersect++;
    //             }
    //             std::cout << "  * Hole " << hi << " intersection summary: YES=" << num_intersect
    //                       << ", NO=" << num_not_intersect << std::endl;
    //         }

    //         if (!seeds_to_export.empty()) {
    //             std::string hole_file = base + "_hole_" + std::to_string(hi) + "_cells.obj";
    //             std::cout << "  * Exporting Voronoi cells for hole " << hi << " to " << hole_file << " ..." << std::endl;
    //             export_voronoi_cells_obj(hole_file, seeds_to_export);
    //         }
    //     }
    // }
    
    for (const auto& v : unique_vertices) {
        obj_file << "v " << std::setprecision(16) 
                 << CGAL::to_double(v.x()) << " " 
                 << CGAL::to_double(v.y()) << " " 
                 << CGAL::to_double(v.z()) << std::endl;
    }
    
    for (const auto& face : face_indices) {
        obj_file << "f " << (face[0] + 1) << " " << (face[1] + 1) << " " << (face[2] + 1) << std::endl;
    }
    
    obj_file.close();
    
    std::cout << "  * Surface mesh saved to " << output_filename << std::endl;
    std::cout << "  * Total vertices: " << unique_vertices.size() << ", triangles: " << face_indices.size() << std::endl;
}

// #include "Generator.h"
// #include <map>
// #include <fstream>
// #include <iomanip>

// #include <CGAL/Exact_predicates_exact_constructions_kernel.h>
// #include <CGAL/Delaunay_triangulation_3.h>
// #include <CGAL/Triangulation_vertex_base_with_info_3.h>
// #include <CGAL/Delaunay_triangulation_cell_base_3.h> // 必须包含

// typedef CGAL::Exact_predicates_exact_constructions_kernel K;
// typedef CGAL::Triangulation_vertex_base_with_info_3<size_t, K> Vb;
// typedef CGAL::Delaunay_triangulation_cell_base_3<K> Cb; // 定义 Cell Base
// typedef CGAL::Triangulation_data_structure_3<Vb, Cb> Tds; // 传入 Vb 和 Cb
// typedef CGAL::Delaunay_triangulation_3<K, Tds> Delaunay;
// typedef Delaunay::Point Point_3;
// typedef Delaunay::Vertex_handle Vertex_handle;
// typedef Delaunay::Cell_handle Cell_handle;
// typedef Delaunay::Edge Edge;
// typedef K::Point_3 Point;
// void Generator::generate_surface_mesh(MeshingTree* seeds, const char* output_filename)
// {
//     std::cout << "Generating surface mesh using CGAL Voronoi..." << std::endl;
    
//     size_t num_seeds = seeds->get_num_tree_points();
//     if (num_seeds == 0) {
//         std::cerr << "Error: No seeds provided." << std::endl;
//         return;
//     }

//     // 1. Build Delaunay triangulation with seed index info
//     Delaunay dt;
//     std::vector<Vertex_handle> vertex_handles(num_seeds);
    
//     for (size_t i = 0; i < num_seeds; i++) {
//         if (!seeds->tree_point_is_active(i)) continue;
//         double* pt = seeds->get_tree_point(i);
//         Vertex_handle vh = dt.insert(Point_3(pt[0], pt[1], pt[2]));
//         vh->info() = i;
//         vertex_handles[i] = vh;
//     }
    
//     std::cout << "  * Delaunay triangulation built with " << dt.number_of_vertices() << " vertices" << std::endl;




//     {
//         std::cout << "  * Exporting inside/outside seed points to OBJ files..." << std::endl;

//         std::ofstream inside_seeds("inside_seeds.obj");
//         std::ofstream outside_seeds("outside_seeds.obj");

//         if (!inside_seeds.is_open() || !outside_seeds.is_open()) {
//             std::cerr << "Error: Cannot create seed point output files" << std::endl;
//         }
//         else {
//             inside_seeds << std::fixed << std::setprecision(16);
//             outside_seeds << std::fixed << std::setprecision(16);

//             inside_seeds << "# Inside Seed Points\n";
//             outside_seeds << "# Outside Seed Points\n";

//             size_t inside_count = 0;
//             size_t outside_count = 0;

//             for (size_t i = 0; i < num_seeds; i++) {
//                 if (!seeds->tree_point_is_active(i)) continue;

//                 double* pt = seeds->get_tree_point(i);
//                 size_t* attrib = seeds->get_tree_point_attrib(i);

//                 // 根据区域ID判断内外
//                 // attrib[5] 是区域ID: 0 = 外部, 非0 = 内部
//                 if (attrib[5] == 1) {
//                     // 外部种子点
//                     outside_seeds << "v " << pt[0] << " " << pt[1] << " " << pt[2] << "\n";
//                     outside_count++;
//                 }
//                 else {
//                     // 内部种子点
//                     inside_seeds << "v " << pt[0] << " " << pt[1] << " " << pt[2] << "\n";
//                     inside_count++;
//                 }
//             }

//             inside_seeds.close();
//             outside_seeds.close();

//             std::cout << "  * Inside seeds: " << inside_count << " (saved to inside_seeds.obj)" << std::endl;
//             std::cout << "  * Outside seeds: " << outside_count << " (saved to outside_seeds.obj)" << std::endl;
//         }
//     }

//     // ==================================================================================
//     // [新增部分] 输出全局维诺图 (Global Voronoi Diagram Wireframe)
//     // ==================================================================================
//     {
//         std::cout << "  * Exporting global Voronoi diagram to global_voronoi.obj ..." << std::endl;
//         std::ofstream vor_out("global_voronoi.obj");
//         if (vor_out.is_open()) {
//             vor_out << std::fixed << std::setprecision(16);
//             vor_out << "# Global Voronoi Diagram (Finite Edges Only)\n";
            
//             size_t edge_v_count = 1;
            
//             // 遍历所有有限的 Delaunay 面 (Finite Facets)
//             // 在 3D 中，Delaunay 面 <---> Voronoi 边
//             for (auto fit = dt.finite_facets_begin(); fit != dt.finite_facets_end(); ++fit) {
//                 // 计算对偶对象
//                 CGAL::Object o = dt.dual(*fit);
                
//                 // 我们只输出有限的线段 (Segment)，忽略射线 (Ray)
//                 if (const K::Segment_3* s = CGAL::object_cast<K::Segment_3>(&o)) {
                    
//                     // (可选) 过滤掉过长的线段，防止某些退化四面体导致的无穷远点干扰视图
//                     if (s->squared_length() > 1e12) continue; 

//                     vor_out << "v " << s->source().x() << " " << s->source().y() << " " << s->source().z() << "\n";
//                     vor_out << "v " << s->target().x() << " " << s->target().y() << " " << s->target().z() << "\n";
//                     vor_out << "l " << edge_v_count << " " << edge_v_count + 1 << "\n";
//                     edge_v_count += 2;
//                 }
//             }
//             vor_out.close();
//             std::cout << "  * Done." << std::endl;
//         } else {
//             std::cerr << "Error: Cannot write to global_voronoi.obj" << std::endl;
//         }
//     }
//     // ==================================================================================
//     {
//     std::cout << "  * Exporting seed pair connections to seed_pairs.obj ..." << std::endl;
    
//     std::ofstream pairs_out("seed_pairs.obj");
//     if (pairs_out.is_open()) {
//         pairs_out << std::fixed << std::setprecision(16);
//         pairs_out << "# Seed Pair Connections\n";
        
//         // 首先输出所有种子点作为顶点
//         for (size_t i = 0; i < num_seeds; i++) {
//             if (!seeds->tree_point_is_active(i)) continue;
//             double* pt = seeds->get_tree_point(i);
//             pairs_out << "v " << pt[0] << " " << pt[1] << " " << pt[2] << "\n";
//         }
        
//         // 然后输出配对连接线
//         size_t pair_count = 0;
//         std::set<std::pair<size_t, size_t>> processed_pairs;
        
//         for (size_t i = 0; i < num_seeds; i++) {
//             if (!seeds->tree_point_is_active(i)) continue;
            
//             size_t* attrib = seeds->get_tree_point_attrib(i);
//             size_t pair_idx = attrib[1]; // attrib[1] 是配对种子索引
            
//             if (pair_idx < num_seeds && seeds->tree_point_is_active(pair_idx)) {
//                 // 避免重复输出同一对
//                 auto pair_key = std::minmax(i, pair_idx);
//                 if (processed_pairs.find(pair_key) == processed_pairs.end()) {
//                     // OBJ使用1-based索引
//                     pairs_out << "l " << (i + 1) << " " << (pair_idx + 1) << "\n";
//                     processed_pairs.insert(pair_key);
//                     pair_count++;
//                 }
//             }
//         }
        
//         pairs_out.close();
//         std::cout << "  * Seed pairs: " << pair_count << " connections (saved to seed_pairs.obj)" << std::endl;
//     } else {
//         std::cerr << "Error: Cannot write to seed_pairs.obj" << std::endl;
//     }
// }
//     // 2. Collect Voronoi facets for inside/outside seed pairs
//     std::vector<std::vector<Point_3>> voronoi_facets;
    
//     for (auto eit = dt.finite_edges_begin(); eit != dt.finite_edges_end(); ++eit) {
//         Cell_handle c = eit->first;
//         int i1 = eit->second;
//         int i2 = eit->third;
        
//         Vertex_handle v1 = c->vertex(i1);
//         Vertex_handle v2 = c->vertex(i2);
        
//         size_t seed_idx1 = v1->info();
//         size_t seed_idx2 = v2->info();
        
//         // Check if this edge connects a seed pair (inside/outside)
//         size_t* attrib1 = seeds->get_tree_point_attrib(seed_idx1);
//         size_t* attrib2 = seeds->get_tree_point_attrib(seed_idx2);
        
//         // attrib[1] is the pair seed index
//         // 这里使用之前生成的配对信息，或者可以改用区域 ID (attrib[5]) 判断: 
//          bool is_interface = (attrib1[5] != attrib2[5]) && (attrib1[5] != 0) && (attrib2[5] != 0);
//         //bool is_pair = (attrib1[1] == seed_idx2) || (attrib2[1] == seed_idx1);
//          if (!is_interface) continue;
//       //  if (!is_pair) continue;
//         int label1 = seed_idx1;
//         int label2 = seed_idx2;
// 		std::cout << "seed pair labels: " << attrib1[5] << " " << attrib2[5] << std::endl;
//         // This edge connects a seed pair - extract the dual Voronoi facet
//         // The Voronoi facet is the convex polygon formed by circumcenters of 
//         // cells incident to this edge
//         std::vector<Point_3> facet_vertices;
        
//         Delaunay::Cell_circulator cc = dt.incident_cells(*eit);
//         Delaunay::Cell_circulator done = cc;
        
//         if (cc == nullptr) continue;
        
//         do {
//             if (!dt.is_infinite(cc)) {
//                 // Compute circumcenter manually or use CGAL dual
//                 // Using dual() is safer and usually cached/optimized in some kernels, 
//                 // but manual construction as you did is also fine for Exact kernel.
//                 Point_3 center = dt.dual(cc); 
//                 facet_vertices.push_back(center);
//             }
//             ++cc;
//         } while (cc != done);
        
//         if (facet_vertices.size() >= 3) {
//             voronoi_facets.push_back(facet_vertices);
//         }
//     }
    
//     std::cout << "  * Found " << voronoi_facets.size() << " Voronoi facets for seed pairs" << std::endl;

//     // 3. Write Surface Mesh OBJ file
//     std::ofstream obj_file(output_filename);
//     if (!obj_file.is_open()) {
//         std::cerr << "Error: Cannot open output file " << output_filename << std::endl;
//         return;
//     }
    
//     obj_file << "# Voronoi surface mesh generated by VoroCrust" << std::endl;
//     obj_file << "# Number of facets: " << voronoi_facets.size() << std::endl;
    
//     // Use a map to avoid duplicate vertices
//     // Note: Using Exact kernel points as map keys might be slow or tricky with rounding.
//     // Your rounding strategy is good for merging close points.
//     std::map<std::tuple<double, double, double>, size_t> vertex_map;
//     std::vector<Point_3> unique_vertices;
//     std::vector<std::vector<size_t>> face_indices;
    
//     auto get_vertex_index = [&](const Point_3& p) -> size_t {
//         // Round to avoid floating point precision issues
//         double x = std::round(CGAL::to_double(p.x()) * 1e10) / 1e10;
//         double y = std::round(CGAL::to_double(p.y()) * 1e10) / 1e10;
//         double z = std::round(CGAL::to_double(p.z()) * 1e10) / 1e10;
//         auto key = std::make_tuple(x, y, z);
        
//         auto it = vertex_map.find(key);
//         if (it != vertex_map.end()) {
//             return it->second;
//         }
//         size_t idx = unique_vertices.size();
//         vertex_map[key] = idx;
//         unique_vertices.push_back(p);
//         return idx;
//     };
    
//     // Process facets and triangulate polygons
//     for (const auto& facet : voronoi_facets) {
//         if (facet.size() < 3) continue;
        
//         std::vector<size_t> indices;
//         for (const auto& pt : facet) {
//             indices.push_back(get_vertex_index(pt));
//         }
        
//         // Fan triangulation for polygon with more than 3 vertices
//         // Improve triangulation: use centroid fan to create better triangles
//         Point_3 centroid = CGAL::ORIGIN;
//         for(const auto& pt : facet) centroid = centroid + (pt - CGAL::ORIGIN);
//         double s = static_cast<double>(facet.size());
// 		centroid = Point_3(centroid.x() / s, centroid.y() / s, centroid.z() / s);
//         size_t centroid_idx = get_vertex_index(centroid);

//         for (size_t i = 0; i < indices.size(); i++) {
//             size_t idx0 = indices[i];
//             size_t idx1 = indices[(i + 1) % indices.size()];
            
//             // 跳过退化三角形
//             if (idx0 == idx1 || idx0 == centroid_idx || idx1 == centroid_idx) continue;

//             face_indices.push_back({centroid_idx, idx0, idx1});
//         }
//     }
    
//     // Write vertices
//     for (const auto& v : unique_vertices) {
//         obj_file << "v " << std::setprecision(16) 
//                  << CGAL::to_double(v.x()) << " " 
//                  << CGAL::to_double(v.y()) << " " 
//                  << CGAL::to_double(v.z()) << std::endl;
//     }
    
//     // Write faces (OBJ uses 1-based indexing)
//     for (const auto& face : face_indices) {
//         obj_file << "f " << (face[0] + 1) << " " << (face[1] + 1) << " " << (face[2] + 1) << std::endl;
//     }
    
//     obj_file.close();
    
//     std::cout << "  * Surface mesh saved to " << output_filename << std::endl;
//     std::cout << "  * Total vertices: " << unique_vertices.size() << ", triangles: " << face_indices.size() << std::endl;
// }struct LabeledPoint {
//     Point point;
//     int label;

//     LabeledPoint(double x, double y, double z, int l)
//         : point(x, y, z), label(l) {
//     }
// };


// // Get Voronoi face vertices for a Delaunay edge
// std::vector<Point> getVoronoiFace(const Delaunay& dt,
//     Vertex_handle v1,
//     Vertex_handle v2) {
//     std::vector<Point> face_vertices;

//     // Get all cells incident to both vertices
//     std::vector<Cell_handle> cells;
//     dt.incident_cells(v1, std::back_inserter(cells));

//     std::vector<Cell_handle> common_cells;
//     for (const auto& cell : cells) {
//         if (!dt.is_infinite(cell)) {
//             for (int i = 0; i < 4; ++i) {
//                 if (cell->vertex(i) == v2) {
//                     common_cells.push_back(cell);
//                     break;
//                 }
//             }
//         }
//     }

//     // Get circumcenters of common cells (Voronoi vertices)
//     for (const auto& cell : common_cells) {
//         face_vertices.push_back(dt.dual(cell));
//     }

//     return face_vertices;
// }

// // Export boundary surface as triangulated mesh
// void exportBoundarySurface(const Delaunay& dt,
//     const std::map<Vertex_handle, int>& vertex_labels,
//     const std::string& filename) {
//     std::ofstream out(filename);
//     if (!out) {
//         std::cerr << "Error: Cannot create file " << filename << std::endl;
//         return;
//     }

//     out << "# Voronoi boundary surface mesh" << std::endl;
//     out << "# Triangulated surface between different labels" << std::endl;

//     std::cout << "\nExtracting boundary surface..." << std::endl;

//     // Find all edges between vertices with different labels
//     std::set<std::pair<Vertex_handle, Vertex_handle>> boundary_edges;

//     for (auto eit = dt.finite_edges_begin(); eit != dt.finite_edges_end(); ++eit) {
//         Cell_handle c = eit->first;
//         int i = eit->second;
//         int j = eit->third;

//         Vertex_handle v1 = c->vertex(i);
//         Vertex_handle v2 = c->vertex(j);

//         if (vertex_labels.find(v1) != vertex_labels.end() &&
//             vertex_labels.find(v2) != vertex_labels.end()) {
//             if (vertex_labels.at(v1) != vertex_labels.at(v2)) {
//                 // Normalize edge (smaller pointer first)
//                 if (v1 < v2) {
//                     boundary_edges.insert({ v1, v2 });
//                 }
//                 else {
//                     boundary_edges.insert({ v2, v1 });
//                 }
//             }
//         }
//     }

//     std::cout << "Found " << boundary_edges.size() << " boundary edges" << std::endl;

//     // Collect all Voronoi vertices and faces
//     std::map<Point, int> point_indices;
//     int vertex_count = 0;
//     std::vector<std::vector<int>> faces;

//     for (const auto& edge_pair : boundary_edges) {
//         Vertex_handle v1 = edge_pair.first;
//         Vertex_handle v2 = edge_pair.second;

//         std::vector<Point> face_points = getVoronoiFace(dt, v1, v2);

//         if (face_points.size() >= 3) {
//             std::vector<int> face_indices;
//             for (const auto& p : face_points) {
//                 if (point_indices.find(p) == point_indices.end()) {
//                     point_indices[p] = ++vertex_count;
//                 }
//                 face_indices.push_back(point_indices[p]);
//             }
//             faces.push_back(face_indices);
//         }
//     }

//     // Write vertices
//     std::vector<Point> ordered_points(vertex_count);
//     for (const auto& pair : point_indices) {
//         ordered_points[pair.second - 1] = pair.first;
//     }

//     for (const auto& p : ordered_points) {
//         out << "v " << p.x() << " " << p.y() << " " << p.z() << std::endl;
//     }

//     // Triangulate and write faces
//     int triangle_count = 0;
//     for (const auto& face : faces) {
//         if (face.size() >= 3) {
//             // Simple fan triangulation from first vertex
//             for (size_t i = 1; i < face.size() - 1; ++i) {
//                 out << "f " << face[0] << " " << face[i] << " " << face[i + 1] << std::endl;
//                 triangle_count++;
//             }
//         }
//     }

//     out.close();

//     std::cout << "Surface mesh exported:" << std::endl;
//     std::cout << "  Vertices: " << vertex_count << std::endl;
//     std::cout << "  Voronoi faces: " << faces.size() << std::endl;
//     std::cout << "  Triangles: " << triangle_count << std::endl;
//     std::cout << "  Saved to: " << filename << std::endl;
// }

// // Export Delaunay facets as triangulated surface
// void exportDelaunayBoundarySurface(const Delaunay& dt,
//     const std::map<Vertex_handle, int>& vertex_labels,
//     const std::string& filename) {
//     std::ofstream out(filename);
//     if (!out) {
//         std::cerr << "Error: Cannot create file " << filename << std::endl;
//         return;
//     }

//     out << "# Delaunay-based boundary surface" << std::endl;

//     std::cout << "\nExtracting Delaunay boundary surface..." << std::endl;

//     // Find boundary facets (Delaunay triangles with different labels on each side)
//     std::set<Point> unique_points;
//     std::vector<std::array<Point, 3>> triangles;

//     for (auto fit = dt.finite_facets_begin(); fit != dt.finite_facets_end(); ++fit) {
//         Cell_handle c1 = fit->first;
//         int idx = fit->second;
//         Cell_handle c2 = c1->neighbor(idx);

//         if (dt.is_infinite(c1) || dt.is_infinite(c2)) continue;

//         // Check if this facet separates different labels
//         std::set<int> labels;
//         for (int i = 0; i < 4; ++i) {
//             Vertex_handle v = c1->vertex(i);
//             if (vertex_labels.find(v) != vertex_labels.end()) {
//                 labels.insert(vertex_labels.at(v));
//             }
//         }
//         for (int i = 0; i < 4; ++i) {
//             Vertex_handle v = c2->vertex(i);
//             if (vertex_labels.find(v) != vertex_labels.end()) {
//                 labels.insert(vertex_labels.at(v));
//             }
//         }

//         if (labels.size() >= 2) {
//             // This is a boundary facet - get its 3 vertices
//             std::array<Point, 3> triangle;
//             int tri_idx = 0;
//             for (int i = 0; i < 4; ++i) {
//                 if (i != idx) {
//                     triangle[tri_idx++] = c1->vertex(i)->point();
//                 }
//             }

//             triangles.push_back(triangle);
//             for (const auto& p : triangle) {
//                 unique_points.insert(p);
//             }
//         }
//     }

//     std::cout << "Found " << triangles.size() << " boundary triangles" << std::endl;
//     std::cout << "With " << unique_points.size() << " unique vertices" << std::endl;

//     // Map points to indices
//     std::map<Point, int> point_to_index;
//     int idx = 1;
//     for (const auto& p : unique_points) {
//         point_to_index[p] = idx++;
//     }

//     // Write vertices
//     for (const auto& p : unique_points) {
//         out << "v " << p.x() << " " << p.y() << " " << p.z() << std::endl;
//     }

//     // Write faces
//     for (const auto& tri : triangles) {
//         out << "f " << point_to_index[tri[0]]
//             << " " << point_to_index[tri[1]]
//             << " " << point_to_index[tri[2]] << std::endl;
//     }

//     out.close();
//     std::cout << "Delaunay surface saved to: " << filename << std::endl;
// }
// // 改进版：带详细统计和验证的边界面提取
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

//     // 收集围绕该边的所有四面体的外接球心
//     std::vector<Point_3> polygon_vertices;

//     // 需要从边构造edge
//     // 在CGAL中，edge是通过(cell, i, j)三元组表示的
//     // 这里我们用另一种方法：找到包含v1和v2的所有cell

//     std::vector<Cell_handle> cells_v1;
//     dt.incident_cells(v1, std::back_inserter(cells_v1));

//     for (auto cell : cells_v1) {
//         if (dt.is_infinite(cell)) continue;

//         // 检查这个cell是否也包含v2
//         bool contains_v2 = false;
//         for (int i = 0; i < 4; i++) {
//             if (cell->vertex(i) == v2) {
//                 contains_v2 = true;
//                 break;
//             }
//         }

//         if (contains_v2) {
//             // 这个cell同时包含v1和v2，其外接球心是多边形的一个顶点
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

//     // 输出多边形顶点
//     for (size_t i = 0; i < polygon_vertices.size(); i++) {
//         const auto& v = polygon_vertices[i];
//         obj_file << "v " << CGAL::to_double(v.x()) << " "
//             << CGAL::to_double(v.y()) << " "
//             << CGAL::to_double(v.z()) << "\n";
//     }

//     // 输出多边形的边（线段）
//     obj_file << "\n# Polygon edges\n";
//     for (size_t i = 0; i < polygon_vertices.size(); i++) {
//         size_t next_i = (i + 1) % polygon_vertices.size();
//         obj_file << "l " << (i + 1) << " " << (next_i + 1) << "\n";
//     }

//     // 可选：也输出多边形的面（如果顶点共面）
//     obj_file << "\n# Polygon face\n";
//     obj_file << "f";
//     for (size_t i = 0; i < polygon_vertices.size(); i++) {
//         obj_file << " " << (i + 1);
//     }
//     obj_file << "\n";

//     // 输出对应的种子点位置（用于参考）
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
// void write_voronoi_facets_to_obj_dedup(
//     const std::string& filename,
//     const std::vector<std::vector<Point_3>>& voronoi_facets
// ) {
//     std::ofstream out(filename);
//     if (!out) {
//         std::cerr << "Cannot open file " << filename << std::endl;
//         return;
//     }

//     out << "# Voronoi facets with deduplicated vertices\n";

//     std::map<Point_3, int> vertex_index;
//     std::vector<Point_3> vertices;
//     std::vector<std::vector<int>> faces;

//     // ---------- build vertex list + face indices ----------
//     for (const auto& facet : voronoi_facets) {
//         if (facet.size() < 3) continue;

//         std::vector<int> face;
//         for (const auto& p : facet) {
//             auto it = vertex_index.find(p);
//             if (it == vertex_index.end()) {
//                 int idx = static_cast<int>(vertices.size()) + 1; // OBJ: 1-based
//                 vertex_index[p] = idx;
//                 vertices.push_back(p);
//                 face.push_back(idx);
//             }
//             else {
//                 face.push_back(it->second);
//             }
//         }
//         faces.push_back(face);
//     }

//     // ---------- write vertices ----------
//     for (const auto& p : vertices) {
//         out << "v "
//             << p.x() << " "
//             << p.y() << " "
//             << p.z() << "\n";
//     }

//     // ---------- write faces ----------
//     for (const auto& f : faces) {
//         out << "f";
//         for (int idx : f) {
//             out << " " << idx;
//         }
//         out << "\n";
//     }

//     out.close();
// }



// void write_voronoi_facets_triangulated_obj(
//     const std::string& filename,
//     const std::vector<std::vector<Point_3>>& voronoi_facets
// ) {
//     std::ofstream out(filename);
//     if (!out) {
//         std::cerr << "Cannot open file " << filename << std::endl;
//         return;
//     }

//     out << "# Triangulated Voronoi facets\n";

//     std::map<Point_3, int> vertex_index;
//     std::vector<Point_3> vertices;
//     std::vector<std::array<int, 3>> triangles;

//     auto get_index = [&](const Point_3& p) {
//         auto it = vertex_index.find(p);
//         if (it != vertex_index.end())
//             return it->second;

//         int idx = static_cast<int>(vertices.size()) + 1; // OBJ 1-based
//         vertex_index[p] = idx;
//         vertices.push_back(p);
//         return idx;
//         };

//     // ---------- triangulate each facet ----------
//     for (const auto& facet : voronoi_facets) {
//         if (facet.size() < 3) continue;

//         int v0 = get_index(facet[0]);

//         for (size_t i = 1; i + 1 < facet.size(); ++i) {
//             int v1 = get_index(facet[i]);
//             int v2 = get_index(facet[i + 1]);

//             triangles.push_back({ v0, v1, v2 });
//         }
//     }

//     // ---------- write vertices ----------
//     for (const auto& p : vertices) {
//         out << "v "
//             << p.x() << " "
//             << p.y() << " "
//             << p.z() << "\n";
//     }

//     // ---------- write triangles ----------
//     for (const auto& t : triangles) {
//         out << "f "
//             << t[0] << " "
//             << t[1] << " "
//             << t[2] << "\n";
//     }

//     out.close();
// }
// void write_voronoi_facets_to_obj(
//     const std::string& filename,
//     const std::vector<std::vector<Point_3>>& voronoi_facets
// ) {
//     std::ofstream out(filename);
//     if (!out) {
//         std::cerr << "Cannot open file " << filename << std::endl;
//         return;
//     }

//     out << "# Voronoi facets OBJ\n";

//     size_t vertex_offset = 1; // OBJ index starts from 1

//     for (const auto& facet : voronoi_facets) {
//         if (facet.size() < 3) continue;

//         // 1. write vertices
//         for (const auto& p : facet) {
//             out << "v "
//                 << p.x() << " "
//                 << p.y() << " "
//                 << p.z() << "\n";
//         }

//         // 2. write face
//         out << "f";
//         for (size_t i = 0; i < facet.size(); ++i) {
//             out << " " << (vertex_offset + i);
//         }
//         out << "\n";

//         vertex_offset += facet.size();
//     }

//     out.close();
// }
// void exportBoundarySurfaceImproved(const Delaunay& dt,
//     const std::map<Vertex_handle, int>& vertex_labels,
//     const std::string& filename) {
//     std::ofstream out(filename);
//     if (!out) {
//         std::cerr << "Error: Cannot create file " << filename << std::endl;
//         return;
//     }

//     out << "# Voronoi boundary surface (improved with validation)" << std::endl;
//     std::cout << "\n=== Extracting Boundary Surface (Improved) ===" << std::endl;

//     // ========================================
//     // 统计：帮助理解哪些面被保留/丢弃
//     // ========================================
//     int total_edges = 0;
//     int boundary_edges_count = 0;  // label 1-2
//     int inner_edges = 0;            // label 1-1
//     int outer_edges = 0;            // label 2-2
//     int no_label_edges = 0;         // 没有标签的边

//     std::set<std::pair<Vertex_handle, Vertex_handle>> boundary_edges;

//     // ========================================
//     // 步骤1：遍历所有Delaunay边
//     // 每条边对应一个Voronoi面
//     // ========================================
//     std::vector<std::vector<Point_3>> voronoi_facets;
//     int fla = 0;
//     for (auto eit = dt.finite_edges_begin(); eit != dt.finite_edges_end(); ++eit) {
//         total_edges++;

//         Cell_handle c = eit->first;
//         int i = eit->second;
//         int j = eit->third;

//         // 获取边的两个端点（站点）
//         Vertex_handle v1 = c->vertex(i);
//         Vertex_handle v2 = c->vertex(j);

//         // 检查是否都有标签
//         if (vertex_labels.find(v1) == vertex_labels.end() ||
//             vertex_labels.find(v2) == vertex_labels.end()) {
//             no_label_edges++;
//             continue;
//         }

//         int label1 = vertex_labels.at(v1);
//         int label2 = vertex_labels.at(v2);

//         // ========================================
//         // 核心判断：保留条件
//         // 
//         // Voronoi面由2个站点生成
//         // 保留条件：这2个站点的label不同
//         // ========================================
//         if (label1 != label2) {
        
//             std::vector<Point_3> facet_vertices;

//             Delaunay::Cell_circulator cc = dt.incident_cells(*eit);
//             Delaunay::Cell_circulator done = cc;

//             if (cc == nullptr) continue;

//             do {
//                 if (!dt.is_infinite(cc)) {
//                     // Compute circumcenter manually or use CGAL dual
//                     // Using dual() is safer and usually cached/optimized in some kernels, 
//                     // but manual construction as you did is also fine for Exact kernel.
//                     Point_3 center = dt.dual(cc);
//                     facet_vertices.push_back(center);
//                 }
//                 ++cc;
//             } while (cc != done);

//             if (facet_vertices.size() >= 3) {
//               /*  if (fla == 0)
//                 {
//                     export_single_voronoi_polygon(dt, v1, v2, "single_voronoi_polygon.obj");
//                     fla = 1;
//                 }*/
//                 voronoi_facets.push_back(facet_vertices);
//             }
       
//             // ✅ 边界面：一个内点 + 一个外点
//             boundary_edges_count++;

//             // 规范化边（保证唯一性）
//             if (v1 < v2) {
//                 boundary_edges.insert({ v1, v2 });
//             }
//             else {
//                 boundary_edges.insert({ v2, v1 });
//             }
//         }
//         else {
//             // ❌ 非边界面
//             if (label1 == 1) {
//                 inner_edges++;  // 内-内
//             }
//             else if (label1 == 2) {
//                 outer_edges++;  // 外-外
//             }
//         }
//     }

//     write_voronoi_facets_to_obj("voronoi.obj", voronoi_facets);
//     write_voronoi_facets_to_obj_dedup("voronoi_dedup.obj", voronoi_facets);
//     write_voronoi_facets_triangulated_obj("voronoi_dedup_triangulated.obj", voronoi_facets);
// std::cout << "  * Found " << voronoi_facets.size() << " Voronoi facets for seed pairs" << std::endl;

// // 3. Write Surface Mesh OBJ file
// std::ofstream obj_file("output_filename.obj");
// if (!obj_file.is_open()) {
//     //  std::cerr << "Error: Cannot open output file " << output_filename << std::endl;
//     return;
// }

// obj_file << "# Voronoi surface mesh generated by VoroCrust" << std::endl;
// obj_file << "# Number of facets: " << voronoi_facets.size() << std::endl;

// // Use a map to avoid duplicate vertices
// // Note: Using Exact kernel points as map keys might be slow or tricky with rounding.
// // Your rounding strategy is good for merging close points.
// std::map<std::tuple<double, double, double>, size_t> vertex_map;
// std::vector<Point_3> unique_vertices;
// std::vector<std::vector<size_t>> face_indices;

// auto get_vertex_index = [&](const Point_3& p) -> size_t {
//     // Round to avoid floating point precision issues
//     double x = std::round(CGAL::to_double(p.x()) * 1e10) / 1e10;
//     double y = std::round(CGAL::to_double(p.y()) * 1e10) / 1e10;
//     double z = std::round(CGAL::to_double(p.z()) * 1e10) / 1e10;
//     auto key = std::make_tuple(x, y, z);

//     auto it = vertex_map.find(key);
//     if (it != vertex_map.end()) {
//         return it->second;
//     }
//     size_t idx = unique_vertices.size();
//     vertex_map[key] = idx;
//     unique_vertices.push_back(p);
//     return idx;
//     };

// // Process facets and triangulate polygons
// for (const auto& facet : voronoi_facets) {
//     if (facet.size() < 3) continue;

//     std::vector<size_t> indices;
//     for (const auto& pt : facet) {
//         indices.push_back(get_vertex_index(pt));
//     }

//     // Fan triangulation for polygon with more than 3 vertices
//     // Improve triangulation: use centroid fan to create better triangles
//     Point_3 centroid = CGAL::ORIGIN;
//     for (const auto& pt : facet) centroid = centroid + (pt - CGAL::ORIGIN);
//     double s = static_cast<double>(facet.size());
//     centroid = Point_3(centroid.x() / s, centroid.y() / s, centroid.z() / s);
//     size_t centroid_idx = get_vertex_index(centroid);

//     for (size_t i = 0; i < indices.size(); i++) {
//         size_t idx0 = indices[i];
//         size_t idx1 = indices[(i + 1) % indices.size()];

//         // 跳过退化三角形
//         if (idx0 == idx1 || idx0 == centroid_idx || idx1 == centroid_idx) continue;

//         face_indices.push_back({ centroid_idx, idx0, idx1 });
//     }
// }

// // Write vertices
// for (const auto& v : unique_vertices) {
//     obj_file << "v " << std::setprecision(16)
//         << CGAL::to_double(v.x()) << " "
//         << CGAL::to_double(v.y()) << " "
//         << CGAL::to_double(v.z()) << std::endl;
// }

// // Write faces (OBJ uses 1-based indexing)
// for (const auto& face : face_indices) {
//     obj_file << "f " << (face[0] + 1) << " " << (face[1] + 1) << " " << (face[2] + 1) << std::endl;
// }

// obj_file.close();

// //  std::cout << "  * Surface mesh saved to " << output_filename << std::endl;
// std::cout << "  * Total vertices: " << unique_vertices.size() << ", triangles: " << face_indices.size() << std::endl;
//     // ========================================
//     // 打印统计：验证结果是否合理
//     // ========================================
//     std::cout << "\n--- Edge Classification ---" << std::endl;
//     std::cout << "Total Delaunay edges: " << total_edges << std::endl;
//     std::cout << "├─ Boundary edges (1-2): " << boundary_edges_count
//         << " (" << (100.0 * boundary_edges_count / total_edges) << "%)" << std::endl;
//     std::cout << "├─ Inner edges (1-1): " << inner_edges
//         << " (" << (100.0 * inner_edges / total_edges) << "%)" << std::endl;
//     std::cout << "├─ Outer edges (2-2): " << outer_edges
//         << " (" << (100.0 * outer_edges / total_edges) << "%)" << std::endl;
//     std::cout << "└─ No label edges: " << no_label_edges << std::endl;





//     // ========================================
//     // 步骤2：为每个边界边提取Voronoi面
//     // ========================================
//     std::map<Point, int> point_indices;
//     int vertex_count = 0;
//     std::vector<std::vector<int>> faces;

//     int valid_faces = 0;
//     int degenerate_faces = 0;  // 退化的面（顶点<3或重复）

//     std::cout << "\n--- Extracting Voronoi Faces ---" << std::endl;

//     for (const auto& edge_pair : boundary_edges) {
//         Vertex_handle v1 = edge_pair.first;
//         Vertex_handle v2 = edge_pair.second;

//         // 获取这条边对应的Voronoi面的顶点
//         std::vector<Point> face_points = getVoronoiFace(dt, v1, v2);

//         // 过滤退化的面
//         if (face_points.size() < 3) {
//             degenerate_faces++;
//             continue;
//         }

//         // 检查是否有重复顶点
//         std::set<Point> unique_points(face_points.begin(), face_points.end());
//         if (unique_points.size() < 3) {
//             degenerate_faces++;
//             continue;
//         }

//         valid_faces++;

//         // 收集顶点索引
//         std::vector<int> face_indices;
//         for (const auto& p : face_points) {
//             if (point_indices.find(p) == point_indices.end()) {
//                 point_indices[p] = ++vertex_count;
//             }
//             face_indices.push_back(point_indices[p]);
//         }
//         faces.push_back(face_indices);
//     }

//     std::cout << "Valid faces: " << valid_faces << std::endl;
//     std::cout << "Degenerate faces (filtered out): " << degenerate_faces << std::endl;

//     // ========================================
//     // 步骤3：写入OBJ文件
//     // ========================================

//     // 写入顶点
//     std::vector<Point> ordered_points(vertex_count);
//     for (const auto& pair : point_indices) {
//         ordered_points[pair.second - 1] = pair.first;
//     }

//     for (const auto& p : ordered_points) {
//         out << "v " << p.x() << " " << p.y() << " " << p.z() << std::endl;
//     }

//     // 三角化并写入面
//     int triangle_count = 0;
//     std::map<int, int> face_vertex_distribution;

//     for (const auto& face : faces) {
//         int num_vertices = face.size();
//         face_vertex_distribution[num_vertices]++;

//         if (face.size() >= 3) {
//             // 扇形三角化从第一个顶点开始
//             for (size_t i = 1; i < face.size() - 1; ++i) {
//                 out << "f " << face[0] << " "
//                     << face[i] << " " << face[i + 1] << std::endl;
//                 triangle_count++;
//             }
//         }
//     }

//     out.close();
//     // ========================================
//     // 最终报告
//     // ========================================
//     std::cout << "\n--- Final Statistics ---" << std::endl;
//     std::cout << "Output file: " << filename << std::endl;
//     std::cout << "├─ Vertices: " << vertex_count << std::endl;
//     std::cout << "├─ Voronoi faces: " << valid_faces << std::endl;
//     std::cout << "└─ Triangles: " << triangle_count << std::endl;

//     std::cout << "\nFace complexity (vertices per face):" << std::endl;
//     for (const auto& pair : face_vertex_distribution) {
//         std::cout << "  " << pair.first << " vertices: "
//             << pair.second << " faces" << std::endl;
//     }

//     // ========================================
//     // 验证结果
//     // ========================================
//     std::cout << "\n--- Validation ---" << std::endl;

//     double boundary_percentage = 100.0 * boundary_edges_count / total_edges;

//     std::cout << "Boundary edges percentage: " << boundary_percentage << "%" << std::endl;

//     // 期望值：如果内外点数量接近，边界边应该占50-70%
//     if (boundary_percentage < 10.0) {
//         std::cout << "⚠️  WARNING: Very few boundary edges!" << std::endl;
//         std::cout << "   → Check if labels are correct" << std::endl;
//         std::cout << "   → Both labels should be present in the data" << std::endl;
//     }
//     else if (boundary_percentage > 90.0) {
//         std::cout << "⚠️  WARNING: Almost all edges are boundary!" << std::endl;
//         std::cout << "   → This is unusual, check your labels" << std::endl;
//     }
//     else if (boundary_percentage >= 40.0 && boundary_percentage <= 70.0) {
//         std::cout << "✅ Boundary percentage looks reasonable" << std::endl;
//         std::cout << "   → Labels appear to be well-distributed" << std::endl;
//     }

//     // 检查面的有效性
//     double valid_ratio = 100.0 * valid_faces / boundary_edges_count;
//     std::cout << "\nValid face ratio: " << valid_ratio << "%" << std::endl;

//     if (valid_ratio < 80.0) {
//         std::cout << "⚠️  WARNING: Many faces were filtered as degenerate" << std::endl;
//         std::cout << "   → This might indicate issues with the triangulation" << std::endl;
//     }
//     else {
//         std::cout << "✅ Most faces are valid (good!)" << std::endl;
//     }

//     std::cout << "\n=== Extraction Complete ===" << std::endl;
// }
// void Generator::generate_surface_mesh1(MeshingTree* seeds, const char* output_filename)
// {
//     std::cout << "Generating surface mesh using CGAL Voronoi1..." << std::endl;

//     size_t num_seeds = seeds->get_num_tree_points();
//     if (num_seeds == 0) {
//         std::cerr << "Error: No seeds provided." << std::endl;
//         return;
//     }

//     // 1. Build Delaunay triangulation with seed index info
//   //  Delaunay dt;
//     std::vector<Vertex_handle> vertex_handles(num_seeds);
//     std::vector<LabeledPoint> labeled_points;
//     for (size_t i = 0; i < num_seeds; i++) {
//         if (!seeds->tree_point_is_active(i)) continue;
//         double* pt = seeds->get_tree_point(i);
//       /*  Vertex_handle vh = dt.insert(Point_3(pt[0], pt[1], pt[2]));
//         vh->info() = i;
//         vertex_handles[i] = vh;*/ 
//         size_t* attrib = seeds->get_tree_point_attrib(i);

//         // 根据区域ID判断内外
//         // attrib[5] 是区域ID: 0 = 外部, 非0 = 内部
//        // if (attrib[5] == 1)
//         labeled_points.push_back(LabeledPoint(
//            pt[0],
//             pt[1],
//             pt[2],
//             attrib[5]
//         ));
//     }

//     if (labeled_points.empty()) {
//         std::cerr << "Error: No points loaded!" << std::endl;
//         return ;
//     }

//     std::cout << "Loaded " << labeled_points.size() << " points" << std::endl;

//     // Count labels
//     std::map<int, int> label_counts;
//     for (const auto& lp : labeled_points) {
//         label_counts[lp.label]++;
//     }

//     std::cout << "\nLabel distribution:" << std::endl;
//     for (const auto& pair : label_counts) {
//         std::cout << "  Label " << pair.first << ": " << pair.second << " points" << std::endl;
//     }

//     // Build Delaunay triangulation
//     std::cout << "\nBuilding Delaunay triangulation..." << std::endl;
//     Delaunay dt;
//     std::map<Vertex_handle, int> vertex_labels;

//     for (const auto& lp : labeled_points) {
//         Vertex_handle vh = dt.insert(lp.point);
//         vertex_labels[vh] = lp.label;
//     }

//     std::cout << "Number of vertices: " << dt.number_of_vertices() << std::endl;
//     std::cout << "Number of cells: " << dt.number_of_cells() << std::endl;
//     std::cout << "Number of facets: " << dt.number_of_facets() << std::endl;

//     // Export results
//     std::cout << "\n=== Exporting Surface Meshes ===" << std::endl;

//     // Method 1: Voronoi-based surface (from Voronoi faces)
//     exportBoundarySurfaceImproved(dt, vertex_labels, "voronoi_surface.obj");

//     // Method 2: Delaunay-based surface (from Delaunay facets)
//     exportDelaunayBoundarySurface(dt, vertex_labels, "delaunay_surface.obj");

//     // Export input points for reference
//     std::ofstream points_out("input_points.obj");
//     std::ofstream points_out1("output_points.obj");
//     if (points_out) {
//         points_out << "# Input points" << std::endl;
//         for (const auto& lp : labeled_points) {
//             if (lp.label != 1)  // Skip label 0 if needed
//                 points_out << "v " << lp.point.x() << " "
//                 << lp.point.y() << " " << lp.point.z() << std::endl;

//             else
//                 points_out1 << "v " << lp.point.x() << " "
//                 << lp.point.y() << " " << lp.point.z() << std::endl;
//         }
//         points_out.close();
//         points_out1.close();
//         std::cout << "\nInput points exported to input_points.obj" << std::endl;
//     }

//     std::cout << "\n=== Complete ===" << std::endl;
//     std::cout << "\nGenerated surface files:" << std::endl;
//     std::cout << "  - voronoi_surface.obj    : Voronoi-based surface mesh (triangulated)" << std::endl;
//     std::cout << "  - delaunay_surface.obj   : Delaunay-based surface mesh (faster)" << std::endl;
//     std::cout << "  - input_points.obj       : Original input points" << std::endl;
//     std::cout << "\nRecommendation:" << std::endl;
//     std::cout << "  - Use delaunay_surface.obj for faster results" << std::endl;
//     std::cout << "  - Use voronoi_surface.obj for true Voronoi surface" << std::endl;
// }