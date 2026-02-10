#include "Generator.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <fstream>
#include <iomanip>
#include <set>
#include <string>
#include <vector>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Delaunay_triangulation_cell_base_3.h> // 必须包含
#define M_PI 3.14159265358979323846
typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_with_info_3<size_t, K> Vb;
typedef CGAL::Delaunay_triangulation_cell_base_3<K> Cb; // 定义 Cell Base
typedef CGAL::Triangulation_data_structure_3<Vb, Cb> Tds; // 传入 Vb 和 Cb
typedef CGAL::Delaunay_triangulation_3<K, Tds> Delaunay;
typedef Delaunay::Point Point_3;
typedef Delaunay::Vertex_handle Vertex_handle;
typedef Delaunay::Cell_handle Cell_handle;
typedef Delaunay::Edge Edge;
typedef K::Point_3 Point;
void Generator::generate_surface_mesh(MeshingTree* seeds, const char* output_filename)
{
    std::cout << "Generating surface mesh using CGAL Voronoi..." << "\n";
    
    size_t num_seeds = seeds->get_num_tree_points();
    if (num_seeds == 0) {
        std::cerr << "Error: No seeds provided." << "\n";
        return;
    }

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
    
    std::cout << "  * Delaunay triangulation built with " << dt.number_of_vertices() << " vertices" << "\n";




    {
        std::cout << "  * Exporting inside/outside seed points to OBJ files..." << "\n";

        std::ofstream inside_seeds("inside_seeds.obj");
        std::ofstream outside_seeds("outside_seeds.obj");

        if (!inside_seeds.is_open() || !outside_seeds.is_open()) {
            std::cerr << "Error: Cannot create seed point output files" << "\n";
        }
        else {
            inside_seeds << std::fixed << std::setprecision(16);
            outside_seeds << std::fixed << std::setprecision(16);

            inside_seeds << "# Inside Seed Points\n";
            outside_seeds << "# Outside Seed Points\n";

            size_t inside_count = 0;
            size_t outside_count = 0;

            for (size_t i = 0; i < num_seeds; i++) {
                if (!seeds->tree_point_is_active(i)) continue;

                double* pt = seeds->get_tree_point(i);
                size_t* attrib = seeds->get_tree_point_attrib(i);

                // 根据区域ID判断内外
                // attrib[5] 是区域ID: 0 = 外部, 非0 = 内部
                if (attrib[5] == 1) {
                    // 外部种子点
                    outside_seeds << "v " << pt[0] << " " << pt[1] << " " << pt[2] << "\n";
                    outside_count++;
                }
                else {
                    // 内部种子点
                    inside_seeds << "v " << pt[0] << " " << pt[1] << " " << pt[2] << "\n";
                    inside_count++;
                }
            }

            inside_seeds.close();
            outside_seeds.close();

            std::cout << "  * Inside seeds: " << inside_count << " (saved to inside_seeds.obj)" << "\n";
            std::cout << "  * Outside seeds: " << outside_count << " (saved to outside_seeds.obj)" << "\n";
        }
    }

    {
    std::cout << "  * Exporting seed pair connections to seed_pairs.obj ..." << "\n";
    
    std::ofstream pairs_out("seed_pairs.obj");
    if (pairs_out.is_open()) {
        pairs_out << std::fixed << std::setprecision(16);
        pairs_out << "# Seed Pair Connections\n";
        
        // 首先输出所有种子点作为顶点
        for (size_t i = 0; i < num_seeds; i++) {
            if (!seeds->tree_point_is_active(i)) continue;
            double* pt = seeds->get_tree_point(i);
            pairs_out << "v " << pt[0] << " " << pt[1] << " " << pt[2] << "\n";
        }
        
        // 然后输出配对连接线
        size_t pair_count = 0;
        std::set<std::pair<size_t, size_t>> processed_pairs;
        
        for (size_t i = 0; i < num_seeds; i++) {
            if (!seeds->tree_point_is_active(i)) continue;
            
            size_t* attrib = seeds->get_tree_point_attrib(i);
            size_t pair_idx = attrib[1]; // attrib[1] 是配对种子索引
            
            if (pair_idx < num_seeds && seeds->tree_point_is_active(pair_idx)) {
                // 避免重复输出同一对
                auto pair_key = std::minmax(i, pair_idx);
                if (processed_pairs.find(pair_key) == processed_pairs.end()) {
                    // OBJ使用1-based索引
                    pairs_out << "l " << (i + 1) << " " << (pair_idx + 1) << "\n";
                    processed_pairs.insert(pair_key);
                    pair_count++;
                }
            }
        }
        
        pairs_out.close();
        std::cout << "  * Seed pairs: " << pair_count << " connections (saved to seed_pairs.obj)" << "\n";
    } else {
        std::cerr << "Error: Cannot write to seed_pairs.obj" << "\n";
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
        // 这里使用之前生成的配对信息，或者可以改用区域 ID (attrib[5]) 判断: 
        //  bool is_interface = (attrib1[5] != attrib2[5]);
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
    
    std::cout << "  * Found " << voronoi_facets.size() << " Voronoi facets for seed pairs" << "\n";

    // 3. Write Surface Mesh OBJ file
    std::ofstream obj_file(output_filename);
    if (!obj_file.is_open()) {
        std::cerr << "Error: Cannot open output file " << output_filename << "\n";
        return;
    }
    
    obj_file << "# Voronoi surface mesh generated by VoroCrust" << "\n";
    obj_file << "# Number of facets: " << voronoi_facets.size() << "\n";
    
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
            
            // 跳过退化三角形
            if (idx0 == idx1 || idx0 == centroid_idx || idx1 == centroid_idx) continue;

            face_indices.push_back({centroid_idx, idx0, idx1});
        }
    }
    
    // Write vertices
    for (const auto& v : unique_vertices) {
        obj_file << "v " << std::setprecision(16) 
                 << CGAL::to_double(v.x()) << " " 
                 << CGAL::to_double(v.y()) << " " 
                 << CGAL::to_double(v.z()) << "\n";
    }
    
    // Write faces (OBJ uses 1-based indexing)
    for (const auto& face : face_indices) {
        obj_file << "f " << (face[0] + 1) << " " << (face[1] + 1) << " " << (face[2] + 1) << "\n";
    }
    
    obj_file.close();
    
    std::cout << "  * Surface mesh saved to " << output_filename << "\n";
    std::cout << "  * Total vertices: " << unique_vertices.size() << ", triangles: " << face_indices.size() << "\n";

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
                std::cerr << "Cannot open file " << filename << "\n";
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
                std::cerr << "Cannot open file " << filename << "\n";
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
            double epsilon = 1e-1,
            double angle_threshold_deg = 170.0,
            bool debug = false
        ) {
            std::cout << "  * Starting optimized deduplication export to " << filename << "..." << "\n";
            std::ofstream out(filename);
            if (!out) {
                std::cerr << "Cannot open file " << filename << "\n";
                return;
            }
            out << "# Voronoi facets (epsilon=" << epsilon
                << ", angle_threshold=" << angle_threshold_deg << "°)\n";

            std::vector<Point_3> vertices;
            std::vector<std::vector<int>> faces;
            
            // --- 优化开始：使用 map 加速查找 ---
            // 使用 tuple<long, long, long> 作为 key，通过 epsilon 缩放将坐标离散化
            using VertexKey = std::tuple<long long, long long, long long>;
            std::map<VertexKey, int> vertex_map; 
            double scale_factor = 1.0 / (epsilon > 0 ? epsilon : 1e-6);

            int total_input_facets = 0;
            int total_output_facets = 0;

            // 优化后的查找函数：O(log N)
            auto find_or_add_vertex = [&](const Point_3& p) -> int {
                double x = CGAL::to_double(p.x());
                double y = CGAL::to_double(p.y());
                double z = CGAL::to_double(p.z());
                
                // 坐标离散化
                long long ix = std::llround(x * scale_factor);
                long long iy = std::llround(y * scale_factor);
                long long iz = std::llround(z * scale_factor);
                
                VertexKey key = std::make_tuple(ix, iy, iz);
                
                auto it = vertex_map.find(key);
                if (it != vertex_map.end()) {
                    return it->second;
                }
                
                vertices.push_back(p);
                int idx = static_cast<int>(vertices.size());
                vertex_map[key] = idx;
                return idx;
            };
            // --- 优化结束 ---

            auto compute_angle = [&](const Point_3& p1, const Point_3& p2, const Point_3& p3) -> double {
                typedef CGAL::Vector_3<K> Vector_3;
                Vector_3 v1 = p1 - p2;
                Vector_3 v2 = p3 - p2;

                double len1_sq = CGAL::to_double(v1.squared_length());
                double len2_sq = CGAL::to_double(v2.squared_length());
                if (len1_sq < 1e-20 || len2_sq < 1e-20) return 0.0;

                double dot = CGAL::to_double(v1 * v2);
                double cos_angle = dot / std::sqrt(len1_sq * len2_sq);
                cos_angle = std::max(-1.0, std::min(1.0, cos_angle));

                return std::acos(cos_angle) * 180.0 / M_PI;
            };

            auto extract_corner_vertices = [&](const std::vector<Point_3>& polygon) -> std::vector<Point_3> {
                if (polygon.size() <= 3) return polygon;
                std::vector<Point_3> corners;
                std::vector<bool> is_corner(polygon.size(), false);
                int n = static_cast<int>(polygon.size());

                for (int i = 0; i < n; ++i) {
                    const Point_3& prev = polygon[(i - 1 + n) % n];
                    const Point_3& curr = polygon[i];
                    const Point_3& next = polygon[(i + 1) % n];
                    if (compute_angle(prev, curr, next) < angle_threshold_deg) is_corner[i] = true;
                }

                for (int i = 0; i < n; ++i) {
                    if (is_corner[i]) corners.push_back(polygon[i]);
                }
                return corners.size() < 3 ? polygon : corners;
            };

            // 处理所有面片
            for (const auto& facet : voronoi_facets) {
                total_input_facets++;
                if (facet.size() < 3) continue;

                std::vector<Point_3> corner_points = extract_corner_vertices(facet);
                if (corner_points.size() < 3) continue;

                std::vector<int> face;
                for (const auto& p : corner_points) {
                    face.push_back(find_or_add_vertex(p));
                }

                // 清理连续重复索引
                std::vector<int> cleaned_face;
                for (size_t i = 0; i < face.size(); ++i) {
                    if (i == 0 || face[i] != face[i - 1]) cleaned_face.push_back(face[i]);
                }
                if (cleaned_face.size() > 1 && cleaned_face.front() == cleaned_face.back()) {
                    cleaned_face.pop_back();
                }

                if (cleaned_face.size() >= 3) {
                    faces.push_back(cleaned_face);
                    total_output_facets++;
                }
            }

            // 写入文件
            out << std::fixed << std::setprecision(16);
            for (const auto& p : vertices) {
                out << "v " << CGAL::to_double(p.x()) << " " << CGAL::to_double(p.y()) << " " << CGAL::to_double(p.z()) << "\n";
            }
            for (const auto& f : faces) {
                out << "f";
                for (int idx : f) out << " " << idx;
                out << "\n";
            }

            std::cout << "  * Exported " << total_output_facets << " facets (" 
                      << vertices.size() << " vertices) to " << filename << "\n";
        };

        auto export_single_voronoi_polygon = [&](Vertex_handle v1, Vertex_handle v2, const std::string& filename) {
            std::cout << "=== Exporting single Voronoi polygon ===" << "\n";
            std::cout << "  * Edge: seed " << v1->info() << " <-> seed " << v2->info() << "\n";

            std::ofstream obj_file(filename);
            if (!obj_file.is_open()) {
                std::cerr << "Error: Cannot open " << filename << "\n";
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

            std::cout << "  * Polygon has " << polygon_vertices.size() << " vertices" << "\n";

            if (polygon_vertices.size() < 3) {
                std::cout << "  * Warning: Not enough vertices to form a polygon" << "\n";
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
            std::cout << "  * Saved to " << filename << "\n";
        };

        auto exportBoundarySurfaceImproved = [&](const std::string& filename) {
            std::ofstream out(filename);
            if (!out) {
                std::cerr << "Error: Cannot create file " << filename << "\n";
                return;
            }

            out << "# Voronoi boundary surface (improved with validation)" << "\n";
            std::cout << "\n=== Extracting Boundary Surface (Improved) ===" << "\n";

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
                        if (!did_export_single) {
                            export_single_voronoi_polygon(v1, v2, "single_voronoi_polygon.obj");
                            did_export_single = true;
                        }
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
            std::cout << "  * Found " << facets.size() << " Voronoi facets for seed pairs" << "\n";

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
                        << CGAL::to_double(v.z()) << "\n";
                }
                for (const auto& face : face_indices) {
                    out << "f " << (face[0] + 1) << " " << (face[1] + 1) << " " << (face[2] + 1) << "\n";
                }

                std::cout << "  * Total vertices: " << unique_vertices.size() << ", triangles: " << face_indices.size() << "\n";
            }

            std::cout << "\n--- Edge Classification ---" << "\n";
            std::cout << "Total Delaunay edges: " << total_edges << "\n";
            std::cout << "├─ Boundary edges (1-2): " << boundary_edges_count
                << " (" << (100.0 * boundary_edges_count / std::max(1, total_edges)) << "%)" << "\n";
            std::cout << "├─ Inner edges (1-1): " << inner_edges
                << " (" << (100.0 * inner_edges / std::max(1, total_edges)) << "%)" << "\n";
            std::cout << "├─ Outer edges (2-2): " << outer_edges
                << " (" << (100.0 * outer_edges / std::max(1, total_edges)) << "%)" << "\n";
            std::cout << "└─ No label edges: " << no_label_edges << "\n";
        };

        auto exportDelaunayBoundarySurface = [&](const std::string& filename) {
            std::ofstream out(filename);
            if (!out) {
                std::cerr << "Error: Cannot create file " << filename << "\n";
                return;
            }

            out << "# Delaunay-based boundary surface" << "\n";
            std::cout << "\nExtracting Delaunay boundary surface..." << "\n";

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
                out << "v " << p.x() << " " << p.y() << " " << p.z() << "\n";
            }
            for (const auto& tri : triangles) {
                out << "f " << point_to_index[tri[0]]
                    << " " << point_to_index[tri[1]]
                    << " " << point_to_index[tri[2]] << "\n";
            }
            std::cout << "Delaunay surface saved to: " << filename << "\n";
        };

        exportBoundarySurfaceImproved("voronoi_surface.obj");
        exportDelaunayBoundarySurface("delaunay_surface.obj");
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
        std::cerr << "Error: Cannot create file " << filename << "\n";
        return;
    }

    out << "# Voronoi boundary surface mesh" << "\n";
    out << "# Triangulated surface between different labels" << "\n";

    std::cout << "\nExtracting boundary surface..." << "\n";

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

    std::cout << "Found " << boundary_edges.size() << " boundary edges" << "\n";

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
        out << "v " << p.x() << " " << p.y() << " " << p.z() << "\n";
    }

    // Triangulate and write faces
    int triangle_count = 0;
    for (const auto& face : faces) {
        if (face.size() >= 3) {
            // Simple fan triangulation from first vertex
            for (size_t i = 1; i < face.size() - 1; ++i) {
                out << "f " << face[0] << " " << face[i] << " " << face[i + 1] << "\n";
                triangle_count++;
            }
        }
    }

    out.close();

    std::cout << "Surface mesh exported:" << "\n";
    std::cout << "  Vertices: " << vertex_count << "\n";
    std::cout << "  Voronoi faces: " << faces.size() << "\n";
    std::cout << "  Triangles: " << triangle_count << "\n";
    std::cout << "  Saved to: " << filename << "\n";
}

// Export Delaunay facets as triangulated surface
void exportDelaunayBoundarySurface(const Delaunay& dt,
    const std::map<Vertex_handle, int>& vertex_labels,
    const std::string& filename) {
    std::ofstream out(filename);
    if (!out) {
        std::cerr << "Error: Cannot create file " << filename << "\n";
        return;
    }

    out << "# Delaunay-based boundary surface" << "\n";

    std::cout << "\nExtracting Delaunay boundary surface..." << "\n";

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

    std::cout << "Found " << triangles.size() << " boundary triangles" << "\n";
    std::cout << "With " << unique_points.size() << " unique vertices" << "\n";

    // Map points to indices
    std::map<Point, int> point_to_index;
    int idx = 1;
    for (const auto& p : unique_points) {
        point_to_index[p] = idx++;
    }

    // Write vertices
    for (const auto& p : unique_points) {
        out << "v " << p.x() << " " << p.y() << " " << p.z() << "\n";
    }

    // Write faces
    for (const auto& tri : triangles) {
        out << "f " << point_to_index[tri[0]]
            << " " << point_to_index[tri[1]]
            << " " << point_to_index[tri[2]] << "\n";
    }

    out.close();
    std::cout << "Delaunay surface saved to: " << filename << "\n";
}
// 改进版：带详细统计和验证的边界面提取
void export_single_voronoi_polygon(
    const Delaunay& dt,
    Vertex_handle v1,
    Vertex_handle v2,
    const std::string& filename)
{
    std::cout << "=== Exporting single Voronoi polygon ===" << "\n";
    std::cout << "  * Edge: seed " << v1->info() << " <-> seed " << v2->info() << "\n";

    std::ofstream obj_file(filename);
    if (!obj_file.is_open()) {
        std::cerr << "Error: Cannot open " << filename << "\n";
        return;
    }

    obj_file << std::fixed << std::setprecision(16);
    obj_file << "# Voronoi polygon for Delaunay edge\n";
    obj_file << "# Seed 1 (index " << v1->info() << ")\n";
    obj_file << "# Seed 2 (index " << v2->info() << ")\n\n";

    // 收集围绕该边的所有四面体的外接球心
    std::vector<Point_3> polygon_vertices;

    // 需要从边构造edge
    // 在CGAL中，edge是通过(cell, i, j)三元组表示的
    // 这里我们用另一种方法：找到包含v1和v2的所有cell

    std::vector<Cell_handle> cells_v1;
    dt.incident_cells(v1, std::back_inserter(cells_v1));

    for (auto cell : cells_v1) {
        if (dt.is_infinite(cell)) continue;

        // 检查这个cell是否也包含v2
        bool contains_v2 = false;
        for (int i = 0; i < 4; i++) {
            if (cell->vertex(i) == v2) {
                contains_v2 = true;
                break;
            }
        }

        if (contains_v2) {
            // 这个cell同时包含v1和v2，其外接球心是多边形的一个顶点
            Point_3 center = dt.dual(cell);
            polygon_vertices.push_back(center);
        }
    }

    std::cout << "  * Polygon has " << polygon_vertices.size() << " vertices" << "\n";

    if (polygon_vertices.size() < 3) {
        std::cout << "  * Warning: Not enough vertices to form a polygon" << "\n";
        obj_file.close();
        return;
    }

    // 输出多边形顶点
    for (size_t i = 0; i < polygon_vertices.size(); i++) {
        const auto& v = polygon_vertices[i];
        obj_file << "v " << CGAL::to_double(v.x()) << " "
            << CGAL::to_double(v.y()) << " "
            << CGAL::to_double(v.z()) << "\n";
    }

    // 输出多边形的边（线段）
    obj_file << "\n# Polygon edges\n";
    for (size_t i = 0; i < polygon_vertices.size(); i++) {
        size_t next_i = (i + 1) % polygon_vertices.size();
        obj_file << "l " << (i + 1) << " " << (next_i + 1) << "\n";
    }

    // 可选：也输出多边形的面（如果顶点共面）
    obj_file << "\n# Polygon face\n";
    obj_file << "f";
    for (size_t i = 0; i < polygon_vertices.size(); i++) {
        obj_file << " " << (i + 1);
    }
    obj_file << "\n";

    // 输出对应的种子点位置（用于参考）
    obj_file << "\n# Corresponding seed points\n";
    obj_file << "# v " << CGAL::to_double(v1->point().x()) << " "
        << CGAL::to_double(v1->point().y()) << " "
        << CGAL::to_double(v1->point().z()) << " # Seed 1\n";
    obj_file << "# v " << CGAL::to_double(v2->point().x()) << " "
        << CGAL::to_double(v2->point().y()) << " "
        << CGAL::to_double(v2->point().z()) << " # Seed 2\n";

    obj_file.close();
    std::cout << "  * Saved to " << filename << "\n";
}

void write_voronoi_facets_to_obj_dedup(
    const std::string& filename,
    const std::vector<std::vector<Point_3>>& voronoi_facets,
    double epsilon = 1e-1,
    double angle_threshold_deg = 170.0,  // 角度阈值（度），大于此值视为共线
    bool debug = false  // 是否输出调试信息
) {
    std::ofstream out(filename);
    if (!out) {
        std::cerr << "Cannot open file " << filename << "\n";
        return;
    }
    out << "# Voronoi facets (epsilon=" << epsilon
        << ", angle_threshold=" << angle_threshold_deg << "°)\n";

    std::vector<Point_3> vertices;
    std::vector<std::vector<int>> faces;
    double cos_threshold = std::cos(angle_threshold_deg * M_PI / 180.0);

    int total_input_facets = 0;
    int total_output_facets = 0;

    // ---------- 辅助函数：查找或添加顶点（网格化加速 O(log N)）----------
    using VertexKey_dedup = std::tuple<long long, long long, long long>;
    std::map<VertexKey_dedup, int> vertex_grid;
    double scale_factor = 1.0 / (epsilon > 0 ? epsilon : 1e-6);

    auto find_or_add_vertex = [&](const Point_3& p) -> int {
        double x = CGAL::to_double(p.x());
        double y = CGAL::to_double(p.y());
        double z = CGAL::to_double(p.z());
        long long ix = std::llround(x * scale_factor);
        long long iy = std::llround(y * scale_factor);
        long long iz = std::llround(z * scale_factor);
        VertexKey_dedup key = std::make_tuple(ix, iy, iz);
        auto it = vertex_grid.find(key);
        if (it != vertex_grid.end()) return it->second;
        vertices.push_back(p);
        int idx = static_cast<int>(vertices.size());
        vertex_grid[key] = idx;
        return idx;
        };

    // ---------- 辅助函数：计算三点形成的角度（度数）----------
    auto compute_angle = [&](const Point_3& p1, const Point_3& p2, const Point_3& p3) -> double {
        typedef CGAL::Vector_3<K> Vector_3;
        Vector_3 v1 = p1 - p2;  // 从p2指向p1
        Vector_3 v2 = p3 - p2;  // 从p2指向p3

        double len1_sq = CGAL::to_double(v1.squared_length());
        double len2_sq = CGAL::to_double(v2.squared_length());

        if (len1_sq < 1e-20 || len2_sq < 1e-20) return 0.0;

        double dot = CGAL::to_double(v1 * v2);
        double cos_angle = dot / std::sqrt(len1_sq * len2_sq);

        // 限制cos值在[-1, 1]范围内，避免浮点误差
        cos_angle = std::max(-1.0, std::min(1.0, cos_angle));

        double angle_rad = std::acos(cos_angle);
        return angle_rad * 180.0 / M_PI;
        };

    // ---------- 改进的角点提取函数 ----------
    auto extract_corner_vertices = [&](const std::vector<Point_3>& polygon) -> std::vector<Point_3> {
        if (polygon.size() <= 3) return polygon;

        std::vector<Point_3> corners;
        std::vector<bool> is_corner(polygon.size(), false);
        int n = polygon.size();

        // 第一遍：标记所有角点
        for (int i = 0; i < n; ++i) {
            const Point_3& prev = polygon[(i - 1 + n) % n];
            const Point_3& curr = polygon[i];
            const Point_3& next = polygon[(i + 1) % n];

            double angle = compute_angle(prev, curr, next);

            // 如果角度小于阈值（不接近180度），则是角点
            if (angle < angle_threshold_deg) {
                is_corner[i] = true;
            }
        }

        // 收集角点
        for (int i = 0; i < n; ++i) {
            if (is_corner[i]) {
                corners.push_back(polygon[i]);
            }
        }

        // 如果检测到的角点太少，直接返回原始多边形
        if (corners.size() < 3) {
            if (debug) {
                std::cout << "Warning: Only found " << corners.size()
                    << " corners, using original polygon with "
                    << polygon.size() << " vertices" << "\n";
            }
            return polygon;
        }

        return corners;
        };

    // ---------- 构建顶点列表和面索引 ----------
    for (const auto& facet : voronoi_facets) {
        total_input_facets++;

        if (facet.size() < 3) {
            if (debug) {
                std::cout << "Skipping facet with < 3 vertices" << "\n";
            }
            continue;
        }

        // 提取角点
        std::vector<Point_3> corner_points = extract_corner_vertices(facet);

        if (debug && corner_points.size() != facet.size()) {
            std::cout << "Facet " << total_input_facets << ": "
                << facet.size() << " vertices -> "
                << corner_points.size() << " corners" << "\n";
        }

        if (corner_points.size() < 3) continue;

        std::vector<int> face;
        for (const auto& p : corner_points) {
            int idx = find_or_add_vertex(p);
            face.push_back(idx);
        }

        // 移除重复的连续索引
        std::vector<int> cleaned_face;
        for (size_t i = 0; i < face.size(); ++i) {
            if (i == 0 || face[i] != face[i - 1]) {
                cleaned_face.push_back(face[i]);
            }
        }

        // 检查首尾是否相同
        if (cleaned_face.size() > 1 && cleaned_face.front() == cleaned_face.back()) {
            cleaned_face.pop_back();
        }

        if (cleaned_face.size() >= 3) {
            faces.push_back(cleaned_face);
            total_output_facets++;
        }
    }

    // ---------- 输出顶点 ----------
    for (const auto& p : vertices) {
        out << "v "
            << p.x() << " "
            << p.y() << " "
            << p.z() << "\n";
    }

    // ---------- 输出面 ----------
    for (const auto& f : faces) {
        out << "f";
        for (int idx : f) {
            out << " " << idx;
        }
        out << "\n";
    }

    out.close();

    std::cout << "=== Summary ===" << "\n";
    std::cout << "Input facets: " << total_input_facets << "\n";
    std::cout << "Output facets: " << total_output_facets << "\n";
    std::cout << "Unique vertices: " << vertices.size() << "\n";
    std::cout << "Wrote to: " << filename << "\n";
}

void write_voronoi_facets_triangulated_obj(
    const std::string& filename,
    const std::vector<std::vector<Point_3>>& voronoi_facets
) {
    std::ofstream out(filename);
    if (!out) {
        std::cerr << "Cannot open file " << filename << "\n";
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
        std::cerr << "Cannot open file " << filename << "\n";
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
        std::cerr << "Error: Cannot create file " << filename << "\n";
        return;
    }

    out << "# Voronoi boundary surface (improved with validation)" << "\n";
    std::cout << "\n=== Extracting Boundary Surface (Improved) ===" << "\n";

    // ========================================
    // 统计：帮助理解哪些面被保留/丢弃
    // ========================================
    int total_edges = 0;
    int boundary_edges_count = 0;  // label 1-2
    int inner_edges = 0;            // label 1-1
    int outer_edges = 0;            // label 2-2
    int no_label_edges = 0;         // 没有标签的边

    std::set<std::pair<Vertex_handle, Vertex_handle>> boundary_edges;

    // ========================================
    // 步骤1：遍历所有Delaunay边
    // 每条边对应一个Voronoi面
    // ========================================
    std::vector<std::vector<Point_3>> voronoi_facets;
    int fla = 0;
    for (auto eit = dt.finite_edges_begin(); eit != dt.finite_edges_end(); ++eit) {
        total_edges++;

        Cell_handle c = eit->first;
        int i = eit->second;
        int j = eit->third;

        // 获取边的两个端点（站点）
        Vertex_handle v1 = c->vertex(i);
        Vertex_handle v2 = c->vertex(j);

        // 检查是否都有标签
        if (vertex_labels.find(v1) == vertex_labels.end() ||
            vertex_labels.find(v2) == vertex_labels.end()) {
            no_label_edges++;
            continue;
        }

        int label1 = vertex_labels.at(v1);
        int label2 = vertex_labels.at(v2);

        // ========================================
        // 核心判断：保留条件
        // 
        // Voronoi面由2个站点生成
        // 保留条件：这2个站点的label不同
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
       
            // ✅ 边界面：一个内点 + 一个外点
            boundary_edges_count++;

            // 规范化边（保证唯一性）
            if (v1 < v2) {
                boundary_edges.insert({ v1, v2 });
            }
            else {
                boundary_edges.insert({ v2, v1 });
            }
        }
        else {
            // ❌ 非边界面
            if (label1 == 1) {
                inner_edges++;  // 内-内
            }
            else if (label1 == 2) {
                outer_edges++;  // 外-外
            }
        }
    }

    write_voronoi_facets_to_obj("voronoi.obj", voronoi_facets);
    write_voronoi_facets_to_obj_dedup("voronoi_dedup.obj", voronoi_facets);
    write_voronoi_facets_triangulated_obj("voronoi_dedup_triangulated.obj", voronoi_facets);
std::cout << "  * Found " << voronoi_facets.size() << " Voronoi facets for seed pairs" << "\n";

// 3. Write Surface Mesh OBJ file
std::ofstream obj_file("output_filename.obj");
if (!obj_file.is_open()) {
    //  std::cerr << "Error: Cannot open output file " << output_filename << "\n";
    return;
}

obj_file << "# Voronoi surface mesh generated by VoroCrust" << "\n";
obj_file << "# Number of facets: " << voronoi_facets.size() << "\n";

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

        // 跳过退化三角形
        if (idx0 == idx1 || idx0 == centroid_idx || idx1 == centroid_idx) continue;

        face_indices.push_back({ centroid_idx, idx0, idx1 });
    }
}

// Write vertices
for (const auto& v : unique_vertices) {
    obj_file << "v " << std::setprecision(16)
        << CGAL::to_double(v.x()) << " "
        << CGAL::to_double(v.y()) << " "
        << CGAL::to_double(v.z()) << "\n";
}

// Write faces (OBJ uses 1-based indexing)
for (const auto& face : face_indices) {
    obj_file << "f " << (face[0] + 1) << " " << (face[1] + 1) << " " << (face[2] + 1) << "\n";
}

obj_file.close();

//  std::cout << "  * Surface mesh saved to " << output_filename << "\n";
std::cout << "  * Total vertices: " << unique_vertices.size() << ", triangles: " << face_indices.size() << "\n";
    // ========================================
    // 打印统计：验证结果是否合理
    // ========================================
    std::cout << "\n--- Edge Classification ---" << "\n";
    std::cout << "Total Delaunay edges: " << total_edges << "\n";
    std::cout << "├─ Boundary edges (1-2): " << boundary_edges_count
        << " (" << (100.0 * boundary_edges_count / total_edges) << "%)" << "\n";
    std::cout << "├─ Inner edges (1-1): " << inner_edges
        << " (" << (100.0 * inner_edges / total_edges) << "%)" << "\n";
    std::cout << "├─ Outer edges (2-2): " << outer_edges
        << " (" << (100.0 * outer_edges / total_edges) << "%)" << "\n";
    std::cout << "└─ No label edges: " << no_label_edges << "\n";





    // ========================================
    // 步骤2：复用已收集的Voronoi面（避免重复调用getVoronoiFace）
    // ========================================
    std::map<Point, int> point_indices;
    int vertex_count = 0;
    std::vector<std::vector<int>> faces;

    int valid_faces = 0;
    int degenerate_faces = 0;

    std::cout << "\n--- Extracting Voronoi Faces ---" << "\n";

    for (const auto& facet : voronoi_facets) {
        if (facet.size() < 3) {
            degenerate_faces++;
            continue;
        }

        valid_faces++;

        std::vector<int> face_indices;
        for (const auto& p : facet) {
            if (point_indices.find(p) == point_indices.end()) {
                point_indices[p] = ++vertex_count;
            }
            face_indices.push_back(point_indices[p]);
        }
        faces.push_back(face_indices);
    }

    std::cout << "Valid faces: " << valid_faces << "\n";
    std::cout << "Degenerate faces (filtered out): " << degenerate_faces << "\n";

    // ========================================
    // 步骤3：写入OBJ文件
    // ========================================

    // 写入顶点
    std::vector<Point> ordered_points(vertex_count);
    for (const auto& pair : point_indices) {
        ordered_points[pair.second - 1] = pair.first;
    }

    for (const auto& p : ordered_points) {
        out << "v " << p.x() << " " << p.y() << " " << p.z() << "\n";
    }

    // 三角化并写入面
    int triangle_count = 0;
    std::map<int, int> face_vertex_distribution;

    for (const auto& face : faces) {
        int num_vertices = face.size();
        face_vertex_distribution[num_vertices]++;

        if (face.size() >= 3) {
            // 扇形三角化从第一个顶点开始
            for (size_t i = 1; i < face.size() - 1; ++i) {
                out << "f " << face[0] << " "
                    << face[i] << " " << face[i + 1] << "\n";
                triangle_count++;
            }
        }
    }

    out.close();

    // ========================================
    // 最终报告
    // ========================================
    std::cout << "\n--- Final Statistics ---" << "\n";
    std::cout << "Output file: " << filename << "\n";
    std::cout << "├─ Vertices: " << vertex_count << "\n";
    std::cout << "├─ Voronoi faces: " << valid_faces << "\n";
    std::cout << "└─ Triangles: " << triangle_count << "\n";

    std::cout << "\nFace complexity (vertices per face):" << "\n";
    for (const auto& pair : face_vertex_distribution) {
        std::cout << "  " << pair.first << " vertices: "
            << pair.second << " faces" << "\n";
    }

    // ========================================
    // 验证结果
    // ========================================
    std::cout << "\n--- Validation ---" << "\n";

    double boundary_percentage = 100.0 * boundary_edges_count / total_edges;

    std::cout << "Boundary edges percentage: " << boundary_percentage << "%" << "\n";

    // 期望值：如果内外点数量接近，边界边应该占50-70%
    if (boundary_percentage < 10.0) {
        std::cout << "⚠️  WARNING: Very few boundary edges!" << "\n";
        std::cout << "   → Check if labels are correct" << "\n";
        std::cout << "   → Both labels should be present in the data" << "\n";
    }
    else if (boundary_percentage > 90.0) {
        std::cout << "⚠️  WARNING: Almost all edges are boundary!" << "\n";
        std::cout << "   → This is unusual, check your labels" << "\n";
    }
    else if (boundary_percentage >= 40.0 && boundary_percentage <= 70.0) {
        std::cout << "✅ Boundary percentage looks reasonable" << "\n";
        std::cout << "   → Labels appear to be well-distributed" << "\n";
    }

    // 检查面的有效性
    double valid_ratio = 100.0 * valid_faces / boundary_edges_count;
    std::cout << "\nValid face ratio: " << valid_ratio << "%" << "\n";

    if (valid_ratio < 80.0) {
        std::cout << "⚠️  WARNING: Many faces were filtered as degenerate" << "\n";
        std::cout << "   → This might indicate issues with the triangulation" << "\n";
    }
    else {
        std::cout << "✅ Most faces are valid (good!)" << "\n";
    }

    std::cout << "\n=== Extraction Complete ===" << "\n";
}
void Generator::generate_surface_mesh1(MeshingTree* seeds, const char* output_filename)
{
    std::cout << "Generating surface mesh using CGAL Voronoi1..." << "\n";

    size_t num_seeds = seeds->get_num_tree_points();
    if (num_seeds == 0) {
        std::cerr << "Error: No seeds provided." << "\n";
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

        // 根据区域ID判断内外
        // attrib[5] 是区域ID: 0 = 外部, 非0 = 内部
       // if (attrib[5] == 1)
        labeled_points.push_back(LabeledPoint(
           pt[0],
            pt[1],
            pt[2],
            attrib[5]
        ));
    }

    if (labeled_points.empty()) {
        std::cerr << "Error: No points loaded!" << "\n";
        return ;
    }

    std::cout << "Loaded " << labeled_points.size() << " points" << "\n";

    // Count labels
    std::map<int, int> label_counts;
    for (const auto& lp : labeled_points) {
        label_counts[lp.label]++;
    }

    std::cout << "\nLabel distribution:" << "\n";
    for (const auto& pair : label_counts) {
        std::cout << "  Label " << pair.first << ": " << pair.second << " points" << "\n";
    }

    // Build Delaunay triangulation
    std::cout << "\nBuilding Delaunay triangulation..." << "\n";
    Delaunay dt;
    std::map<Vertex_handle, int> vertex_labels;

    for (const auto& lp : labeled_points) {
        Vertex_handle vh = dt.insert(lp.point);
        vertex_labels[vh] = lp.label;
    }

    std::cout << "Number of vertices: " << dt.number_of_vertices() << "\n";
    std::cout << "Number of cells: " << dt.number_of_cells() << "\n";
    std::cout << "Number of facets: " << dt.number_of_facets() << "\n";

    // Export results
    std::cout << "\n=== Exporting Surface Meshes ===" << "\n";

    // Method 1: Voronoi-based surface (from Voronoi faces)
    exportBoundarySurfaceImproved(dt, vertex_labels, "voronoi_surface.obj");

    // Method 2: Delaunay-based surface (from Delaunay facets)
    exportDelaunayBoundarySurface(dt, vertex_labels, "delaunay_surface.obj");

    // Export input points for reference
    std::ofstream points_out("input_points.obj");
    std::ofstream points_out1("output_points.obj");
    if (points_out) {
        points_out << "# Input points" << "\n";
        for (const auto& lp : labeled_points) {
            if (lp.label != 1)  // Skip label 0 if needed
                points_out << "v " << lp.point.x() << " "
                << lp.point.y() << " " << lp.point.z() << "\n";

            else
                points_out1 << "v " << lp.point.x() << " "
                << lp.point.y() << " " << lp.point.z() << "\n";
        }
        points_out.close();
        points_out1.close();
        std::cout << "\nInput points exported to input_points.obj" << "\n";
    }

    std::cout << "\n=== Complete ===" << "\n";
    std::cout << "\nGenerated surface files:" << "\n";
    std::cout << "  - voronoi_surface.obj    : Voronoi-based surface mesh (triangulated)" << "\n";
    std::cout << "  - delaunay_surface.obj   : Delaunay-based surface mesh (faster)" << "\n";
    std::cout << "  - input_points.obj       : Original input points" << "\n";
    std::cout << "\nRecommendation:" << "\n";
    std::cout << "  - Use delaunay_surface.obj for faster results" << "\n";
    std::cout << "  - Use voronoi_surface.obj for true Voronoi surface" << "\n";
}