#include "Generator.h"
#include <map>
#include <fstream>
#include <iomanip>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Delaunay_triangulation_cell_base_3.h> // 必须包含

typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_with_info_3<size_t, K> Vb;
typedef CGAL::Delaunay_triangulation_cell_base_3<K> Cb; // 定义 Cell Base
typedef CGAL::Triangulation_data_structure_3<Vb, Cb> Tds; // 传入 Vb 和 Cb
typedef CGAL::Delaunay_triangulation_3<K, Tds> Delaunay;
typedef Delaunay::Point Point_3;
typedef Delaunay::Vertex_handle Vertex_handle;
typedef Delaunay::Cell_handle Cell_handle;
typedef Delaunay::Edge Edge;

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
        
        bool is_pair = (attrib1[1] == seed_idx2) || (attrib2[1] == seed_idx1);
        
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
    
    std::map<std::tuple<double, double, double>, size_t> vertex_map;
    std::vector<Point_3> unique_vertices;
    std::vector<std::vector<size_t>> face_indices;
    
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
        if (facet.size() < 3) continue;
        
        std::vector<size_t> indices;
        for (const auto& pt : facet) {
            indices.push_back(get_vertex_index(pt));
        }
        
        Point_3 centroid = CGAL::ORIGIN;
        for(const auto& pt : facet) centroid = centroid + (pt - CGAL::ORIGIN);
        double s = static_cast<double>(facet.size());
        centroid = Point_3(centroid.x() / s, centroid.y() / s, centroid.z() / s);
        size_t centroid_idx = get_vertex_index(centroid);

        for (size_t i = 0; i < indices.size(); i++) {
            size_t idx0 = indices[i];
            size_t idx1 = indices[(i + 1) % indices.size()];
            if (idx0 == idx1 || idx0 == centroid_idx || idx1 == centroid_idx) continue;
            face_indices.push_back({centroid_idx, idx0, idx1});
        }
    }
    
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