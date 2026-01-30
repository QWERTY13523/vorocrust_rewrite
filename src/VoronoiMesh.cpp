#include "Generator.h"
#include <map>
#include <fstream>
#include <iomanip>

#ifdef USE_CGAL
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

    // ==================================================================================
    // [新增部分] 输出全局维诺图 (Global Voronoi Diagram Wireframe)
    // ==================================================================================
    {
        std::cout << "  * Exporting global Voronoi diagram to global_voronoi.obj ..." << std::endl;
        std::ofstream vor_out("global_voronoi.obj");
        if (vor_out.is_open()) {
            vor_out << std::fixed << std::setprecision(16);
            vor_out << "# Global Voronoi Diagram (Finite Edges Only)\n";
            
            size_t edge_v_count = 1;
            
            // 遍历所有有限的 Delaunay 面 (Finite Facets)
            // 在 3D 中，Delaunay 面 <---> Voronoi 边
            for (auto fit = dt.finite_facets_begin(); fit != dt.finite_facets_end(); ++fit) {
                // 计算对偶对象
                CGAL::Object o = dt.dual(*fit);
                
                // 我们只输出有限的线段 (Segment)，忽略射线 (Ray)
                if (const K::Segment_3* s = CGAL::object_cast<K::Segment_3>(&o)) {
                    
                    // (可选) 过滤掉过长的线段，防止某些退化四面体导致的无穷远点干扰视图
                    if (s->squared_length() > 1e12) continue; 

                    vor_out << "v " << s->source().x() << " " << s->source().y() << " " << s->source().z() << "\n";
                    vor_out << "v " << s->target().x() << " " << s->target().y() << " " << s->target().z() << "\n";
                    vor_out << "l " << edge_v_count << " " << edge_v_count + 1 << "\n";
                    edge_v_count += 2;
                }
            }
            vor_out.close();
            std::cout << "  * Done." << std::endl;
        } else {
            std::cerr << "Error: Cannot write to global_voronoi.obj" << std::endl;
        }
    }
    // ==================================================================================

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
        // bool is_interface = (attrib1[5] != attrib2[5]) && (attrib1[5] != 0) && (attrib2[5] != 0);
        bool is_pair = (attrib1[1] == seed_idx2) || (attrib2[1] == seed_idx1);
        
        if (!is_pair) continue;
        
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
                 << CGAL::to_double(v.z()) << std::endl;
    }
    
    // Write faces (OBJ uses 1-based indexing)
    for (const auto& face : face_indices) {
        obj_file << "f " << (face[0] + 1) << " " << (face[1] + 1) << " " << (face[2] + 1) << std::endl;
    }
    
    obj_file.close();
    
    std::cout << "  * Surface mesh saved to " << output_filename << std::endl;
    std::cout << "  * Total vertices: " << unique_vertices.size() << ", triangles: " << face_indices.size() << std::endl;
}

#else
// Fallback when CGAL is not available
void Generator::generate_surface_mesh(MeshingTree* seeds, const char* output_filename)
{
    std::cerr << "Error: generate_surface_mesh requires CGAL. Please rebuild with -DUSE_CGAL=ON" << std::endl;
}
#endif