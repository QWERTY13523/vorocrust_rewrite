//构建球树：长度，[x,y,z,r]，添加顶点
//查询：最近邻查询（球，半径），顶点查询，顶点数，
//树里构建顶点之间连接关系

#pragma once

#include<iostream>
#include<cstring>
#include<cmath>
#include<vector>
#include<stack>
#include<queue>
#include<cstdlib>
#include<climits>
#include "Geometry.h"

class MeshingTree
{

public:

	MeshingTree();

	~MeshingTree();

    int disable_auto_balance() { _auto_balance = false; return 0;}

    int reset_graph() { _graph.clear(); return 0;}

    int graph_get_neighbors(size_t ipoint, std::vector<size_t>& neighbors) { 
        if (ipoint < _graph.size()) neighbors = _graph[ipoint]; 
        else neighbors.clear();
        return 0;
    };

    int graph_connect_nodes(size_t ipoint, size_t jpoint);
    
    int graph_connect(size_t ipoint, size_t jpoint);

    bool graph_connected(size_t ipoint, size_t jpoint);

    int get_num_tree_points() { return _points.size(); };

    int add_tree_point(size_t num_dim, double* x, double* normal, size_t* attrib);

    int get_tree_point(size_t point_index, double* x);

	int get_tree_point(size_t point_index, size_t num_dim, double* x);

	bool get_tree_point_attrib(size_t point_index, size_t attrib_index, size_t &point_attrib);

	double get_tree_point_attrib(size_t point_index, size_t attrib_index);

    void clear_memory();

    double* get_tree_point(size_t point_index){
        return _points[point_index].data();
    }

    double* get_tree_point_normal(size_t point_index){
        return _points_normal[point_index].data();
    }

    size_t* get_tree_point_attrib(size_t point_index){
        return _points_attrib[point_index].data();
    }

    int set_tree_point_attrib(size_t point_index, size_t attrib_index, size_t attrib);

	int set_tree_point_attrib(size_t point_index, size_t attrib_index, double attrib);

    int get_closest_tree_point(size_t tree_point_index, size_t &closest_tree_point, double &closest_distance);
    int get_closest_tree_point(double* x, size_t &closest_tree_point, double &closest_distance);
    int get_closest_tree_point(double* x, size_t num_exculded_tree_points, size_t* exculded_tree_points, size_t &closest_tree_point, double &closest_distance);

	int get_tree_points_in_sphere(double* x, double r, size_t &num_points_in_sphere, size_t* &points_in_sphere, size_t &capacity);

	int lazy_delete_tree_point(const std::vector<double>& tree_point);
	int lazy_delete_tree_point(size_t point_index);
	bool tree_point_is_active(size_t point_index) const;

private:

	int lazy_delete_tree_point_internal(size_t point_index);
	int lazy_delete_tree_rebuild_compact();
	int lazy_delete_tree_rebuild_if_needed();

    int init_global_variables();

    int get_normal_component(size_t num_dim, size_t num_basis, double** basis, double* vect, double &norm);

	int kd_tree_balance_quicksort(size_t target_pos, size_t left, size_t right, size_t active_dim, size_t* tree_nodes_sorted);

	int kd_tree_quicksort_adjust_target_position(size_t target_pos, size_t left, size_t right, size_t active_dim, size_t* tree_nodes_sorted);

	int kd_tree_add_point(size_t point_index);
	int kd_tree_build_balanced();

	int kd_tree_get_nodes_order(size_t d_index, size_t node_index,                                // indices to traverse the kd-tree
		                        size_t &num_traversed, size_t* ordered_indices
	                            );


	int kd_tree_get_seeds_in_sphere(double* x,                                                    // Sphere center
		                            double r,                                                     // Sphere radius
		                            size_t d_index, size_t node_index,                            // indices to traverse the kd-tree
		                            size_t &num_points_in_sphere, size_t* &points_in_sphere,      // number of points in sphere and their indices
		                            size_t &capacity                                              // Size of points in sphere array
	                               );

	int kd_tree_get_seeds_in_sphere(double* x,                                                    // Sphere center
		                            double r,                                                     // neighborhood radius
		                            size_t max_index,                                             // maximum index to be considered
		                            size_t d_index, size_t node_index,                            // indices to traverse the kd-tree
		                            size_t &num_points_in_sphere, size_t* &points_in_sphere,      // number of points in sphere and their indices
		                            size_t &capacity                                              // Size of points in sphere array
	                               );

	int kd_tree_get_closest_seed(size_t seed_index,                                  // seed index
		                         size_t d_index, size_t node_index,                  // indices to traverse the kd-tree
		                         size_t &closest_seed, double &closest_distance      // index of closest seed and distance from it
	                            );

	int kd_tree_get_closest_seed(double* x,                                            // point coordinates
		                         size_t d_index, size_t node_index,                  // indices to traverse the kd-tree
		                         size_t &closest_seed, double &closest_distance,     // index of closest seed and distance from it
		                         size_t &num_nodes_visited                            // nodes visited
	                            );

	int kd_tree_get_closest_seed(double* x,                                            // point coordinates
		                         size_t d_index, size_t node_index,                  // indices to traverse the kd-tree
		                         size_t num_exculded_tree_points, size_t* exculded_tree_points,
		                         size_t &closest_seed, double &closest_distance,     // index of closest seed and distance from it
		                         size_t &num_nodes_visited                            // nodes visited
	                            );

	int add_entry(size_t entry, size_t &num_points_in_sphere, size_t* &points_in_sphere, size_t &capacity);


	bool find_brute(size_t entry, size_t* I, size_t num_entries);

	int quicksort(double* x, double* y1, double* y2, size_t* I, size_t left, size_t right);
    size_t _num_dim;
    size_t _capacity;
    std::vector<std::vector<double>> _points;
    std::vector<std::vector<double>> _points_normal;

    std::vector<double> _xmin;
    std::vector<double> _xmax;

    std::vector<std::vector<size_t>> _points_attrib;

    std::vector<std::vector<size_t>> _graph;

    size_t _tree_root;
    size_t _tree_max_height;
    std::vector<size_t> _tree_left;
    std::vector<size_t> _tree_right;
    
    bool _auto_balance;
    
    bool _marked_only;
    bool _is_redundant;
    std::vector<bool> _marked;

	size_t _lazy_deleted_since_rebuild = 0;
	size_t _lazy_delete_rebuild_threshold = 512;
	
	Geometry _geom;

};