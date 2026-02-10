#include "MeshingTree.h"

const double DST_TOL = 1E-6;

MeshingTree::MeshingTree()
{
    _xmin.assign(3,DBL_MAX);
    _xmax.assign(3,-DBL_MAX);
    _num_dim = 0;
    _capacity = 0;
    _tree_root = 0;
    _tree_max_height = 0;
    _auto_balance = true;
    _marked_only = false;
    _is_redundant = false;
}

int MeshingTree::graph_connect_nodes(size_t ipoint, size_t jpoint)
{
    graph_connect(ipoint, jpoint);   
    graph_connect(jpoint, ipoint);   
    return 0;
}

int MeshingTree::graph_connect(size_t ipoint, size_t jpoint)
{
    size_t max_index = std::max(ipoint, jpoint);
    if (_graph.size() <= max_index)
        _graph.resize(max_index + 1);
    if(graph_connected(ipoint, jpoint)) return 0;
    _graph[ipoint].push_back(jpoint);
    return 0;
}

bool MeshingTree::graph_connected(size_t ipoint, size_t jpoint)
{
    if (ipoint >= _graph.size() || jpoint >= _graph.size()) return false;
    
    // 检查 ipoint 的邻接列表中是否包含 jpoint
    for (size_t neighbor : _graph[ipoint])
    {
        if (neighbor == jpoint) return true;
    }
    return false;
}

int MeshingTree::add_tree_point(size_t num_dim, double* x, double* normal, size_t* attrib)
{
	
    _num_dim = static_cast<int>(num_dim);
    if (_xmin.size() < num_dim)
        _xmin.assign(num_dim, DBL_MAX);
    if (_xmax.size() < num_dim)
        _xmax.assign(num_dim, -DBL_MAX);

    for (size_t idim = 0; idim < num_dim; idim++)
    {
        _xmin[idim] = std::min(_xmin[idim], x[idim]);
        _xmax[idim] = std::max(_xmax[idim], x[idim]);
    }

    size_t new_index = _points.size();
    
    // 先添加点到 _points
    _points.emplace_back(x, x + num_dim);
    if (normal)
        _points_normal.emplace_back(normal, normal + num_dim);
    else
        _points_normal.emplace_back(num_dim, 0.0);

    if (attrib)
    {
        size_t num_attrib = attrib[0];
        _points_attrib.emplace_back(attrib, attrib + num_attrib);
    }
    else
        _points_attrib.emplace_back();

    // 再更新 kd-tree
    if (new_index == 0)
    {
        _tree_left.push_back(0);
        _tree_right.push_back(0);
        _tree_root = 0;
        _tree_max_height = 1;
    }
    else
    {
        _tree_left.push_back(new_index);
        _tree_right.push_back(new_index);
        kd_tree_add_point(new_index);
    }

    if (_graph.size() < _points.size())
        _graph.resize(_points.size());
    if (_marked.size() < _points.size())
        _marked.resize(_points.size(), true);
    else
        _marked.push_back(true);

    if (_auto_balance && _tree_max_height > log2(1 + _points.size()) * 10)
		kd_tree_build_balanced();
    return 0;
}

int MeshingTree::get_tree_point(size_t point_index, double* x)
{
    if (point_index >= _points.size()) return -1;
    for (size_t idim = 0; idim < _num_dim; idim++) x[idim] = _points[point_index][idim];
    return 0;
}

int MeshingTree::get_tree_point(size_t point_index, size_t num_dim, double* x)
{
    if (point_index >= _points.size()) return -1;
    for (size_t idim = 0; idim < num_dim; idim++) x[idim] = _points[point_index][idim];
    return 0;
}

bool MeshingTree::get_tree_point_attrib(size_t point_index, size_t attrib_index, size_t &point_attrib)
{
    if (point_index >= _points.size()) return false;
    if (_points_attrib[point_index].empty()) return false;
    point_attrib = _points_attrib[point_index][1 + attrib_index];
    return true;
}

double MeshingTree::get_tree_point_attrib(size_t point_index, size_t attrib_index)
{
	return _points[point_index][_num_dim + attrib_index];
}

int MeshingTree::set_tree_point_attrib(size_t point_index, size_t attrib_index, size_t attrib)
{
	_points_attrib[point_index][1 + attrib_index] = attrib;
	return 0;
}

// int MeshingTree::set_tree_point_attrib(size_t point_index, size_t attrib_index, double attrib)
// {
// 	_points[point_index][_num_dim + attrib_index] = attrib;
// 	return 0;
// }

void MeshingTree::clear_memory()
{
    _points.clear();
    _points_attrib.clear();
    _points_normal.clear();
    _xmax.clear();
    _xmin.clear();
    _marked.clear();
    _graph.clear();
    _tree_left.clear();
    _tree_right.clear();
    _tree_root = 0;
    _tree_max_height = 0;
    _marked_only = false;
    _auto_balance = true;
}

int MeshingTree::kd_tree_build_balanced()
{
    #pragma region Build Balanced kd-tree:
    std::vector<size_t> tree_nodes_sorted(_points.size());
    for (size_t i = 0; i < _points.size(); i++) tree_nodes_sorted[i] = i;
    for(size_t iseed = 0; iseed < _points.size(); iseed++)
    {
        _tree_left[iseed] = iseed;
        _tree_right[iseed] = iseed;
    }
    _tree_root = _points.size() - 1;

    size_t target_pos = _points.size() / 2;
    _tree_max_height = 0;
    kd_tree_balance_quicksort(target_pos, 0, _points.size() - 1, 0, tree_nodes_sorted.data());
    tree_nodes_sorted.clear();
    tree_nodes_sorted.shrink_to_fit();
    return 0;
    #pragma endregion
}

int MeshingTree::get_closest_tree_point(size_t tree_point_index, size_t &closest_tree_point, double &closest_distance)
{
	#pragma region Closest Neighbor Search using tree:

	closest_tree_point = _points.size();
	kd_tree_get_closest_seed(tree_point_index, 0, _tree_root, closest_tree_point, closest_distance);

	return 0;
	#pragma endregion
}

int MeshingTree::kd_tree_get_closest_seed(double* x,                                            // point coordinates
	                                     size_t d_index, size_t node_index,                  // indices to traverse the kd-tree
	                                     size_t &closest_seed, double &closest_distance,     // index of closest seed and distance from it
	                                     size_t &num_nodes_visited                            // nodes visited
	                                    )
{
	(void)d_index;
	(void)node_index;
	num_nodes_visited = 0;
	closest_seed = _points.size();
	closest_distance = DBL_MAX;
	for (size_t i = 0; i < _points.size(); i++)
	{
		if (_marked_only && i < _marked.size() && !_marked[i]) continue;
		double dst_sq = 0.0;
		for (size_t idim = 0; idim < _num_dim; idim++)
		{
			double dx = _points[i][idim] - x[idim];
			dst_sq += dx * dx;
		}
		if (dst_sq < closest_distance * closest_distance)
		{
			closest_seed = i;
			closest_distance = sqrt(dst_sq);
		}
		num_nodes_visited++;
	}
	return 0;
}

int MeshingTree::kd_tree_get_closest_seed(double* x,                                            // point coordinates
	                                     size_t d_index, size_t node_index,                  // indices to traverse the kd-tree
	                                     size_t num_exculded_tree_points, size_t* exculded_tree_points,
	                                     size_t &closest_seed, double &closest_distance,     // index of closest seed and distance from it
	                                     size_t &num_nodes_visited                            // nodes visited
	                                    )
{
	(void)d_index;
	(void)node_index;
	num_nodes_visited = 0;
	closest_seed = _points.size();
	closest_distance = DBL_MAX;
	for (size_t i = 0; i < _points.size(); i++)
	{
		if (_marked_only && i < _marked.size() && !_marked[i]) continue;
		if (find_brute(i, exculded_tree_points, num_exculded_tree_points)) continue;
		double dst_sq = 0.0;
		for (size_t idim = 0; idim < _num_dim; idim++)
		{
			double dx = _points[i][idim] - x[idim];
			dst_sq += dx * dx;
		}
		if (dst_sq < closest_distance * closest_distance)
		{
			closest_seed = i;
			closest_distance = sqrt(dst_sq);
		}
		num_nodes_visited++;
	}
	return 0;
}

int MeshingTree::get_closest_tree_point(double* x, size_t &closest_tree_point, double &closest_distance)
{
	#pragma region Closest Neighbor Search using tree:
	closest_tree_point = _points.size();
	size_t num_nodes_visited = 0;
	if (_points.size() == 0) return 1;
	kd_tree_get_closest_seed(x, 0, _tree_root, closest_tree_point, closest_distance, num_nodes_visited);
	return 0;
	#pragma endregion
}

int MeshingTree::get_closest_tree_point(double* x, size_t num_exculded_tree_points, size_t* exculded_tree_points, size_t &closest_tree_point, double &closest_distance)
{
	#pragma region Closest Neighbor Search using a tree:

	closest_tree_point = _points.size();
	size_t num_nodes_visited = 0;
	if (_points.size() == 0) return 1;
	kd_tree_get_closest_seed(x, 0, _tree_root, num_exculded_tree_points, exculded_tree_points, closest_tree_point, closest_distance, num_nodes_visited);
	return 0;
	#pragma endregion
}


int MeshingTree::get_tree_points_in_sphere(double* x, double r, size_t &num_points_in_sphere, size_t* &points_in_sphere, size_t &capacity)
{
	#pragma region tree neighbor search:
	num_points_in_sphere = 0;
	kd_tree_get_seeds_in_sphere(x, r, 0, _tree_root, num_points_in_sphere, points_in_sphere, capacity);
	return 0;
	#pragma endregion
}

int MeshingTree::get_normal_component(size_t num_dim, size_t num_basis, double** basis, double* vect, double &norm)
{
	#pragma region Get Normal component to some basis:

	double* comp = new double[num_basis];

	// project point to current basis
	for (size_t ibasis = 0; ibasis < num_basis; ibasis++)
	{
		comp[ibasis] = 0.0;
		for (size_t idim = 0; idim < num_dim; idim++) comp[ibasis] += vect[idim] * basis[ibasis][idim];
	}

	// get vector component orthogonal to current basis
	for (size_t ibasis = 0; ibasis < num_basis; ibasis++)
	{
		for (size_t idim = 0; idim < num_dim; idim++) vect[idim] -= comp[ibasis] * basis[ibasis][idim];
	}

	delete[] comp;

	norm = 0.0;
	for (size_t idim = 0; idim < num_dim; idim++) norm += vect[idim] * vect[idim];

	if (fabs(norm) < 1E-10) return 1;

	norm = 1.0 / sqrt(norm);
	for (size_t idim = 0; idim < num_dim; idim++) vect[idim] *= norm;

	return 0;
	#pragma endregion
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// kd-tree  Methods  /////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


int MeshingTree::kd_tree_balance_quicksort(size_t target_pos, size_t left, size_t right, size_t active_dim, size_t* tree_nodes_sorted)
{
	#pragma region kd tree balance:
	kd_tree_quicksort_adjust_target_position(target_pos, left, right, active_dim, tree_nodes_sorted);

	// target position is correct .. add to tree
	if (_tree_root == _capacity)    _tree_root = tree_nodes_sorted[target_pos];
	else                              kd_tree_add_point(tree_nodes_sorted[target_pos]);

	/* recursion */
	active_dim++;
	if (active_dim == _num_dim) active_dim = 0;

	if (target_pos > left + 1)  kd_tree_balance_quicksort((left + target_pos - 1) / 2, left, target_pos - 1, active_dim, tree_nodes_sorted);
	else if (left < target_pos) kd_tree_add_point(tree_nodes_sorted[left]);

	if (target_pos + 1 < right)  kd_tree_balance_quicksort((target_pos + 1 + right) / 2, target_pos + 1, right, active_dim, tree_nodes_sorted);
	else if (right > target_pos) kd_tree_add_point(tree_nodes_sorted[right]);

	return 0;
	#pragma endregion
}

int MeshingTree::kd_tree_quicksort_adjust_target_position(size_t target_pos, size_t left, size_t right, size_t active_dim, size_t* tree_nodes_sorted)
{
	#pragma region kd tree Quick sort pivot:
	size_t i = left, j = right;

	size_t pivot_seed = tree_nodes_sorted[(left + right) / 2];
	double pivot = _points[pivot_seed][active_dim];

	/* partition */
	while (i <= j)
	{
		while (_points[tree_nodes_sorted[i]][active_dim] < pivot)
			i++;
		while (_points[tree_nodes_sorted[j]][active_dim] > pivot)
			j--;
		if (i <= j)
		{
			size_t tmp_index = tree_nodes_sorted[i];
			tree_nodes_sorted[i] = tree_nodes_sorted[j];
			tree_nodes_sorted[j] = tmp_index;

			i++;
			if (j > 0) j--;
		}
	};

	/* recursion */
	if (j > 0 && left < j && left <= target_pos && j >= target_pos)
		kd_tree_quicksort_adjust_target_position(target_pos, left, j, active_dim, tree_nodes_sorted);
	if (i < right && i <= target_pos && right >= target_pos)
		kd_tree_quicksort_adjust_target_position(target_pos, i, right, active_dim, tree_nodes_sorted);

	return 0;
	#pragma endregion
}

int MeshingTree::kd_tree_add_point(size_t seed_index)
{
	#pragma region kd tree add point:
	// insert sphere into tree
	if (seed_index >= _points.size()) return 1;
	if (_points.empty()) return 1;
	if (_num_dim <= 0) return 1;
	size_t parent_index(_tree_root); size_t d_index(0);
	if (parent_index >= _points.size()) parent_index = 0;
	size_t branch_height(1);
	while (true)
	{
		if (_points[seed_index][d_index] >  _points[parent_index][d_index])
		{
			if (_tree_right[parent_index] == parent_index)
			{
				_tree_right[parent_index] = seed_index;
				branch_height++;
				break;
			}
			else
			{
				parent_index = _tree_right[parent_index];
				branch_height++;
			}
		}
		else
		{
			if (_tree_left[parent_index] == parent_index)
			{
				_tree_left[parent_index] = seed_index;
				branch_height++;
				break;
			}
			else
			{
				parent_index = _tree_left[parent_index];
				branch_height++;
			}
		}
		d_index++;
		if (d_index == _num_dim) d_index = 0;
	}
	if (branch_height > _tree_max_height) _tree_max_height = branch_height;
	return 0;
	#pragma endregion
}

int MeshingTree::kd_tree_get_seeds_in_sphere(double* x,                                                    // Sphere center
	                                       double r,                                                     // neighborhood radius
	                                       size_t d_index, size_t node_index,                            // indices to traverse the kd-tree
	                                       size_t &num_points_in_sphere, size_t* &points_in_sphere,      // number of points in sphere and their indices
	                                       size_t &capacity                                              // Size of points in sphere array
	                                      )
{
	#pragma region kd tree neighbor search:
	if (d_index == _num_dim) d_index = 0;

	if (!_marked_only || _marked[node_index])
	{
		double dst_sq(0.0);
		for (size_t idim = 0; idim < _num_dim; idim++)
		{
			double dx = _points[node_index][idim] - x[idim];
			dst_sq += dx * dx;
		}

		if (dst_sq < (1 + 2.0E-6) * r * r) add_entry(node_index, num_points_in_sphere, points_in_sphere, capacity);
	}

	bool check_right(false), check_left(false);
	double neighbor_min(x[d_index] - r), neighbor_max(x[d_index] + r);

	if (_tree_right[node_index] != node_index && neighbor_max >  _points[node_index][d_index]) check_right = true;
	if (_tree_left[node_index] != node_index && neighbor_min < _points[node_index][d_index]) check_left = true;

	if (check_right) kd_tree_get_seeds_in_sphere(x, r, d_index + 1, _tree_right[node_index], num_points_in_sphere, points_in_sphere, capacity);
	if (check_left)  kd_tree_get_seeds_in_sphere(x, r, d_index + 1, _tree_left[node_index], num_points_in_sphere, points_in_sphere, capacity);

	return 0;
	#pragma endregion
}

int MeshingTree::kd_tree_get_seeds_in_sphere(double* x,                                                    // Sphere center
	                                       double r,                                                     // neighborhood radius
	                                       size_t max_index,                                             // maximum index to be considered
	                                       size_t d_index, size_t node_index,                            // indices to traverse the kd-tree
	                                       size_t &num_points_in_sphere, size_t* &points_in_sphere,      // number of points in sphere and their indices
	                                       size_t &capacity                                              // Size of points in sphere array
)
{
	#pragma region kd tree neighbor search:
	if (d_index == _num_dim) d_index = 0;

	if (node_index <= max_index)
	{
		double dst_sq(0.0);
		for (size_t idim = 0; idim < _num_dim; idim++)
		{
			double dx = _points[node_index][idim] - x[idim];
			dst_sq += dx * dx;
		}
		if (dst_sq < r * r + DST_TOL) add_entry(node_index, num_points_in_sphere, points_in_sphere, capacity);
	}

	bool check_right(false), check_left(false);
	double neighbor_min(x[d_index] - r), neighbor_max(x[d_index] + r);

	if (_tree_right[node_index] != node_index && neighbor_max >  _points[node_index][d_index]) check_right = true;
	if (_tree_left[node_index] != node_index && neighbor_min < _points[node_index][d_index]) check_left = true;

	if (check_right) kd_tree_get_seeds_in_sphere(x, r, max_index,  d_index + 1, _tree_right[node_index], num_points_in_sphere, points_in_sphere, capacity);
	if (check_left)  kd_tree_get_seeds_in_sphere(x, r, max_index, d_index + 1, _tree_left[node_index], num_points_in_sphere, points_in_sphere, capacity);

	return 0;
	#pragma endregion
}

int MeshingTree::kd_tree_get_closest_seed(size_t seed_index,                                  // seed index
	                                    size_t d_index, size_t node_index,                  // indices to traverse the kd-tree
	                                    size_t &closest_seed, double &closest_distance      // index of closest seed and distance from it
	                                   )
{
	#pragma region kd tree closest neighbor search:
	if (d_index == _num_dim) d_index = 0;

	if (!_marked_only || _marked[node_index])
	{
		if (seed_index != node_index)
		{
			double dst_sq(0.0);
			for (size_t idim = 0; idim < _num_dim; idim++)
			{
				double dx = _points[seed_index][idim] - _points[node_index][idim];
				dst_sq += dx * dx;
			}

			if (dst_sq < closest_distance * closest_distance)
			{
				closest_seed = node_index;
				closest_distance = sqrt(dst_sq);
			}
		}
	}

	bool check_right(false), check_left(false);

	double neighbor_min, neighbor_max;
	if (closest_distance == DBL_MAX)
	{
		if (_tree_right[node_index] != node_index) check_right = true;
		if (_tree_left[node_index] != node_index) check_left = true;
	}
	else
	{
		neighbor_min = _points[seed_index][d_index] - closest_distance; neighbor_max = _points[seed_index][d_index] + closest_distance;
		if (_tree_right[node_index] != node_index && neighbor_max >  _points[node_index][d_index]) check_right = true;
		if (_tree_left[node_index] != node_index && neighbor_min < _points[node_index][d_index]) check_left = true;
	}


	if (check_right && check_left)
	{
		if (_points[seed_index][d_index] > _points[node_index][d_index])
		{
			kd_tree_get_closest_seed(seed_index, d_index + 1, _tree_right[node_index], closest_seed, closest_distance);
			neighbor_min = _points[seed_index][d_index] - closest_distance;
			if (neighbor_min < _points[node_index][d_index])
			{
				kd_tree_get_closest_seed(seed_index, d_index + 1, _tree_left[node_index], closest_seed, closest_distance);
			}
		}
		else
		{
			kd_tree_get_closest_seed(seed_index, d_index + 1, _tree_left[node_index], closest_seed, closest_distance);
			neighbor_max = _points[seed_index][d_index] + closest_distance;
			if (neighbor_max > _points[node_index][d_index])
			{
				kd_tree_get_closest_seed(seed_index, d_index + 1, _tree_right[node_index], closest_seed, closest_distance);
			}
		}
	}
	else if (check_right) kd_tree_get_closest_seed(seed_index, d_index + 1, _tree_right[node_index], closest_seed, closest_distance);
	else if (check_left)  kd_tree_get_closest_seed(seed_index, d_index + 1, _tree_left[node_index], closest_seed, closest_distance);

	return 0;
	#pragma endregion
}

bool MeshingTree::find_brute(size_t entry, size_t* I, size_t num_entries)
{
	#pragma region find using Brutal search:
	for (size_t i = 0; i < num_entries; i++) if (I[i] == entry) return true;
	return false;
	#pragma endregion
}

int MeshingTree::add_entry(size_t entry, size_t &num_points_in_sphere, size_t* &points_in_sphere, size_t &capacity)
{
	if (num_points_in_sphere >= capacity)
	{
		size_t new_capacity = capacity == 0 ? 8 : capacity * 2;
		size_t* new_points = new size_t[new_capacity];
		for (size_t i = 0; i < num_points_in_sphere; i++) new_points[i] = points_in_sphere[i];
		delete[] points_in_sphere;
		points_in_sphere = new_points;
		capacity = new_capacity;
	}
	points_in_sphere[num_points_in_sphere++] = entry;
	return 0;
}

int MeshingTree::quicksort(double* x, double* y1, double* y2, size_t* I, size_t left, size_t right)
{
	#pragma region Quick Sort:
	size_t i = left, j = right;
	double pivot = x[(left + right) / 2];

	/* partition */
	while (i <= j)
	{
		while (x[i] < pivot)
			i++;
		while (x[j] > pivot)
			j--;
		if (i <= j)
		{
			double tmp = x[i]; x[i] = x[j]; x[j] = tmp;

			tmp = y1[i]; y1[i] = y1[j]; y1[j] = tmp;

			tmp = y2[i]; y2[i] = y2[j]; y2[j] = tmp;

			size_t tmpi = I[i]; I[i] = I[j]; I[j] = tmpi;

			i++;
			if (j > 0) j--;
		}
	};

	/* recursion */

	if (j > 0 && left < j)
		quicksort(x, y1, y2, I, left, j);
	if (i < right)
		quicksort(x, y1, y2, I, i, right);

	return 0;
	#pragma endregion
}

int MeshingTree::lazy_delete_tree_point(const std::vector<double>& tree_point)
{
    if (tree_point.empty()) return -1;
    if (_points.empty()) return -1;
 
    for (size_t i = 0; i < _points.size(); i++)
    {
        if (_points[i].size() != tree_point.size()) continue;
 
        bool same = true;
        for (size_t j = 0; j < tree_point.size(); j++)
        {
            if (_points[i][j] != tree_point[j]) { same = false; break; }
        }
        if (same) return lazy_delete_tree_point_internal(i);
    }
    return -1;
}

int MeshingTree::lazy_delete_tree_point(size_t point_index)
{
	return lazy_delete_tree_point_internal(point_index);
}

bool MeshingTree::tree_point_is_active(size_t point_index) const
{
	if (point_index >= _points.size()) return false;
	if (!_marked_only) return true;
	if (point_index >= _marked.size()) return true;
	return _marked[point_index];
}
 
int MeshingTree::lazy_delete_tree_point_internal(size_t point_index)
{
    if (point_index >= _points.size()) return -1;
 
    // 第一次进入惰性删除模式：开启过滤，并把所有点先标记为“有效”
    if (!_marked_only)
    {
        _marked.assign(_points.size(), true);
        _marked_only = true;
        _lazy_deleted_since_rebuild = 0;
    }
    if (_marked.size() < _points.size()) _marked.resize(_points.size(), true);
 
    // 已经删除过
    if (!_marked[point_index]) return 1;
 
    _marked[point_index] = false;
    _lazy_deleted_since_rebuild++;
    return lazy_delete_tree_rebuild_if_needed();
}
 
int MeshingTree::lazy_delete_tree_rebuild_if_needed()
{
    // if (!_marked_only) return 0;
    // if (_lazy_deleted_since_rebuild < _lazy_delete_rebuild_threshold) return 0;
    // return lazy_delete_tree_rebuild_compact();
	return 0;
}
 
int MeshingTree::lazy_delete_tree_rebuild_compact()
{
    if (_points.empty())
    {
        _lazy_deleted_since_rebuild = 0;
        return 0;
    }
 
    // 1) compact 数据：只保留 _marked==true 的点
    const size_t old_n = _points.size();
    std::vector<size_t> old_to_new(old_n, static_cast<size_t>(-1));
 
    std::vector<std::vector<double>> new_points;
    std::vector<std::vector<double>> new_points_normal;
    std::vector<std::vector<size_t>> new_points_attrib;
    std::vector<std::vector<size_t>> new_graph;
 
    new_points.reserve(old_n);
    new_points_normal.reserve(_points_normal.size());
    new_points_attrib.reserve(_points_attrib.size());
 
    for (size_t i = 0; i < old_n; i++)
    {
        if (i < _marked.size() && _marked[i])
        {
            old_to_new[i] = new_points.size();
            new_points.push_back(_points[i]);
            if (i < _points_normal.size()) new_points_normal.push_back(_points_normal[i]);
            if (i < _points_attrib.size()) new_points_attrib.push_back(_points_attrib[i]);
        }
    }
 
    // 2) compact graph（可选但更完整）：重映射索引
    new_graph.resize(new_points.size());
    if (_graph.size() == old_n)
    {
        for (size_t i = 0; i < old_n; i++)
        {
            if (old_to_new[i] == static_cast<size_t>(-1)) continue;
            const size_t ni = old_to_new[i];
 
            for (size_t k = 0; k < _graph[i].size(); k++)
            {
                const size_t uj = _graph[i][k];
                if (uj >= old_n) continue;
                if (old_to_new[uj] == static_cast<size_t>(-1)) continue;

                new_graph[ni].push_back(old_to_new[uj]);
            }
        }
    }
 
    _points.swap(new_points);
    _points_normal.swap(new_points_normal);
    _points_attrib.swap(new_points_attrib);
    _graph.swap(new_graph);
 
    // 3) 重置 marked 状态（compact 后全部有效）
    _marked.assign(_points.size(), true);
    _marked_only = true;
    _lazy_deleted_since_rebuild = 0;
 
    // 4) 复用现有逻辑重构树：重置 tree 数组并 rebuild
    _tree_left.clear();
    _tree_right.clear();
    _tree_root = 0;
    _tree_max_height = 0;
 
    if (!_points.empty())
    {
        const size_t n = _points.size();
        _tree_left.resize(n);
        _tree_right.resize(n);
        for (size_t i = 0; i < n; i++)
        {
            _tree_left[i] = i;
            _tree_right[i] = i;
        }
 
        _tree_root = 0;
        _tree_max_height = 1;
        for (size_t i = 1; i < n; i++) kd_tree_add_point(i);
        // 如果你更想复用“平衡重构”，这里也可以直接调用 kd_tree_build_balanced();
    }
 
    // 5) 维护 xmin/xmax/_num_dim/_capacity
    if (!_points.empty())
    {
        const size_t dims = _points[0].size();
        _xmin.assign(dims, DBL_MAX);
        _xmax.assign(dims, -DBL_MAX);
 
        for (size_t i = 0; i < _points.size(); i++)
        {
            for (size_t d = 0; d < dims; d++)
            {
                _xmin[d] = std::min(_xmin[d], _points[i][d]);
                _xmax[d] = std::max(_xmax[d], _points[i][d]);
            }
        }
 
        if (_num_dim == 0) _num_dim = static_cast<int>(dims);
        _capacity = static_cast<int>(_points.size());
    }
 
    return 0;
}