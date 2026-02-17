#include "Generator.h"
#include <unordered_set>
#include <array>
#include <algorithm>
#include <sstream>
#include <iostream>
#include <string>
#include <iomanip>
#include <cfloat>

Generator::Generator()
{
}

Generator::~Generator()
{
}

int Generator::read_input_obj_file(std::string filename, size_t &num_points, double** &points, size_t &num_faces, size_t** &faces)
{
	#pragma region Reading Obj File:

	//open file
	std::ifstream tmpfile(filename.c_str());

	// count vertices and faces
	num_points = 0; num_faces = 0;
	if (tmpfile.is_open() && tmpfile.good())
	{
		std::string line = "";
		while (getline(tmpfile, line))
		{
			std::istringstream iss(line);
			std::vector<std::string> tokens;
			copy(std::istream_iterator<std::string>(iss), std::istream_iterator<std::string>(),
				std::back_inserter<std::vector<std::string> >(tokens));

			if (tokens.size() == 0) continue;

			if (tokens[0] == "v") num_points++;
			if (tokens[0] == "f") num_faces++;
		}
	}

	points = new double*[num_points];
	faces = new size_t*[num_faces];

	std::ifstream myfile(filename.c_str());

	num_points = 0; num_faces = 0;
	if (myfile.is_open() && myfile.good())
	{
		std::string line = "";
		while (getline(myfile, line))
		{
			std::istringstream iss(line);
			std::vector<std::string> tokens;
			copy(std::istream_iterator<std::string>(iss), std::istream_iterator<std::string>(),
				std::back_inserter<std::vector<std::string> >(tokens));

			if (tokens.size() == 0) continue;

			if (tokens[0] == "v")
			{
				points[num_points] = new double[3];
				for (size_t idim = 0; idim < 3; idim++) points[num_points][idim] = std::stod(tokens[1 + idim]);
				num_points++;
			}
			if (tokens[0] == "f")
			{

				for (size_t i = 1; i <= 3; i++)
				{
					std::size_t pos = tokens[i].find("/");
					tokens[i] = tokens[i].substr(0, pos);
				}

				faces[num_faces] = new size_t[4];
				faces[num_faces][0] = 3;
				for (size_t i = 1; i <= 3; i++) faces[num_faces][i] = static_cast<size_t>(std::stod(tokens[i]) - 1);

				bool redundant(false);
				if (false)
				{
					// check if redundant facet
					for (size_t jface = 0; jface < num_faces; jface++)
					{
						size_t num_common(0);
						for (size_t ii = 1; ii <= 3; ii++)
						{
							for (size_t jj = 1; jj <= 3; jj++)
							{
								if (faces[num_faces][ii] == faces[jface][jj]) num_common++;
							}
						}
						if (num_common == 3)
						{
							redundant = true;
							break;
						}
					}
				}

				if (redundant)
				{
					delete faces[num_faces];
					continue;
				}
				num_faces++;
			}
		}
	}

	// weld_nearby_points(num_points, points, num_faces, faces, 1E-8);
	// save_mesh_obj("clean.obj", num_points, points, num_faces, faces);

	std::cout << "  * Number of Input mesh points = " << num_points << std::endl;
	std::cout << "  * Number of Input mesh faces = " << num_faces << std::endl;

	

	#pragma endregion
	return 0;
}


void Generator::generate_spheres(const char* filename, MeshingTree *&spheres)
{
    std::ifstream fin(filename);
    if (!fin) {
        std::cerr << "open failed: " << filename << "\n";
        return;
    }

    std::string line;
    if (std::getline(fin, line)) { 
        if (line.find("x") == std::string::npos) {
            fin.clear();
            fin.seekg(0);
        }
    }

    size_t count = 0;
    while (std::getline(fin, line)) {
        if (line.empty()) {
            continue;
        }
        for (char &c : line) {
            if (c == ',') {
                c = ' ';
            }
        }
        std::istringstream ss(line);

        double x = 0.0, y = 0.0, z = 0.0, r = 0.0;
        int fid = 0;
        if (!(ss >> x >> y >> z >> r >> fid)) {
            continue;
        }
        double nx = 0.0, ny = 0.0, nz = 0.0;
        if (!(ss >> nx >> ny >> nz)) {
            nx = 0.0;
            ny = 0.0;
            nz = 0.0;
        }
        double s[4] = {x, y, z, r + 1e-4};
        double n[4] = {nx, ny, nz, 0.0};
        size_t attrib[3] = {2, 0, count + 1}; 
        attrib[1] = static_cast<size_t>(fid);
        spheres->add_tree_point(4, s, n, attrib);
        ++count;
    }
}

void Generator::read_obj_faces(const char* filename, std::vector<int>& faces_flat, size_t& num_faces)
{
	std::ifstream fin(filename);
	if (!fin)
	{
		std::cerr << "[Error] Cannot open file: " << filename << std::endl;
		num_faces = 0;
		faces_flat.clear();
		return;
	}

	faces_flat.clear();
	std::string line;
	while (std::getline(fin, line))
	{
		if (line.empty() || line[0] != 'f') continue;
		std::istringstream ss(line);
		std::string token;
		ss >> token; // skip 'f'
		std::array<int, 3> f{};
		int cnt = 0;
		while (cnt < 3 && (ss >> token))
		{
			size_t pos = token.find('/');
			int vid = 0;
			try {
				if (pos == std::string::npos) vid = std::stoi(token);
				else vid = std::stoi(token.substr(0, pos));
			} catch (...) {
				continue;
			}
			
			if (vid <= 0) continue; // OBJ indices are 1-based
			f[cnt] = vid - 1;
			cnt++;
		}
		if (cnt == 3)
		{
			faces_flat.push_back(f[0]);
			faces_flat.push_back(f[1]);
			faces_flat.push_back(f[2]);
		}
	}

	num_faces = static_cast<int>(faces_flat.size() / 3);
	std::cout << "[Info] Read " << num_faces << " faces from " << filename << std::endl;
}

void Generator::generate_surface_seeds(size_t num_points, double **points, size_t num_faces, size_t **faces,
         MeshingTree *spheres,MeshingTree *upper_seeds, MeshingTree *lower_seeds, size_t number_of_facets, std::vector<int> faces_flat)
{
    // Memory allocation
    double* sphere_i = new double[4];
    double* sphere_j = new double[4];
    double* sphere_k = new double[4];

    double* c_ijk = new double[3];
    double** centers = new double*[3];
    double* radii = new double[3];
    double* triplet_normal = new double[3];
    double* upper_seed = new double[4];
    double* lower_seed = new double[4];
    
    std::vector<size_t> attrib;
    attrib.resize(6);
    attrib[0] = 6, attrib[5] = 0;
    for(int i = 0; i < number_of_facets; i++)
    {
        int i_sphere = faces_flat[i * 3];
        int j_sphere = faces_flat[i * 3 + 1];
        int k_sphere = faces_flat[i * 3 + 2];

        // Bounds check
        size_t num_spheres = spheres->get_num_tree_points();
        if (i_sphere < 0 || (size_t)i_sphere >= num_spheres ||
            j_sphere < 0 || (size_t)j_sphere >= num_spheres ||
            k_sphere < 0 || (size_t)k_sphere >= num_spheres) {
            std::cerr << "[Error] Sphere index out of bounds in face " << i 
                      << ": " << i_sphere << ", " << j_sphere << ", " << k_sphere 
                      << " (total spheres: " << num_spheres << ")" << std::endl;
            continue;
        }

        bool failed = false;

        // 1. Load Spheres
        spheres->get_tree_point(i_sphere, 4, sphere_i);
        spheres->get_tree_point(j_sphere, 4, sphere_j);
        spheres->get_tree_point(k_sphere, 4, sphere_k);

        // 2. Pairwise overlapping checks
        // Relaxed check: Allow touching (tangent) spheres. Only fail if strictly separated.
        // Using +1E-10 allows for slight numerical error where dist is slightly > sum_radii
        double dst_ij = _geom.distance(3, sphere_i, sphere_j);
        if (dst_ij > sphere_i[3] + sphere_j[3] + 1E-10) failed = true;

        if (!failed) {
            double dst_ik = _geom.distance(3, sphere_i, sphere_k);
            if (dst_ik > sphere_i[3] + sphere_k[3] + 1E-10) failed = true;
        }
        if (!failed) {
            double dst_jk = _geom.distance(3, sphere_j, sphere_k);
            if (dst_jk > sphere_j[3] + sphere_k[3] + 1E-10) failed = true;
        }

        double vi = 0.0;
        
        // 3. Compute Radical Center (Power Vertex)
        if (!failed) {
            centers[0] = sphere_i; centers[1] = sphere_j; centers[2] = sphere_k;
            radii[0] = sphere_i[3]; radii[1] = sphere_j[3]; radii[2] = sphere_k[3];
            
            bool pv_valid = _geom.get_power_vertex(3, 3, centers, radii, c_ijk);
            
            if (!pv_valid) {
                failed = true;
            } else {
                double hi = _geom.distance(3, c_ijk, sphere_i);
                
                if (hi > sphere_i[3] - 1e-10) { 
                    failed = true;
                } else {
                    // Safe sqrt: if hi > radius (but within tolerance), clamp to 0 (tangent point)
                    vi = sqrt(fmax(0.0, sphere_i[3] * sphere_i[3] - hi * hi));
                    
                    if (vi < 1e-4) vi = 1e-4;
                }
            }
        }

        // 4. Check Collinearity via Area
        if (!failed) {
            double area = _geom.get_3d_triangle_area(centers);
            if (area < 1E-10) failed = true;
        }

        // --- SUCCESS: GENERATE SEEDS ---
        _geom.get_3d_triangle_normal(centers, triplet_normal);

        size_t sphere_i_face;
        spheres->get_tree_point_attrib(i_sphere, 0, sphere_i_face);
        size_t fi_i1 = faces[sphere_i_face][1];
        size_t fi_i2 = faces[sphere_i_face][2];
        size_t fi_i3 = faces[sphere_i_face][3];
        
        double* fi_corners[3] = { points[fi_i1], points[fi_i2], points[fi_i3] };
        double fi_normal[3];
        _geom.get_3d_triangle_normal(fi_corners, fi_normal);

        double dot = _geom.dot_product(3, fi_normal, triplet_normal);
        
        if (dot < 0.0)
        {
            attrib[2] = i_sphere; attrib[3] = k_sphere; attrib[4] = j_sphere;
            for (size_t idim = 0; idim < 3; idim++) triplet_normal[idim] = -triplet_normal[idim];
        }
        else
        {
            attrib[2] = i_sphere; attrib[3] = j_sphere; attrib[4] = k_sphere;
        }

        // Even if vi is very small (close points), we generate both.
        // Upper and Lower seeds are stored in separate trees, so they won't overwrite each other.
        for (size_t idim = 0; idim < 3; idim++)
        {
            upper_seed[idim] = c_ijk[idim] + vi * triplet_normal[idim];
            lower_seed[idim] = c_ijk[idim] - vi * triplet_normal[idim];
        }

        #pragma region Valid Intersection Pairs Storage
        auto add_current_seed_pair = [&](size_t& upper_seed_index,
                                         size_t& lower_seed_index,
                                         bool& upper_added,
                                         bool& lower_added)
        {
            size_t iclosest;
            double hclosest(DBL_MAX);
            upper_seed_index = upper_seeds->get_num_tree_points();
            lower_seed_index = lower_seeds->get_num_tree_points();
            upper_added = false;
            lower_added = false;

            // Keep duplicate suppression for stability with neighboring faces.
            lower_seeds->get_closest_tree_point(lower_seed, iclosest, hclosest);
            if (hclosest < 1E-10) lower_seed_index = iclosest;

            hclosest = DBL_MAX;
            upper_seeds->get_closest_tree_point(upper_seed, iclosest, hclosest);
            if (hclosest < 1E-10) upper_seed_index = iclosest;

            if (lower_seed_index == lower_seeds->get_num_tree_points())
            {
                lower_seed[3] = fmin(sphere_i[3], fmin(sphere_j[3], sphere_k[3]));
                attrib[1] = upper_seed_index;
                lower_seeds->add_tree_point(4, lower_seed, triplet_normal, attrib.data());
                lower_added = true;
            }

            if (upper_seed_index == upper_seeds->get_num_tree_points())
            {
                upper_seed[3] = fmin(sphere_i[3], fmin(sphere_j[3], sphere_k[3]));
                attrib[1] = lower_seed_index;
                upper_seeds->add_tree_point(4, upper_seed, triplet_normal, attrib.data());
                upper_added = true;
            }
        };

        // First add the currently computed seed pair.
        size_t upper_seed_index = 0, lower_seed_index = 0;
        bool upper_added = false, lower_added = false;
        add_current_seed_pair(upper_seed_index, lower_seed_index, upper_added, lower_added);

        // If the pair is too close, delete old seeds, enlarge radius, recompute, then add new seeds.
        double seed_pair_distance = _geom.distance(3, upper_seed, lower_seed);
        if (seed_pair_distance < 1e-4)
        {
            double original_upper_seed[4] = {upper_seed[0], upper_seed[1], upper_seed[2], upper_seed[3]};
            double original_lower_seed[4] = {lower_seed[0], lower_seed[1], lower_seed[2], lower_seed[3]};
            double original_triplet_normal[3] = {triplet_normal[0], triplet_normal[1], triplet_normal[2]};
            size_t original_attr2 = attrib[2], original_attr3 = attrib[3], original_attr4 = attrib[4];

            if (lower_added) lower_seeds->lazy_delete_tree_point(lower_seed_index);
            if (upper_added) upper_seeds->lazy_delete_tree_point(upper_seed_index);

            double* sphere_i_ptr = spheres->get_tree_point((size_t)i_sphere);
            double original_radius = sphere_i_ptr[3];
            const int max_expand_iters = 20;
            const double expand_factor = 1.05;
            bool fixed_pair = false;

            for (int expand_iter = 0; expand_iter < max_expand_iters; expand_iter++)
            {
                sphere_i_ptr[3] *= expand_factor;

                bool recompute_failed = false;
                spheres->get_tree_point(i_sphere, 4, sphere_i);
                spheres->get_tree_point(j_sphere, 4, sphere_j);
                spheres->get_tree_point(k_sphere, 4, sphere_k);

                double dst_ij = _geom.distance(3, sphere_i, sphere_j);
                if (dst_ij > sphere_i[3] + sphere_j[3] + 1E-10) recompute_failed = true;

                if (!recompute_failed)
                {
                    double dst_ik = _geom.distance(3, sphere_i, sphere_k);
                    if (dst_ik > sphere_i[3] + sphere_k[3] + 1E-10) recompute_failed = true;
                }

                if (!recompute_failed)
                {
                    double dst_jk = _geom.distance(3, sphere_j, sphere_k);
                    if (dst_jk > sphere_j[3] + sphere_k[3] + 1E-10) recompute_failed = true;
                }

                double vi_new = 0.0;
                if (!recompute_failed)
                {
                    centers[0] = sphere_i; centers[1] = sphere_j; centers[2] = sphere_k;
                    radii[0] = sphere_i[3]; radii[1] = sphere_j[3]; radii[2] = sphere_k[3];

                    bool pv_valid = _geom.get_power_vertex(3, 3, centers, radii, c_ijk);
                    if (!pv_valid)
                    {
                        recompute_failed = true;
                    }
                    else
                    {
                        double hi = _geom.distance(3, c_ijk, sphere_i);
                        if (hi > sphere_i[3] - 1e-10)
                        {
                            recompute_failed = true;
                        }
                        else
                        {
                            vi_new = sqrt(fmax(0.0, sphere_i[3] * sphere_i[3] - hi * hi));
                            if (vi_new < 1e-4) vi_new = 1e-4;
                        }
                    }
                }

                if (!recompute_failed)
                {
                    double area = _geom.get_3d_triangle_area(centers);
                    if (area < 1E-10) recompute_failed = true;
                }

                if (recompute_failed) continue;
                if (!_geom.get_3d_triangle_normal(centers, triplet_normal)) continue;

                size_t sphere_i_face;
                if (!spheres->get_tree_point_attrib(i_sphere, 0, sphere_i_face)) continue;
                if (sphere_i_face >= num_faces) continue;

                size_t fi_i1 = faces[sphere_i_face][1];
                size_t fi_i2 = faces[sphere_i_face][2];
                size_t fi_i3 = faces[sphere_i_face][3];
                double* fi_corners[3] = { points[fi_i1], points[fi_i2], points[fi_i3] };
                double fi_normal[3];
                if (!_geom.get_3d_triangle_normal(fi_corners, fi_normal)) continue;

                double dot = _geom.dot_product(3, fi_normal, triplet_normal);
                if (dot < 0.0)
                {
                    attrib[2] = i_sphere; attrib[3] = k_sphere; attrib[4] = j_sphere;
                    for (size_t idim = 0; idim < 3; idim++) triplet_normal[idim] = -triplet_normal[idim];
                }
                else
                {
                    attrib[2] = i_sphere; attrib[3] = j_sphere; attrib[4] = k_sphere;
                }

                for (size_t idim = 0; idim < 3; idim++)
                {
                    upper_seed[idim] = c_ijk[idim] + vi_new * triplet_normal[idim];
                    lower_seed[idim] = c_ijk[idim] - vi_new * triplet_normal[idim];
                }

                seed_pair_distance = _geom.distance(3, upper_seed, lower_seed);
                if (seed_pair_distance >= 1e-4)
                {
                    fixed_pair = true;
                    break;
                }
            }

            if (!fixed_pair)
            {
                sphere_i_ptr[3] = original_radius;

                for (size_t idim = 0; idim < 4; idim++) {
                    upper_seed[idim] = original_upper_seed[idim];
                    lower_seed[idim] = original_lower_seed[idim];
                }
                for (size_t idim = 0; idim < 3; idim++) {
                    triplet_normal[idim] = original_triplet_normal[idim];
                }
                attrib[2] = original_attr2;
                attrib[3] = original_attr3;
                attrib[4] = original_attr4;

                add_current_seed_pair(upper_seed_index, lower_seed_index, upper_added, lower_added);
                std::cerr << "[Warning] Face " << i
                          << " seed pair distance still < 1e-4 after radius expansion." << std::endl;
                continue;
            }

            // Add recomputed replacement seeds.
            add_current_seed_pair(upper_seed_index, lower_seed_index, upper_added, lower_added);
        }
        #pragma endregion
    }

    delete[] sphere_i; delete[] sphere_j; delete[] sphere_k;
    delete[] c_ijk; delete[] centers; delete[] radii;
    delete[] triplet_normal; delete[] upper_seed; delete[] lower_seed;

}

void Generator::generate_surface_seeds(size_t num_points, double **points, size_t num_faces, size_t **faces,
         MeshingTree *spheres,MeshingTree *upper_seeds, MeshingTree *lower_seeds)
{
    // Initialize sliver points OBJ file
    std::ofstream sliver_file("sliver_points.obj", std::ios::out);
    size_t sliver_vertex_offset = 1;
    if (sliver_file.is_open()) {
        sliver_file << std::fixed << std::setprecision(16);
        sliver_file << "# Sliver points detected by VoroCrust\n";
        sliver_file << "# Vertices include: triangle face corners, sphere centers, optional covering sphere, seeds\n";
        sliver_file << "# Format: v x y z r g b (vertex colors optional)\n";
    }

    size_t num_spheres = spheres->get_num_tree_points();
    double* sphere = new double[4];

    double* sphere_i = new double[4];
    double* sphere_j = new double[4];
    double* sphere_k = new double[4];

    double* c_ijk = new double[3];
    double** centers = new double*[3];
    double* radii = new double[3];
    double* triplet_normal = new double[3];
    double* upper_seed = new double[4];
    double* lower_seed = new double[4];
    std::vector<size_t> attrib;
    attrib.resize(6);
    attrib[0] = 6, attrib[5] = 0;
    for(size_t isphere = 0; isphere < num_spheres; isphere++)
    {
        size_t sphere_index_i(isphere);
        spheres->get_tree_point(isphere, 4, sphere_i);
        size_t sphere_i_face;
        spheres->get_tree_point_attrib(isphere, 0, sphere_i_face);
        size_t num_near_by_spheres(0);
        size_t *near_by_spheres(nullptr);

        _methods.get_overlapping_spheres(sphere_i, spheres, num_near_by_spheres, near_by_spheres);
        for(size_t jsphere = 0; jsphere < num_near_by_spheres; jsphere++)
        {
            size_t sphere_index_j(near_by_spheres[jsphere]);
            if(sphere_index_j <= sphere_index_i) continue;
            
            spheres->get_tree_point(sphere_index_j, 4, sphere_j);
            double dst_ij = _geom.distance(3, sphere_i, sphere_j);
            if (dst_ij > sphere_i[3] + sphere_j[3] - 1E-10)
                continue;
            
            for(size_t ksphere = 0; ksphere < num_near_by_spheres; ksphere++)
            {
                size_t sphere_index_k = near_by_spheres[ksphere];
                if (sphere_index_k <= sphere_index_i)
                  continue;
                if (sphere_index_k <= sphere_index_j)
                  continue;

                spheres->get_tree_point(sphere_index_k, 4, sphere_k);

                double dst_ik = _geom.distance(3, sphere_i, sphere_k);
                if (dst_ik > sphere_i[3] + sphere_k[3] - 1E-10)
                  continue;

                double dst_jk = _geom.distance(3, sphere_j, sphere_k);
                if (dst_jk > sphere_j[3] + sphere_k[3] - 1E-10)
                  continue;

                centers[0] = sphere_i; centers[1] = sphere_j; centers[2] = sphere_k;
                radii[0] = sphere_i[3]; radii[1] = sphere_j[3]; radii[2] = sphere_k[3];
                
                // Debug: check power vertex calculation
                bool pv_valid = _geom.get_power_vertex(3, 3, centers, radii, c_ijk);
                if (!pv_valid) {
                    std::cout << "Power vertex calculation failed for spheres " 
                              << sphere_index_i << ", " << sphere_index_j << ", " << sphere_index_k << std::endl;
                    continue;
                }
                
                // Debug: verify distances from power vertex to each sphere center
                double hi = _geom.distance(3, c_ijk, sphere_i);
                double hj = _geom.distance(3, c_ijk, sphere_j);
                double hk = _geom.distance(3, c_ijk, sphere_k);
                if (hi > sphere_i[3] - 1E-10) continue;

                double vi = sqrt(sphere_i[3] * sphere_i[3] - hi * hi);
                if (vi < 1e-4) vi = 1e-4;

                // Three overlapping spheres with an intersection pair
                double area = _geom.get_3d_triangle_area(centers);

                if (area < 1E-10)
                  continue; // spheres are colinear

                _geom.get_3d_triangle_normal(centers, triplet_normal);

                // adjust triplet normal

                size_t fi_i1 = faces[sphere_i_face][1];
                size_t fi_i2 = faces[sphere_i_face][2];
                size_t fi_i3 = faces[sphere_i_face][3];
                double** fi_corners = new double*[3];
                fi_corners[0] = points[fi_i1]; fi_corners[1] = points[fi_i2], fi_corners[2] = points[fi_i3];
                double* fi_normal = new double[3];
                _geom.get_3d_triangle_normal(fi_corners, fi_normal);

                double dot = _geom.dot_product(3, fi_normal, triplet_normal);
                delete[] fi_normal; delete[] fi_corners;
                
                // Check for sliver detection
               
                if (dot < 0.0)
                {
                  attrib[2] = sphere_index_i;
                  attrib[3] = sphere_index_k;
                  attrib[4] = sphere_index_j;
                  for (size_t idim = 0; idim < 3; idim++) triplet_normal[idim] = -triplet_normal[idim];
                }
                else
                {
                  attrib[2] = sphere_index_i;
                  attrib[3] = sphere_index_j;
                  attrib[4] = sphere_index_k;
                }

                for (size_t idim = 0; idim < 3; idim++)
                {
                  upper_seed[idim] = c_ijk[idim] + vi * triplet_normal[idim];
                  lower_seed[idim] = c_ijk[idim] - vi * triplet_normal[idim];
                }

				size_t num1 = 0, num2 = 0, cap = 100;
                size_t* covering_spheres1 = new size_t[cap],*covering_spheres2 = new size_t[cap];
                bool upper_covered = _methods.point_covered(upper_seed, spheres, 0.0, sphere_index_i, sphere_index_j, sphere_index_k, num1, cap, covering_spheres1);
                bool lower_covered = _methods.point_covered(lower_seed, spheres, 0.0, sphere_index_i, sphere_index_j, sphere_index_k, num2, cap, covering_spheres2);
                if(upper_covered || lower_covered){
                    if (sliver_file.is_open()) {
                        size_t cover_idx = SIZE_MAX;
                        if (upper_covered && num1 > 0) cover_idx = covering_spheres1[0];
                        else if (lower_covered && num2 > 0) cover_idx = covering_spheres2[0];

                        double* p1 = points[fi_i1];
                        double* p2 = points[fi_i2];
                        double* p3 = points[fi_i3];

                        size_t tri_base = sliver_vertex_offset;
                        sliver_file << "v " << p1[0] << " " << p1[1] << " " << p1[2] << " 0.8 0.8 0.8\n";
                        sliver_file << "v " << p2[0] << " " << p2[1] << " " << p2[2] << " 0.8 0.8 0.8\n";
                        sliver_file << "v " << p3[0] << " " << p3[1] << " " << p3[2] << " 0.8 0.8 0.8\n";
                        sliver_file << "f " << tri_base << " " << (tri_base + 1) << " " << (tri_base + 2) << "\n";
                        sliver_vertex_offset += 3;

                        sliver_file << "v " << sphere_i[0] << " " << sphere_i[1] << " " << sphere_i[2] << " 0.2 1.0 0.2\n";
                        sliver_file << "v " << sphere_j[0] << " " << sphere_j[1] << " " << sphere_j[2] << " 0.2 1.0 0.2\n";
                        sliver_file << "v " << sphere_k[0] << " " << sphere_k[1] << " " << sphere_k[2] << " 0.2 1.0 0.2\n";
                        sliver_vertex_offset += 3;

                        if (cover_idx != SIZE_MAX && cover_idx < num_spheres) {
                            double cover_sp[4];
                            spheres->get_tree_point(cover_idx, 4, cover_sp);
                            sliver_file << "v " << cover_sp[0] << " " << cover_sp[1] << " " << cover_sp[2] << " 1.0 0.6 0.1\n";
                            sliver_vertex_offset += 1;
                        }

                        sliver_file << "v " << upper_seed[0] << " " << upper_seed[1] << " " << upper_seed[2] << " 1.0 0.0 0.0\n";
                        sliver_file << "v " << lower_seed[0] << " " << lower_seed[1] << " " << lower_seed[2] << " 0.0 0.0 1.0\n";
                        sliver_vertex_offset += 2;
                    }
                    delete[] covering_spheres1; delete[] covering_spheres2; continue;
                }

                // Sliver resolution: shrink least-impactful sphere to minimum viable radius
                // if(upper_covered != lower_covered)
                // {
                //     bool sliver_triplet_invalid = false;
                //     const int max_sliver_iters = 5;
                //     for (int sliver_iter = 0; sliver_iter < max_sliver_iters; sliver_iter++)
                //     {
                //         // Identify the 4th sphere (the one covering a seed)
                //         size_t fourth_sphere = SIZE_MAX;
                //         if (upper_covered && num1 > 0) fourth_sphere = covering_spheres1[0];
                //         else if (lower_covered && num2 > 0) fourth_sphere = covering_spheres2[0];
                //         if (fourth_sphere == SIZE_MAX) break;

                //         // The four candidate spheres involved in the sliver
                //         size_t candidates[4] = {sphere_index_i, sphere_index_j, sphere_index_k, fourth_sphere};

                //         // Find the sphere with fewest overlapping neighbors (least impact)
                //         int best_idx = -1;
                //         size_t min_neighbor_count = SIZE_MAX;
                //         for (int c = 0; c < 4; c++)
                //         {
                //             double cs[4];
                //             spheres->get_tree_point(candidates[c], 4, cs);
                //             size_t n_overlap = 0;
                //             size_t* overlap_list = nullptr;
                //             _methods.get_overlapping_spheres(cs, spheres, n_overlap, overlap_list);
                //             if (n_overlap < min_neighbor_count) {
                //                 min_neighbor_count = n_overlap;
                //                 best_idx = c;
                //             }
                //             delete[] overlap_list;
                //         }
                //         if (best_idx < 0) break;

                //         // Compute the minimum radius that still intersects all neighbors
                //         size_t shrink_idx = candidates[best_idx];
                //         double shrink_sp[4];
                //         spheres->get_tree_point(shrink_idx, 4, shrink_sp);
                //         size_t n_overlap = 0;
                //         size_t* overlap_list = nullptr;
                //         _methods.get_overlapping_spheres(shrink_sp, spheres, n_overlap, overlap_list);

                //         double min_viable_radius = 1E-10;
                //         for (size_t n = 0; n < n_overlap; n++)
                //         {
                //             if (overlap_list[n] == shrink_idx) continue;
                //             double ns[4];
                //             spheres->get_tree_point(overlap_list[n], 4, ns);
                //             double d = _geom.distance(3, shrink_sp, ns);
                //             double required = d - ns[3] + 1E-10;
                //             if (required > min_viable_radius) min_viable_radius = required;
                //         }
                //         delete[] overlap_list;

                //         // Shrink the sphere radius
                //         spheres->get_tree_point(shrink_idx)[3] = min_viable_radius;

                //         // Reload triplet sphere data after shrinking
                //         spheres->get_tree_point(sphere_index_i, 4, sphere_i);
                //         spheres->get_tree_point(sphere_index_j, 4, sphere_j);
                //         spheres->get_tree_point(sphere_index_k, 4, sphere_k);

                //         // Check if triplet still overlaps after shrinking
                //         if (_geom.distance(3, sphere_i, sphere_j) > sphere_i[3] + sphere_j[3] - 1E-10 ||
                //             _geom.distance(3, sphere_i, sphere_k) > sphere_i[3] + sphere_k[3] - 1E-10 ||
                //             _geom.distance(3, sphere_j, sphere_k) > sphere_j[3] + sphere_k[3] - 1E-10)
                //         { sliver_triplet_invalid = true; break; }

                //         // Recompute power vertex
                //         centers[0] = sphere_i; centers[1] = sphere_j; centers[2] = sphere_k;
                //         radii[0] = sphere_i[3]; radii[1] = sphere_j[3]; radii[2] = sphere_k[3];
                //         if (!_geom.get_power_vertex(3, 3, centers, radii, c_ijk))
                //         { sliver_triplet_invalid = true; break; }

                //         double new_hi = _geom.distance(3, c_ijk, sphere_i);
                //         if (new_hi > sphere_i[3] - 1E-10)
                //         { sliver_triplet_invalid = true; break; }
                //         double new_vi = sqrt(sphere_i[3] * sphere_i[3] - new_hi * new_hi);

                //         if (_geom.get_3d_triangle_area(centers) < 1E-10)
                //         { sliver_triplet_invalid = true; break; }

                //         _geom.get_3d_triangle_normal(centers, triplet_normal);

                //         // Re-adjust triplet normal direction
                //         double** fc = new double*[3];
                //         fc[0] = points[fi_i1]; fc[1] = points[fi_i2]; fc[2] = points[fi_i3];
                //         double* fn = new double[3];
                //         _geom.get_3d_triangle_normal(fc, fn);
                //         double nd = _geom.dot_product(3, fn, triplet_normal);
                //         delete[] fn; delete[] fc;

                //         if (nd < 0.0) {
                //             attrib[2] = sphere_index_i; attrib[3] = sphere_index_k; attrib[4] = sphere_index_j;
                //             for (size_t idim = 0; idim < 3; idim++) triplet_normal[idim] = -triplet_normal[idim];
                //         } else {
                //             attrib[2] = sphere_index_i; attrib[3] = sphere_index_j; attrib[4] = sphere_index_k;
                //         }

                //         for (size_t idim = 0; idim < 3; idim++)
                //         {
                //             upper_seed[idim] = c_ijk[idim] + new_vi * triplet_normal[idim];
                //             lower_seed[idim] = c_ijk[idim] - new_vi * triplet_normal[idim];
                //         }

                //         // Re-check coverage after recomputation
                //         delete[] covering_spheres1; delete[] covering_spheres2;
                //         num1 = 0; num2 = 0;
                //         covering_spheres1 = new size_t[cap]; covering_spheres2 = new size_t[cap];
                //         upper_covered = _methods.point_covered(upper_seed, spheres, 0.0, sphere_index_i, sphere_index_j, sphere_index_k, num1, cap, covering_spheres1);
                //         lower_covered = _methods.point_covered(lower_seed, spheres, 0.0, sphere_index_i, sphere_index_j, sphere_index_k, num2, cap, covering_spheres2);

                //         if (upper_covered && lower_covered) break;
                //         if (upper_covered == lower_covered) break; // Both uncovered - sliver resolved
                //     }

                //     if (sliver_triplet_invalid || (upper_covered && lower_covered))
                //     { delete[] covering_spheres1; delete[] covering_spheres2; continue; }
                // }
                delete[] covering_spheres1; delete[] covering_spheres2;
                #pragma region A valid Intersection Pairs:
                size_t iclosest; double hclosest(DBL_MAX);
                size_t upper_seed_index(upper_seeds->get_num_tree_points());
                size_t lower_seed_index(lower_seeds->get_num_tree_points());

                lower_seeds->get_closest_tree_point(lower_seed, iclosest, hclosest);
                if (hclosest < 1E-10) lower_seed_index = iclosest;

                hclosest = DBL_MAX;
                upper_seeds->get_closest_tree_point(upper_seed, iclosest, hclosest);
                if (hclosest < 1E-10) upper_seed_index = iclosest;

                if (lower_seed_index == lower_seeds->get_num_tree_points())
                {
                  lower_seed[3] = fmin(sphere_i[3], sphere_j[3]);
                  lower_seed[3] = fmin(lower_seed[3], sphere_k[3]);

                  attrib[1] = upper_seed_index;
                  lower_seeds->add_tree_point(4, lower_seed, triplet_normal, attrib.data());
                }

                if (upper_seed_index == upper_seeds->get_num_tree_points())
                {
                  upper_seed[3] = fmin(sphere_i[3], sphere_j[3]);
                  upper_seed[3] = fmin(upper_seed[3], sphere_k[3]);

                  attrib[1] = lower_seed_index;
                  upper_seeds->add_tree_point(4, upper_seed, triplet_normal, attrib.data());
                }
                #pragma endregion
            }
        }
    }

    delete[] sphere;
    delete[] sphere_i; delete[] sphere_j; delete[] sphere_k;
    delete[] c_ijk; delete[] centers; delete[] radii;
    delete[] triplet_normal; delete[] upper_seed; delete[] lower_seed;
    attrib.clear();
    attrib.shrink_to_fit();

    double min_dst_uu(DBL_MAX), min_dst_ll(DBL_MAX), min_dst_ul(DBL_MAX);
		size_t num_upper_seeds = upper_seeds->get_num_tree_points();
		size_t num_lower_seeds = lower_seeds->get_num_tree_points();
		size_t num_spheres_s = spheres->get_num_tree_points();
		for (size_t iseed = 0; iseed < num_upper_seeds; iseed++)
		{
        size_t* attrib = upper_seeds->get_tree_point_attrib(iseed);
			size_t si(attrib[2]), sj(attrib[3]), sk(attrib[4]);
			double min_radius(DBL_MAX);
			double sphere_buf[4];
			if (si < num_spheres_s)
			{
				spheres->get_tree_point(si, 4, sphere_buf);
				if (sphere_buf[3] < min_radius) min_radius = sphere_buf[3];
			}
			if (sj < num_spheres_s)
			{
				spheres->get_tree_point(sj, 4, sphere_buf);
				if (sphere_buf[3] < min_radius) min_radius = sphere_buf[3];
			}
			if (sk < num_spheres_s)
			{
				spheres->get_tree_point(sk, 4, sphere_buf);
				if (sphere_buf[3] < min_radius) min_radius = sphere_buf[3];
			}

			size_t iclosest(iseed); double hclosest(DBL_MAX);
			upper_seeds->get_closest_tree_point(iseed, iclosest, hclosest);
			hclosest /= min_radius;

			if (hclosest < min_dst_uu) min_dst_uu = hclosest;

			double seed_buf[4];
			upper_seeds->get_tree_point(iseed, 4, seed_buf);
			size_t jclosest = num_lower_seeds; double vclosest = DBL_MAX;
			lower_seeds->get_closest_tree_point(seed_buf, jclosest, vclosest);
			vclosest /= min_radius;
			if (vclosest < min_dst_ul) min_dst_ul = vclosest;
		}

		for (size_t iseed = 0; iseed < num_lower_seeds; iseed++)
		{
        size_t* attrib = lower_seeds->get_tree_point_attrib(iseed);
			size_t si(attrib[2]), sj(attrib[3]), sk(attrib[4]);
			double min_radius(DBL_MAX);
			double sphere_buf[4];
			if (si < num_spheres_s)
			{
				spheres->get_tree_point(si, 4, sphere_buf);
				if (sphere_buf[3] < min_radius) min_radius = sphere_buf[3];
			}
			if (sj < num_spheres_s)
			{
				spheres->get_tree_point(sj, 4, sphere_buf);
				if (sphere_buf[3] < min_radius) min_radius = sphere_buf[3];
			}
			if (sk < num_spheres_s)
			{
				spheres->get_tree_point(sk, 4, sphere_buf);
				if (sphere_buf[3] < min_radius) min_radius = sphere_buf[3];
			}

			size_t iclosest(iseed); double hclosest(DBL_MAX);
			lower_seeds->get_closest_tree_point(iseed, iclosest, hclosest);
			hclosest /= min_radius;
			if (hclosest < min_dst_ll) min_dst_ll = hclosest;

			double seed_buf[4];
			lower_seeds->get_tree_point(iseed, 4, seed_buf);
			size_t jclosest = num_lower_seeds; double vclosest = DBL_MAX;
			upper_seeds->get_closest_tree_point(seed_buf, jclosest, vclosest);
			vclosest /= min_radius;
			if (vclosest < min_dst_ul) min_dst_ul = vclosest;
		}
		
		std::cout << "  * Min. relative distance between upper seeds = " << min_dst_uu << std::endl;
		std::cout << "  * Min. relative distance between lower seeds = " << min_dst_ll << std::endl;
		std::cout << "  * Min. relative distance between lower and upper seeds = " << min_dst_ul << std::endl;

    // Export all seeds with links to their three generating sphere centers
    {
        std::ofstream seed_links_file("seeds_with_links.obj");
        if (seed_links_file.is_open()) {
            seed_links_file << std::fixed << std::setprecision(16);
            seed_links_file << "# All seed points with links to generating sphere centers\n";

            size_t vertex_offset = 1;

            auto write_seeds_with_links = [&](MeshingTree* seed_tree, const char* label) {
                size_t n = seed_tree->get_num_tree_points();
                for (size_t i = 0; i < n; ++i) {
                    if (!seed_tree->tree_point_is_active(i)) continue;

                    double* pt = seed_tree->get_tree_point(i);
                    size_t* att = seed_tree->get_tree_point_attrib(i);
                    size_t seed_vi = vertex_offset;

                    seed_links_file << "v " << pt[0] << " " << pt[1] << " " << pt[2] << "\n";
                    ++vertex_offset;

                    size_t sids[3] = { att[2], att[3], att[4] };
                    for (size_t si : sids) {
                        if (si >= num_spheres_s) continue;
                        double sp[4];
                        spheres->get_tree_point(si, 4, sp);
                        size_t sp_vi = vertex_offset;
                        seed_links_file << "v " << sp[0] << " " << sp[1] << " " << sp[2] << "\n";
                        seed_links_file << "l " << seed_vi << " " << sp_vi << "\n";
                        ++vertex_offset;
                    }
                }
            };

            write_seeds_with_links(upper_seeds, "upper");
            write_seeds_with_links(lower_seeds, "lower");

            seed_links_file.close();
            std::cout << "  * Exported seeds with sphere links to seeds_with_links.obj" << std::endl;
        }
    }
}


void Generator::color_surface_seeds(int num_faces, MeshingTree *surface_spheres, MeshingTree *upper_seeds, MeshingTree *lower_seeds, MeshingTree *seeds,
    std::vector<int> face, double** points, double* &seedes, size_t* &seeds_region_id, double* &seeds_sizing)
{
    std::cout << "[Info] Starting robust color_surface_seeds..." << std::endl;
    size_t num_upper_seeds(upper_seeds->get_num_tree_points());
    size_t num_lower_seeds(lower_seeds->get_num_tree_points());

    // Adjust id of seed pair for upper seeds
    for (size_t iseed = 0; iseed < num_upper_seeds; iseed++)
    {
        double x[4];
        upper_seeds->get_tree_point(iseed, 4, x);
        double* normal = upper_seeds->get_tree_point_normal(iseed);
        size_t* attrib = upper_seeds->get_tree_point_attrib(iseed);
        attrib[1] += num_upper_seeds; // Update pair index
        seeds->add_tree_point(4, x, normal, attrib);
    }
    upper_seeds->clear_memory();

    for (size_t iseed = 0; iseed < num_lower_seeds; iseed++)
    {
        double x[4];
        lower_seeds->get_tree_point(iseed, 4, x);
        double* normal = lower_seeds->get_tree_point_normal(iseed);
        size_t* attrib = lower_seeds->get_tree_point_attrib(iseed);
        seeds->add_tree_point(4, x, normal, attrib);
    }
    lower_seeds->clear_memory();

    // 2. Filter seeds based on OBJ faces (if provided)
    size_t num_seeds(seeds->get_num_tree_points());
    if (num_faces > 0 && !face.empty())
    {
        struct FaceKeyHash {
            size_t operator()(const std::array<size_t, 3>& a) const noexcept {
                size_t h = 0;
                h ^= std::hash<size_t>{}(a[0]) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
                h ^= std::hash<size_t>{}(a[1]) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
                h ^= std::hash<size_t>{}(a[2]) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
                return h;
            }
        };

        std::unordered_set<std::array<size_t, 3>, FaceKeyHash> face_set;
        face_set.reserve(static_cast<size_t>(num_faces) * 2);
        for (int i = 0; i < num_faces; i++)
        {
            size_t a = static_cast<size_t>(face[3 * i + 0]);
            size_t b = static_cast<size_t>(face[3 * i + 1]);
            size_t c = static_cast<size_t>(face[3 * i + 2]);
            std::array<size_t, 3> key{ a, b, c };
            std::sort(key.begin(), key.end());
            face_set.insert(key);
        }

        for (size_t iseed = 0; iseed < num_seeds; iseed++)
        {
            if (!seeds->tree_point_is_active(iseed)) continue;
            size_t* attrib = seeds->get_tree_point_attrib(iseed);
            const size_t jseed = attrib[1];
            if (jseed >= num_seeds) continue;
            if (iseed > jseed) continue; // Only process one of the pair
            if (!seeds->tree_point_is_active(jseed)) continue;

            std::array<size_t, 3> key{ attrib[2], attrib[3], attrib[4] };
            std::sort(key.begin(), key.end());
            if(face_set.find(key) == face_set.end())
            {
                seeds->lazy_delete_tree_point(iseed);
                seeds->lazy_delete_tree_point(jseed);
            }
        }
    }

    // 3. Color seeds using surface sphere normals (outside=0, inside=1)
    size_t num_spheres = surface_spheres->get_num_tree_points();
    auto compute_seed_dot = [&](size_t seed_index) -> double {
        double seed[4];
        seeds->get_tree_point(seed_index, 4, seed);
        size_t* attrib = seeds->get_tree_point_attrib(seed_index);
        size_t sphere_indices[3] = {attrib[2], attrib[3], attrib[4]};

        double dot_sum = 0.0;
        int count = 0;
        for (size_t i = 0; i < 3; ++i) {
            size_t si = sphere_indices[i];
            if (si >= num_spheres) continue;
            double center[4];
            surface_spheres->get_tree_point(si, 4, center);
            double* normal = surface_spheres->get_tree_point_normal(si);
            double nx = normal[0];
            double ny = normal[1];
            double nz = normal[2];
            double nlen_sq = nx * nx + ny * ny + nz * nz;
            if (nlen_sq < 1e-20) continue;
            double dx = seed[0] - center[0];
            double dy = seed[1] - center[1];
            double dz = seed[2] - center[2];
            dot_sum += dx * nx + dy * ny + dz * nz;
            count++;
        }
        if (count == 0) return 0.0;
        return dot_sum / static_cast<double>(count);
    };

    for (size_t iseed = 0; iseed < num_seeds; ++iseed)
    {
        if (!seeds->tree_point_is_active(iseed)) continue;
        size_t* attrib = seeds->get_tree_point_attrib(iseed);
        size_t jseed = attrib[1];
        if (jseed >= num_seeds || !seeds->tree_point_is_active(jseed))
        {
            double dot = compute_seed_dot(iseed);
            attrib[5] = (dot < 0.0) ? 1 : 0;
            continue;
        }
        if (iseed > jseed) continue;

        double dot_i = compute_seed_dot(iseed);
        double dot_j = compute_seed_dot(jseed);
        size_t outside_seed = (dot_i >= dot_j) ? iseed : jseed;
        size_t inside_seed = (outside_seed == iseed) ? jseed : iseed;
        seeds->get_tree_point_attrib(outside_seed)[5] = 0;
        seeds->get_tree_point_attrib(inside_seed)[5] = 1;
    }
    
}
// void Generator::color_surface_seeds(int num_faces, MeshingTree *surface_spheres, MeshingTree *upper_seeds, MeshingTree *lower_seeds, MeshingTree *seeds,
//     std::vector<int> face, double** points, double* &seedes, size_t* &seeds_region_id, double* &seeds_sizing)
// {
//     std::cout << "[Info] Starting robust color_surface_seeds..." << std::endl;

//     // 1. Merge Upper and Lower seeds into 'seeds' tree
//     size_t num_upper_seeds(upper_seeds->get_num_tree_points());
//     size_t num_lower_seeds(lower_seeds->get_num_tree_points());

//     // Adjust id of seed pair for upper seeds
//     for (size_t iseed = 0; iseed < num_upper_seeds; iseed++)
//     {
//         double x[4];
//         upper_seeds->get_tree_point(iseed, 4, x);
//         double* normal = upper_seeds->get_tree_point_normal(iseed);
//         size_t* attrib = upper_seeds->get_tree_point_attrib(iseed);
//         attrib[1] += num_upper_seeds; // Update pair index
//         seeds->add_tree_point(4, x, normal, attrib);
//     }
//     upper_seeds->clear_memory();

//     for (size_t iseed = 0; iseed < num_lower_seeds; iseed++)
//     {
//         double x[4];
//         lower_seeds->get_tree_point(iseed, 4, x);
//         double* normal = lower_seeds->get_tree_point_normal(iseed);
//         size_t* attrib = lower_seeds->get_tree_point_attrib(iseed);
//         seeds->add_tree_point(4, x, normal, attrib);
//     }
//     lower_seeds->clear_memory();

//     // 2. Filter seeds based on OBJ faces (if provided)
//     size_t num_seeds(seeds->get_num_tree_points());
//     if (num_faces > 0 && !face.empty())
//     {
//         struct FaceKeyHash {
//             size_t operator()(const std::array<size_t, 3>& a) const noexcept {
//                 size_t h = 0;
//                 h ^= std::hash<size_t>{}(a[0]) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
//                 h ^= std::hash<size_t>{}(a[1]) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
//                 h ^= std::hash<size_t>{}(a[2]) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
//                 return h;
//             }
//         };

//         std::unordered_set<std::array<size_t, 3>, FaceKeyHash> face_set;
//         face_set.reserve(static_cast<size_t>(num_faces) * 2);
//         for (int i = 0; i < num_faces; i++)
//         {
//             size_t a = static_cast<size_t>(face[3 * i + 0]);
//             size_t b = static_cast<size_t>(face[3 * i + 1]);
//             size_t c = static_cast<size_t>(face[3 * i + 2]);
//             std::array<size_t, 3> key{ a, b, c };
//             std::sort(key.begin(), key.end());
//             face_set.insert(key);
//         }

//         for (size_t iseed = 0; iseed < num_seeds; iseed++)
//         {
//             if (!seeds->tree_point_is_active(iseed)) continue;
//             size_t* attrib = seeds->get_tree_point_attrib(iseed);
//             const size_t jseed = attrib[1];
//             if (jseed >= num_seeds) continue;
//             if (iseed > jseed) continue; // Only process one of the pair
//             if (!seeds->tree_point_is_active(jseed)) continue;

//             std::array<size_t, 3> key{ attrib[2], attrib[3], attrib[4] };
//             std::sort(key.begin(), key.end());
//             // if(key[0]==723 && key[1]==778 && key[2]==1989){
//             //     std::cout<<"Found sliver!"<<std::endl;
//             // }    
//             if(face_set.find(key) == face_set.end())
//             {
//                 // std::cout<<key[0]<<" "<<key[1]<<" "<<key[2]<<std::endl;
//                 seeds->lazy_delete_tree_point(iseed);
//                 seeds->lazy_delete_tree_point(jseed);
//             }
//         }
//     }

//     // 3. Build connectivity graph
//     for (size_t iseed = 0; iseed < num_seeds; iseed++)
//     {
//         if (!seeds->tree_point_is_active(iseed)) continue;
//         size_t* attrib = seeds->get_tree_point_attrib(iseed);
//         size_t si(attrib[2]), sj(attrib[3]), sk(attrib[4]);
//         surface_spheres->graph_connect(si, iseed);
//         surface_spheres->graph_connect(sj, iseed);
//         surface_spheres->graph_connect(sk, iseed);
//     }

//     double* face_normal = new double[3];
//     double** face_corners = new double*[3];
//     size_t num_spheres(surface_spheres->get_num_tree_points());

//     // 4. Propagate Orientation (The Vorocrust Method)
//     for (size_t isphere = 0; isphere < num_spheres; isphere++)
//     {
//         std::vector<size_t> sphere_seeds;
//         surface_spheres->graph_get_neighbors(isphere, sphere_seeds);
//         if (sphere_seeds.empty()) continue;

//         size_t num_sphere_seeds = sphere_seeds.size();

//         // [Step A] Reset markers for this sphere
//         for (size_t i = 0; i < num_sphere_seeds; i++)
//         {
//             size_t seed_i(sphere_seeds[i]);
//             size_t* attrib = seeds->get_tree_point_attrib(seed_i);
//             attrib[5] = 0; // 0 means unvisited/undefined for this local loop
//         }

//         size_t first_neighbor(isphere), jsphere(isphere);

//         // [Step B] Pick the first seed to dictate orientation
//         if (true)
//         {
//             size_t seed_1(sphere_seeds[0]); // Pick the first valid one
//             size_t* attrib = seeds->get_tree_point_attrib(seed_1);

//             size_t seed_2(attrib[1]);
//             size_t si(attrib[2]), sj(attrib[3]), sk(attrib[4]);

//             // [CRITICAL] Rotate indices so 'si' is the current sphere
//             while (si != isphere)
//             {
//                 size_t tmp = si; si = sj; sj = sk; sk = tmp;
//             }

//             double* x_1 = seeds->get_tree_point(seed_1);
//             double* x_2 = seeds->get_tree_point(seed_2);

//             face_corners[0] = surface_spheres->get_tree_point(si);
//             face_corners[1] = surface_spheres->get_tree_point(sj);
//             face_corners[2] = surface_spheres->get_tree_point(sk);

//             _geom.get_3d_triangle_normal(face_corners, face_normal);

//             double dot(0.0);
//             for (size_t idim = 0; idim < 3; idim++) dot += (x_2[idim] - x_1[idim]) * face_normal[idim];

//             if (dot > 0.0)
//             {
//                 attrib[5] = 1; // 1: Inside
//                 attrib = seeds->get_tree_point_attrib(seed_2);
//                 attrib[5] = 2; // 2: Outside
//             }
//             else
//             {
//                 attrib[5] = 2; // Outside
//                 attrib = seeds->get_tree_point_attrib(seed_2);
//                 attrib[5] = 1; // Inside
//             }
//             first_neighbor = sj; // Record where we started
//             jsphere = sk;        // The vertex we are moving towards
//         }

//         // [Step C] Walk the loop (Topological Propagation)
//         while (true)
//         {
//             bool done(true);

//             // Look for the next neighbor in the chain
//             for (size_t i = 0; i < num_sphere_seeds; i++)
//             {
//                 size_t seed_1(sphere_seeds[i]);
//                 size_t* attrib = seeds->get_tree_point_attrib(seed_1);

//                 if (attrib[5] != 0) continue; // Already processed

//                 size_t seed_2(attrib[1]);
//                 size_t si(attrib[2]), sj(attrib[3]), sk(attrib[4]);

//                 // [CRITICAL] Rotate indices again
//                 while (si != isphere)
//                 {
//                     size_t tmp = si; si = sj; sj = sk; sk = tmp;
//                 }

//                 // We are looking for a triangle that shares the edge (si, jsphere)
//                 // In a valid traversal, one of sj or sk MUST be jsphere.
//                 if (sj != jsphere && sk != jsphere) continue; // Not the neighbor we want

//                 // [CRITICAL FIX] Enforce winding order
//                 // If sj is NOT the shared vertex, it means the triangle is flipped relative to our traversal.
//                 // We MUST swap sj and sk to maintain consistent normal direction.
//                 if (sj != jsphere)
//                 {
//                     size_t tmp = sj; sj = sk; sk = tmp;
//                 }

//                 // Now calculate normal with the CONSISTENT vertex order (si, sj, sk)
//                 double* x_1 = seeds->get_tree_point(seed_1);
//                 double* x_2 = seeds->get_tree_point(seed_2);

//                 face_corners[0] = surface_spheres->get_tree_point(si);
//                 face_corners[1] = surface_spheres->get_tree_point(sj);
//                 face_corners[2] = surface_spheres->get_tree_point(sk);

//                 _geom.get_3d_triangle_normal(face_corners, face_normal);

//                 double dot(0.0);
//                 for (size_t idim = 0; idim < 3; idim++) dot += (x_2[idim] - x_1[idim]) * face_normal[idim];

//                 if (dot > 0.0)
//                 {
//                     attrib[5] = 1; // Inside
//                     attrib = seeds->get_tree_point_attrib(seed_2);
//                     attrib[5] = 2; // Outside
//                 }
//                 else
//                 {
//                     attrib[5] = 2; // Outside
//                     attrib = seeds->get_tree_point_attrib(seed_2);
//                     attrib[5] = 1; // Inside
//                 }

//                 jsphere = sk; // Advance to the next vertex
//                 if (jsphere == first_neighbor) done = true; // Loop closed
                
//                 done = false; // Found a valid step, continue the while loop
//                 break; 
//             }
//             if (done) break; // No more neighbors found or loop closed
//         }

//         // [Step D] Connect compatible seeds
//         // Now that seeds on this sphere are locally consistent, connect them in the graph
//         for (size_t i = 0; i < num_sphere_seeds; i++)
//         {
//             size_t seed_i(sphere_seeds[i]);
//             size_t* attrib_i = seeds->get_tree_point_attrib(seed_i);
            
//             // Only connect processed seeds
//             if(attrib_i[5] == 0) continue;

//             for (size_t j = i + 1; j < num_sphere_seeds; j++)
//             {
//                 size_t seed_j(sphere_seeds[j]);
//                 size_t* attrib_j = seeds->get_tree_point_attrib(seed_j);

//                 if (attrib_i[5] == attrib_j[5])
//                 {
//                     seeds->graph_connect_nodes(seed_i, seed_j);
//                 }
//             }
//         }

//         // [Step E] Reset for global coloring
//         // We wipe the local 1/2 markers because we will re-flood globally
//         for (size_t i = 0; i < num_sphere_seeds; i++)
//         {
//             size_t seed_i(sphere_seeds[i]);
//             size_t* attrib = seeds->get_tree_point_attrib(seed_i);
//             attrib[5] = 0;
//         }
//     }

//     delete[] face_normal;
//     delete[] face_corners;

//     // 5. Global Coloring (Flooding Disjoint Subgraphs)
//     size_t region_id(0);
//     while (true)
//     {
//         bool done = true;
//         // Find a seed not yet colored
//         for (size_t iseed = 0; iseed < num_seeds; iseed++)
//         {
//             if (!seeds->tree_point_is_active(iseed)) continue;
//             size_t* attrib = seeds->get_tree_point_attrib(iseed);
//             if (attrib[5] != 0) continue;
//             region_id++; attrib[5] = region_id;
//             done = false;
//             break;
//         }
//         if (done) break;

//         // Flood fill from that seed
//         while (true)
//         {
//             bool subregion_done(true);
//             for (size_t iseed = 0; iseed < num_seeds; iseed++)
//             {
//                 if (!seeds->tree_point_is_active(iseed)) continue;
//                 size_t* attrib = seeds->get_tree_point_attrib(iseed);
//                 if (attrib[5] != region_id) continue;

//                 std::vector<size_t> neighbor_seeds;
//                 seeds->graph_get_neighbors(iseed, neighbor_seeds);
//                 if (neighbor_seeds.empty()) continue;

//                 for (size_t i = 0; i < neighbor_seeds.size(); i++)
//                 {
//                     size_t neighbor_seed = neighbor_seeds[i];
//                     if (!seeds->tree_point_is_active(neighbor_seed)) continue;

//                     size_t* neighbor_attrib = seeds->get_tree_point_attrib(neighbor_seed);

//                     if (neighbor_attrib[5] != 0)
//                     {
//                         if (neighbor_attrib[5] != region_id)
//                         {
//                             // This indicates a topology error or inconsistency in local orientation
//                             // But usually safe to ignore in simple reconstruction
//                         }
//                         continue;
//                     }
//                     neighbor_attrib[5] = region_id;
//                     subregion_done = false;
//                 }
//             }
//             if (subregion_done) break;
//         }
//     }

//     size_t num_subregions = region_id;
//     std::cout << "  * Number of subregions detected: " << num_subregions << std::endl;

//     // 6. Export Data
//     size_t num_active_seeds = 0;
//     for (size_t iseed = 0; iseed < num_seeds; iseed++)
//     {
//         if (seeds->tree_point_is_active(iseed)) {
//             num_active_seeds++;
//         }
//     }

//     std::cout << "  * Compacting seeds: Total " << num_seeds << " -> Active " << num_active_seeds << std::endl;

//     // [2] 按照活跃数量分配内存
//     // 注意：如果 num_active_seeds 为 0，这里可能需要特殊处理，防止 new 0
//     if (num_active_seeds > 0) {
//         seedes = new double[num_active_seeds * 3];
//         seeds_region_id = new size_t[num_active_seeds];
//         seeds_sizing = new double[num_active_seeds];

//         // [3] 数据压实：只拷贝活跃点
//         size_t current_idx = 0;
//         for (size_t iseed = 0; iseed < num_seeds; iseed++)
//         {
//             if (!seeds->tree_point_is_active(iseed)) continue; // 跳过已删除的点

//             double x[4];
//             seeds->get_tree_point(iseed, 4, x);
//             size_t* attrib = seeds->get_tree_point_attrib(iseed);

//             // 调试输出（保留你原来的逻辑）
//            // if(!attrib[5]) std::cout << "[Warning] Seed " << iseed << " has region ID 0!" << std::endl;

//             // 填充数据到紧凑的数组中
//             for (size_t idim = 0; idim < 3; idim++) {
//                 seedes[current_idx * 3 + idim] = x[idim];
//             }
            
//             seeds_region_id[current_idx] = attrib[5]; 
//             seeds_sizing[current_idx] = x[3];

//             current_idx++; // 只有写入了数据才移动指针
//         }

//         // [4] 生成 CSV，传入真实的活跃点数量
//         generate_seed_csv("seeds.csv", 3, num_active_seeds, seedes, seeds_sizing, seeds_region_id);
//     }
//     else {
//         // 防止空指针崩溃
//         seedes = nullptr;
//         seeds_region_id = nullptr;
//         seeds_sizing = nullptr;
//         std::cout << "[Warning] No active seeds to export." << std::endl;
//     }

//     generate_seed_csv("seeds.csv", 3, num_seeds, seedes, seeds_sizing, seeds_region_id);
// }
	

void Generator::generate_seed_csv(const char* filename, int num_dim, size_t num_seeds, double* spheres, double* spheres_sizing, size_t* spheres_region_id)
{
	#pragma region Save tree to CSV file:
	std::fstream file(filename, std::ios::out);
	// Spheres
	file << "x1coord";
	for (size_t idim = 1; idim < num_dim; idim++) file << ", x" << idim + 1 << "coord";
	file << ", radius" << std::endl;

	for (size_t i = 0; i < num_seeds; i++)
	{
		if(spheres_region_id[i]==0)continue;
		if (spheres_sizing != 0)
		{
			//if (fabs(spheres[i * num_dim] - 0.5) > spheres_sizing[i]) continue;
			file << std::setprecision(16) << spheres[i * num_dim];
			for (size_t idim = 1; idim < num_dim; idim++) file << ", " << spheres[i * num_dim + idim];
			file << ", " << spheres_sizing[i];
			if (spheres_region_id != 0) file << ", " << spheres_region_id[i];
		}
		else
		{
			file << std::setprecision(16) << spheres[i * (num_dim + 1)];
			for (size_t idim = 1; idim <= num_dim; idim++) file << ", " << spheres[i * (num_dim + 1) + idim];
			if (spheres_region_id != 0) file << ", " << spheres_region_id[i];
		}
		file << std::endl;
	}
	#pragma endregion
}
