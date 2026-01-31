#include "Generator.h"
#include <unordered_set>
#include <array>
#include <algorithm>
#include <sstream>
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
  if (std::getline(fin, line)) { // 跳表头（若无表头会回退）
    if (line.find("x") == std::string::npos) {
      fin.clear();
      fin.seekg(0);
    }
  }

  size_t count = 0;
  while (std::getline(fin, line)) {
    if (line.empty())
      continue;
    for (char &c : line)
      if (c == ',')
        c = ' ';
    std::istringstream ss(line);

    double x, y, z, r;
    int fid;
    if (!(ss >> x >> y >> z >> r >> fid))
      continue;
    double s[4] = {x, y, z, r};
    size_t attrib[3] = {2, 0, count + 1}; // {长度, 来源面索引(占位或真实面号)}
    attrib[1] = fid;
    spheres->add_tree_point(4, s, /*normal=*/0, attrib);
    ++count;
  }
}

void Generator::read_obj_faces(const char* filename, std::vector<int>& faces_flat, size_t& num_faces)
{
	std::ifstream fin(filename);
	if (!fin)
	{
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
		ss >> token;
		std::array<int, 3> f{};
		int cnt = 0;
		while (cnt < 3 && (ss >> token))
		{
			size_t pos = token.find('/');
			int vid = 0;
			if (pos == std::string::npos) vid = std::stoi(token);
			else vid = std::stoi(token.substr(0, pos));
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
}

// void Generator::generate_connections(const char* filename, std::vector<std::vector<std::tuple<size_t,size_t,size_t>>> &connections)
// {
//     std::ifstream fin(filename);
//     if (!fin) {
//         std::cerr << "Failed to open OBJ file: " << filename << "\n";
//         return;
//     }

//     std::vector<std::vector<size_t>> faces;
//     std::string line;
//     while (std::getline(fin, line)) {
//         if (line.empty() || line[0] != 'f') continue;
//         std::istringstream ss(line);
//         std::string token;
//         ss >> token; // skip 'f'
//         std::vector<size_t> face;
//         while (ss >> token) {
//             size_t pos = token.find('/');
//             size_t vid = std::stoul(token.substr(0, pos));
//             face.push_back(vid - 1); // convert to 0-based
//         }
//         if (face.size() == 3) faces.push_back(face);
//     }

//     // Build vertex -> faces map
//     std::vector<std::vector<size_t>> vertex_to_faces;
//     for (const auto& f : faces) {
//         for (size_t v : f) {
//             if (vertex_to_faces.size() <= v) vertex_to_faces.resize(v + 1);
//             vertex_to_faces[v].push_back(&f - &faces[0]);
//         }
//     }

//     // Build face adjacency via shared vertices, store full vertex triples as tuples
//     connections.resize(faces.size());
//     for (size_t i = 0; i < faces.size(); ++i) {
//         std::unordered_set<size_t> neighbor_indices;
//         for (size_t v : faces[i]) {
//             for (size_t adj : vertex_to_faces[v]) {
//                 if (adj != i) neighbor_indices.insert(adj);
//             }
//         }
//         // Convert neighbor indices to full vertex triples as tuples
//         for (size_t adj : neighbor_indices) {
//             connections[i].emplace_back(faces[adj][0], faces[adj][1], faces[adj][2]);
//         }
//     }
// }

void Generator::generate_surface_seeds(size_t num_points, double **points, size_t num_faces, size_t **faces,
         MeshingTree *spheres,MeshingTree *upper_seeds, MeshingTree *lower_seeds)
{
    // Initialize sliver points OBJ file
    std::ofstream sliver_file("sliver_points.obj", std::ios::out);
    if (sliver_file.is_open()) {
        sliver_file << "# Sliver points detected by VoroCrust" << std::endl;
        sliver_file << "# Format: v x y z r g b (r=1.0 for upper seeds, b=1.0 for lower seeds)" << std::endl;
        sliver_file.close();
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
				delete[] covering_spheres1;delete[] covering_spheres2;
                if(upper_covered && lower_covered)continue;
                if(upper_covered != lower_covered)
                {
                    std::cout << "Sliver! caused by spheres " << sphere_index_i << " " << sphere_index_j << " " << sphere_index_k << std::endl;
                    for(int i=0;i<3;i++)
                    {
                        std::cout<<upper_seed[i]<<" ";
                    }
                    std::cout<<std::endl;
                    for(int i=0;i<3;i++)
                    {
                        std::cout<<lower_seed[i]<<" ";
                    }
                    std::cout<<std::endl;
                    
                    // Get sphere centers for the three spheres causing the sliver
                    double sphere_i_center[4], sphere_j_center[4], sphere_k_center[4];
                    spheres->get_tree_point(sphere_index_i, 4, sphere_i_center);
                    spheres->get_tree_point(sphere_index_j, 4, sphere_j_center);
                    spheres->get_tree_point(sphere_index_k, 4, sphere_k_center);
                    
                    // Output sliver points and sphere centers to OBJ file
                    std::ofstream sliver_file("sliver_points.obj", std::ios::app);
                    if (sliver_file.is_open()) {
                        // Write upper seed point (red)
                        sliver_file << "v " << upper_seed[0] << " " << upper_seed[1] << " " << upper_seed[2] << " 1.0 0.0 0.0" << std::endl;
                        // Write lower seed point (blue)
                        sliver_file << "v " << lower_seed[0] << " " << lower_seed[1] << " " << lower_seed[2] << " 0.0 0.0 1.0" << std::endl;
                        // Write sphere centers (green)
                        sliver_file << "v " << sphere_i_center[0] << " " << sphere_i_center[1] << " " << sphere_i_center[2] << " 0.0 1.0 0.0" << std::endl;
                        sliver_file << "v " << sphere_j_center[0] << " " << sphere_j_center[1] << " " << sphere_j_center[2] << " 0.0 1.0 0.0" << std::endl;
                        sliver_file << "v " << sphere_k_center[0] << " " << sphere_k_center[1] << " " << sphere_k_center[2] << " 0.0 1.0 0.0" << std::endl;
                        
                        // Write triangle formed by the three sphere centers (yellow)
                        // Get current vertex count for face indexing
                        static int vertex_count = 0;
                        if (vertex_count == 0) {
                            // Count existing vertices in file
                            std::ifstream count_file("sliver_points.obj");
                            std::string line;
                            while (std::getline(count_file, line)) {
                                if (line.substr(0, 2) == "v ") {
                                    vertex_count++;
                                }
                            }
                            count_file.close();
                        }
                        
                        // Write triangle face using the last 3 vertices (sphere centers)
                        sliver_file << "f " << (vertex_count - 2) << " " << (vertex_count - 1) << " " << vertex_count << std::endl;
                        vertex_count += 5; // 2 seed points + 3 sphere centers
                        
                        sliver_file.close();
                    }
                }
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
}

void Generator::color_surface_seeds(int num_faces, MeshingTree *surface_spheres, MeshingTree *upper_seeds, MeshingTree *lower_seeds, MeshingTree *seeds,
    std::vector<int> face, double** points, double* &seedes, size_t* &seeds_region_id, double* &seeds_sizing)
{
    std::cout << "[Info] Starting robust color_surface_seeds..." << std::endl;

    // 1. Merge Upper and Lower seeds into 'seeds' tree
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
            if(key[0]==723 && key[1]==778 && key[2]==1989){
                std::cout<<"Found sliver!"<<std::endl;
            }
            if(face_set.find(key) == face_set.end())
            {
                std::cout<<key[0]<<" "<<key[1]<<" "<<key[2]<<std::endl;
                seeds->lazy_delete_tree_point(iseed);
                seeds->lazy_delete_tree_point(jseed);
            }
        }
    }

    // 3. Build connectivity graph
    for (size_t iseed = 0; iseed < num_seeds; iseed++)
    {
        if (!seeds->tree_point_is_active(iseed)) continue;
        size_t* attrib = seeds->get_tree_point_attrib(iseed);
        size_t si(attrib[2]), sj(attrib[3]), sk(attrib[4]);
        surface_spheres->graph_connect(si, iseed);
        surface_spheres->graph_connect(sj, iseed);
        surface_spheres->graph_connect(sk, iseed);
    }

    double* face_normal = new double[3];
    double** face_corners = new double*[3];
    size_t num_spheres(surface_spheres->get_num_tree_points());

    // 4. Propagate Orientation (The Vorocrust Method)
    for (size_t isphere = 0; isphere < num_spheres; isphere++)
    {
        std::vector<size_t> sphere_seeds;
        surface_spheres->graph_get_neighbors(isphere, sphere_seeds);
        if (sphere_seeds.empty()) continue;

        size_t num_sphere_seeds = sphere_seeds.size();

        // [Step A] Reset markers for this sphere
        for (size_t i = 0; i < num_sphere_seeds; i++)
        {
            size_t seed_i(sphere_seeds[i]);
            size_t* attrib = seeds->get_tree_point_attrib(seed_i);
            attrib[5] = 0; // 0 means unvisited/undefined for this local loop
        }

        size_t first_neighbor(isphere), jsphere(isphere);

        // [Step B] Pick the first seed to dictate orientation
        if (true)
        {
            size_t seed_1(sphere_seeds[0]); // Pick the first valid one
            size_t* attrib = seeds->get_tree_point_attrib(seed_1);

            size_t seed_2(attrib[1]);
            size_t si(attrib[2]), sj(attrib[3]), sk(attrib[4]);

            // [CRITICAL] Rotate indices so 'si' is the current sphere
            while (si != isphere)
            {
                size_t tmp = si; si = sj; sj = sk; sk = tmp;
            }

            double* x_1 = seeds->get_tree_point(seed_1);
            double* x_2 = seeds->get_tree_point(seed_2);

            face_corners[0] = surface_spheres->get_tree_point(si);
            face_corners[1] = surface_spheres->get_tree_point(sj);
            face_corners[2] = surface_spheres->get_tree_point(sk);

            _geom.get_3d_triangle_normal(face_corners, face_normal);

            double dot(0.0);
            for (size_t idim = 0; idim < 3; idim++) dot += (x_2[idim] - x_1[idim]) * face_normal[idim];

            if (dot > 0.0)
            {
                attrib[5] = 1; // 1: Inside
                attrib = seeds->get_tree_point_attrib(seed_2);
                attrib[5] = 2; // 2: Outside
            }
            else
            {
                attrib[5] = 2; // Outside
                attrib = seeds->get_tree_point_attrib(seed_2);
                attrib[5] = 1; // Inside
            }
            first_neighbor = sj; // Record where we started
            jsphere = sk;        // The vertex we are moving towards
        }

        // [Step C] Walk the loop (Topological Propagation)
        while (true)
        {
            bool done(true);

            // Look for the next neighbor in the chain
            for (size_t i = 0; i < num_sphere_seeds; i++)
            {
                size_t seed_1(sphere_seeds[i]);
                size_t* attrib = seeds->get_tree_point_attrib(seed_1);

                if (attrib[5] != 0) continue; // Already processed

                size_t seed_2(attrib[1]);
                size_t si(attrib[2]), sj(attrib[3]), sk(attrib[4]);

                // [CRITICAL] Rotate indices again
                while (si != isphere)
                {
                    size_t tmp = si; si = sj; sj = sk; sk = tmp;
                }

                // We are looking for a triangle that shares the edge (si, jsphere)
                // In a valid traversal, one of sj or sk MUST be jsphere.
                if (sj != jsphere && sk != jsphere) continue; // Not the neighbor we want

                // [CRITICAL FIX] Enforce winding order
                // If sj is NOT the shared vertex, it means the triangle is flipped relative to our traversal.
                // We MUST swap sj and sk to maintain consistent normal direction.
                if (sj != jsphere)
                {
                    size_t tmp = sj; sj = sk; sk = tmp;
                }

                // Now calculate normal with the CONSISTENT vertex order (si, sj, sk)
                double* x_1 = seeds->get_tree_point(seed_1);
                double* x_2 = seeds->get_tree_point(seed_2);

                face_corners[0] = surface_spheres->get_tree_point(si);
                face_corners[1] = surface_spheres->get_tree_point(sj);
                face_corners[2] = surface_spheres->get_tree_point(sk);

                _geom.get_3d_triangle_normal(face_corners, face_normal);

                double dot(0.0);
                for (size_t idim = 0; idim < 3; idim++) dot += (x_2[idim] - x_1[idim]) * face_normal[idim];

                if (dot > 0.0)
                {
                    attrib[5] = 1; // Inside
                    attrib = seeds->get_tree_point_attrib(seed_2);
                    attrib[5] = 2; // Outside
                }
                else
                {
                    attrib[5] = 2; // Outside
                    attrib = seeds->get_tree_point_attrib(seed_2);
                    attrib[5] = 1; // Inside
                }

                jsphere = sk; // Advance to the next vertex
                if (jsphere == first_neighbor) done = true; // Loop closed
                
                done = false; // Found a valid step, continue the while loop
                break; 
            }
            if (done) break; // No more neighbors found or loop closed
        }

        // [Step D] Connect compatible seeds
        // Now that seeds on this sphere are locally consistent, connect them in the graph
        for (size_t i = 0; i < num_sphere_seeds; i++)
        {
            size_t seed_i(sphere_seeds[i]);
            size_t* attrib_i = seeds->get_tree_point_attrib(seed_i);
            
            // Only connect processed seeds
            if(attrib_i[5] == 0) continue;

            for (size_t j = i + 1; j < num_sphere_seeds; j++)
            {
                size_t seed_j(sphere_seeds[j]);
                size_t* attrib_j = seeds->get_tree_point_attrib(seed_j);

                if (attrib_i[5] == attrib_j[5])
                {
                    seeds->graph_connect_nodes(seed_i, seed_j);
                }
            }
        }

        // [Step E] Reset for global coloring
        // We wipe the local 1/2 markers because we will re-flood globally
        for (size_t i = 0; i < num_sphere_seeds; i++)
        {
            size_t seed_i(sphere_seeds[i]);
            size_t* attrib = seeds->get_tree_point_attrib(seed_i);
            attrib[5] = 0;
        }
    }

    delete[] face_normal;
    delete[] face_corners;

    // 5. Global Coloring (Flooding Disjoint Subgraphs)
    size_t region_id(0);
    while (true)
    {
        bool done = true;
        // Find a seed not yet colored
        for (size_t iseed = 0; iseed < num_seeds; iseed++)
        {
            if (!seeds->tree_point_is_active(iseed)) continue;
            size_t* attrib = seeds->get_tree_point_attrib(iseed);
            if (attrib[5] != 0) continue;
            region_id++; attrib[5] = region_id;
            done = false;
            break;
        }
        if (done) break;

        // Flood fill from that seed
        while (true)
        {
            bool subregion_done(true);
            for (size_t iseed = 0; iseed < num_seeds; iseed++)
            {
                if (!seeds->tree_point_is_active(iseed)) continue;
                size_t* attrib = seeds->get_tree_point_attrib(iseed);
                if (attrib[5] != region_id) continue;

                std::vector<size_t> neighbor_seeds;
                seeds->graph_get_neighbors(iseed, neighbor_seeds);
                if (neighbor_seeds.empty()) continue;

                for (size_t i = 0; i < neighbor_seeds.size(); i++)
                {
                    size_t neighbor_seed = neighbor_seeds[i];
                    if (!seeds->tree_point_is_active(neighbor_seed)) continue;

                    size_t* neighbor_attrib = seeds->get_tree_point_attrib(neighbor_seed);

                    if (neighbor_attrib[5] != 0)
                    {
                        if (neighbor_attrib[5] != region_id)
                        {
                            // This indicates a topology error or inconsistency in local orientation
                            // But usually safe to ignore in simple reconstruction
                        }
                        continue;
                    }
                    neighbor_attrib[5] = region_id;
                    subregion_done = false;
                }
            }
            if (subregion_done) break;
        }
    }

    size_t num_subregions = region_id;
    std::cout << "  * Number of subregions detected: " << num_subregions << std::endl;

    // 6. Export Data
    size_t num_active_seeds = 0;
    for (size_t iseed = 0; iseed < num_seeds; iseed++)
    {
        if (seeds->tree_point_is_active(iseed)) {
            num_active_seeds++;
        }
    }

    std::cout << "  * Compacting seeds: Total " << num_seeds << " -> Active " << num_active_seeds << std::endl;

    // [2] 按照活跃数量分配内存
    // 注意：如果 num_active_seeds 为 0，这里可能需要特殊处理，防止 new 0
    if (num_active_seeds > 0) {
        seedes = new double[num_active_seeds * 3];
        seeds_region_id = new size_t[num_active_seeds];
        seeds_sizing = new double[num_active_seeds];

        // [3] 数据压实：只拷贝活跃点
        size_t current_idx = 0;
        for (size_t iseed = 0; iseed < num_seeds; iseed++)
        {
            if (!seeds->tree_point_is_active(iseed)) continue; // 跳过已删除的点

            double x[4];
            seeds->get_tree_point(iseed, 4, x);
            size_t* attrib = seeds->get_tree_point_attrib(iseed);

            // 调试输出（保留你原来的逻辑）
            if(!attrib[5]) std::cout << "[Warning] Seed " << iseed << " has region ID 0!" << std::endl;

            // 填充数据到紧凑的数组中
            for (size_t idim = 0; idim < 3; idim++) {
                seedes[current_idx * 3 + idim] = x[idim];
            }
            
            seeds_region_id[current_idx] = attrib[5]; 
            seeds_sizing[current_idx] = x[3];

            current_idx++; // 只有写入了数据才移动指针
        }

        // [4] 生成 CSV，传入真实的活跃点数量
        generate_seed_csv("seeds.csv", 3, num_active_seeds, seedes, seeds_sizing, seeds_region_id);
    }
    else {
        // 防止空指针崩溃
        seedes = nullptr;
        seeds_region_id = nullptr;
        seeds_sizing = nullptr;
        std::cout << "[Warning] No active seeds to export." << std::endl;
    }

    // Rebuild seeds tree with only active points so downstream (e.g. generate_surface_mesh)
    // sees the filtered/deleted set and pairing indices remain consistent.
    if (num_active_seeds > 0 && num_active_seeds < num_seeds)
    {
        std::vector<size_t> old_to_new(num_seeds, static_cast<size_t>(-1));
        size_t new_idx = 0;
        for (size_t iseed = 0; iseed < num_seeds; iseed++)
        {
            if (!seeds->tree_point_is_active(iseed)) continue;
            old_to_new[iseed] = new_idx;
            new_idx++;
        }

        struct SeedSnapshot
        {
            double x[4];
            double normal[4];
            std::vector<size_t> attrib;
        };
        std::vector<SeedSnapshot> snapshots;
        snapshots.reserve(num_active_seeds);

        for (size_t iseed = 0; iseed < num_seeds; iseed++)
        {
            if (!seeds->tree_point_is_active(iseed)) continue;

            SeedSnapshot s;
            seeds->get_tree_point(iseed, 4, s.x);
            double* n = seeds->get_tree_point_normal(iseed);
            for (size_t idim = 0; idim < 4; idim++) s.normal[idim] = n[idim];

            size_t* attrib = seeds->get_tree_point_attrib(iseed);
            const size_t num_attrib = attrib ? attrib[0] : 0;
            if (num_attrib > 0)
            {
                s.attrib.assign(attrib, attrib + num_attrib);
                const size_t old_pair = (s.attrib.size() > 1) ? s.attrib[1] : static_cast<size_t>(-1);
                if (old_pair < num_seeds && old_to_new[old_pair] != static_cast<size_t>(-1))
                {
                    s.attrib[1] = old_to_new[old_pair];
                }
                else if (s.attrib.size() > 1)
                {
                    s.attrib[1] = old_to_new[iseed];
                }
            }
            snapshots.push_back(std::move(s));
        }

        seeds->clear_memory();
        for (size_t i = 0; i < snapshots.size(); i++)
        {
            size_t* attrib_ptr = snapshots[i].attrib.empty() ? nullptr : snapshots[i].attrib.data();
            seeds->add_tree_point(4, snapshots[i].x, snapshots[i].normal, attrib_ptr);
        }
    }
}
	

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