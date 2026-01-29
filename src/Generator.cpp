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

void Generator::color_surface_seeds(int num_faces, MeshingTree *surface_spheres,MeshingTree *upper_seeds, MeshingTree *lower_seeds,MeshingTree *seeds,
        std::vector<int> face, double* &seedes, size_t* &seeds_region_id, double* &seeds_sizing)
{

	size_t actual_num_faces = face.size() / 3;
    std::cout << "[Debug] In color_surface_seeds:" << std::endl;
    std::cout << "  - Face vector size: " << face.size() << " (implies " << actual_num_faces << " faces)" << std::endl;
    
    size_t num_upper_seeds(upper_seeds->get_num_tree_points());
    size_t num_lower_seeds(lower_seeds->get_num_tree_points());
    std::cout << "  - Upper seeds: " << num_upper_seeds << ", Lower seeds: " << num_lower_seeds << std::endl;

    // adjust id of seed pair for upper seeds
    for (size_t iseed = 0; iseed < num_upper_seeds; iseed++)
    {
      // copy point location, normal and attributes
      double x[4];
      upper_seeds->get_tree_point(iseed, 4, x);
      double* normal = upper_seeds->get_tree_point_normal(iseed);
      size_t* attrib = upper_seeds->get_tree_point_attrib(iseed);
      attrib[1] += num_upper_seeds;
      seeds->add_tree_point(4, x, normal, attrib);
    }
    upper_seeds->clear_memory();

    for (size_t iseed = 0; iseed < num_lower_seeds; iseed++)
    {
      // copy point location, normal and attributes
      double x[4];
      lower_seeds->get_tree_point(iseed, 4, x);
      double* normal = lower_seeds->get_tree_point_normal(iseed);
      size_t* attrib = lower_seeds->get_tree_point_attrib(iseed);
      seeds->add_tree_point(4, x, normal, attrib);
    }
    lower_seeds->clear_memory();

    size_t num_seeds(seeds->get_num_tree_points());
		if (num_faces > 0 && !face.empty())
		{
			struct FaceKeyHash
			{
				size_t operator()(const std::array<size_t, 3>& a) const noexcept
				{
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
				std::array<size_t, 3> key{a, b, c};
				std::sort(key.begin(), key.end());
				face_set.insert(key);
			}

			int max_face_idx = -1;
        for(int idx : face) {
            if(idx > max_face_idx) max_face_idx = idx;
        }
        std::cout << "[Debug] Max Vertex Index in OBJ: " << max_face_idx << std::endl;
        std::cout << "[Debug] Total Spheres (Seeds' source): " << surface_spheres->get_num_tree_points() << std::endl;
        
        if (max_face_idx >= surface_spheres->get_num_tree_points()) {
            std::cout << "[ERROR] MISMATCH DETECTED: OBJ references vertices up to index " << max_face_idx 
                      << ", but there are only " << surface_spheres->get_num_tree_points() << " spheres!" << std::endl;
        }

			for (size_t iseed = 0; iseed < num_seeds; iseed++)
			{
				if (!seeds->tree_point_is_active(iseed)) continue;
				size_t* attrib = seeds->get_tree_point_attrib(iseed);
				const size_t jseed = attrib[1];
				if (jseed >= num_seeds) continue;
				if (iseed > jseed) continue;
				if (!seeds->tree_point_is_active(jseed)) continue;

				std::array<size_t, 3> key{attrib[2], attrib[3], attrib[4]};
				std::sort(key.begin(), key.end());
				if (face_set.find(key) == face_set.end())
				{
					seeds->lazy_delete_tree_point(iseed);
					seeds->lazy_delete_tree_point(jseed);
				}
			}
		}
		
		
		for (size_t iseed = 0; iseed < num_seeds; iseed++)
		{
			// Establish Sphere -> Seeds directional graph
			if (!seeds->tree_point_is_active(iseed)) continue;
			size_t* attrib = seeds->get_tree_point_attrib(iseed);
			size_t si(attrib[2]), sj(attrib[3]), sk(attrib[4]);
			surface_spheres->graph_connect(si, iseed); surface_spheres->graph_connect(sj, iseed); surface_spheres->graph_connect(sk, iseed);
		}

		double* face_normal = new double[3];
		double** face_corners = new double*[3];
		size_t num_spheres(surface_spheres->get_num_tree_points());
		int cand = -1;
		for (size_t isphere = 0; isphere < num_spheres; isphere++)
		{
			#pragma region Connect Seeds of Surface Spheres:
			std::vector<size_t> sphere_seeds;
			surface_spheres->graph_get_neighbors(isphere, sphere_seeds);
			if (sphere_seeds.empty()) continue;

			size_t num_sphere_seeds = sphere_seeds.size();
			for (size_t i = 0; i < num_sphere_seeds; i++)
			{
				size_t seed_i(sphere_seeds[i]);
				if (!seeds->tree_point_is_active(seed_i)) continue;
				size_t attrib_pair = 0;
				seeds->get_tree_point_attrib(seed_i, 1, attrib_pair);
				if (attrib_pair < num_seeds && seeds->tree_point_is_active(attrib_pair))
				{
					size_t* attrib = seeds->get_tree_point_attrib(seed_i);
					attrib[5] = 0;
				}
			}

			size_t first_neighbor(isphere), jsphere(isphere);
			
			#pragma region first face dictates orientation:
			size_t seed_1 = static_cast<size_t>(-1);
			for (size_t i = 0; i < num_sphere_seeds; i++)
			{
				const size_t candidate = sphere_seeds[i];
				if (seeds->tree_point_is_active(candidate)) { seed_1 = candidate; cand = i; break; }
			}
			if (seed_1 == static_cast<size_t>(-1)) continue;
			size_t* attrib = seeds->get_tree_point_attrib(seed_1);

			size_t seed_2(attrib[1]);
			if (seed_2 >= num_seeds) continue;
			if (!seeds->tree_point_is_active(seed_2)) continue;
			size_t si(attrib[2]), sj(attrib[3]), sk(attrib[4]);

			while (si != isphere)
			{
				size_t tmp = si; si = sj; sj = sk; sk = tmp; // rotate indices
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
				attrib[5] = 1; // inside
				attrib = seeds->get_tree_point_attrib(seed_2);
				attrib[5] = 2; // outside
			}
			else
			{
				attrib[5] = 2; // outside
				attrib = seeds->get_tree_point_attrib(seed_2);
				attrib[5] = 1; // inside
			}
			first_neighbor = sj;
			jsphere = sk;
			#pragma endregion
			
			//std::cout<<seed_1<<" "<<seeds->get_tree_point_attrib(seed_1,6)<<" "<<seed_2<<" "<<seeds->get_tree_point_attrib(seed_2,6)<<std::endl;
			while (true)
			{
				#pragma region Mark loop seeds:
				bool done(true);

				for (size_t i = 0; i < num_sphere_seeds; i++)
				{
					size_t seed_1(sphere_seeds[cand + 1 + i]);
					if (!seeds->tree_point_is_active(seed_1)) continue;
					size_t* attrib = seeds->get_tree_point_attrib(seed_1);

					size_t seed_2(attrib[1]);
					if (seed_2 >= num_seeds) continue;
					if (!seeds->tree_point_is_active(seed_2)) continue;
					size_t si(attrib[2]), sj(attrib[3]), sk(attrib[4]);

					while (si != isphere)
					{
						size_t tmp = si; si = sj; sj = sk; sk = tmp; // rotate indices
					}

					if (attrib[5] != 0) continue; // pair is already marked
					if (sj != jsphere && sk != jsphere) continue; // wrong triangle

					done = false;
					break;
				}
				if (done) break;

				for (size_t i = 0; i < num_sphere_seeds; i++)
				{
					size_t seed_1(sphere_seeds[cand + 1 + i]);
					if (!seeds->tree_point_is_active(seed_1)) continue;
					size_t* attrib = seeds->get_tree_point_attrib(seed_1);

					size_t seed_2(attrib[1]);
					if (seed_2 >= num_seeds || !seeds->tree_point_is_active(seed_2)) continue;
					size_t si(attrib[2]), sj(attrib[3]), sk(attrib[4]);

					while (si != isphere)
					{
						size_t tmp = si; si = sj; sj = sk; sk = tmp; // rotate indices
					}

					if (attrib[5] != 0) continue; // pair is already marked
					if (sj != jsphere && sk != jsphere) continue; // wrong triangle

					if (sj != jsphere)
					{
						size_t tmp = sj; sj = sk; sk = tmp; // flip triangle
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
						attrib[5] = 1; // inside
						attrib = seeds->get_tree_point_attrib(seed_2);
						attrib[5] = 2; // outside
					}
					else
					{
						attrib[5] = 2; // outside
						attrib = seeds->get_tree_point_attrib(seed_2);
						attrib[5] = 1; // inside
					}
					jsphere = sk;
					if (jsphere == first_neighbor) done = true;
					break;
				}
				if (done) break;
				#pragma endregion
			}

			for (size_t i = 0; i < num_sphere_seeds; i++)
			{
				#pragma region connect Sphere seeds:
				size_t seed_i(sphere_seeds[i]);
				if (!seeds->tree_point_is_active(seed_i)) continue;
				size_t* attrib_i = seeds->get_tree_point_attrib(seed_i);
				for (size_t j = i + 1; j < num_sphere_seeds; j++)
				{
					size_t seed_j(sphere_seeds[j]);
					if (!seeds->tree_point_is_active(seed_j)) continue;
					size_t* attrib_j = seeds->get_tree_point_attrib(seed_j);

					if (attrib_i[5] == attrib_j[5])
					{
						seeds->graph_connect_nodes(seed_i, seed_j);
					}
				}
				#pragma endregion
			}

			for (size_t i = 0; i < num_sphere_seeds; i++)
			{
				size_t seed_i(sphere_seeds[i]);
				if (!seeds->tree_point_is_active(seed_i)) continue;
				seeds->set_tree_point_attrib(seed_i, 5, static_cast<size_t>(0));
			}
			#pragma endregion
		}

	delete[] face_normal;
	delete[] face_corners;

	// Coloring Disjoint Subgraphs
	size_t region_id(0);
	while (true)
	{
		bool done = true;
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

		while (true)
		{
			// Flooding Subregion
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
							std::cout << "Error::Invalid coloring!!!" << std::endl;
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

	size_t num_subregions = region_id > 0 ? region_id - 1 : 0;

	// Marking ghost seeds
	if (num_seeds > 0)
	{
		double xmin[3], xmax[3];
		for (size_t idim = 0; idim < 3; idim++) xmin[idim] = DBL_MAX;
		for (size_t idim = 0; idim < 3; idim++) xmax[idim] = -DBL_MAX;

		for (size_t iseed = 0; iseed < num_seeds; iseed++)
		{
			if (!seeds->tree_point_is_active(iseed)) continue;
			double* xseed = seeds->get_tree_point(iseed);
			for (size_t idim = 0; idim < 3; idim++)
			{
				if (xseed[idim] < xmin[idim]) xmin[idim] = xseed[idim];
				if (xseed[idim] > xmax[idim]) xmax[idim] = xseed[idim];
			}
		}
		double xfar[3];
		for (size_t idim = 0; idim < 3; idim++) xfar[idim] = xmax[idim] + 2.0 * (xmax[idim] - xmin[idim]);

		size_t iclosest(0); double hclosest(DBL_MAX);
		seeds->get_closest_tree_point(xfar, iclosest, hclosest);

		size_t* attrib = seeds->get_tree_point_attrib(iclosest);
		size_t ghost_region_id = attrib[5];

		// Mark ghost region as 0, shift other region IDs down
		for (size_t iseed = 0; iseed < num_seeds; iseed++)
		{
			if (!seeds->tree_point_is_active(iseed)) continue;
			size_t* attrib = seeds->get_tree_point_attrib(iseed);
			if (attrib[5] == ghost_region_id) attrib[5] = 0;
		}
		for (size_t iseed = 0; iseed < num_seeds; iseed++)
		{
			if (!seeds->tree_point_is_active(iseed)) continue;
			size_t* attrib = seeds->get_tree_point_attrib(iseed);
			if (attrib[5] > ghost_region_id) attrib[5]--;
		}
	}

	std::cout << "  * Number of subregions is " << num_subregions << std::endl;

	num_seeds = seeds->get_num_tree_points();
	seedes = new double[num_seeds * 3];
	seeds_region_id = new size_t[num_seeds];
	seeds_sizing = new double[num_seeds];

	for (size_t iseed = 0; iseed < num_seeds; iseed++)
	{
		double x[4];
		seeds->get_tree_point(iseed, 4, x);
		size_t* attrib = seeds->get_tree_point_attrib(iseed);
		for (size_t idim = 0; idim < 3; idim++) seedes[iseed * 3 + idim] = x[idim];
		seeds_region_id[iseed] = attrib[5];
		seeds_sizing[iseed] = x[3];
	}

   generate_seed_csv("seeds.csv", 3, num_seeds, seedes, seeds_sizing, seeds_region_id);
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