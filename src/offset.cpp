#include "Offset.h"
#include "Geometry.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#ifdef USE_CGAL
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Regular_triangulation_cell_base_3.h>
#endif

namespace {

struct SeedRecord {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    double r = 0.0;
    size_t region = 0; // 0 outside, 1 inside
    size_t seed_id = static_cast<size_t>(-1);
    size_t pair_id = static_cast<size_t>(-1);
    bool has_face_indices = false;
    size_t face_idx[3] = {static_cast<size_t>(-1), static_cast<size_t>(-1), static_cast<size_t>(-1)};
    bool has_face_vertices = false;
    double face_v[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
};

struct WeightedSite {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    double w = 0.0;
    size_t pair_id = static_cast<size_t>(-1);
    size_t seed_id = static_cast<size_t>(-1);
    size_t owner_seed_id = static_cast<size_t>(-1); // basepoint/offset pair key
    int type = 0;   // 0 = seed, 1 = basepoint
    int region = 0; // 0 outside, 1 inside
};

struct BaseOffsetPairRecord {
    size_t pair_id = static_cast<size_t>(-1);
    size_t owner_seed_id = static_cast<size_t>(-1);
    int region = 0; // 0 outside, 1 inside
    double base[3] = {0.0, 0.0, 0.0};
    double offset[3] = {0.0, 0.0, 0.0};
    double normal[3] = {0.0, 0.0, 0.0}; // unit direction from base -> offset
};

#ifdef USE_CGAL
using K = CGAL::Exact_predicates_exact_constructions_kernel;
using Point_3 = K::Point_3;
using Weighted_point = K::Weighted_point_3;
using RVb = CGAL::Regular_triangulation_vertex_base_3<K>;
using Vb = CGAL::Triangulation_vertex_base_with_info_3<size_t, K, RVb>;
using Cb = CGAL::Regular_triangulation_cell_base_3<K>;
using Tds = CGAL::Triangulation_data_structure_3<Vb, Cb>;
using Regular = CGAL::Regular_triangulation_3<K, Tds>;
#endif

static inline void trim_in_place(std::string& s)
{
    size_t start = 0;
    while (start < s.size() && std::isspace(static_cast<unsigned char>(s[start]))) {
        start++;
    }
    size_t end = s.size();
    while (end > start && std::isspace(static_cast<unsigned char>(s[end - 1]))) {
        end--;
    }
    if (start == 0 && end == s.size()) return;
    s = s.substr(start, end - start);
}

static bool parse_seed_csv_line(const std::string& line, SeedRecord& out)
{
    std::string s = line;
    trim_in_place(s);
    if (s.empty()) return false;
    if (s[0] == '#') return false;

    std::vector<std::string> tokens;
    std::stringstream ss(s);
    std::string token;
    while (std::getline(ss, token, ',')) {
        trim_in_place(token);
        if (!token.empty()) tokens.push_back(token);
    }

    if (tokens.size() < 4) return false;
    if (!tokens[0].empty() && std::isalpha(static_cast<unsigned char>(tokens[0][0]))) {
        return false;
    }

    try {
        out.x = std::stod(tokens[0]);
        out.y = std::stod(tokens[1]);
        out.z = std::stod(tokens[2]);
        out.r = std::stod(tokens[3]);
        if (tokens.size() >= 5) out.region = static_cast<size_t>(std::stoull(tokens[4]));
        if (tokens.size() >= 6) out.seed_id = static_cast<size_t>(std::stoull(tokens[5]));
        if (tokens.size() >= 7) out.pair_id = static_cast<size_t>(std::stoull(tokens[6]));
        if (tokens.size() >= 10) {
            out.has_face_indices = true;
            out.face_idx[0] = static_cast<size_t>(std::stoull(tokens[7]));
            out.face_idx[1] = static_cast<size_t>(std::stoull(tokens[8]));
            out.face_idx[2] = static_cast<size_t>(std::stoull(tokens[9]));
        }
        if (tokens.size() >= 19) {
            out.has_face_vertices = true;
            out.face_v[0][0] = std::stod(tokens[10]);
            out.face_v[0][1] = std::stod(tokens[11]);
            out.face_v[0][2] = std::stod(tokens[12]);
            out.face_v[1][0] = std::stod(tokens[13]);
            out.face_v[1][1] = std::stod(tokens[14]);
            out.face_v[1][2] = std::stod(tokens[15]);
            out.face_v[2][0] = std::stod(tokens[16]);
            out.face_v[2][1] = std::stod(tokens[17]);
            out.face_v[2][2] = std::stod(tokens[18]);
        }
    } catch (...) {
        return false;
    }

    return true;
}

static inline double distance3(const double* a, const double* b)
{
    const double dx = a[0] - b[0];
    const double dy = a[1] - b[1];
    const double dz = a[2] - b[2];
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

static void compute_unit_normal_base_to_offset(
    const double base[3],
    const double offset[3],
    double normal_out[3]
)
{
    const double vx = offset[0] - base[0];
    const double vy = offset[1] - base[1];
    const double vz = offset[2] - base[2];
    const double len = std::sqrt(vx * vx + vy * vy + vz * vz);
    if (len < 1e-20) {
        normal_out[0] = 0.0;
        normal_out[1] = 0.0;
        normal_out[2] = 0.0;
        return;
    }
    normal_out[0] = vx / len;
    normal_out[1] = vy / len;
    normal_out[2] = vz / len;
}

static bool compute_basepoint(
    const SeedRecord& s0,
    const SeedRecord& s1,
    double base_out[3],
    double& d0,
    double& d1
)
{
    const double p0[3] = {s0.x, s0.y, s0.z};
    const double p1[3] = {s1.x, s1.y, s1.z};

    if (!s0.has_face_vertices && !s1.has_face_vertices) {
        base_out[0] = 0.5 * (p0[0] + p1[0]);
        base_out[1] = 0.5 * (p0[1] + p1[1]);
        base_out[2] = 0.5 * (p0[2] + p1[2]);
        d0 = distance3(p0, base_out);
        d1 = distance3(p1, base_out);
        return true;
    }

    const double (*face)[3] = s0.has_face_vertices ? s0.face_v : s1.face_v;
    const double a[3] = {face[0][0], face[0][1], face[0][2]};
    const double b[3] = {face[1][0], face[1][1], face[1][2]};
    const double c[3] = {face[2][0], face[2][1], face[2][2]};

    double ab[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]};
    double ac[3] = {c[0] - a[0], c[1] - a[1], c[2] - a[2]};
    double n[3] = {
        ab[1] * ac[2] - ab[2] * ac[1],
        ab[2] * ac[0] - ab[0] * ac[2],
        ab[0] * ac[1] - ab[1] * ac[0]
    };
    const double nlen2 = n[0]*n[0] + n[1]*n[1] + n[2]*n[2];
    if (nlen2 < 1e-20) {
        base_out[0] = 0.5 * (p0[0] + p1[0]);
        base_out[1] = 0.5 * (p0[1] + p1[1]);
        base_out[2] = 0.5 * (p0[2] + p1[2]);
        d0 = distance3(p0, base_out);
        d1 = distance3(p1, base_out);
        return true;
    }

    const double dir[3] = {p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]};
    const double denom = n[0]*dir[0] + n[1]*dir[1] + n[2]*dir[2];
    if (std::abs(denom) < 1e-20) {
        base_out[0] = 0.5 * (p0[0] + p1[0]);
        base_out[1] = 0.5 * (p0[1] + p1[1]);
        base_out[2] = 0.5 * (p0[2] + p1[2]);
        d0 = distance3(p0, base_out);
        d1 = distance3(p1, base_out);
        return true;
    }

    const double t = (n[0] * (a[0] - p0[0]) + n[1] * (a[1] - p0[1]) + n[2] * (a[2] - p0[2])) / denom;
    double base[3] = {p0[0] + t * dir[0], p0[1] + t * dir[1], p0[2] + t * dir[2]};

    Geometry geom;
    double* corners[3];
    double ca[3] = {a[0], a[1], a[2]};
    double cb[3] = {b[0], b[1], b[2]};
    double cc[3] = {c[0], c[1], c[2]};
    corners[0] = ca; corners[1] = cb; corners[2] = cc;

    double q[3];
    double proj_dist = 0.0;
    geom.project_to_3d_triangle(base, corners, q, proj_dist);

    base_out[0] = q[0];
    base_out[1] = q[1];
    base_out[2] = q[2];

    d0 = distance3(p0, base_out);
    d1 = distance3(p1, base_out);
    return true;
}

static double compute_R_value(const Offset::PowerDiagramParams& params, const double base[3], int region)
{
    (void)base;
    double r = params.r_value;
    if (params.r_signed_by_region && region == 1) r = -r;
    return r;
}

static bool export_power_diagram_obj(
    const std::vector<WeightedSite>& sites,
    const std::string& out_obj,
    double merge_epsilon,
    std::vector<std::vector<Point_3>>* facets_out,
    std::string* error
)
{
    if (sites.empty()) {
        if (error) *error = "No weighted sites to export.";
        return false;
    }

    Regular rt;
    for (size_t i = 0; i < sites.size(); ++i) {
        const auto& s = sites[i];
        Weighted_point wp(Point_3(s.x, s.y, s.z), s.w);
        Regular::Vertex_handle vh = rt.insert(wp);
        if (vh != Regular::Vertex_handle()) {
            vh->info() = i;
        }
    }

    double min_x = std::numeric_limits<double>::infinity();
    double min_y = std::numeric_limits<double>::infinity();
    double min_z = std::numeric_limits<double>::infinity();
    double max_x = -std::numeric_limits<double>::infinity();
    double max_y = -std::numeric_limits<double>::infinity();
    double max_z = -std::numeric_limits<double>::infinity();
    for (const auto& s : sites) {
        min_x = std::min(min_x, s.x);
        min_y = std::min(min_y, s.y);
        min_z = std::min(min_z, s.z);
        max_x = std::max(max_x, s.x);
        max_y = std::max(max_y, s.y);
        max_z = std::max(max_z, s.z);
    }
    const double global_dx = max_x - min_x;
    const double global_dy = max_y - min_y;
    const double global_dz = max_z - min_z;
    const double global_diag = std::sqrt(global_dx * global_dx + global_dy * global_dy + global_dz * global_dz);

    std::vector<std::vector<Point_3>> facets;
    size_t far_vertices_filtered = 0;
    for (auto eit = rt.finite_edges_begin(); eit != rt.finite_edges_end(); ++eit) {
        Regular::Cell_handle c = eit->first;
        int i1 = eit->second;
        int i2 = eit->third;
        Regular::Vertex_handle v1 = c->vertex(i1);
        Regular::Vertex_handle v2 = c->vertex(i2);
        if (rt.is_infinite(v1) || rt.is_infinite(v2)) continue;

        size_t idx1 = v1->info();
        size_t idx2 = v2->info();
        if (idx1 >= sites.size() || idx2 >= sites.size()) continue;
        const auto& s1 = sites[idx1];
        const auto& s2 = sites[idx2];
        const double mx = 0.5 * (s1.x + s2.x);
        const double my = 0.5 * (s1.y + s2.y);
        const double mz = 0.5 * (s1.z + s2.z);
        const double pdx = s1.x - s2.x;
        const double pdy = s1.y - s2.y;
        const double pdz = s1.z - s2.z;
        const double pair_dist = std::sqrt(pdx * pdx + pdy * pdy + pdz * pdz);
        const double far_limit = 8.0 * pair_dist + 0.25 * global_diag;

        std::vector<Point_3> facet_vertices;
        Regular::Cell_circulator cc = rt.incident_cells(*eit);
        Regular::Cell_circulator done = cc;
        if (cc == nullptr) continue;
        do {
            if (!rt.is_infinite(cc)) {
                Point_3 center = rt.dual(cc);
                const double cx = CGAL::to_double(center.x());
                const double cy = CGAL::to_double(center.y());
                const double cz = CGAL::to_double(center.z());
                const double dmx = cx - mx;
                const double dmy = cy - my;
                const double dmz = cz - mz;
                const double dist_mid = std::sqrt(dmx * dmx + dmy * dmy + dmz * dmz);
                if (dist_mid <= far_limit) {
                    facet_vertices.push_back(center);
                } else {
                    far_vertices_filtered++;
                }
            }
            ++cc;
        } while (cc != done);

        if (facet_vertices.size() >= 3) {
            facets.push_back(std::move(facet_vertices));
        }
    }
    std::cout << "  * Power facets (" << out_obj << ") filtered far dual vertices: "
              << far_vertices_filtered << std::endl;

    if (facets_out != nullptr) {
        *facets_out = facets;
    }

    std::ofstream out(out_obj);
    if (!out.is_open()) {
        if (error) *error = "Cannot open output OBJ: " + out_obj;
        return false;
    }

    out << "# Power diagram facets\n";
    out << "# sites=" << sites.size() << " facets=" << facets.size() << "\n";

    std::map<std::tuple<double, double, double>, size_t> vertex_map;
    std::vector<Point_3> unique_vertices;

    const double eps = (merge_epsilon > 0.0) ? merge_epsilon : 1e-10;
    auto round_key = [&](const Point_3& p) -> std::tuple<double, double, double> {
        const double x = std::round(CGAL::to_double(p.x()) / eps) * eps;
        const double y = std::round(CGAL::to_double(p.y()) / eps) * eps;
        const double z = std::round(CGAL::to_double(p.z()) / eps) * eps;
        return std::make_tuple(x, y, z);
    };

    auto get_vertex_index = [&](const Point_3& p) -> size_t {
        auto key = round_key(p);
        auto it = vertex_map.find(key);
        if (it != vertex_map.end()) return it->second;
        size_t idx = unique_vertices.size();
        vertex_map.emplace(key, idx);
        unique_vertices.push_back(p);
        return idx;
    };

    std::vector<std::vector<size_t>> face_indices;
    face_indices.reserve(facets.size());
    for (const auto& facet : facets) {
        if (facet.size() < 3) continue;
        std::vector<size_t> inds;
        inds.reserve(facet.size());
        for (const auto& p : facet) {
            inds.push_back(get_vertex_index(p));
        }
        face_indices.push_back(std::move(inds));
    }

    out << std::fixed << std::setprecision(16);
    for (const auto& p : unique_vertices) {
        out << "v " << CGAL::to_double(p.x()) << " "
            << CGAL::to_double(p.y()) << " "
            << CGAL::to_double(p.z()) << "\n";
    }

    for (const auto& inds : face_indices) {
        if (inds.size() < 3) continue;
        out << "f";
        for (size_t idx : inds) {
            out << " " << (idx + 1);
        }
        out << "\n";
    }

    return true;
}

static bool write_facets_dedup_obj(
    const std::string& filename,
    const std::vector<std::vector<Point_3>>& facets,
    double epsilon,
    double angle_threshold_deg,
    std::string* error
)
{
    std::ofstream out(filename);
    if (!out.is_open()) {
        if (error) *error = "Cannot open dedup OBJ: " + filename;
        return false;
    }
    out << "# Deduplicated facets (epsilon=" << epsilon
        << ", angle_threshold=" << angle_threshold_deg << "deg)\n";

    struct VertexKey {
        long long x, y, z;
        bool operator==(const VertexKey& other) const {
            return x == other.x && y == other.y && z == other.z;
        }
    };
    struct VertexKeyHash {
        std::size_t operator()(const VertexKey& k) const {
            size_t h1 = std::hash<long long>{}(k.x);
            size_t h2 = std::hash<long long>{}(k.y);
            size_t h3 = std::hash<long long>{}(k.z);
            return h1 ^ (h2 << 1) ^ (h3 << 2);
        }
    };

    const double eps = (epsilon > 0.0) ? epsilon : 1e-10;
    const double scale_factor = 1.0 / eps;
    const double rad_to_deg = 180.0 / 3.14159265358979323846;

    std::vector<Point_3> vertices;
    vertices.reserve(facets.size() * 2);
    std::vector<std::vector<int>> faces;
    faces.reserve(facets.size());
    std::unordered_map<VertexKey, int, VertexKeyHash> vertex_map;
    vertex_map.reserve(facets.size() * 3);

    auto find_or_add_vertex = [&](const Point_3& p) -> int {
        double x = CGAL::to_double(p.x());
        double y = CGAL::to_double(p.y());
        double z = CGAL::to_double(p.z());
        long long ix = std::llround(x * scale_factor);
        long long iy = std::llround(y * scale_factor);
        long long iz = std::llround(z * scale_factor);
        VertexKey key{ix, iy, iz};
        auto it = vertex_map.find(key);
        if (it != vertex_map.end()) return it->second;
        vertices.push_back(p);
        int idx = static_cast<int>(vertices.size()); // 1-based
        vertex_map.emplace(key, idx);
        return idx;
    };

    auto compute_angle = [&](const Point_3& p1, const Point_3& p2, const Point_3& p3) -> double {
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
        if (cos_angle > 1.0) cos_angle = 1.0;
        else if (cos_angle < -1.0) cos_angle = -1.0;
        return std::acos(cos_angle) * rad_to_deg;
    };

    std::vector<Point_3> corner_points;
    corner_points.reserve(16);
    std::vector<int> face_indices;
    face_indices.reserve(16);
    std::vector<int> cleaned_face;
    cleaned_face.reserve(16);

    for (const auto& facet : facets) {
        if (facet.size() < 3) continue;
        corner_points.clear();
        int n = static_cast<int>(facet.size());
        if (n <= 3) {
            corner_points = facet;
        } else {
            for (int i = 0; i < n; ++i) {
                const Point_3& prev = facet[(i - 1 + n) % n];
                const Point_3& curr = facet[i];
                const Point_3& next = facet[(i + 1) % n];
                if (compute_angle(prev, curr, next) < angle_threshold_deg) {
                    corner_points.push_back(curr);
                }
            }
            if (corner_points.size() < 3) corner_points = facet;
        }

        if (corner_points.size() < 3) continue;
        face_indices.clear();
        for (const auto& p : corner_points) {
            face_indices.push_back(find_or_add_vertex(p));
        }

        cleaned_face.clear();
        if (!face_indices.empty()) {
            cleaned_face.push_back(face_indices[0]);
            for (size_t i = 1; i < face_indices.size(); ++i) {
                if (face_indices[i] != face_indices[i - 1]) {
                    cleaned_face.push_back(face_indices[i]);
                }
            }
            if (cleaned_face.size() > 1 && cleaned_face.front() == cleaned_face.back()) {
                cleaned_face.pop_back();
            }
        }

        if (cleaned_face.size() >= 3) {
            bool has_duplicate = false;
            for (size_t i = 0; i < cleaned_face.size() && !has_duplicate; ++i) {
                for (size_t j = i + 1; j < cleaned_face.size() && !has_duplicate; ++j) {
                    if (cleaned_face[i] == cleaned_face[j]) has_duplicate = true;
                }
            }
            if (!has_duplicate) {
                faces.push_back(cleaned_face);
            }
        }
    }

    out << std::fixed << std::setprecision(16);
    for (const auto& p : vertices) {
        out << "v " << CGAL::to_double(p.x()) << " "
            << CGAL::to_double(p.y()) << " "
            << CGAL::to_double(p.z()) << "\n";
    }
    for (const auto& f : faces) {
        out << "f";
        for (int idx : f) out << " " << idx;
        out << "\n";
    }

    return true;
}

static bool write_sites_points_obj(
    const std::string& filename,
    const std::vector<WeightedSite>& sites,
    int target_type,
    std::string* error
)
{
    if (filename.empty()) return true;

    std::ofstream out(filename);
    if (!out.is_open()) {
        if (error) *error = "Cannot open sites OBJ: " + filename;
        return false;
    }

    out << "# Sites points OBJ\n";
    out << "# type=" << target_type << " (0=offset,1=basepoint)\n";
    out << std::fixed << std::setprecision(16);

    size_t count = 0;
    for (const auto& s : sites) {
        if (s.type != target_type) continue;
        out << "v " << s.x << " " << s.y << " " << s.z << "\n";
        count++;
    }
    std::cout << "  * Wrote " << count << " points to " << filename << std::endl;
    return true;
}

static bool write_base_offset_pairs_csv(
    const std::string& filename,
    const std::vector<BaseOffsetPairRecord>& records,
    std::string* error
)
{
    std::ofstream out(filename);
    if (!out.is_open()) {
        if (error) *error = "Cannot open pair-normal CSV: " + filename;
        return false;
    }

    out << "pair_id,owner_seed_id,region,base_x,base_y,base_z,offset_x,offset_y,offset_z,normal_x,normal_y,normal_z\n";
    out << std::fixed << std::setprecision(16);
    for (const auto& r : records) {
        out << r.pair_id << ","
            << r.owner_seed_id << ","
            << r.region << ","
            << r.base[0] << "," << r.base[1] << "," << r.base[2] << ","
            << r.offset[0] << "," << r.offset[1] << "," << r.offset[2] << ","
            << r.normal[0] << "," << r.normal[1] << "," << r.normal[2] << "\n";
    }
    std::cout << "  * Wrote " << records.size()
              << " base-offset pairs with normals to " << filename << std::endl;
    return true;
}

static bool collect_power_base_seed_interface_facets(
    const std::vector<WeightedSite>& sites,
    bool require_same_region,
    bool require_owner_match,
    std::vector<std::vector<Point_3>>& facets,
    std::string* error
)
{
    facets.clear();
    if (sites.empty()) {
        if (error) *error = "No weighted sites for interface facets.";
        return false;
    }

    Regular rt;
    size_t inserted_vertices = 0;
    for (size_t site_idx = 0; site_idx < sites.size(); ++site_idx) {
        const auto& s = sites[site_idx];
        Weighted_point wp(Point_3(s.x, s.y, s.z), s.w);
        Regular::Vertex_handle vh = rt.insert(wp);
        if (vh != Regular::Vertex_handle()) {
            if (vh->info() != site_idx) {
                // same vertex can be reused; info gets overwritten by latest site index
            }
            vh->info() = site_idx;
            inserted_vertices++;
        }
    }

    std::cout << "  * Interface sites total: " << sites.size()
              << ", finite vertices in RT: " << rt.number_of_vertices()
              << ", insert handles: " << inserted_vertices << std::endl;

    double min_x = std::numeric_limits<double>::infinity();
    double min_y = std::numeric_limits<double>::infinity();
    double min_z = std::numeric_limits<double>::infinity();
    double max_x = -std::numeric_limits<double>::infinity();
    double max_y = -std::numeric_limits<double>::infinity();
    double max_z = -std::numeric_limits<double>::infinity();
    for (const auto& s : sites) {
        min_x = std::min(min_x, s.x);
        min_y = std::min(min_y, s.y);
        min_z = std::min(min_z, s.z);
        max_x = std::max(max_x, s.x);
        max_y = std::max(max_y, s.y);
        max_z = std::max(max_z, s.z);
    }
    const double global_dx = max_x - min_x;
    const double global_dy = max_y - min_y;
    const double global_dz = max_z - min_z;
    const double global_diag = std::sqrt(global_dx * global_dx + global_dy * global_dy + global_dz * global_dz);

    size_t edge_total = 0;
    size_t edge_base_seed = 0;
    size_t edge_same_region = 0;
    size_t edge_owner_match = 0;
    size_t edge_kept = 0;
    size_t far_vertices_filtered = 0;
    for (auto eit = rt.finite_edges_begin(); eit != rt.finite_edges_end(); ++eit) {
        edge_total++;
        Regular::Cell_handle c = eit->first;
        int i1 = eit->second;
        int i2 = eit->third;
        Regular::Vertex_handle v1 = c->vertex(i1);
        Regular::Vertex_handle v2 = c->vertex(i2);
        if (rt.is_infinite(v1) || rt.is_infinite(v2)) continue;

        size_t idx1 = v1->info();
        size_t idx2 = v2->info();
        if (idx1 >= sites.size() || idx2 >= sites.size()) continue;
        const auto& s1 = sites[idx1];
        const auto& s2 = sites[idx2];
        if (!((s1.type == 1 && s2.type == 0) || (s1.type == 0 && s2.type == 1))) continue;
        edge_base_seed++;
        if (require_same_region) {
            if (s1.region != s2.region) continue;
            edge_same_region++;
        } else {
            edge_same_region++;
        }
        if (require_owner_match) {
            if (s1.owner_seed_id != s2.owner_seed_id) continue;
            edge_owner_match++;
        } else {
            edge_owner_match++;
        }
        edge_kept++;

        std::vector<Point_3> facet_vertices;
        const double mx = 0.5 * (s1.x + s2.x);
        const double my = 0.5 * (s1.y + s2.y);
        const double mz = 0.5 * (s1.z + s2.z);
        const double pdx = s1.x - s2.x;
        const double pdy = s1.y - s2.y;
        const double pdz = s1.z - s2.z;
        const double pair_dist = std::sqrt(pdx * pdx + pdy * pdy + pdz * pdz);
        const double far_limit = 8.0 * pair_dist + 0.25 * global_diag;
        Regular::Cell_circulator cc = rt.incident_cells(*eit);
        Regular::Cell_circulator done = cc;
        if (cc == nullptr) continue;
        do {
            if (!rt.is_infinite(cc)) {
                Point_3 center = rt.dual(cc);
                const double cx = CGAL::to_double(center.x());
                const double cy = CGAL::to_double(center.y());
                const double cz = CGAL::to_double(center.z());
                const double dmx = cx - mx;
                const double dmy = cy - my;
                const double dmz = cz - mz;
                const double dist_mid = std::sqrt(dmx * dmx + dmy * dmy + dmz * dmz);
                if (dist_mid <= far_limit) {
                    facet_vertices.push_back(center);
                } else {
                    far_vertices_filtered++;
                }
            }
            ++cc;
        } while (cc != done);

        if (facet_vertices.size() >= 3) {
            facets.push_back(std::move(facet_vertices));
        }
    }

    std::cout << "  * Interface edge stats: total=" << edge_total
              << ", base-seed=" << edge_base_seed
              << ", same-region=" << edge_same_region
              << ", owner-match=" << edge_owner_match
              << ", kept=" << edge_kept << std::endl;
    std::cout << "  * Interface far dual vertices filtered: " << far_vertices_filtered << std::endl;
    std::cout << "  * Interface facets collected: " << facets.size() << std::endl;
    return true;
}

static bool export_power_interface_obj(
    const std::vector<WeightedSite>& sites,
    const std::string& out_obj,
    double merge_epsilon,
    std::string* error
)
{
    if (out_obj.empty()) return true;

    std::vector<std::vector<Point_3>> facets;
    if (!collect_power_base_seed_interface_facets(
            sites,
            /*require_same_region=*/true,
            /*require_owner_match=*/true,
            facets,
            error)) {
        return false;
    }

    std::ofstream out(out_obj);
    if (!out.is_open()) {
        if (error) *error = "Cannot open interface OBJ: " + out_obj;
        return false;
    }
    out << "# Power interface facets (raw facets with far dual filtering)\n";
    out << "# sites=" << sites.size() << " facets=" << facets.size() << "\n";

    std::map<std::tuple<double, double, double>, size_t> vertex_map;
    std::vector<Point_3> unique_vertices;
    const double eps = (merge_epsilon > 0.0) ? merge_epsilon : 1e-10;
    auto round_key = [&](const Point_3& p) -> std::tuple<double, double, double> {
        const double x = std::round(CGAL::to_double(p.x()) / eps) * eps;
        const double y = std::round(CGAL::to_double(p.y()) / eps) * eps;
        const double z = std::round(CGAL::to_double(p.z()) / eps) * eps;
        return std::make_tuple(x, y, z);
    };
    auto get_vertex_index = [&](const Point_3& p) -> size_t {
        auto key = round_key(p);
        auto it = vertex_map.find(key);
        if (it != vertex_map.end()) return it->second;
        size_t idx = unique_vertices.size();
        vertex_map.emplace(key, idx);
        unique_vertices.push_back(p);
        return idx;
    };

    std::vector<std::vector<size_t>> face_indices;
    face_indices.reserve(facets.size());
    for (const auto& facet : facets) {
        if (facet.size() < 3) continue;
        std::vector<size_t> inds;
        inds.reserve(facet.size());
        for (const auto& p : facet) {
            inds.push_back(get_vertex_index(p));
        }
        face_indices.push_back(std::move(inds));
    }

    out << std::fixed << std::setprecision(16);
    for (const auto& p : unique_vertices) {
        out << "v " << CGAL::to_double(p.x()) << " "
            << CGAL::to_double(p.y()) << " "
            << CGAL::to_double(p.z()) << "\n";
    }
    for (const auto& inds : face_indices) {
        if (inds.size() < 3) continue;
        out << "f";
        for (size_t idx : inds) {
            out << " " << (idx + 1);
        }
        out << "\n";
    }
    std::cout << "  * Interface OBJ written (raw-filtered): " << out_obj << std::endl;

    return true;
}

} // namespace

namespace Offset {

bool build_power_diagram_from_seed_csv(
    const std::string& seed_csv,
    const std::string& outside_obj,
    const std::string& inside_obj,
    const std::string& interface_obj,
    const std::string& basepoints_obj,
    const std::string& offsets_obj,
    const PowerDiagramParams& params,
    std::string* error
)
{
    std::cout << "=== Build power diagram from seed CSV ===" << std::endl;
    std::cout << "  * Input CSV: " << seed_csv << std::endl;

    std::ifstream file(seed_csv);
    if (!file.is_open()) {
        if (error) *error = "Cannot open seed CSV: " + seed_csv;
        return false;
    }

    std::vector<SeedRecord> seeds;
    std::string line;
    while (std::getline(file, line)) {
        SeedRecord rec;
        if (parse_seed_csv_line(line, rec)) {
            seeds.push_back(rec);
        }
    }

    if (seeds.empty()) {
        if (error) *error = "No valid seed rows parsed.";
        return false;
    }
    std::cout << "  * Parsed seeds: " << seeds.size() << std::endl;

    std::unordered_map<size_t, size_t> id_to_index;
    id_to_index.reserve(seeds.size());
    for (size_t i = 0; i < seeds.size(); ++i) {
        size_t id = seeds[i].seed_id;
        if (id == static_cast<size_t>(-1)) id = i;
        id_to_index[id] = i;
    }

    std::vector<WeightedSite> outside_sites;
    std::vector<WeightedSite> inside_sites;
    std::vector<WeightedSite> combined_sites;
    std::vector<BaseOffsetPairRecord> outside_pair_records;
    std::vector<BaseOffsetPairRecord> inside_pair_records;
    outside_sites.reserve(seeds.size() * 2);
    inside_sites.reserve(seeds.size() * 2);
    combined_sites.reserve(seeds.size() * 4);
    outside_pair_records.reserve(seeds.size());
    inside_pair_records.reserve(seeds.size());

    size_t pair_count = 0;
    std::vector<bool> visited(seeds.size(), false);
    for (size_t i = 0; i < seeds.size(); ++i) {
        if (visited[i]) continue;
        const SeedRecord& s0 = seeds[i];
        if (s0.pair_id == static_cast<size_t>(-1)) continue;
        auto it = id_to_index.find(s0.pair_id);
        if (it == id_to_index.end()) continue;
        size_t j = it->second;
        if (j >= seeds.size()) continue;
        if (j == i) continue;

        visited[i] = true;
        visited[j] = true;

        const SeedRecord& s1 = seeds[j];

        const SeedRecord* outside_seed = nullptr;
        const SeedRecord* inside_seed = nullptr;
        if (s0.region == 0 && s1.region == 1) {
            outside_seed = &s0;
            inside_seed = &s1;
        } else if (s0.region == 1 && s1.region == 0) {
            outside_seed = &s1;
            inside_seed = &s0;
        } else {
            continue;
        }

        double base[3];
        double d_out = 0.0, d_in = 0.0;
        if (!compute_basepoint(*outside_seed, *inside_seed, base, d_out, d_in)) continue;

        const double p_out[3] = {outside_seed->x, outside_seed->y, outside_seed->z};
        const double p_in[3] = {inside_seed->x, inside_seed->y, inside_seed->z};

        const double r_out = compute_R_value(params, base, 0);
        const double r_in = compute_R_value(params, base, 1);

        const double w_out = 2.0 * d_out * r_out - d_out * d_out;
        const double w_in = 2.0 * d_in * r_in - d_in * d_in;

        WeightedSite base_out;
        base_out.x = base[0];
        base_out.y = base[1];
        base_out.z = base[2];
        base_out.w = w_out;
        base_out.pair_id = outside_seed->pair_id;
        base_out.seed_id = static_cast<size_t>(-1);
        base_out.owner_seed_id = outside_seed->seed_id;
        base_out.type = 1;
        base_out.region = 0;

        WeightedSite seed_out;
        seed_out.x = p_out[0];
        seed_out.y = p_out[1];
        seed_out.z = p_out[2];
        seed_out.w = 0.0;
        seed_out.pair_id = outside_seed->pair_id;
        seed_out.seed_id = outside_seed->seed_id;
        seed_out.owner_seed_id = outside_seed->seed_id;
        seed_out.type = 0;
        seed_out.region = 0;

        WeightedSite base_in = base_out;
        base_in.w = w_in;
        base_in.seed_id = static_cast<size_t>(-1);
        base_in.owner_seed_id = inside_seed->seed_id;
        base_in.region = 1;

        WeightedSite seed_in;
        seed_in.x = p_in[0];
        seed_in.y = p_in[1];
        seed_in.z = p_in[2];
        seed_in.w = 0.0;
        seed_in.pair_id = inside_seed->pair_id;
        seed_in.seed_id = inside_seed->seed_id;
        seed_in.owner_seed_id = inside_seed->seed_id;
        seed_in.type = 0;
        seed_in.region = 1;

        outside_sites.push_back(base_out);
        outside_sites.push_back(seed_out);
        inside_sites.push_back(base_in);
        inside_sites.push_back(seed_in);
        combined_sites.push_back(base_out);
        combined_sites.push_back(seed_out);
        combined_sites.push_back(base_in);
        combined_sites.push_back(seed_in);

        BaseOffsetPairRecord out_rec;
        out_rec.pair_id = seed_out.pair_id;
        out_rec.owner_seed_id = seed_out.owner_seed_id;
        out_rec.region = 0;
        out_rec.base[0] = base_out.x;
        out_rec.base[1] = base_out.y;
        out_rec.base[2] = base_out.z;
        out_rec.offset[0] = seed_out.x;
        out_rec.offset[1] = seed_out.y;
        out_rec.offset[2] = seed_out.z;
        compute_unit_normal_base_to_offset(out_rec.base, out_rec.offset, out_rec.normal);
        outside_pair_records.push_back(out_rec);

        BaseOffsetPairRecord in_rec;
        in_rec.pair_id = seed_in.pair_id;
        in_rec.owner_seed_id = seed_in.owner_seed_id;
        in_rec.region = 1;
        in_rec.base[0] = base_in.x;
        in_rec.base[1] = base_in.y;
        in_rec.base[2] = base_in.z;
        in_rec.offset[0] = seed_in.x;
        in_rec.offset[1] = seed_in.y;
        in_rec.offset[2] = seed_in.z;
        compute_unit_normal_base_to_offset(in_rec.base, in_rec.offset, in_rec.normal);
        inside_pair_records.push_back(in_rec);

        pair_count++;
    }

    if (outside_sites.empty() || inside_sites.empty()) {
        if (error) *error = "No valid seed pairs found for power diagram.";
        return false;
    }

    std::cout << "  * Paired seeds: " << pair_count << std::endl;
    std::cout << "  * Outside sites: " << outside_sites.size()
              << ", inside sites: " << inside_sites.size()
              << ", combined sites: " << combined_sites.size() << std::endl;

    if (!write_sites_points_obj(basepoints_obj, combined_sites, 1, error)) {
        return false;
    }
    if (!write_sites_points_obj(offsets_obj, combined_sites, 0, error)) {
        return false;
    }
    if (!write_base_offset_pairs_csv("offset_pairs_outside.csv", outside_pair_records, error)) {
        return false;
    }
    if (!write_base_offset_pairs_csv("offset_pairs_inside.csv", inside_pair_records, error)) {
        return false;
    }

    std::vector<std::vector<Point_3>> outside_interface_facets;
    std::vector<std::vector<Point_3>> inside_interface_facets;
    std::string local_error;
    if (!export_power_diagram_obj(outside_sites, outside_obj, params.vertex_merge_epsilon, nullptr, &local_error)) {
        if (error) *error = local_error;
        return false;
    }
    if (!export_power_diagram_obj(inside_sites, inside_obj, params.vertex_merge_epsilon, nullptr, &local_error)) {
        if (error) *error = local_error;
        return false;
    }

    if (!collect_power_base_seed_interface_facets(
            outside_sites,
            /*require_same_region=*/true,
            /*require_owner_match=*/false,
            outside_interface_facets,
            &local_error)) {
        if (error) *error = local_error;
        return false;
    }
    if (!collect_power_base_seed_interface_facets(
            inside_sites,
            /*require_same_region=*/true,
            /*require_owner_match=*/false,
            inside_interface_facets,
            &local_error)) {
        if (error) *error = local_error;
        return false;
    }

    if (!write_facets_dedup_obj("interface_dedup.obj", outside_interface_facets, params.vertex_merge_epsilon, 170.0, &local_error)) {
        if (error) *error = local_error;
        return false;
    }
    if (!write_facets_dedup_obj("interface_dedup_inside.obj", inside_interface_facets, params.vertex_merge_epsilon, 170.0, &local_error)) {
        if (error) *error = local_error;
        return false;
    }

    if (!interface_obj.empty()) {
        std::cout << "  * Exporting interface to " << interface_obj << std::endl;
        if (!export_power_interface_obj(combined_sites, interface_obj, params.vertex_merge_epsilon, &local_error)) {
            if (error) *error = local_error;
            return false;
        }
    }

    return true;
}

} // namespace Offset
