#include <algorithm>
#include <array>
#include <chrono>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <Eigen/Core>
// libigl: orientation repair + FWN baseline
#include <igl/bfs_orient.h>
#include <igl/fast_winding_number.h>

// CGAL exact-predicate / inexact-construction kernel
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/helpers.h>

#ifdef USE_EMBREE
#include <embree3/rtcore.h>
#endif

namespace fs = std::filesystem;
using Clock = std::chrono::steady_clock;

// ============================================================
// User paths
// ============================================================
static const std::string MODEL_OBJ = "inout/obj/10_u_frame_box_voronoi.obj";
static const std::string SEED_CSV  = "inout/seeds/10_u_frame_box_seed_pairs_with_faces.csv";
static const std::string OUT_CSV   = "result_inout_exact_epick_surface_perturb.csv";
static const std::string ORIENTED_MODEL_OBJ =
    (fs::path("inout") / (fs::path(MODEL_OBJ).stem().string() + "_oriented.obj")).string();

// Query settings
static constexpr std::size_t NUM_TESTS = 1000000;

// Query points are sampled on the reconstructed mesh surface, then displaced
// by a small signed amount along the oriented face normal.
// With the current Armadillo normalized bbox, max bbox span is about 1.0,
// so 2.5e-3 means an absolute perturbation of roughly 0.0025.
static constexpr double SURFACE_PERTURB_REL = 2.5e-3;
static constexpr double SURFACE_PERTURB_MIN_RATIO = 0.10; // avoid exact-boundary samples
static constexpr unsigned long long QUERY_RANDOM_SEED = 20260430ull;

// ============================================================
// Basic data structures
// ============================================================
struct SeedPoint {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    unsigned char inside = 0; // 1 = inside seed, 0 = outside seed
};

struct BBox3d {
    double min_x = std::numeric_limits<double>::infinity();
    double min_y = std::numeric_limits<double>::infinity();
    double min_z = std::numeric_limits<double>::infinity();
    double max_x = -std::numeric_limits<double>::infinity();
    double max_y = -std::numeric_limits<double>::infinity();
    double max_z = -std::numeric_limits<double>::infinity();

    void include(double x, double y, double z) noexcept
    {
        min_x = std::min(min_x, x);
        min_y = std::min(min_y, y);
        min_z = std::min(min_z, z);
        max_x = std::max(max_x, x);
        max_y = std::max(max_y, y);
        max_z = std::max(max_z, z);
    }

    bool valid() const noexcept
    {
        return min_x <= max_x && min_y <= max_y && min_z <= max_z;
    }

    bool outside(double x, double y, double z) const noexcept
    {
        return x < min_x || x > max_x ||
               y < min_y || y > max_y ||
               z < min_z || z > max_z;
    }

    Eigen::Vector3d center() const noexcept
    {
        return Eigen::Vector3d(
            0.5 * (min_x + max_x),
            0.5 * (min_y + max_y),
            0.5 * (min_z + max_z)
        );
    }

    Eigen::Vector3d half_extent() const noexcept
    {
        Eigen::Vector3d h(
            0.5 * (max_x - min_x),
            0.5 * (max_y - min_y),
            0.5 * (max_z - min_z)
        );
        constexpr double eps = 1e-12;
        h.x() = std::max(h.x(), eps);
        h.y() = std::max(h.y(), eps);
        h.z() = std::max(h.z(), eps);
        return h;
    }
};

static double elapsed_ms(const Clock::time_point& a, const Clock::time_point& b)
{
    return std::chrono::duration<double, std::milli>(b - a).count();
}

static std::string trim_copy(std::string s)
{
    std::size_t begin = 0;
    while (begin < s.size() && std::isspace(static_cast<unsigned char>(s[begin]))) {
        ++begin;
    }
    std::size_t end = s.size();
    while (end > begin && std::isspace(static_cast<unsigned char>(s[end - 1]))) {
        --end;
    }
    return s.substr(begin, end - begin);
}

static bool parse_double_strict(const std::string& token, double& value)
{
    try {
        std::size_t pos = 0;
        const std::string t = trim_copy(token);
        if (t.empty()) return false;
        value = std::stod(t, &pos);
        return pos == t.size();
    } catch (...) {
        return false;
    }
}

static std::vector<std::string> split_csv_line(const std::string& line)
{
    std::vector<std::string> tokens;
    std::stringstream ss(line);
    std::string token;
    while (std::getline(ss, token, ',')) {
        tokens.push_back(trim_copy(token));
    }
    return tokens;
}

// ============================================================
// OBJ reader: supports v and triangular / polygonal f lines.
// Polygons are fan-triangulated.
// ============================================================
static bool parse_obj_vertex_index(const std::string& token, int vertex_count, int& out_index)
{
    const std::size_t slash = token.find('/');
    const std::string head = token.substr(0, slash);
    if (head.empty()) return false;

    int raw = 0;
    try {
        raw = std::stoi(head);
    } catch (...) {
        return false;
    }
    if (raw == 0) return false;

    const int idx = raw > 0 ? raw - 1 : vertex_count + raw;
    if (idx < 0 || idx >= vertex_count) return false;
    out_index = idx;
    return true;
}

static bool load_obj_as_triangles(
    const std::string& filename,
    Eigen::MatrixXd& vertices,
    Eigen::MatrixXi& faces
)
{
    std::ifstream in(filename);
    if (!in.is_open()) return false;

    std::vector<Eigen::Vector3d> vertex_list;
    std::vector<Eigen::Vector3i> face_list;

    std::string line;
    while (std::getline(in, line)) {
        line = trim_copy(line);
        if (line.empty() || line[0] == '#') continue;

        std::stringstream ss(line);
        std::string tag;
        ss >> tag;

        if (tag == "v") {
            Eigen::Vector3d v;
            if (!(ss >> v.x() >> v.y() >> v.z())) {
                return false;
            }
            vertex_list.push_back(v);
        } else if (tag == "f") {
            std::vector<int> polygon;
            std::string tok;
            while (ss >> tok) {
                int idx = -1;
                if (!parse_obj_vertex_index(tok, static_cast<int>(vertex_list.size()), idx)) {
                    return false;
                }
                polygon.push_back(idx);
            }
            if (polygon.size() < 3) continue;
            for (std::size_t i = 1; i + 1 < polygon.size(); ++i) {
                face_list.emplace_back(polygon[0], polygon[i], polygon[i + 1]);
            }
        }
    }

    if (vertex_list.empty() || face_list.empty()) return false;

    vertices.resize(static_cast<Eigen::Index>(vertex_list.size()), 3);
    for (std::size_t i = 0; i < vertex_list.size(); ++i) {
        vertices.row(static_cast<Eigen::Index>(i)) = vertex_list[i].transpose();
    }

    faces.resize(static_cast<Eigen::Index>(face_list.size()), 3);
    for (std::size_t i = 0; i < face_list.size(); ++i) {
        faces.row(static_cast<Eigen::Index>(i)) = face_list[i].transpose();
    }
    return true;
}

static BBox3d compute_bbox(const Eigen::MatrixXd& vertices)
{
    BBox3d bbox;
    for (Eigen::Index i = 0; i < vertices.rows(); ++i) {
        bbox.include(vertices(i, 0), vertices(i, 1), vertices(i, 2));
    }
    return bbox;
}

static double bbox_max_span(const BBox3d& bbox) noexcept
{
    if (!bbox.valid()) return 1e-12;
    return std::max({
        bbox.max_x - bbox.min_x,
        bbox.max_y - bbox.min_y,
        bbox.max_z - bbox.min_z,
        1e-12
    });
}

struct SurfaceSampleFace {
    int a = -1;
    int b = -1;
    int c = -1;
    Eigen::Vector3d normal = Eigen::Vector3d::Zero();
    double cumulative_area = 0.0;
};

static bool generate_surface_perturbed_queries(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& faces,
    std::size_t num_tests,
    double perturb_abs,
    unsigned long long random_seed,
    Eigen::MatrixXd& query_points,
    std::vector<double>& query_x,
    std::vector<double>& query_y,
    std::vector<double>& query_z,
    std::size_t& expected_inside_by_offset,
    std::size_t& skipped_degenerate_faces,
    double& total_surface_area
)
{
    std::vector<SurfaceSampleFace> sample_faces;
    sample_faces.reserve(static_cast<std::size_t>(faces.rows()));

    expected_inside_by_offset = 0;
    skipped_degenerate_faces = 0;
    total_surface_area = 0.0;

    for (Eigen::Index i = 0; i < faces.rows(); ++i) {
        const int ia = faces(i, 0);
        const int ib = faces(i, 1);
        const int ic = faces(i, 2);
        if (ia < 0 || ib < 0 || ic < 0 ||
            ia >= vertices.rows() || ib >= vertices.rows() || ic >= vertices.rows()) {
            ++skipped_degenerate_faces;
            continue;
        }

        const Eigen::Vector3d A = vertices.row(ia).transpose();
        const Eigen::Vector3d B = vertices.row(ib).transpose();
        const Eigen::Vector3d C = vertices.row(ic).transpose();
        const Eigen::Vector3d cross = (B - A).cross(C - A);
        const double double_area = cross.norm();
        if (double_area <= 1e-20) {
            ++skipped_degenerate_faces;
            continue;
        }

        total_surface_area += 0.5 * double_area;
        SurfaceSampleFace rec;
        rec.a = ia;
        rec.b = ib;
        rec.c = ic;
        rec.normal = cross / double_area; // outward after orient_faces_outward()
        rec.cumulative_area = total_surface_area;
        sample_faces.push_back(rec);
    }

    if (sample_faces.empty() || total_surface_area <= 0.0 || perturb_abs <= 0.0) {
        return false;
    }

    std::mt19937_64 rng(random_seed);
    std::uniform_real_distribution<double> area_dist(0.0, total_surface_area);
    std::uniform_real_distribution<double> unit_dist(0.0, 1.0);
    std::uniform_real_distribution<double> offset_mag_dist(
        std::max(0.0, SURFACE_PERTURB_MIN_RATIO) * perturb_abs,
        perturb_abs
    );
    std::bernoulli_distribution outside_side_dist(0.5);

    query_points.resize(static_cast<Eigen::Index>(num_tests), 3);
    query_x.resize(num_tests);
    query_y.resize(num_tests);
    query_z.resize(num_tests);

    for (std::size_t i = 0; i < num_tests; ++i) {
        const double target_area = area_dist(rng);
        auto it = std::lower_bound(
            sample_faces.begin(),
            sample_faces.end(),
            target_area,
            [](const SurfaceSampleFace& f, double value) {
                return f.cumulative_area < value;
            }
        );
        if (it == sample_faces.end()) {
            it = sample_faces.end() - 1;
        }

        const Eigen::Vector3d A = vertices.row(it->a).transpose();
        const Eigen::Vector3d B = vertices.row(it->b).transpose();
        const Eigen::Vector3d C = vertices.row(it->c).transpose();

        // Uniform point on triangle: barycentric sampling via sqrt trick.
        const double u = unit_dist(rng);
        const double v = unit_dist(rng);
        const double su = std::sqrt(u);
        const double w0 = 1.0 - su;
        const double w1 = su * (1.0 - v);
        const double w2 = su * v;
        Eigen::Vector3d q = w0 * A + w1 * B + w2 * C;

        // Positive displacement goes along outward normal, negative goes inward.
        const bool sample_outside = outside_side_dist(rng);
        const double mag = offset_mag_dist(rng);
        const double signed_offset = sample_outside ? mag : -mag;
        if (!sample_outside) {
            ++expected_inside_by_offset;
        }
        q += signed_offset * it->normal;

        query_x[i] = q.x();
        query_y[i] = q.y();
        query_z[i] = q.z();
        query_points(static_cast<Eigen::Index>(i), 0) = q.x();
        query_points(static_cast<Eigen::Index>(i), 1) = q.y();
        query_points(static_cast<Eigen::Index>(i), 2) = q.z();
    }

    return true;
}


static double signed_mesh_volume(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& faces
)
{
    if (vertices.rows() == 0 || faces.rows() == 0) {
        return 0.0;
    }

    // Translate by centroid to reduce numerical cancellation.
    const Eigen::Vector3d center = vertices.colwise().mean();
    double volume6 = 0.0;

    for (Eigen::Index i = 0; i < faces.rows(); ++i) {
        const Eigen::Vector3d a = vertices.row(faces(i, 0)).transpose() - center;
        const Eigen::Vector3d b = vertices.row(faces(i, 1)).transpose() - center;
        const Eigen::Vector3d c = vertices.row(faces(i, 2)).transpose() - center;
        volume6 += a.dot(b.cross(c));
    }

    return volume6 / 6.0;
}

static void orient_faces_outward(
    const Eigen::MatrixXd& vertices,
    Eigen::MatrixXi& faces
)
{
    if (vertices.rows() == 0 || faces.rows() == 0) {
        return;
    }

    Eigen::MatrixXi oriented_faces;
    Eigen::VectorXi components;

    // bfs_orient fixes local orientation consistency inside each connected component.
    igl::bfs_orient(faces, oriented_faces, components);
    faces = std::move(oriented_faces);

    // For a closed single-component mesh, positive signed volume usually means outward orientation.
    // If the global orientation is inward, flip all triangles.
    const double vol = signed_mesh_volume(vertices, faces);
    if (vol < 0.0) {
        faces = faces.rowwise().reverse().eval();
    }
}

static bool write_obj_triangles(
    const std::string& filename,
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& faces
)
{
    const fs::path output_path(filename);
    if (!output_path.parent_path().empty()) {
        std::error_code ec;
        fs::create_directories(output_path.parent_path(), ec);
        if (ec) {
            std::cerr << "Warning: cannot create output directory "
                      << output_path.parent_path().string() << ": " << ec.message() << "\n";
            return false;
        }
    }

    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "Warning: cannot write oriented OBJ: " << filename << "\n";
        return false;
    }

    out << std::fixed << std::setprecision(16);
    out << "# Oriented mesh exported by VoroCrust in/out benchmark\n";
    out << "# Source OBJ: " << MODEL_OBJ << "\n";
    out << "# Vertices: " << vertices.rows() << "\n";
    out << "# Faces: " << faces.rows() << "\n";

    for (Eigen::Index i = 0; i < vertices.rows(); ++i) {
        out << "v " << vertices(i, 0) << ' ' << vertices(i, 1) << ' ' << vertices(i, 2) << "\n";
    }
    for (Eigen::Index i = 0; i < faces.rows(); ++i) {
        out << "f " << (faces(i, 0) + 1) << ' '
            << (faces(i, 1) + 1) << ' '
            << (faces(i, 2) + 1) << "\n";
    }

    return true;
}

// ============================================================
// Seed CSV reader
// Expected format:
// inside_x, inside_y, inside_z, outside_x, outside_y, outside_z,
// face_v0_x, face_v0_y, face_v0_z, face_v1_x, face_v1_y, face_v1_z,
// face_v2_x, face_v2_y, face_v2_z
// Only first 6 columns are used for nearest-seed in/out classification.
// ============================================================
static bool parse_first_6_doubles_csv_line(const std::string& line, double v[6])
{
    const char* p = line.c_str();
    while (*p == ' ' || *p == '\t' || *p == '\r' || *p == '\n') ++p;

    // Empty/comment/header lines.
    if (*p == '\0' || *p == '#') return false;
    if ((*p >= 'a' && *p <= 'z') || (*p >= 'A' && *p <= 'Z')) return false;

    for (int i = 0; i < 6; ++i) {
        char* end_ptr = nullptr;
        v[i] = std::strtod(p, &end_ptr);
        if (end_ptr == p) return false;
        p = end_ptr;

        while (*p == ' ' || *p == '\t') ++p;
        if (i < 5) {
            if (*p == ',') {
                ++p;
            }
            while (*p == ' ' || *p == '\t') ++p;
        }
    }
    return true;
}

static bool load_seed_pairs_csv(const std::string& filename, std::vector<SeedPoint>& seeds)
{
    std::ifstream in(filename);
    if (!in.is_open()) return false;

    seeds.clear();
    // Armadillo has 16008 seed pairs -> 32016 seeds.  Reserving avoids repeated
    // reallocations without requiring a first pass over the CSV.
    seeds.reserve(40000);

    std::string line;
    double v[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    while (std::getline(in, line)) {
        if (!parse_first_6_doubles_csv_line(line, v)) {
            continue;
        }
        seeds.push_back(SeedPoint{v[0], v[1], v[2], 1});
        seeds.push_back(SeedPoint{v[3], v[4], v[5], 0});
    }

    return !seeds.empty();
}

// ============================================================
// Fast single-thread nearest seed index
// Exact 1-NN over labeled seeds.
// Build-time optimized version:
// - midpoint partition first, fallback to nth_element only for bad splits
// - leaf size increased to reduce node count and recursive build work
// - query layout stores only coordinates + label; original-index map is removed
// - build-only buffers are released after constructing the query layout
// ============================================================
class FastSeedIndex {
public:
    explicit FastSeedIndex(const std::vector<SeedPoint>& seeds)
    {
        reset(seeds);
    }

    void reset(const std::vector<SeedPoint>& seeds)
    {
        const std::size_t n = seeds.size();
        if (n == 0) {
            throw std::runtime_error("FastSeedIndex: empty seed set.");
        }
        if (n > static_cast<std::size_t>(invalid_index() - 1)) {
            throw std::runtime_error("FastSeedIndex: too many seeds for uint32 index.");
        }

        n_points_ = n;
        const Index nidx = static_cast<Index>(n);

        xs_.resize(n);
        ys_.resize(n);
        zs_.resize(n);
        labels_.resize(n);
        order_.resize(n);
        tmp_.resize(n);

        Node root_bbox;
        reset_bbox(root_bbox, 0, nidx);

        for (Index i = 0; i < nidx; ++i) {
            const std::size_t si = static_cast<std::size_t>(i);
            xs_[si] = seeds[si].x;
            ys_[si] = seeds[si].y;
            zs_[si] = seeds[si].z;
            labels_[si] = static_cast<unsigned char>(seeds[si].inside != 0);
            order_[si] = i;
            include_original_index(root_bbox, i);
        }

        nodes_.clear();
        split_axis_.clear();
        split_value_.clear();

        const std::size_t approx_leaves = (n + leaf_size_ - 1) / leaf_size_;
        nodes_.reserve(2 * approx_leaves + 64);
        split_axis_.reserve(2 * approx_leaves + 64);
        split_value_.reserve(2 * approx_leaves + 64);

        if (n <= linear_scan_limit_) {
            linear_only_ = true;
            root_ = invalid_index();
        } else {
            linear_only_ = false;
            root_ = build_recursive(0, nidx, root_bbox);
        }

        build_query_layout();

        // These arrays are needed only while building the tree / layout.
        // Dropping them keeps memory pressure lower for later benchmark sections.
        std::vector<double>().swap(xs_);
        std::vector<double>().swap(ys_);
        std::vector<double>().swap(zs_);
        std::vector<unsigned char>().swap(labels_);
        std::vector<Index>().swap(order_);
        std::vector<Index>().swap(tmp_);
    }

    std::size_t size() const noexcept { return n_points_; }

    inline bool inside_fast(
        double qx,
        double qy,
        double qz,
        double* out_dist2 = nullptr,
        std::size_t* out_index = nullptr
    ) const noexcept
    {
        Index best_index = invalid_index();
        unsigned char best_label = 0;
        double best_dist2 = std::numeric_limits<double>::infinity();

        if (linear_only_ || root_ == invalid_index()) {
            linear_query(qx, qy, qz, best_dist2, best_index, best_label);
        } else {
            kd_query(qx, qy, qz, best_dist2, best_index, best_label);
        }

        if (out_dist2) *out_dist2 = best_dist2;
        if (out_index) *out_index = static_cast<std::size_t>(best_index);
        return best_label != 0;
    }

    // Hot path for benchmark classification.
    // It computes the exact nearest labeled seed but does not return distance/index,
    // avoiding extra writes in the timed query loop.  Traversal is intentionally
    // the stable two-child bbox-pruned variant, not the previous near-only variant.
    inline bool inside_label_only(double qx, double qy, double qz) const noexcept
    {
        unsigned char best_label = 0;
        double best_dist2 = std::numeric_limits<double>::infinity();

        if (linear_only_ || root_ == invalid_index()) {
            linear_query_label_only(qx, qy, qz, best_dist2, best_label);
        } else {
            kd_query_label_only(qx, qy, qz, best_dist2, best_label);
        }

        return best_label != 0;
    }

private:
    using Index = std::uint32_t;

    struct Node {
        double min_x = 0.0;
        double min_y = 0.0;
        double min_z = 0.0;
        double max_x = 0.0;
        double max_y = 0.0;
        double max_z = 0.0;
        Index left = std::numeric_limits<Index>::max();
        Index right = std::numeric_limits<Index>::max();
        Index begin = 0;
        Index end = 0;

        bool is_leaf() const noexcept { return left == invalid_index(); }
    };

    // Larger leaves reduce node count and build time. 64 is a good default for
    // ~3e4 labeled seeds: build is much cheaper, while per-query time usually
    // changes only slightly because leaf scans are contiguous and unrolled.
    static constexpr Index leaf_size_ = 64;
    static constexpr std::size_t linear_scan_limit_ = 1024;
    static constexpr Index fallback_split_divisor_ = 8;
    static constexpr std::size_t stack_capacity_ = 128;
    static constexpr unsigned char invalid_axis_ = 255;

    static constexpr Index invalid_index() noexcept
    {
        return std::numeric_limits<Index>::max();
    }

    static inline void reset_bbox(Node& node, Index begin, Index end) noexcept
    {
        node.min_x = node.min_y = node.min_z = std::numeric_limits<double>::infinity();
        node.max_x = node.max_y = node.max_z = -std::numeric_limits<double>::infinity();
        node.left = invalid_index();
        node.right = invalid_index();
        node.begin = begin;
        node.end = end;
    }

    inline void include_original_index(Node& node, Index idx) const noexcept
    {
        const std::size_t i = static_cast<std::size_t>(idx);
        node.min_x = std::min(node.min_x, xs_[i]);
        node.min_y = std::min(node.min_y, ys_[i]);
        node.min_z = std::min(node.min_z, zs_[i]);
        node.max_x = std::max(node.max_x, xs_[i]);
        node.max_y = std::max(node.max_y, ys_[i]);
        node.max_z = std::max(node.max_z, zs_[i]);
    }

    static inline void include_point(Node& node, double x, double y, double z) noexcept
    {
        node.min_x = std::min(node.min_x, x);
        node.min_y = std::min(node.min_y, y);
        node.min_z = std::min(node.min_z, z);
        node.max_x = std::max(node.max_x, x);
        node.max_y = std::max(node.max_y, y);
        node.max_z = std::max(node.max_z, z);
    }

    static inline double outside_axis_distance(double q, double lo, double hi) noexcept
    {
        return q < lo ? (lo - q) : (q > hi ? (q - hi) : 0.0);
    }

    static inline double bbox_distance2(const Node& n, double qx, double qy, double qz) noexcept
    {
        const double dx = outside_axis_distance(qx, n.min_x, n.max_x);
        const double dy = outside_axis_distance(qy, n.min_y, n.max_y);
        const double dz = outside_axis_distance(qz, n.min_z, n.max_z);
        return dx * dx + dy * dy + dz * dz;
    }

    static inline double axis_value(double x, double y, double z, unsigned char axis) noexcept
    {
        return axis == 0 ? x : (axis == 1 ? y : z);
    }

    static int longest_axis(const Node& n) noexcept
    {
        const double sx = n.max_x - n.min_x;
        const double sy = n.max_y - n.min_y;
        const double sz = n.max_z - n.min_z;
        if (sy > sx && sy >= sz) return 1;
        if (sz > sx && sz > sy) return 2;
        return 0;
    }

    static double axis_midpoint(const Node& n, int axis) noexcept
    {
        if (axis == 0) return 0.5 * (n.min_x + n.max_x);
        if (axis == 1) return 0.5 * (n.min_y + n.max_y);
        return 0.5 * (n.min_z + n.max_z);
    }

    double coord(Index idx, int axis) const noexcept
    {
        const std::size_t i = static_cast<std::size_t>(idx);
        return axis == 0 ? xs_[i] : (axis == 1 ? ys_[i] : zs_[i]);
    }

    void compute_bbox(Index begin, Index end, Node& out) const noexcept
    {
        reset_bbox(out, begin, end);
        for (Index k = begin; k < end; ++k) {
            include_original_index(out, order_[static_cast<std::size_t>(k)]);
        }
    }

    bool midpoint_partition_with_bboxes(
        Index begin,
        Index end,
        int axis,
        double split,
        Index& mid,
        Node& left_bbox,
        Node& right_bbox
    )
    {
        reset_bbox(left_bbox, begin, begin);
        reset_bbox(right_bbox, begin, end);

        const Index* const order = order_.data();
        Index* const tmp = tmp_.data();
        const double* const xs = xs_.data();
        const double* const ys = ys_.data();
        const double* const zs = zs_.data();

        Index left = begin;
        Index right = end;

        for (Index k = begin; k < end; ++k) {
            const Index idx = order[static_cast<std::size_t>(k)];
            const std::size_t si = static_cast<std::size_t>(idx);
            const double x = xs[si];
            const double y = ys[si];
            const double z = zs[si];
            const double v = axis == 0 ? x : (axis == 1 ? y : z);

            if (v < split) {
                tmp[static_cast<std::size_t>(left++)] = idx;
                include_point(left_bbox, x, y, z);
            } else {
                tmp[static_cast<std::size_t>(--right)] = idx;
                include_point(right_bbox, x, y, z);
            }
        }

        mid = left;
        if (left != right || mid == begin || mid == end) {
            return false;
        }

        left_bbox.begin = begin;
        left_bbox.end = mid;
        right_bbox.begin = mid;
        right_bbox.end = end;

        std::copy(
            tmp_.begin() + static_cast<std::ptrdiff_t>(begin),
            tmp_.begin() + static_cast<std::ptrdiff_t>(end),
            order_.begin() + static_cast<std::ptrdiff_t>(begin)
        );
        return true;
    }

    Index build_recursive(Index begin, Index end, const Node& bbox)
    {
        Node node = bbox;
        node.begin = begin;
        node.end = end;
        node.left = invalid_index();
        node.right = invalid_index();

        const Index node_id = static_cast<Index>(nodes_.size());
        nodes_.push_back(node);
        split_axis_.push_back(invalid_axis_);
        split_value_.push_back(0.0);

        const Index count = end - begin;
        if (count <= leaf_size_) {
            return node_id;
        }

        const int axis = longest_axis(bbox);
        double split = axis_midpoint(bbox, axis);

        Index mid = begin;
        Node left_bbox;
        Node right_bbox;
        bool ok = midpoint_partition_with_bboxes(begin, end, axis, split, mid, left_bbox, right_bbox);

        // Very unbalanced midpoint splits produce slow queries.  Fallback to median
        // only when necessary; this keeps the common build path linear per level.
        const Index min_side = std::max<Index>(leaf_size_, count / fallback_split_divisor_);
        const bool bad_split = (!ok) || (mid <= begin + min_side) || (mid >= end - min_side);
        if (bad_split) {
            mid = begin + count / 2;
            auto first = order_.begin() + static_cast<std::ptrdiff_t>(begin);
            auto nth = order_.begin() + static_cast<std::ptrdiff_t>(mid);
            auto last = order_.begin() + static_cast<std::ptrdiff_t>(end);
            std::nth_element(first, nth, last, [&](Index a, Index b) {
                return coord(a, axis) < coord(b, axis);
            });
            split = coord(order_[static_cast<std::size_t>(mid)], axis);
            compute_bbox(begin, mid, left_bbox);
            compute_bbox(mid, end, right_bbox);
        }

        split_axis_[static_cast<std::size_t>(node_id)] = static_cast<unsigned char>(axis);
        split_value_[static_cast<std::size_t>(node_id)] = split;

        nodes_[static_cast<std::size_t>(node_id)].left = build_recursive(begin, mid, left_bbox);
        nodes_[static_cast<std::size_t>(node_id)].right = build_recursive(mid, end, right_bbox);
        return node_id;
    }

    void build_query_layout()
    {
        const std::size_t n = n_points_;
        qxs_.resize(n);
        qys_.resize(n);
        qzs_.resize(n);
        qlabels_.resize(n);

        if (linear_only_) {
            for (std::size_t i = 0; i < n; ++i) {
                qxs_[i] = xs_[i];
                qys_[i] = ys_[i];
                qzs_[i] = zs_[i];
                qlabels_[i] = labels_[i];
            }
            return;
        }

        for (std::size_t k = 0; k < n; ++k) {
            const Index idx = order_[k];
            const std::size_t i = static_cast<std::size_t>(idx);
            qxs_[k] = xs_[i];
            qys_[k] = ys_[i];
            qzs_[k] = zs_[i];
            qlabels_[k] = labels_[i];
        }
    }

    inline void test_candidate(
        std::size_t k,
        double qx,
        double qy,
        double qz,
        double& best_dist2,
        Index& best_index,
        unsigned char& best_label
    ) const noexcept
    {
        const double dx = qx - qxs_[k];
        const double dy = qy - qys_[k];
        const double dz = qz - qzs_[k];
        const double d2 = dx * dx + dy * dy + dz * dz;
        if (d2 < best_dist2) {
            best_dist2 = d2;
            // This is the compact query-layout index, not the original CSV index.
            // It is used only for optional debug/checksum output.
            best_index = static_cast<Index>(k);
            best_label = qlabels_[k];
        }
    }

    inline void scan_leaf(
        const Node& node,
        double qx,
        double qy,
        double qz,
        double& best_dist2,
        Index& best_index,
        unsigned char& best_label
    ) const noexcept
    {
        Index k = node.begin;
        const Index end = node.end;
        const Index end8 = k + ((end - k) & ~static_cast<Index>(7));

        for (; k < end8; k += 8) {
            const std::size_t kk = static_cast<std::size_t>(k);
            test_candidate(kk + 0, qx, qy, qz, best_dist2, best_index, best_label);
            test_candidate(kk + 1, qx, qy, qz, best_dist2, best_index, best_label);
            test_candidate(kk + 2, qx, qy, qz, best_dist2, best_index, best_label);
            test_candidate(kk + 3, qx, qy, qz, best_dist2, best_index, best_label);
            test_candidate(kk + 4, qx, qy, qz, best_dist2, best_index, best_label);
            test_candidate(kk + 5, qx, qy, qz, best_dist2, best_index, best_label);
            test_candidate(kk + 6, qx, qy, qz, best_dist2, best_index, best_label);
            test_candidate(kk + 7, qx, qy, qz, best_dist2, best_index, best_label);
        }
        for (; k < end; ++k) {
            test_candidate(static_cast<std::size_t>(k), qx, qy, qz, best_dist2, best_index, best_label);
        }
    }

    void linear_query(
        double qx,
        double qy,
        double qz,
        double& best_dist2,
        Index& best_index,
        unsigned char& best_label
    ) const noexcept
    {
        Node fake;
        fake.begin = 0;
        fake.end = static_cast<Index>(n_points_);
        scan_leaf(fake, qx, qy, qz, best_dist2, best_index, best_label);
    }

    inline void scan_leaf_label_only(
        const Node& node,
        double qx,
        double qy,
        double qz,
        double& best_dist2,
        unsigned char& best_label
    ) const noexcept
    {
        const double* const xs = qxs_.data();
        const double* const ys = qys_.data();
        const double* const zs = qzs_.data();
        const unsigned char* const labels = qlabels_.data();

        Index k = node.begin;
        const Index end = node.end;
        const Index end8 = k + ((end - k) & ~static_cast<Index>(7));

#define TEST_LABEL_ONLY(KK) \
        do { \
            const std::size_t kk_ = static_cast<std::size_t>(KK); \
            const double dx = qx - xs[kk_]; \
            const double dy = qy - ys[kk_]; \
            const double dz = qz - zs[kk_]; \
            const double d2 = dx * dx + dy * dy + dz * dz; \
            if (d2 < best_dist2) { \
                best_dist2 = d2; \
                best_label = labels[kk_]; \
            } \
        } while (false)

        for (; k < end8; k += 8) {
            TEST_LABEL_ONLY(k + 0);
            TEST_LABEL_ONLY(k + 1);
            TEST_LABEL_ONLY(k + 2);
            TEST_LABEL_ONLY(k + 3);
            TEST_LABEL_ONLY(k + 4);
            TEST_LABEL_ONLY(k + 5);
            TEST_LABEL_ONLY(k + 6);
            TEST_LABEL_ONLY(k + 7);
        }
        for (; k < end; ++k) {
            TEST_LABEL_ONLY(k);
        }

#undef TEST_LABEL_ONLY
    }

    void linear_query_label_only(
        double qx,
        double qy,
        double qz,
        double& best_dist2,
        unsigned char& best_label
    ) const noexcept
    {
        Node fake;
        fake.begin = 0;
        fake.end = static_cast<Index>(n_points_);
        scan_leaf_label_only(fake, qx, qy, qz, best_dist2, best_label);
    }

    void kd_query_label_only(
        double qx,
        double qy,
        double qz,
        double& best_dist2,
        unsigned char& best_label
    ) const noexcept
    {
        Index stack_nodes[stack_capacity_];
        double stack_dists[stack_capacity_];
        std::size_t stack_size = 0;

        const auto push = [&](Index node_id, double dist2) noexcept {
            if (stack_size < stack_capacity_) {
                stack_nodes[stack_size] = node_id;
                stack_dists[stack_size] = dist2;
                ++stack_size;
            }
        };

        Index node_id = root_;
        while (node_id != invalid_index()) {
            const Node& node = nodes_[static_cast<std::size_t>(node_id)];
            if (node.is_leaf()) {
                scan_leaf_label_only(node, qx, qy, qz, best_dist2, best_label);
                break;
            }

            const std::size_t ni = static_cast<std::size_t>(node_id);
            const unsigned char axis = split_axis_[ni];
            const double q_axis = axis_value(qx, qy, qz, axis);
            const bool go_left = q_axis < split_value_[ni];
            const Index near_child = go_left ? node.left : node.right;
            const Index far_child = go_left ? node.right : node.left;
            const double far_d2 = bbox_distance2(nodes_[static_cast<std::size_t>(far_child)], qx, qy, qz);
            push(far_child, far_d2);
            node_id = near_child;
        }

        while (stack_size > 0) {
            --stack_size;
            const double entry_d2 = stack_dists[stack_size];
            if (entry_d2 >= best_dist2) continue;

            const Index entry_id = stack_nodes[stack_size];
            const Node& node = nodes_[static_cast<std::size_t>(entry_id)];
            if (node.is_leaf()) {
                scan_leaf_label_only(node, qx, qy, qz, best_dist2, best_label);
                continue;
            }

            const Node& left = nodes_[static_cast<std::size_t>(node.left)];
            const Node& right = nodes_[static_cast<std::size_t>(node.right)];
            const double ld2 = bbox_distance2(left, qx, qy, qz);
            const double rd2 = bbox_distance2(right, qx, qy, qz);

            if (ld2 < rd2) {
                if (rd2 < best_dist2) push(node.right, rd2);
                if (ld2 < best_dist2) push(node.left, ld2);
            } else {
                if (ld2 < best_dist2) push(node.left, ld2);
                if (rd2 < best_dist2) push(node.right, rd2);
            }
        }
    }

    void kd_query(
        double qx,
        double qy,
        double qz,
        double& best_dist2,
        Index& best_index,
        unsigned char& best_label
    ) const noexcept
    {
        Index stack_nodes[stack_capacity_];
        double stack_dists[stack_capacity_];
        std::size_t stack_size = 0;

        const auto push = [&](Index node_id, double dist2) noexcept {
            if (stack_size < stack_capacity_) {
                stack_nodes[stack_size] = node_id;
                stack_dists[stack_size] = dist2;
                ++stack_size;
            }
        };

        Index node_id = root_;
        while (node_id != invalid_index()) {
            const Node& node = nodes_[static_cast<std::size_t>(node_id)];
            if (node.is_leaf()) {
                scan_leaf(node, qx, qy, qz, best_dist2, best_index, best_label);
                break;
            }

            const std::size_t ni = static_cast<std::size_t>(node_id);
            const unsigned char axis = split_axis_[ni];
            const double q_axis = axis_value(qx, qy, qz, axis);
            const bool go_left = q_axis < split_value_[ni];
            const Index near_child = go_left ? node.left : node.right;
            const Index far_child = go_left ? node.right : node.left;
            const double far_d2 = bbox_distance2(nodes_[static_cast<std::size_t>(far_child)], qx, qy, qz);
            push(far_child, far_d2);
            node_id = near_child;
        }

        while (stack_size > 0) {
            --stack_size;
            const double entry_d2 = stack_dists[stack_size];
            if (entry_d2 >= best_dist2) continue;

            const Index entry_id = stack_nodes[stack_size];
            const Node& node = nodes_[static_cast<std::size_t>(entry_id)];
            if (node.is_leaf()) {
                scan_leaf(node, qx, qy, qz, best_dist2, best_index, best_label);
                continue;
            }

            const Node& left = nodes_[static_cast<std::size_t>(node.left)];
            const Node& right = nodes_[static_cast<std::size_t>(node.right)];
            const double ld2 = bbox_distance2(left, qx, qy, qz);
            const double rd2 = bbox_distance2(right, qx, qy, qz);

            if (ld2 < rd2) {
                if (rd2 < best_dist2) push(node.right, rd2);
                if (ld2 < best_dist2) push(node.left, ld2);
            } else {
                if (ld2 < best_dist2) push(node.left, ld2);
                if (rd2 < best_dist2) push(node.right, rd2);
            }
        }
    }

    std::size_t n_points_ = 0;
    std::vector<double> xs_, ys_, zs_;
    std::vector<unsigned char> labels_;
    std::vector<Index> order_, tmp_;

    std::vector<double> qxs_, qys_, qzs_;
    std::vector<unsigned char> qlabels_;

    std::vector<Node> nodes_;
    std::vector<unsigned char> split_axis_;
    std::vector<double> split_value_;
    Index root_ = invalid_index();
    bool linear_only_ = true;
};

// ============================================================
// CGAL EPICK Side_of_triangle_mesh classifier
// ============================================================
class CgalEpikInOutClassifier {
public:
    using K = CGAL::Exact_predicates_inexact_constructions_kernel;
    using Point = K::Point_3;
    using Mesh = CGAL::Surface_mesh<Point>;
    using SideTester = CGAL::Side_of_triangle_mesh<Mesh, K>;

    CgalEpikInOutClassifier(const Eigen::MatrixXd& vertices, const Eigen::MatrixXi& faces)
    {
        std::vector<Mesh::Vertex_index> v_indices;
        v_indices.reserve(static_cast<std::size_t>(vertices.rows()));

        const auto t0 = Clock::now();
        for (Eigen::Index i = 0; i < vertices.rows(); ++i) {
            v_indices.push_back(mesh_.add_vertex(Point(vertices(i, 0), vertices(i, 1), vertices(i, 2))));
        }

        for (Eigen::Index i = 0; i < faces.rows(); ++i) {
            const int a = faces(i, 0);
            const int b = faces(i, 1);
            const int c = faces(i, 2);
            if (a < 0 || b < 0 || c < 0 ||
                a >= vertices.rows() || b >= vertices.rows() || c >= vertices.rows()) {
                ++skipped_faces_;
                continue;
            }

            const auto f = mesh_.add_face(
                v_indices[static_cast<std::size_t>(a)],
                v_indices[static_cast<std::size_t>(b)],
                v_indices[static_cast<std::size_t>(c)]
            );
            if (f == Mesh::null_face()) {
                ++skipped_faces_;
            }
        }

        closed_ = CGAL::is_closed(mesh_);
        tester_ = std::make_unique<SideTester>(mesh_);
        const auto t1 = Clock::now();
        build_ms_ = elapsed_ms(t0, t1);
    }

    CGAL::Bounded_side side(double x, double y, double z) const
    {
        return (*tester_)(Point(x, y, z));
    }

    double build_ms() const noexcept { return build_ms_; }
    std::size_t skipped_faces() const noexcept { return skipped_faces_; }
    bool closed() const noexcept { return closed_; }

private:
    Mesh mesh_;
    std::unique_ptr<SideTester> tester_;
    double build_ms_ = 0.0;
    std::size_t skipped_faces_ = 0;
    bool closed_ = false;
};

#ifdef USE_EMBREE
// ============================================================
// Embree ray-parity classifier
// ============================================================
struct EmbreeVertex { float x, y, z; };
struct EmbreeTriangle { unsigned int v0, v1, v2; };

class EmbreeParityClassifier {
public:
    EmbreeParityClassifier(const Eigen::MatrixXd& vertices, const Eigen::MatrixXi& faces)
    {
        const auto t0 = Clock::now();
        device_ = rtcNewDevice(nullptr);
        if (!device_) throw std::runtime_error("Embree: rtcNewDevice failed.");

        scene_ = rtcNewScene(device_);
        if (!scene_) throw std::runtime_error("Embree: rtcNewScene failed.");

        RTCGeometry geom = rtcNewGeometry(device_, RTC_GEOMETRY_TYPE_TRIANGLE);
        auto* embree_vertices = static_cast<EmbreeVertex*>(rtcSetNewGeometryBuffer(
            geom, RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3,
            sizeof(EmbreeVertex), static_cast<std::size_t>(vertices.rows())
        ));
        auto* embree_triangles = static_cast<EmbreeTriangle*>(rtcSetNewGeometryBuffer(
            geom, RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3,
            sizeof(EmbreeTriangle), static_cast<std::size_t>(faces.rows())
        ));
        if (!embree_vertices || !embree_triangles) {
            throw std::runtime_error("Embree: rtcSetNewGeometryBuffer failed.");
        }

        for (Eigen::Index i = 0; i < vertices.rows(); ++i) {
            embree_vertices[i] = EmbreeVertex{
                static_cast<float>(vertices(i, 0)),
                static_cast<float>(vertices(i, 1)),
                static_cast<float>(vertices(i, 2))
            };
        }
        for (Eigen::Index i = 0; i < faces.rows(); ++i) {
            embree_triangles[i] = EmbreeTriangle{
                static_cast<unsigned int>(faces(i, 0)),
                static_cast<unsigned int>(faces(i, 1)),
                static_cast<unsigned int>(faces(i, 2))
            };
        }

        rtcCommitGeometry(geom);
        rtcAttachGeometry(scene_, geom);
        rtcReleaseGeometry(geom);
        rtcCommitScene(scene_);

        const auto t1 = Clock::now();
        build_ms_ = elapsed_ms(t0, t1);
    }

    ~EmbreeParityClassifier()
    {
        if (scene_) rtcReleaseScene(scene_);
        if (device_) rtcReleaseDevice(device_);
    }

    EmbreeParityClassifier(const EmbreeParityClassifier&) = delete;
    EmbreeParityClassifier& operator=(const EmbreeParityClassifier&) = delete;

    bool inside(double x, double y, double z) const
    {
        constexpr float dx = 0.86054051f;
        constexpr float dy = 0.43027025f;
        constexpr float dz = 0.27017059f;

        RTCIntersectContext context;
        rtcInitIntersectContext(&context);

        unsigned int hit_count = 0;
        float tnear = 1e-6f;

        while (true) {
            RTCRayHit rayhit{};
            rayhit.ray.org_x = static_cast<float>(x);
            rayhit.ray.org_y = static_cast<float>(y);
            rayhit.ray.org_z = static_cast<float>(z);
            rayhit.ray.dir_x = dx;
            rayhit.ray.dir_y = dy;
            rayhit.ray.dir_z = dz;
            rayhit.ray.tnear = tnear;
            rayhit.ray.tfar = std::numeric_limits<float>::infinity();
            rayhit.ray.mask = 0xFFFFFFFFu;
            rayhit.ray.flags = 0;
            rayhit.hit.geomID = RTC_INVALID_GEOMETRY_ID;
            rayhit.hit.primID = RTC_INVALID_GEOMETRY_ID;

            rtcIntersect1(scene_, &context, &rayhit);
            if (rayhit.hit.geomID == RTC_INVALID_GEOMETRY_ID) break;

            ++hit_count;
            tnear = rayhit.ray.tfar + 1e-5f;
            if (hit_count > 1000000u) break;
        }

        return (hit_count & 1u) != 0u;
    }

    double build_ms() const noexcept { return build_ms_; }

private:
    RTCDevice device_ = nullptr;
    RTCScene scene_ = nullptr;
    double build_ms_ = 0.0;
};
#endif

// ============================================================
// Main benchmark
// ============================================================
int main()
{
    if (!fs::exists(MODEL_OBJ)) {
        std::cerr << "Error: model OBJ does not exist: " << MODEL_OBJ << "\n";
        return 1;
    }
    if (!fs::exists(SEED_CSV)) {
        std::cerr << "Error: seed CSV does not exist: " << SEED_CSV << "\n";
        return 1;
    }

    Eigen::MatrixXd vertices;
    Eigen::MatrixXi faces;
    if (!load_obj_as_triangles(MODEL_OBJ, vertices, faces)) {
        std::cerr << "Error: failed to read triangular OBJ: " << MODEL_OBJ << "\n";
        return 1;
    }

    const double volume_before_orient = signed_mesh_volume(vertices, faces);
    orient_faces_outward(vertices, faces);
    const double volume_after_orient = signed_mesh_volume(vertices, faces);
    if (!write_obj_triangles(ORIENTED_MODEL_OBJ, vertices, faces)) {
        return 1;
    }

    std::vector<SeedPoint> seeds;
    if (!load_seed_pairs_csv(SEED_CSV, seeds)) {
        std::cerr << "Error: failed to read seed pair CSV: " << SEED_CSV << "\n";
        return 1;
    }

    BBox3d mesh_bbox = compute_bbox(vertices);
    const double perturb_abs = SURFACE_PERTURB_REL * bbox_max_span(mesh_bbox);

    Eigen::MatrixXd query_points;
    std::vector<double> query_x, query_y, query_z;
    std::size_t expected_inside_by_offset = 0;
    std::size_t skipped_surface_sample_faces = 0;
    double total_surface_area = 0.0;

    if (!generate_surface_perturbed_queries(
            vertices,
            faces,
            NUM_TESTS,
            perturb_abs,
            QUERY_RANDOM_SEED,
            query_points,
            query_x,
            query_y,
            query_z,
            expected_inside_by_offset,
            skipped_surface_sample_faces,
            total_surface_area)) {
        std::cerr << "Error: failed to generate surface-perturbed query points.\n";
        return 1;
    }

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "[Input]\n"
              << "  mesh       = " << MODEL_OBJ << "\n"
              << "  oriented   = " << ORIENTED_MODEL_OBJ << "\n"
              << "  seeds      = " << SEED_CSV << "\n"
              << "  vertices   = " << vertices.rows() << "\n"
              << "  faces      = " << faces.rows() << "\n"
              << "  volume before orientation fix = " << volume_before_orient << "\n"
              << "  volume after orientation fix  = " << volume_after_orient << "\n"
              << "  seeds      = " << seeds.size() << " (inside/outside points)\n"
              << "  tests      = " << NUM_TESTS << "\n"
              << "  query mode = area-weighted surface sampling + signed normal perturbation\n"
              << "  perturb_abs = " << perturb_abs
              << " (" << SURFACE_PERTURB_REL << " * mesh bbox max span)\n"
              << "  expected_inside_by_offset = " << expected_inside_by_offset << "\n"
              << "  sampled_surface_area = " << total_surface_area << "\n"
              << "  skipped_surface_sample_faces = " << skipped_surface_sample_faces << "\n";

    // Ours: nearest labeled seed
    const auto seed_build_t0 = Clock::now();
    FastSeedIndex seed_classifier(seeds);
    const auto seed_build_t1 = Clock::now();
    const double seed_build_ms = elapsed_ms(seed_build_t0, seed_build_t1);

    std::vector<unsigned char> ours_inside(NUM_TESTS, 0);
    std::size_t ours_inside_count = 0;
    std::size_t ours_bbox_reject_count = 0;

    // Single-thread optimized query path.
    // Important: the timed loop only performs classification.
    // Checksums are computed after timing so they do not pollute average per-query time.
    const auto ours_t0 = Clock::now();
    for (std::size_t i = 0; i < NUM_TESTS; ++i) {
        const double x = query_x[i], y = query_y[i], z = query_z[i];

        // Exact cheap early-out: outside mesh AABB => outside closed mesh.
        if (mesh_bbox.outside(x, y, z)) {
            ours_inside[i] = 0;
            ++ours_bbox_reject_count;
            continue;
        }

        const bool inside = seed_classifier.inside_label_only(x, y, z);
        ours_inside[i] = static_cast<unsigned char>(inside);
        if (inside) ++ours_inside_count;
    }
    const auto ours_t1 = Clock::now();
    const double ours_query_ms = elapsed_ms(ours_t0, ours_t1);

    // Optional post-timing checksum to prevent accidental result neglect in future edits.
    // This is deliberately outside the measured query loop.
    std::size_t seed_index_checksum = 1469598103934665603ull;
    double seed_dist2_checksum = 0.0;
    for (std::size_t i = 0; i < NUM_TESTS; i += 9973) {
        double d2 = 0.0;
        std::size_t idx = 0;
        (void)seed_classifier.inside_fast(query_x[i], query_y[i], query_z[i], &d2, &idx);
        seed_dist2_checksum += d2;
        seed_index_checksum ^= idx + 0x9e3779b97f4a7c15ull + (seed_index_checksum << 6) + (seed_index_checksum >> 2);
    }

    // CGAL EPICK: ground-truth-ish bounded-side test
    CgalEpikInOutClassifier cgal_classifier(vertices, faces);
    if (cgal_classifier.skipped_faces() > 0) {
        std::cerr << "Warning: CGAL skipped " << cgal_classifier.skipped_faces()
                  << " invalid / duplicate / non-manifold faces while building Surface_mesh.\n";
    }
    if (!cgal_classifier.closed()) {
        std::cerr << "Warning: CGAL Surface_mesh is not closed. Inside/outside labels may be unreliable.\n";
    }

    std::vector<unsigned char> cgal_inside(NUM_TESTS, 0);
    std::vector<unsigned char> cgal_valid(NUM_TESTS, 1);
    std::size_t cgal_inside_count = 0;
    std::size_t cgal_boundary_count = 0;
    std::size_t cgal_bbox_reject_count = 0;

    const auto cgal_t0 = Clock::now();
    for (std::size_t i = 0; i < NUM_TESTS; ++i) {
        const double x = query_x[i], y = query_y[i], z = query_z[i];

        // Exact and cheap: outside the vertex AABB implies outside the closed mesh.
        if (mesh_bbox.outside(x, y, z)) {
            cgal_inside[i] = 0;
            ++cgal_bbox_reject_count;
            continue;
        }

        const CGAL::Bounded_side s = cgal_classifier.side(x, y, z);
        if (s == CGAL::ON_BOUNDARY) {
            cgal_valid[i] = 0;
            ++cgal_boundary_count;
            continue;
        }
        const bool inside = (s == CGAL::ON_BOUNDED_SIDE);
        cgal_inside[i] = static_cast<unsigned char>(inside);
        if (inside) ++cgal_inside_count;
    }
    const auto cgal_t1 = Clock::now();
    const double cgal_query_ms = elapsed_ms(cgal_t0, cgal_t1);

    // FWN baseline on reconstruction mesh
    std::vector<unsigned char> fwn_inside(NUM_TESTS, 0);
    std::size_t fwn_inside_count = 0;
    const auto fwn_t0 = Clock::now();
    Eigen::VectorXd winding;
    igl::fast_winding_number(vertices, faces, query_points, winding);
    for (std::size_t i = 0; i < NUM_TESTS; ++i) {
        const bool inside = std::abs(winding(static_cast<Eigen::Index>(i))) > 0.5;
        fwn_inside[i] = static_cast<unsigned char>(inside);
        if (inside) ++fwn_inside_count;
    }
    const auto fwn_t1 = Clock::now();
    const double fwn_total_ms = elapsed_ms(fwn_t0, fwn_t1);

#ifdef USE_EMBREE
    EmbreeParityClassifier embree_classifier(vertices, faces);
    std::vector<unsigned char> embree_inside(NUM_TESTS, 0);
    std::size_t embree_inside_count = 0;
    std::size_t embree_bbox_reject_count = 0;

    const auto embree_t0 = Clock::now();
    for (std::size_t i = 0; i < NUM_TESTS; ++i) {
        const double x = query_x[i], y = query_y[i], z = query_z[i];
        if (mesh_bbox.outside(x, y, z)) {
            embree_inside[i] = 0;
            ++embree_bbox_reject_count;
            continue;
        }
        const bool inside = embree_classifier.inside(x, y, z);
        embree_inside[i] = static_cast<unsigned char>(inside);
        if (inside) ++embree_inside_count;
    }
    const auto embree_t1 = Clock::now();
    const double embree_query_ms = elapsed_ms(embree_t0, embree_t1);
#endif

    // Accuracy against CGAL, excluding CGAL boundary samples.
    std::size_t effective = 0;
    std::size_t ours_mismatch = 0;
    std::size_t fwn_mismatch = 0;
#ifdef USE_EMBREE
    std::size_t embree_mismatch = 0;
#endif

    for (std::size_t i = 0; i < NUM_TESTS; ++i) {
        if (!cgal_valid[i]) continue;
        ++effective;
        if (ours_inside[i] != cgal_inside[i]) ++ours_mismatch;
        if (fwn_inside[i] != cgal_inside[i]) ++fwn_mismatch;
#ifdef USE_EMBREE
        if (embree_inside[i] != cgal_inside[i]) ++embree_mismatch;
#endif
    }

    const auto acc = [&](std::size_t mismatch) -> double {
        if (effective == 0) return 100.0;
        return 100.0 * (1.0 - static_cast<double>(mismatch) / static_cast<double>(effective));
    };

    std::cout << "\n[Timing]\n";
    std::cout << std::left << std::setw(18) << "Method"
              << std::right << std::setw(18) << "Build/Prep(ms)"
              << std::setw(18) << "Query(ms)"
              << std::setw(18) << "Avg(us/query)" << "\n";
    std::cout << std::string(72, '-') << "\n";
    std::cout << std::left << std::setw(18) << "Ours SeedNN-stable"
              << std::right << std::setw(18) << seed_build_ms
              << std::setw(18) << ours_query_ms
              << std::setw(18) << (ours_query_ms * 1000.0 / static_cast<double>(NUM_TESTS)) << "\n";
    std::cout << std::left << std::setw(18) << "CGAL EPICK"
              << std::right << std::setw(18) << cgal_classifier.build_ms()
              << std::setw(18) << cgal_query_ms
              << std::setw(18) << (cgal_query_ms * 1000.0 / static_cast<double>(NUM_TESTS)) << "\n";
    std::cout << std::left << std::setw(18) << "FWN"
              << std::right << std::setw(18) << 0.0
              << std::setw(18) << fwn_total_ms
              << std::setw(18) << (fwn_total_ms * 1000.0 / static_cast<double>(NUM_TESTS)) << "\n";
#ifdef USE_EMBREE
    std::cout << std::left << std::setw(18) << "Embree parity"
              << std::right << std::setw(18) << embree_classifier.build_ms()
              << std::setw(18) << embree_query_ms
              << std::setw(18) << (embree_query_ms * 1000.0 / static_cast<double>(NUM_TESTS)) << "\n";
#endif

    std::cout << "\n[Counts / Accuracy vs CGAL EPICK]\n";
    std::cout << "  effective_tests       = " << effective << " / " << NUM_TESTS << "\n";
    std::cout << "  cgal_boundary_skipped = " << cgal_boundary_count << "\n";
    std::cout << "  cgal_bbox_rejected    = " << cgal_bbox_reject_count << "\n";
    std::cout << "  ours_inside_count     = " << ours_inside_count << "\n";
    std::cout << "  ours_bbox_rejected    = " << ours_bbox_reject_count << "\n";
    std::cout << "  cgal_inside_count     = " << cgal_inside_count << "\n";
    std::cout << "  fwn_inside_count      = " << fwn_inside_count << "\n";
#ifdef USE_EMBREE
    std::cout << "  embree_inside_count   = " << embree_inside_count << "\n";
    std::cout << "  embree_bbox_rejected  = " << embree_bbox_reject_count << "\n";
#endif
    std::cout << "  ours_mismatches       = " << ours_mismatch << ", accuracy = " << acc(ours_mismatch) << "%\n";
    std::cout << "  fwn_mismatches        = " << fwn_mismatch << ", accuracy = " << acc(fwn_mismatch) << "%\n";
#ifdef USE_EMBREE
    std::cout << "  embree_mismatches     = " << embree_mismatch << ", accuracy = " << acc(embree_mismatch) << "%\n";
#endif
    std::cout << "  seed_index_checksum   = " << seed_index_checksum << "\n";
    std::cout << "  seed_dist2_checksum   = " << seed_dist2_checksum << "\n";

    const bool csv_exists = fs::exists(OUT_CSV);
    std::ofstream out(OUT_CSV, std::ios::app);
    if (!out.is_open()) {
        std::cerr << "Warning: cannot write CSV output: " << OUT_CSV << "\n";
        return 0;
    }
    if (!csv_exists) {
        out << "mesh,seed_csv,vertices,faces,seeds,total_tests,effective_tests,"
            << "query_mode,surface_perturb_abs,expected_inside_by_offset,sampled_surface_area,skipped_surface_sample_faces,"
            << "ours_build_ms,ours_query_ms,ours_avg_us,ours_bbox_rejected,"
            << "cgal_build_ms,cgal_query_ms,cgal_avg_us,cgal_boundary_skipped,cgal_bbox_rejected,cgal_skipped_faces,cgal_closed,"
            << "fwn_total_ms,fwn_avg_us,";
#ifdef USE_EMBREE
        out << "embree_build_ms,embree_query_ms,embree_avg_us,embree_bbox_rejected,";
#endif
        out << "ours_inside,cgal_inside,fwn_inside,";
#ifdef USE_EMBREE
        out << "embree_inside,";
#endif
        out << "ours_mismatch,ours_accuracy,fwn_mismatch,fwn_accuracy,";
#ifdef USE_EMBREE
        out << "embree_mismatch,embree_accuracy,";
#endif
        out << "seed_index_checksum,seed_dist2_checksum\n";
    }

    out << std::fixed << std::setprecision(9)
        << MODEL_OBJ << ','
        << SEED_CSV << ','
        << vertices.rows() << ','
        << faces.rows() << ','
        << seeds.size() << ','
        << NUM_TESTS << ','
        << effective << ','
        << "surface_normal_perturb" << ','
        << perturb_abs << ','
        << expected_inside_by_offset << ','
        << total_surface_area << ','
        << skipped_surface_sample_faces << ','
        << seed_build_ms << ','
        << ours_query_ms << ','
        << (ours_query_ms * 1000.0 / static_cast<double>(NUM_TESTS)) << ','
        << ours_bbox_reject_count << ','
        << cgal_classifier.build_ms() << ','
        << cgal_query_ms << ','
        << (cgal_query_ms * 1000.0 / static_cast<double>(NUM_TESTS)) << ','
        << cgal_boundary_count << ','
        << cgal_bbox_reject_count << ','
        << cgal_classifier.skipped_faces() << ','
        << (cgal_classifier.closed() ? 1 : 0) << ','
        << fwn_total_ms << ','
        << (fwn_total_ms * 1000.0 / static_cast<double>(NUM_TESTS)) << ',';
#ifdef USE_EMBREE
    out << embree_classifier.build_ms() << ','
        << embree_query_ms << ','
        << (embree_query_ms * 1000.0 / static_cast<double>(NUM_TESTS)) << ','
        << embree_bbox_reject_count << ',';
#endif
    out << ours_inside_count << ','
        << cgal_inside_count << ','
        << fwn_inside_count << ',';
#ifdef USE_EMBREE
    out << embree_inside_count << ',';
#endif
    out << ours_mismatch << ','
        << acc(ours_mismatch) << ','
        << fwn_mismatch << ','
        << acc(fwn_mismatch) << ',';
#ifdef USE_EMBREE
    out << embree_mismatch << ','
        << acc(embree_mismatch) << ',';
#endif
    out << seed_index_checksum << ','
        << seed_dist2_checksum << '\n';

    std::cout << "\nDone. Results appended to " << OUT_CSV << "\n";
    return 0;
}
