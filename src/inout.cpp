#include "Inout.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <Eigen/Core>
#include <igl/bfs_orient.h>
#include <igl/fast_winding_number.h>
#include <igl/readOBJ.h>
#include <nanoflann.hpp>

#ifdef USE_CGAL
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_traits_3.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/helpers.h>
#endif

#ifdef USE_EMBREE
#include <embree3/rtcore.h>
#endif

namespace {

struct SeedPoint {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    double r = 0.0;
    size_t region = 0;
};

enum class SeedCsvLayout {
    Unknown,
    SingleSeed,
    SeedPair
};

struct SeedPointCloud {
    const std::vector<SeedPoint>& points;

    std::size_t kdtree_get_point_count() const
    {
        return points.size();
    }

    double kdtree_get_pt(const std::size_t idx, const std::size_t dim) const
    {
        if (dim == 0) return points[idx].x;
        if (dim == 1) return points[idx].y;
        return points[idx].z;
    }

    template <class BBOX>
    bool kdtree_get_bbox(BBOX&) const
    {
        return false;
    }
};

using SeedKdTree = nanoflann::KDTreeSingleIndexAdaptor<
    nanoflann::L2_Simple_Adaptor<double, SeedPointCloud>,
    SeedPointCloud,
    3,
    std::size_t
>;


class FastSeedIndex {
public:
    explicit FastSeedIndex(const std::vector<SeedPoint>& points)
    {
        reset(points);
    }

    void reset(const std::vector<SeedPoint>& points)
    {
        const std::size_t n_size = points.size();
        n_points_ = n_size;
        if (n_size > static_cast<std::size_t>(std::numeric_limits<Index>::max() - 1)) {
            throw std::runtime_error("FastSeedIndex: too many seed points for uint32 index.");
        }

        const Index n = static_cast<Index>(n_size);

        xs_.resize(n_size);
        ys_.resize(n_size);
        zs_.resize(n_size);
        regions_.resize(n_size);
        order_.resize(n_size);
        temp_order_.resize(n_size);

        Node root_bbox;
        if (n > 0) {
            reset_bbox(root_bbox, 0, n);
        }

        for (Index i = 0; i < n; ++i) {
            const std::size_t si = static_cast<std::size_t>(i);
            xs_[si] = points[si].x;
            ys_[si] = points[si].y;
            zs_[si] = points[si].z;
            regions_[si] = static_cast<unsigned char>(points[si].region != 0);
            order_[si] = i;
            include_index_in_bbox(root_bbox, i);
        }

        nodes_.clear();
        split_axis_.clear();
        split_value_.clear();
        root_ = invalid_index();
        linear_only_ = (n_size <= linear_scan_limit_);

        if (!linear_only_ && n > 0) {
            nodes_.reserve(2 * ((n_size + leaf_size_ - 1) / leaf_size_) + 64);
            split_axis_.reserve(nodes_.capacity());
            split_value_.reserve(nodes_.capacity());
            root_ = build_recursive(0, n, root_bbox);
        }

        // Query-time layout optimization:
        // reorder coordinates by the final KD-tree order, so each leaf scans
        // qxs_/qys_/qzs_ sequentially instead of xs_[order_[k]] indirectly.
        build_query_layout();

        // Build scratch is not needed after construction.  This reduces the
        // resident memory footprint before the benchmark query loop starts.
        std::vector<Index>().swap(temp_order_);

    }

    std::size_t size() const noexcept
    {
        return n_points_;
    }

    std::size_t region(std::size_t idx) const noexcept
    {
        return static_cast<std::size_t>(regions_[idx]);
    }

    const double* x_data() const noexcept { return xs_.empty() ? qxs_.data() : xs_.data(); }
    const double* y_data() const noexcept { return ys_.empty() ? qys_.data() : ys_.data(); }
    const double* z_data() const noexcept { return zs_.empty() ? qzs_.data() : zs_.data(); }

    std::size_t nearest_index(double qx, double qy, double qz, double* out_dist2 = nullptr) const
    {
        unsigned char dummy_region = 0;
        return nearest_impl<false>(qx, qy, qz, out_dist2, dummy_region);
    }

    std::size_t nearest_index_and_region(
        double qx,
        double qy,
        double qz,
        unsigned char& out_region,
        double* out_dist2 = nullptr
    ) const
    {
        return nearest_impl<true>(qx, qy, qz, out_dist2, out_region);
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
        Index left = invalid_index();
        Index right = invalid_index();
        Index begin = 0;
        Index end = 0;

        bool is_leaf() const noexcept
        {
            return left == invalid_index();
        }
    };

    // 40 is a compromise between traversal depth and leaf scan length for
    // 1e5 queries / <=2e5 seeds. Try 32 if your seeds are highly clustered.
    static constexpr Index leaf_size_ = 40;
    static constexpr std::size_t linear_scan_limit_ = 1024;
    static constexpr Index fallback_split_divisor_ = 8;
    static constexpr std::size_t stack_capacity_ = 128;
    static constexpr unsigned char invalid_axis_ = 255;

    static constexpr Index invalid_index() noexcept
    {
        return std::numeric_limits<Index>::max();
    }

    static constexpr std::size_t invalid_size_t() noexcept
    {
        return std::numeric_limits<std::size_t>::max();
    }

    static inline double outside_axis_distance(double q, double lo, double hi) noexcept
    {
        return q < lo ? (lo - q) : (q > hi ? (q - hi) : 0.0);
    }

    static inline double bbox_distance2(const Node& node, double qx, double qy, double qz) noexcept
    {
        const double dx = outside_axis_distance(qx, node.min_x, node.max_x);
        const double dy = outside_axis_distance(qy, node.min_y, node.max_y);
        const double dz = outside_axis_distance(qz, node.min_z, node.max_z);
        return dx * dx + dy * dy + dz * dz;
    }

    static int longest_axis(const Node& node) noexcept
    {
        const double span_x = node.max_x - node.min_x;
        const double span_y = node.max_y - node.min_y;
        const double span_z = node.max_z - node.min_z;
        int axis = 0;
        if (span_y > span_x && span_y >= span_z) {
            axis = 1;
        } else if (span_z > span_x && span_z > span_y) {
            axis = 2;
        }
        return axis;
    }

    static double axis_midpoint(const Node& node, int axis) noexcept
    {
        if (axis == 0) return 0.5 * (node.min_x + node.max_x);
        if (axis == 1) return 0.5 * (node.min_y + node.max_y);
        return 0.5 * (node.min_z + node.max_z);
    }

    double coord_axis(Index idx, int axis) const noexcept
    {
        const std::size_t si = static_cast<std::size_t>(idx);
        if (axis == 0) return xs_[si];
        if (axis == 1) return ys_[si];
        return zs_[si];
    }

    static inline double query_axis(double qx, double qy, double qz, unsigned char axis) noexcept
    {
        return axis == 0 ? qx : (axis == 1 ? qy : qz);
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

    inline void include_index_in_bbox(Node& node, Index idx) const noexcept
    {
        const std::size_t si = static_cast<std::size_t>(idx);
        const double x = xs_[si];
        const double y = ys_[si];
        const double z = zs_[si];

        node.min_x = std::min(node.min_x, x);
        node.min_y = std::min(node.min_y, y);
        node.min_z = std::min(node.min_z, z);
        node.max_x = std::max(node.max_x, x);
        node.max_y = std::max(node.max_y, y);
        node.max_z = std::max(node.max_z, z);
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
        double split_value = axis_midpoint(bbox, axis);

        Index mid = begin;
        Node left_bbox;
        Node right_bbox;
        bool partition_ok = midpoint_partition_with_bboxes(
            begin,
            end,
            axis,
            split_value,
            mid,
            left_bbox,
            right_bbox
        );

        const Index min_side = std::max<Index>(leaf_size_, count / fallback_split_divisor_);
        const bool bad_split = (!partition_ok) || (mid <= begin + min_side) || (mid >= end - min_side);
        if (bad_split) {
            mid = begin + count / 2;
            nth_element_axis(begin, mid, end, axis);
            split_value = coord_axis(order_[static_cast<std::size_t>(mid)], axis);
            compute_bbox(begin, mid, left_bbox);
            compute_bbox(mid, end, right_bbox);
        }

        split_axis_[static_cast<std::size_t>(node_id)] = static_cast<unsigned char>(axis);
        split_value_[static_cast<std::size_t>(node_id)] = split_value;

        nodes_[static_cast<std::size_t>(node_id)].left = build_recursive(begin, mid, left_bbox);
        nodes_[static_cast<std::size_t>(node_id)].right = build_recursive(mid, end, right_bbox);
        return node_id;
    }

    bool midpoint_partition_with_bboxes(
        Index begin,
        Index end,
        int axis,
        double split_value,
        Index& mid,
        Node& left_bbox,
        Node& right_bbox
    )
    {
        reset_bbox(left_bbox, begin, begin);
        reset_bbox(right_bbox, begin, end);

        const Index* const order = order_.data();
        Index* const scratch = temp_order_.data();
        const double* const xs = xs_.data();
        const double* const ys = ys_.data();
        const double* const zs = zs_.data();

        Index left = begin;
        Index right = end;

        if (axis == 0) {
            for (Index k = begin; k < end; ++k) {
                const Index idx = order[static_cast<std::size_t>(k)];
                const double v = xs[static_cast<std::size_t>(idx)];
                if (v < split_value) {
                    scratch[static_cast<std::size_t>(left++)] = idx;
                    include_index_in_bbox(left_bbox, idx);
                } else {
                    scratch[static_cast<std::size_t>(--right)] = idx;
                    include_index_in_bbox(right_bbox, idx);
                }
            }
        } else if (axis == 1) {
            for (Index k = begin; k < end; ++k) {
                const Index idx = order[static_cast<std::size_t>(k)];
                const double v = ys[static_cast<std::size_t>(idx)];
                if (v < split_value) {
                    scratch[static_cast<std::size_t>(left++)] = idx;
                    include_index_in_bbox(left_bbox, idx);
                } else {
                    scratch[static_cast<std::size_t>(--right)] = idx;
                    include_index_in_bbox(right_bbox, idx);
                }
            }
        } else {
            for (Index k = begin; k < end; ++k) {
                const Index idx = order[static_cast<std::size_t>(k)];
                const double v = zs[static_cast<std::size_t>(idx)];
                if (v < split_value) {
                    scratch[static_cast<std::size_t>(left++)] = idx;
                    include_index_in_bbox(left_bbox, idx);
                } else {
                    scratch[static_cast<std::size_t>(--right)] = idx;
                    include_index_in_bbox(right_bbox, idx);
                }
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
            temp_order_.begin() + static_cast<std::ptrdiff_t>(begin),
            temp_order_.begin() + static_cast<std::ptrdiff_t>(end),
            order_.begin() + static_cast<std::ptrdiff_t>(begin)
        );
        return true;
    }

    void nth_element_axis(Index begin, Index mid, Index end, int axis)
    {
        auto first = order_.begin() + static_cast<std::ptrdiff_t>(begin);
        auto nth = order_.begin() + static_cast<std::ptrdiff_t>(mid);
        auto last = order_.begin() + static_cast<std::ptrdiff_t>(end);

        if (axis == 0) {
            const double* const xs = xs_.data();
            std::nth_element(first, nth, last, [&](Index lhs, Index rhs) {
                return xs[static_cast<std::size_t>(lhs)] < xs[static_cast<std::size_t>(rhs)];
            });
        } else if (axis == 1) {
            const double* const ys = ys_.data();
            std::nth_element(first, nth, last, [&](Index lhs, Index rhs) {
                return ys[static_cast<std::size_t>(lhs)] < ys[static_cast<std::size_t>(rhs)];
            });
        } else {
            const double* const zs = zs_.data();
            std::nth_element(first, nth, last, [&](Index lhs, Index rhs) {
                return zs[static_cast<std::size_t>(lhs)] < zs[static_cast<std::size_t>(rhs)];
            });
        }
    }

    void compute_bbox(Index begin, Index end, Node& node) const
    {
        double min_x = std::numeric_limits<double>::infinity();
        double min_y = std::numeric_limits<double>::infinity();
        double min_z = std::numeric_limits<double>::infinity();
        double max_x = -std::numeric_limits<double>::infinity();
        double max_y = -std::numeric_limits<double>::infinity();
        double max_z = -std::numeric_limits<double>::infinity();

        const double* const xs = xs_.data();
        const double* const ys = ys_.data();
        const double* const zs = zs_.data();
        const Index* const order = order_.data();

        for (Index k = begin; k < end; ++k) {
            const std::size_t si = static_cast<std::size_t>(order[static_cast<std::size_t>(k)]);
            const double x = xs[si];
            const double y = ys[si];
            const double z = zs[si];

            min_x = std::min(min_x, x);
            min_y = std::min(min_y, y);
            min_z = std::min(min_z, z);
            max_x = std::max(max_x, x);
            max_y = std::max(max_y, y);
            max_z = std::max(max_z, z);
        }

        node.min_x = min_x;
        node.min_y = min_y;
        node.min_z = min_z;
        node.max_x = max_x;
        node.max_y = max_y;
        node.max_z = max_z;
        node.left = invalid_index();
        node.right = invalid_index();
        node.begin = begin;
        node.end = end;
    }

    void build_query_layout()
    {
        const std::size_t n = order_.size();
        qxs_.resize(n);
        qys_.resize(n);
        qzs_.resize(n);
        qregions_.resize(n);
        qoriginal_index_.resize(n);

        const double* const xs = xs_.data();
        const double* const ys = ys_.data();
        const double* const zs = zs_.data();
        const unsigned char* const regs = regions_.data();
        const Index* const order = order_.data();

        for (std::size_t k = 0; k < n; ++k) {
            const Index idx = order[k];
            const std::size_t si = static_cast<std::size_t>(idx);
            qxs_[k] = xs[si];
            qys_[k] = ys[si];
            qzs_[k] = zs[si];
            qregions_[k] = regs[si];
            qoriginal_index_[k] = idx;
        }
    }

    template <bool NeedRegion>
    inline void update_candidate(
        std::size_t k,
        double qx,
        double qy,
        double qz,
        double& best_dist2,
        Index& best_original_idx,
        unsigned char& best_region
    ) const noexcept
    {
        const double dx = qx - qxs_[k];
        const double dy = qy - qys_[k];
        const double dz = qz - qzs_[k];
        const double d2 = dx * dx + dy * dy + dz * dz;
        if (d2 < best_dist2) {
            best_dist2 = d2;
            best_original_idx = qoriginal_index_[k];
            if constexpr (NeedRegion) best_region = qregions_[k];
        }
    }

    template <bool NeedRegion>
    inline void scan_leaf(
        const Node& node,
        double qx,
        double qy,
        double qz,
        double& best_dist2,
        Index& best_original_idx,
        unsigned char& best_region
    ) const noexcept
    {
        const double* const xs = qxs_.data();
        const double* const ys = qys_.data();
        const double* const zs = qzs_.data();
        const Index* const orig = qoriginal_index_.data();
        const unsigned char* const regs = qregions_.data();

        auto test_one = [&](std::size_t kk) noexcept {
            const double dx = qx - xs[kk];
            const double dy = qy - ys[kk];
            const double dz = qz - zs[kk];
            const double d2 = dx * dx + dy * dy + dz * dz;
            if (d2 < best_dist2) {
                best_dist2 = d2;
                best_original_idx = orig[kk];
                if constexpr (NeedRegion) best_region = regs[kk];
            }
        };

        Index k = node.begin;
        const Index end = node.end;
        const Index end8 = k + ((end - k) & ~static_cast<Index>(7));

        for (; k < end8; k += 8) {
            const std::size_t k0 = static_cast<std::size_t>(k);
            test_one(k0);
            test_one(k0 + 1);
            test_one(k0 + 2);
            test_one(k0 + 3);
            test_one(k0 + 4);
            test_one(k0 + 5);
            test_one(k0 + 6);
            test_one(k0 + 7);
        }

        for (; k < end; ++k) {
            test_one(static_cast<std::size_t>(k));
        }
    }

    static inline bool push_stack(
        Index* stack_nodes,
        double* stack_dists,
        std::size_t& stack_size,
        Index node_id,
        double bbox_dist2
    ) noexcept
    {
        if (stack_size >= stack_capacity_) return false;
        stack_nodes[stack_size] = node_id;
        stack_dists[stack_size] = bbox_dist2;
        ++stack_size;
        return true;
    }

    template <bool NeedRegion>
    std::size_t nearest_impl(
        double qx,
        double qy,
        double qz,
        double* out_dist2,
        unsigned char& out_region
    ) const
    {
        if (n_points_ == 0) {
            if (out_dist2) *out_dist2 = std::numeric_limits<double>::infinity();
            if constexpr (NeedRegion) out_region = 0;
            return invalid_size_t();
        }

        if (linear_only_ || root_ == invalid_index()) {
            return linear_nearest<NeedRegion>(qx, qy, qz, out_dist2, out_region);
        }

        const Node* const nodes = nodes_.data();
        const unsigned char* const split_axis = split_axis_.data();
        const double* const split_value = split_value_.data();

        Index best_original_idx = invalid_index();
        unsigned char best_region = 0;
        double best_dist2 = std::numeric_limits<double>::infinity();

        Index stack_nodes[stack_capacity_];
        double stack_dists[stack_capacity_];
        std::size_t stack_size = 0;

        // First descent: use split-plane order instead of computing both child
        // bbox distances at every level.  This cheaply obtains an initial best.
        Index node_id = root_;
        while (true) {
            const Node& node = nodes[static_cast<std::size_t>(node_id)];
            if (node.is_leaf()) {
                scan_leaf<NeedRegion>(node, qx, qy, qz, best_dist2, best_original_idx, best_region);
                break;
            }

            const std::size_t ni = static_cast<std::size_t>(node_id);
            const unsigned char axis = split_axis[ni];
            const double q_axis = query_axis(qx, qy, qz, axis);
            const bool go_left = q_axis < split_value[ni];
            const Index near_child = go_left ? node.left : node.right;
            const Index far_child = go_left ? node.right : node.left;
            const double far_dist2 = bbox_distance2(nodes[static_cast<std::size_t>(far_child)], qx, qy, qz);
            if (!push_stack(stack_nodes, stack_dists, stack_size, far_child, far_dist2)) {
                return linear_nearest<NeedRegion>(qx, qy, qz, out_dist2, out_region);
            }
            node_id = near_child;
        }

        while (stack_size > 0) {
            --stack_size;
            const double entry_dist2 = stack_dists[stack_size];
            if (entry_dist2 >= best_dist2) {
                continue;
            }

            const Index entry_node_id = stack_nodes[stack_size];
            const Node& node = nodes[static_cast<std::size_t>(entry_node_id)];
            if (node.is_leaf()) {
                scan_leaf<NeedRegion>(node, qx, qy, qz, best_dist2, best_original_idx, best_region);
                continue;
            }

            const Node& left = nodes[static_cast<std::size_t>(node.left)];
            const Node& right = nodes[static_cast<std::size_t>(node.right)];
            const double left_dist2 = bbox_distance2(left, qx, qy, qz);
            const double right_dist2 = bbox_distance2(right, qx, qy, qz);

            if (left_dist2 < right_dist2) {
                if (right_dist2 < best_dist2 &&
                    !push_stack(stack_nodes, stack_dists, stack_size, node.right, right_dist2)) {
                    return linear_nearest<NeedRegion>(qx, qy, qz, out_dist2, out_region);
                }
                if (left_dist2 < best_dist2 &&
                    !push_stack(stack_nodes, stack_dists, stack_size, node.left, left_dist2)) {
                    return linear_nearest<NeedRegion>(qx, qy, qz, out_dist2, out_region);
                }
            } else {
                if (left_dist2 < best_dist2 &&
                    !push_stack(stack_nodes, stack_dists, stack_size, node.left, left_dist2)) {
                    return linear_nearest<NeedRegion>(qx, qy, qz, out_dist2, out_region);
                }
                if (right_dist2 < best_dist2 &&
                    !push_stack(stack_nodes, stack_dists, stack_size, node.right, right_dist2)) {
                    return linear_nearest<NeedRegion>(qx, qy, qz, out_dist2, out_region);
                }
            }
        }

        if (best_original_idx == invalid_index()) {
            return linear_nearest<NeedRegion>(qx, qy, qz, out_dist2, out_region);
        }

        if (out_dist2) *out_dist2 = best_dist2;
        if constexpr (NeedRegion) out_region = best_region;
        return static_cast<std::size_t>(best_original_idx);
    }

    template <bool NeedRegion>
    std::size_t linear_nearest(
        double qx,
        double qy,
        double qz,
        double* out_dist2,
        unsigned char& out_region
    ) const
    {
        std::size_t best_idx = 0;
        unsigned char best_region = 0;
        double best_dist2 = std::numeric_limits<double>::infinity();

        if (!xs_.empty()) {
            const std::size_t n = xs_.size();
            const double* const xs = xs_.data();
            const double* const ys = ys_.data();
            const double* const zs = zs_.data();
            const unsigned char* const regs = regions_.data();

            for (std::size_t i = 0; i < n; ++i) {
                const double dx = qx - xs[i];
                const double dy = qy - ys[i];
                const double dz = qz - zs[i];
                const double d2 = dx * dx + dy * dy + dz * dz;
                if (d2 < best_dist2) {
                    best_dist2 = d2;
                    best_idx = i;
                    if constexpr (NeedRegion) best_region = regs[i];
                }
            }
        } else {
            // Emergency exact fallback. It should almost never run, but preserves
            // correctness if the explicit DFS stack ever overflows.
            const std::size_t n = qxs_.size();
            const double* const xs = qxs_.data();
            const double* const ys = qys_.data();
            const double* const zs = qzs_.data();
            const Index* const orig = qoriginal_index_.data();
            const unsigned char* const regs = qregions_.data();

            for (std::size_t i = 0; i < n; ++i) {
                const double dx = qx - xs[i];
                const double dy = qy - ys[i];
                const double dz = qz - zs[i];
                const double d2 = dx * dx + dy * dy + dz * dz;
                if (d2 < best_dist2) {
                    best_dist2 = d2;
                    best_idx = static_cast<std::size_t>(orig[i]);
                    if constexpr (NeedRegion) best_region = regs[i];
                }
            }
        }

        if (out_dist2) *out_dist2 = best_dist2;
        if constexpr (NeedRegion) out_region = best_region;
        return best_idx;
    }

    std::size_t n_points_ = 0;

    std::vector<double> xs_;
    std::vector<double> ys_;
    std::vector<double> zs_;
    std::vector<unsigned char> regions_;
    std::vector<Index> order_;
    std::vector<Index> temp_order_;

    std::vector<double> qxs_;
    std::vector<double> qys_;
    std::vector<double> qzs_;
    std::vector<unsigned char> qregions_;
    std::vector<Index> qoriginal_index_;

    std::vector<Node> nodes_;
    std::vector<unsigned char> split_axis_;
    std::vector<double> split_value_;
    Index root_ = invalid_index();
    bool linear_only_ = true;
};



struct FaceRecord {
    Eigen::Vector3d a;
    Eigen::Vector3d b;
    Eigen::Vector3d c;
    Eigen::Vector3d centroid;
    Eigen::Vector3d normal;
};

struct FaceCentroidCloud {
    const std::vector<FaceRecord>& faces;

    std::size_t kdtree_get_point_count() const
    {
        return faces.size();
    }

    double kdtree_get_pt(const std::size_t idx, const std::size_t dim) const
    {
        return faces[idx].centroid(static_cast<Eigen::Index>(dim));
    }

    template <class BBOX>
    bool kdtree_get_bbox(BBOX&) const
    {
        return false;
    }
};

using FaceKdTree = nanoflann::KDTreeSingleIndexAdaptor<
    nanoflann::L2_Simple_Adaptor<double, FaceCentroidCloud>,
    FaceCentroidCloud,
    3,
    std::size_t
>;

static double signed_mesh_volume(const Eigen::MatrixXd& vertices, const Eigen::MatrixXi& faces)
{
    if (vertices.rows() == 0 || faces.rows() == 0) return 0.0;

    const Eigen::Vector3d centroid = vertices.colwise().mean();
    double total = 0.0;
    for (Eigen::Index i = 0; i < faces.rows(); ++i) {
        const Eigen::Vector3d v0 = vertices.row(faces(i, 0)).transpose() - centroid;
        const Eigen::Vector3d v1 = vertices.row(faces(i, 1)).transpose() - centroid;
        const Eigen::Vector3d v2 = vertices.row(faces(i, 2)).transpose() - centroid;
        total += v0.dot(v1.cross(v2));
    }
    return total / 6.0;
}

static void orient_faces_outward(const Eigen::MatrixXd& vertices, Eigen::MatrixXi& faces)
{
    Eigen::MatrixXi oriented_faces;
    Eigen::VectorXi components;
    igl::bfs_orient(faces, oriented_faces, components);
    faces = std::move(oriented_faces);

    if (signed_mesh_volume(vertices, faces) < 0.0) {
        faces = faces.rowwise().reverse().eval();
    }
}

static Eigen::Vector3d closest_point_on_triangle(
    const Eigen::Vector3d& p,
    const Eigen::Vector3d& a,
    const Eigen::Vector3d& b,
    const Eigen::Vector3d& c
)
{
    const Eigen::Vector3d ab = b - a;
    const Eigen::Vector3d ac = c - a;
    const Eigen::Vector3d ap = p - a;
    const double d1 = ab.dot(ap);
    const double d2 = ac.dot(ap);
    if (d1 <= 0.0 && d2 <= 0.0) return a;

    const Eigen::Vector3d bp = p - b;
    const double d3 = ab.dot(bp);
    const double d4 = ac.dot(bp);
    if (d3 >= 0.0 && d4 <= d3) return b;

    const double vc = d1 * d4 - d3 * d2;
    if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) {
        const double v = d1 / (d1 - d3);
        return a + v * ab;
    }

    const Eigen::Vector3d cp = p - c;
    const double d5 = ab.dot(cp);
    const double d6 = ac.dot(cp);
    if (d6 >= 0.0 && d5 <= d6) return c;

    const double vb = d5 * d2 - d1 * d6;
    if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) {
        const double w = d2 / (d2 - d6);
        return a + w * ac;
    }

    const double va = d3 * d6 - d5 * d4;
    if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0) {
        const double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return b + w * (c - b);
    }

    const double denom = 1.0 / (va + vb + vc);
    const double v = vb * denom;
    const double w = vc * denom;
    return a + ab * v + ac * w;
}

static std::vector<FaceRecord> build_face_records(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& faces
)
{
    std::vector<FaceRecord> records;
    records.reserve(static_cast<std::size_t>(faces.rows()));
    for (Eigen::Index i = 0; i < faces.rows(); ++i) {
        FaceRecord face;
        face.a = vertices.row(faces(i, 0)).transpose();
        face.b = vertices.row(faces(i, 1)).transpose();
        face.c = vertices.row(faces(i, 2)).transpose();
        face.centroid = (face.a + face.b + face.c) / 3.0;
        face.normal = (face.b - face.a).cross(face.c - face.a);
        const double norm = face.normal.norm();
        if (norm <= 0.0) continue;
        face.normal /= norm;
        records.push_back(face);
    }
    return records;
}

static bool classify_inside_by_nearest_face(
    const Eigen::Vector3d& query,
    const std::vector<FaceRecord>& faces,
    const FaceKdTree& tree,
    std::size_t candidate_count,
    std::vector<std::size_t>& indices,
    std::vector<double>& centroid_distances
)
{
    candidate_count = std::min(candidate_count, faces.size());
    indices.resize(candidate_count);
    centroid_distances.resize(candidate_count);
    nanoflann::KNNResultSet<double, std::size_t> result_set(candidate_count);
    result_set.init(indices.data(), centroid_distances.data());
    const double q[3] = {query.x(), query.y(), query.z()};
    tree.findNeighbors(result_set, q, nanoflann::SearchParameters());

    double best_distance2 = std::numeric_limits<double>::infinity();
    const FaceRecord* best_face = nullptr;
    Eigen::Vector3d best_closest = Eigen::Vector3d::Zero();

    for (std::size_t idx : indices) {
        if (idx >= faces.size()) continue;
        const FaceRecord& face = faces[idx];
        const Eigen::Vector3d closest = closest_point_on_triangle(query, face.a, face.b, face.c);
        const double distance2 = (query - closest).squaredNorm();
        if (distance2 < best_distance2) {
            best_distance2 = distance2;
            best_face = &face;
            best_closest = closest;
        }
    }

    if (!best_face) return false;
    return (query - best_closest).dot(best_face->normal) < 0.0;
}

static bool point_outside_bbox(
    const Eigen::Vector3d& q,
    double min_x,
    double min_y,
    double min_z,
    double max_x,
    double max_y,
    double max_z
)
{
    return q.x() < min_x || q.x() > max_x ||
           q.y() < min_y || q.y() > max_y ||
           q.z() < min_z || q.z() > max_z;
}

static bool load_reference_mesh_with_fallback(
    const char* reference_mesh_obj,
    Eigen::MatrixXd& vertices,
    Eigen::MatrixXi& faces,
    std::string& used_path
)
{
    if (!reference_mesh_obj || reference_mesh_obj[0] == '\0') {
        return false;
    }

    std::vector<std::string> candidates;
    const auto add_candidate = [&](const std::filesystem::path& p) {
        const std::string s = p.string();
        if (s.empty()) return;
        if (std::find(candidates.begin(), candidates.end(), s) == candidates.end()) {
            candidates.push_back(s);
        }
    };

    const std::filesystem::path ref_path(reference_mesh_obj);
    add_candidate(ref_path);
    if (!ref_path.is_absolute()) {
        add_candidate(std::filesystem::path("build") / ref_path);
        add_candidate(std::filesystem::path("../build") / ref_path);
    }

    for (const auto& candidate : candidates) {
        Eigen::MatrixXd candidate_vertices;
        Eigen::MatrixXi candidate_faces;
        if (!igl::readOBJ(candidate, candidate_vertices, candidate_faces)) {
            continue;
        }
        if (candidate_vertices.rows() == 0 || candidate_faces.rows() == 0 || candidate_faces.cols() != 3) {
            continue;
        }

        vertices = std::move(candidate_vertices);
        faces = std::move(candidate_faces);
        used_path = candidate;
        return true;
    }

    return false;
}

static bool parse_obj_vertex_index(const std::string& token, int& index)
{
    const std::size_t slash = token.find('/');
    const std::string index_token = token.substr(0, slash);
    if (index_token.empty()) return false;
    try {
        index = std::stoi(index_token);
    } catch (...) {
        return false;
    }
    return index != 0;
}

static bool load_obj_as_triangles(
    const std::string& filename,
    Eigen::MatrixXd& vertices,
    Eigen::MatrixXi& faces
)
{
    std::ifstream input(filename);
    if (!input.is_open()) {
        return false;
    }

    std::vector<Eigen::Vector3d> vertex_list;
    std::vector<Eigen::Vector3i> face_list;
    std::string line;
    while (std::getline(input, line)) {
        std::string trimmed = line;
        const std::size_t first = trimmed.find_first_not_of(" \t\r\n");
        if (first == std::string::npos) continue;
        const std::size_t last = trimmed.find_last_not_of(" \t\r\n");
        trimmed = trimmed.substr(first, last - first + 1);
        if (trimmed.empty() || trimmed[0] == '#') continue;

        std::stringstream ss(trimmed);
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
            std::string token;
            while (ss >> token) {
                int raw_index = 0;
                if (!parse_obj_vertex_index(token, raw_index)) {
                    return false;
                }
                int zero_based = raw_index > 0
                    ? raw_index - 1
                    : static_cast<int>(vertex_list.size()) + raw_index;
                if (zero_based < 0 || zero_based >= static_cast<int>(vertex_list.size())) {
                    return false;
                }
                polygon.push_back(zero_based);
            }
            if (polygon.size() < 3) continue;
            for (std::size_t i = 1; i + 1 < polygon.size(); ++i) {
                face_list.emplace_back(polygon[0], polygon[i], polygon[i + 1]);
            }
        }
    }

    if (vertex_list.empty() || face_list.empty()) {
        return false;
    }

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

    void include_vertices(const Eigen::MatrixXd& vertices) noexcept
    {
        for (Eigen::Index i = 0; i < vertices.rows(); ++i) {
            include(vertices(i, 0), vertices(i, 1), vertices(i, 2));
        }
    }

    bool valid() const noexcept
    {
        return min_x <= max_x && min_y <= max_y && min_z <= max_z;
    }

    double center_x() const noexcept { return 0.5 * (min_x + max_x); }
    double center_y() const noexcept { return 0.5 * (min_y + max_y); }
    double center_z() const noexcept { return 0.5 * (min_z + max_z); }

    double max_span() const noexcept
    {
        return std::max({max_x - min_x, max_y - min_y, max_z - min_z, 1e-12});
    }
};

struct SimilarityTransform3d {
    double cx = 0.0;
    double cy = 0.0;
    double cz = 0.0;
    double inv_scale = 1.0;
};

static SimilarityTransform3d make_unit_bbox_transform(const Eigen::MatrixXd& vertices)
{
    BBox3d bbox;
    bbox.include_vertices(vertices);
    SimilarityTransform3d transform;
    if (!bbox.valid()) {
        return transform;
    }
    transform.cx = bbox.center_x();
    transform.cy = bbox.center_y();
    transform.cz = bbox.center_z();
    transform.inv_scale = 1.0 / bbox.max_span();
    return transform;
}

static void apply_transform_to_vertices(Eigen::MatrixXd& vertices, const SimilarityTransform3d& transform)
{
    for (Eigen::Index i = 0; i < vertices.rows(); ++i) {
        vertices(i, 0) = (vertices(i, 0) - transform.cx) * transform.inv_scale;
        vertices(i, 1) = (vertices(i, 1) - transform.cy) * transform.inv_scale;
        vertices(i, 2) = (vertices(i, 2) - transform.cz) * transform.inv_scale;
    }
}

static void translate_vertices(Eigen::MatrixXd& vertices, const Eigen::Vector3d& offset)
{
    for (Eigen::Index i = 0; i < vertices.rows(); ++i) {
        vertices(i, 0) += offset.x();
        vertices(i, 1) += offset.y();
        vertices(i, 2) += offset.z();
    }
}

static void apply_transform_to_seeds(std::vector<SeedPoint>& points, const SimilarityTransform3d& transform)
{
    for (auto& p : points) {
        p.x = (p.x - transform.cx) * transform.inv_scale;
        p.y = (p.y - transform.cy) * transform.inv_scale;
        p.z = (p.z - transform.cz) * transform.inv_scale;
        p.r *= transform.inv_scale;
    }
}

#ifdef USE_EMBREE
struct EmbreeVertex {
    float x;
    float y;
    float z;
};

struct EmbreeTriangle {
    unsigned int v0;
    unsigned int v1;
    unsigned int v2;
};

class EmbreeParityTester {
public:
    EmbreeParityTester(const Eigen::MatrixXd& vertices, const Eigen::MatrixXi& faces)
    {
        device_ = rtcNewDevice(nullptr);
        if (!device_) {
            throw std::runtime_error("rtcNewDevice failed");
        }

        scene_ = rtcNewScene(device_);
        if (!scene_) {
            rtcReleaseDevice(device_);
            device_ = nullptr;
            throw std::runtime_error("rtcNewScene failed");
        }

        RTCGeometry geometry = rtcNewGeometry(device_, RTC_GEOMETRY_TYPE_TRIANGLE);
        auto* embree_vertices = static_cast<EmbreeVertex*>(
            rtcSetNewGeometryBuffer(
                geometry,
                RTC_BUFFER_TYPE_VERTEX,
                0,
                RTC_FORMAT_FLOAT3,
                sizeof(EmbreeVertex),
                static_cast<std::size_t>(vertices.rows())
            )
        );
        auto* embree_triangles = static_cast<EmbreeTriangle*>(
            rtcSetNewGeometryBuffer(
                geometry,
                RTC_BUFFER_TYPE_INDEX,
                0,
                RTC_FORMAT_UINT3,
                sizeof(EmbreeTriangle),
                static_cast<std::size_t>(faces.rows())
            )
        );
        if (!embree_vertices || !embree_triangles) {
            rtcReleaseGeometry(geometry);
            throw std::runtime_error("rtcSetNewGeometryBuffer failed");
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

        rtcCommitGeometry(geometry);
        rtcAttachGeometry(scene_, geometry);
        rtcReleaseGeometry(geometry);
        rtcCommitScene(scene_);
    }

    EmbreeParityTester(const EmbreeParityTester&) = delete;
    EmbreeParityTester& operator=(const EmbreeParityTester&) = delete;

    ~EmbreeParityTester()
    {
        if (scene_) {
            rtcReleaseScene(scene_);
        }
        if (device_) {
            rtcReleaseDevice(device_);
        }
    }

    bool inside(double x, double y, double z) const
    {
        constexpr float direction_x = 0.86054051f;
        constexpr float direction_y = 0.43027025f;
        constexpr float direction_z = 0.27017059f;
        unsigned int hit_count = 0;
        float tnear = 1e-6f;

        RTCIntersectContext context;
        rtcInitIntersectContext(&context);

        while (true) {
            RTCRayHit rayhit{};
            rayhit.ray.org_x = static_cast<float>(x);
            rayhit.ray.org_y = static_cast<float>(y);
            rayhit.ray.org_z = static_cast<float>(z);
            rayhit.ray.tnear = tnear;
            rayhit.ray.dir_x = direction_x;
            rayhit.ray.dir_y = direction_y;
            rayhit.ray.dir_z = direction_z;
            rayhit.ray.time = 0.0f;
            rayhit.ray.tfar = std::numeric_limits<float>::infinity();
            rayhit.ray.mask = 0xFFFFFFFFu;
            rayhit.ray.id = 0;
            rayhit.ray.flags = 0;
            rayhit.hit.geomID = RTC_INVALID_GEOMETRY_ID;
            rayhit.hit.primID = RTC_INVALID_GEOMETRY_ID;

            rtcIntersect1(scene_, &context, &rayhit);
            if (rayhit.hit.geomID == RTC_INVALID_GEOMETRY_ID) {
                break;
            }

            hit_count++;
            tnear = rayhit.ray.tfar + 1e-5f;
            if (hit_count > 1000000u) {
                break;
            }
        }

        return (hit_count % 2u) == 1u;
    }

private:
    RTCDevice device_ = nullptr;
    RTCScene scene_ = nullptr;
};
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

static std::string to_lower_copy(std::string s)
{
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    return s;
}

static bool token_looks_like_integer(const std::string& token)
{
    if (token.empty()) return false;

    size_t start = 0;
    if (token[0] == '+' || token[0] == '-') {
        start = 1;
    }
    if (start >= token.size()) return false;

    for (size_t i = start; i < token.size(); ++i) {
        if (!std::isdigit(static_cast<unsigned char>(token[i]))) {
            return false;
        }
    }

    return true;
}

static bool split_csv_tokens(const std::string& line, std::vector<std::string>& tokens)
{
    std::string s = line;
    trim_in_place(s);
    if (s.empty()) return false;
    if (s[0] == '#') return false;

    // Accept both comma-separated CSV and whitespace/tab-separated files.
    // The expected seed format is: x1coord, x2coord, x3coord, radius.
    if (s.find(',') != std::string::npos) {
        std::stringstream ss(s);
        std::string token;
        while (std::getline(ss, token, ',')) {
            trim_in_place(token);
            if (!token.empty()) tokens.push_back(token);
        }
    } else {
        std::stringstream ss(s);
        std::string token;
        while (ss >> token) {
            trim_in_place(token);
            if (!token.empty()) tokens.push_back(token);
        }
    }

    return !tokens.empty();
}

static SeedCsvLayout infer_seed_csv_layout(const std::vector<std::string>& tokens)
{
    if (tokens.empty()) return SeedCsvLayout::Unknown;

    if (!tokens[0].empty() && std::isalpha(static_cast<unsigned char>(tokens[0][0]))) {
        const std::string header0 = to_lower_copy(tokens[0]);
        if (header0.find("inside") != std::string::npos ||
            header0.find("outside") != std::string::npos) {
            return SeedCsvLayout::SeedPair;
        }
        return SeedCsvLayout::SingleSeed;
    }

    if (tokens.size() >= 15) {
        return SeedCsvLayout::SeedPair;
    }

    if (tokens.size() == 6) {
        if (token_looks_like_integer(tokens[4]) && token_looks_like_integer(tokens[5])) {
            return SeedCsvLayout::SingleSeed;
        }
        return SeedCsvLayout::SeedPair;
    }

    if (tokens.size() >= 4) {
        return SeedCsvLayout::SingleSeed;
    }

    return SeedCsvLayout::Unknown;
}

static bool append_points_from_csv_line(
    const std::string& line,
    SeedCsvLayout& layout,
    std::vector<SeedPoint>& out_points,
    bool& saw_explicit_region
)
{
    std::vector<std::string> tokens;
    if (!split_csv_tokens(line, tokens)) {
        return true;
    }

    if (layout == SeedCsvLayout::Unknown) {
        layout = infer_seed_csv_layout(tokens);
    }

    if (!tokens[0].empty() && std::isalpha(static_cast<unsigned char>(tokens[0][0]))) {
        return true;
    }

    try {
        if (layout == SeedCsvLayout::SeedPair) {
            if (tokens.size() < 6) return false;

            // Required seed-pair CSV convention:
            // inside_x, inside_y, inside_z, outside_x, outside_y, outside_z,
            // face_v0_x, face_v0_y, face_v0_z, face_v1_x, face_v1_y, face_v1_z,
            // face_v2_x, face_v2_y, face_v2_z.
            // The face columns are accepted for compatibility/diagnostics, but the
            // nearest-seed classifier uses only the first six seed-coordinate columns.
            SeedPoint inside_point;
            inside_point.x = std::stod(tokens[0]);
            inside_point.y = std::stod(tokens[1]);
            inside_point.z = std::stod(tokens[2]);
            inside_point.region = 1;

            SeedPoint outside_point;
            outside_point.x = std::stod(tokens[3]);
            outside_point.y = std::stod(tokens[4]);
            outside_point.z = std::stod(tokens[5]);
            outside_point.region = 0;

            saw_explicit_region = true;
            out_points.push_back(inside_point);
            out_points.push_back(outside_point);
        } else {
            if (tokens.size() < 4) return false;

            SeedPoint point;
            point.x = std::stod(tokens[0]);
            point.y = std::stod(tokens[1]);
            point.z = std::stod(tokens[2]);
            point.r = std::stod(tokens[3]);
            if (tokens.size() >= 5) {
                point.region = static_cast<size_t>(std::stoull(tokens[4]));
                saw_explicit_region = true;
            } else {
                point.region = 1;
            }
            out_points.push_back(point);
        }
    } catch (...) {
        return false;
    }

    return true;
}

} // namespace

bool Inout::benchmark_seed_nearest_csv(
    const std::string& filename,
    std::size_t num_tests,
    const std::string& csv_output
)
{
    if (num_tests == 0) {
        std::cerr << "Error: num_tests must be greater than 0." << std::endl;
        return false;
    }

    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: cannot open seed csv: " << filename << std::endl;
        return false;
    }

    std::vector<SeedPoint> points;
    SeedCsvLayout layout = SeedCsvLayout::Unknown;
    bool saw_explicit_region = false;
    std::string line;
    while (std::getline(file, line)) {
        if (!append_points_from_csv_line(line, layout, points, saw_explicit_region)) {
            std::cerr << "Error: failed to parse seed csv line: " << line << std::endl;
            return false;
        }
    }

    if (points.empty()) {
        std::cerr << "Error: no valid seed points parsed from " << filename << std::endl;
        return false;
    }

    double min_x = std::numeric_limits<double>::infinity();
    double min_y = std::numeric_limits<double>::infinity();
    double min_z = std::numeric_limits<double>::infinity();
    double max_x = -std::numeric_limits<double>::infinity();
    double max_y = -std::numeric_limits<double>::infinity();
    double max_z = -std::numeric_limits<double>::infinity();
    for (const auto& p : points) {
        min_x = std::min(min_x, p.x);
        min_y = std::min(min_y, p.y);
        min_z = std::min(min_z, p.z);
        max_x = std::max(max_x, p.x);
        max_y = std::max(max_y, p.y);
        max_z = std::max(max_z, p.z);
    }

    const double span_x = max_x - min_x;
    const double span_y = max_y - min_y;
    const double span_z = max_z - min_z;
    const double max_span = std::max({span_x, span_y, span_z, 1e-9});
    const double padding = 0.02 * max_span;
    min_x -= padding;
    min_y -= padding;
    min_z -= padding;
    max_x += padding;
    max_y += padding;
    max_z += padding;

    const auto build_start = std::chrono::steady_clock::now();
    const FastSeedIndex seed_index(points);
    const auto build_end = std::chrono::steady_clock::now();
    const auto build_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(build_end - build_start);

    std::mt19937_64 gen(20260430);
    std::uniform_real_distribution<double> dist_x(min_x, max_x);
    std::uniform_real_distribution<double> dist_y(min_y, max_y);
    std::uniform_real_distribution<double> dist_z(min_z, max_z);

    std::vector<double> query_x(num_tests);
    std::vector<double> query_y(num_tests);
    std::vector<double> query_z(num_tests);
    for (std::size_t i = 0; i < num_tests; ++i) {
        query_x[i] = dist_x(gen);
        query_y[i] = dist_y(gen);
        query_z[i] = dist_z(gen);
    }

    // Keep results/checksum so the compiler cannot legally erase the nearest-seed loop.
    std::vector<std::size_t> nearest_indices(num_tests);
    double dist2_checksum = 0.0;
    std::size_t index_checksum = 1469598103934665603ull;

    const auto query_start = std::chrono::steady_clock::now();
    for (std::size_t i = 0; i < num_tests; ++i) {
        double dist2 = 0.0;
        const std::size_t nearest_idx = seed_index.nearest_index(query_x[i], query_y[i], query_z[i], &dist2);
        if (nearest_idx >= points.size()) {
            std::cerr << "Error: nearest-neighbor query failed at test " << i << std::endl;
            return false;
        }
        nearest_indices[i] = nearest_idx;
        dist2_checksum += dist2;
        index_checksum ^= nearest_idx + 0x9e3779b97f4a7c15ull + (index_checksum << 6) + (index_checksum >> 2);
    }
    const auto query_end = std::chrono::steady_clock::now();
    const auto query_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(query_end - query_start);

    const double build_ms = std::chrono::duration<double, std::milli>(build_ns).count();
    const double query_total_ms = std::chrono::duration<double, std::milli>(query_ns).count();
    const double query_avg_us = query_total_ms * 1000.0 / static_cast<double>(num_tests);

    std::cout << std::fixed << std::setprecision(6)
              << "[SeedNearest Summary]"
              << " seed_csv=" << filename
              << ", seeds=" << points.size()
              << ", total_tests=" << num_tests
              << ", seed_index_build_ms=" << build_ms
              << ", method_query_total_ms=" << query_total_ms
              << ", method_query_avg_us=" << query_avg_us
              << ", index_checksum=" << index_checksum
              << ", dist2_checksum=" << dist2_checksum
              << ", metric=euclidean_nearest_seed_xyz"
              << ", radius_used=false"
              << std::endl;

    if (!csv_output.empty()) {
        const bool exists = std::filesystem::exists(csv_output);
        std::ofstream out(csv_output, std::ios::app);
        if (!out.is_open()) {
            std::cerr << "Error: cannot write csv output: " << csv_output << std::endl;
            return false;
        }
        if (!exists) {
            out << "seed_csv,seeds,total_tests,seed_index_build_ms,method_query_total_ms,method_query_avg_us,index_checksum,dist2_checksum,metric,radius_used\n";
        }
        out << std::fixed << std::setprecision(9)
            << filename << ','
            << points.size() << ','
            << num_tests << ','
            << build_ms << ','
            << query_total_ms << ','
            << query_avg_us << ','
            << index_checksum << ','
            << dist2_checksum << ','
            << "euclidean_nearest_seed_xyz" << ','
            << "false" << '\n';
    }

    (void)saw_explicit_region;
    return true;
}

bool Inout::analyze_seed_csv(
    const std::string& filename,
    const char* cloud_obj,
    const char* query_obj,
    const char* combined_obj,
    const char* nearest_obj,
    const char* reference_mesh_obj,
    std::size_t num_tests
)
{
    if (num_tests == 0) {
        std::cerr << "Error: num_tests must be greater than 0." << std::endl;
        return false;
    }

    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: cannot open seed csv: " << filename << std::endl;
        return false;
    }

    std::vector<SeedPoint> points;
    SeedCsvLayout layout = SeedCsvLayout::Unknown;
    bool saw_explicit_region = false;
    std::string line;
    while (std::getline(file, line)) {
        if (!append_points_from_csv_line(line, layout, points, saw_explicit_region)) {
            std::cerr << "Error: failed to parse seed csv line: " << line << std::endl;
            return false;
        }
    }

    if (points.empty()) {
        std::cerr << "Error: no valid points parsed from " << filename << std::endl;
        return false;
    }

    if (!reference_mesh_obj || reference_mesh_obj[0] == '\0') {
        std::cerr << "Error: reference_mesh_obj is empty." << std::endl;
        return false;
    }

    double min_x = std::numeric_limits<double>::infinity();
    double min_y = std::numeric_limits<double>::infinity();
    double min_z = std::numeric_limits<double>::infinity();
    double max_x = -std::numeric_limits<double>::infinity();
    double max_y = -std::numeric_limits<double>::infinity();
    double max_z = -std::numeric_limits<double>::infinity();

    for (const auto& p : points) {
        min_x = std::min(min_x, p.x);
        min_y = std::min(min_y, p.y);
        min_z = std::min(min_z, p.z);
        max_x = std::max(max_x, p.x);
        max_y = std::max(max_y, p.y);
        max_z = std::max(max_z, p.z);
    }

    const double span_x = max_x - min_x;
    const double span_y = max_y - min_y;
    const double span_z = max_z - min_z;
    const double max_span = std::max({span_x, span_y, span_z, 1e-9});
    const double padding = 0.02 * max_span;
    min_x -= padding;
    min_y -= padding;
    min_z -= padding;
    max_x += padding;
    max_y += padding;
    max_z += padding;

    const auto seed_index_build_start = std::chrono::steady_clock::now();
    const FastSeedIndex seed_index(points);
    const auto seed_index_build_end = std::chrono::steady_clock::now();
    const std::chrono::nanoseconds seed_index_build_ns =
        std::chrono::duration_cast<std::chrono::nanoseconds>(seed_index_build_end - seed_index_build_start);

    std::cout << "[SeedCSV] parsed " << points.size()
              << " seeds from " << filename
              << "; query metric = Euclidean nearest seed, using columns x1coord/x2coord/x3coord"
              << "; radius is loaded but not used by nearest-point distance"
              << std::endl;
    if (!saw_explicit_region) {
        std::cout << "[SeedCSV] no explicit region column found; validation labels default to region=1."
                  << std::endl;
    }

    Eigen::MatrixXd reference_vertices;
    Eigen::MatrixXi reference_faces;
    std::string loaded_mesh_path;
    if (!load_reference_mesh_with_fallback(reference_mesh_obj, reference_vertices, reference_faces, loaded_mesh_path)) {
        std::cerr << "Error: cannot read a valid triangular reference mesh obj from: "
                  << reference_mesh_obj
                  << " (also tried build/ and ../build/ fallbacks)"
                  << std::endl;
        return false;
    }
    std::cout << "[Validation] reference mesh: " << loaded_mesh_path
              << " (faces=" << reference_faces.rows() << ")" << std::endl;

    std::mt19937_64 gen(std::random_device{}());
    std::uniform_real_distribution<double> dist_x(min_x, max_x);
    std::uniform_real_distribution<double> dist_y(min_y, max_y);
    std::uniform_real_distribution<double> dist_z(min_z, max_z);

    Eigen::MatrixXd query_points(num_tests, 3);
    std::vector<double> query_x(num_tests);
    std::vector<double> query_y(num_tests);
    std::vector<double> query_z(num_tests);
    for (std::size_t i = 0; i < num_tests; i++) {
        query_x[i] = dist_x(gen);
        query_y[i] = dist_y(gen);
        query_z[i] = dist_z(gen);
        query_points(static_cast<Eigen::Index>(i), 0) = query_x[i];
        query_points(static_cast<Eigen::Index>(i), 1) = query_y[i];
        query_points(static_cast<Eigen::Index>(i), 2) = query_z[i];
    }

    const auto winding_start = std::chrono::steady_clock::now();
    Eigen::VectorXd winding_numbers;
    igl::fast_winding_number(reference_vertices, reference_faces, query_points, winding_numbers);
    const auto winding_end = std::chrono::steady_clock::now();
    const std::chrono::nanoseconds winding_time_ns =
        std::chrono::duration_cast<std::chrono::nanoseconds>(winding_end - winding_start);

    struct ErrorPair {
        double qx;
        double qy;
        double qz;
        double nx;
        double ny;
        double nz;
    };

    std::vector<ErrorPair> errors;
    std::size_t mismatch_count = 0;
    std::size_t boundary_count = 0;
    std::size_t valid_count = 0;

    std::vector<std::size_t> valid_indices;
    valid_indices.reserve(num_tests);
    std::vector<unsigned char> reference_inside(num_tests, 0);
    for (std::size_t i = 0; i < num_tests; i++) {
        const double winding = winding_numbers(static_cast<Eigen::Index>(i));
        if (std::abs(std::abs(winding) - 0.5) <= 1e-8) {
            boundary_count++;
            continue;
        }
        reference_inside[i] = static_cast<unsigned char>(std::abs(winding) > 0.5);
        valid_indices.push_back(i);
    }
    valid_count = valid_indices.size();

    std::vector<std::size_t> nearest_indices(valid_count);
    const auto method_start = std::chrono::steady_clock::now();
    for (std::size_t k = 0; k < valid_count; ++k) {
        const std::size_t i = valid_indices[k];
        const std::size_t nearest_idx = seed_index.nearest_index(query_x[i], query_y[i], query_z[i]);
        if (nearest_idx >= points.size()) {
            std::cerr << "Error: nearest-neighbor query failed at test " << i << std::endl;
            return false;
        }
        nearest_indices[k] = nearest_idx;
    }
    const auto method_end = std::chrono::steady_clock::now();
    const std::chrono::nanoseconds method_time_ns =
        std::chrono::duration_cast<std::chrono::nanoseconds>(method_end - method_start);

    for (std::size_t k = 0; k < valid_count; ++k) {
        const std::size_t i = valid_indices[k];
        const std::size_t nearest_idx = nearest_indices[k];
        const bool predicted_inside = points[nearest_idx].region != 0;
        const bool gt_inside = reference_inside[i] != 0;
        if (predicted_inside != gt_inside) {
            mismatch_count++;
            const SeedPoint& nearest = points[nearest_idx];
            errors.push_back({query_x[i], query_y[i], query_z[i], nearest.x, nearest.y, nearest.z});
        }
    }

    const double accuracy = (valid_count > 0)
        ? (1.0 - static_cast<double>(mismatch_count) / static_cast<double>(valid_count))
        : 1.0;

    const double winding_total_ms = std::chrono::duration<double, std::milli>(winding_time_ns).count();
    const double seed_index_build_ms = std::chrono::duration<double, std::milli>(seed_index_build_ns).count();
    const double method_total_ms = std::chrono::duration<double, std::milli>(method_time_ns).count();
    const double winding_avg_us = (num_tests > 0)
        ? (winding_total_ms * 1000.0 / static_cast<double>(num_tests))
        : 0.0;
    const double method_avg_us = (valid_count > 0)
        ? (method_total_ms * 1000.0 / static_cast<double>(valid_count))
        : 0.0;

    std::cout << std::fixed << std::setprecision(6)
              << "[Validation Summary] total_tests=" << num_tests
              << ", effective_tests=" << valid_count
              << ", boundary_skipped=" << boundary_count
              << ", mismatches=" << mismatch_count
              << ", accuracy=" << (accuracy * 100.0) << "%"
              << ", winding_total_ms=" << winding_total_ms
              << ", winding_avg_us=" << winding_avg_us
              << ", seed_index_build_ms=" << seed_index_build_ms
              << ", method_total_ms=" << method_total_ms
              << ", method_avg_us=" << method_avg_us
              << std::endl;

    {
        if (cloud_obj && cloud_obj[0] != '\0') {
            std::ofstream out(cloud_obj);
            if (out.is_open()) {
                out << std::fixed << std::setprecision(16);
                for (const auto& p : points) {
                    out << "v " << p.x << " " << p.y << " " << p.z << "\n";
                }
            } else {
                std::cerr << "Warning: cannot write " << cloud_obj << std::endl;
            }
        }
    }

    {
        if (query_obj && query_obj[0] != '\0') {
            std::ofstream out(query_obj);
            if (out.is_open()) {
                out << std::fixed << std::setprecision(16);
                for (const auto& e : errors) {
                    out << "v " << e.qx << " " << e.qy << " " << e.qz << "\n";
                }
            } else {
                std::cerr << "Warning: cannot write " << query_obj << std::endl;
            }
        }
    }

    {
        if (nearest_obj && nearest_obj[0] != '\0') {
            std::ofstream out(nearest_obj);
            if (out.is_open()) {
                out << std::fixed << std::setprecision(16);
                for (const auto& e : errors) {
                    out << "v " << e.nx << " " << e.ny << " " << e.nz << "\n";
                }
            } else {
                std::cerr << "Warning: cannot write " << nearest_obj << std::endl;
            }
        }
    }

    {
        if (combined_obj && combined_obj[0] != '\0') {
            std::ofstream out(combined_obj);
            if (out.is_open()) {
                const std::size_t num_errors = errors.size();
                out << std::fixed << std::setprecision(16);
                for (const auto& e : errors) {
                    out << "v " << e.qx << " " << e.qy << " " << e.qz << " 1 0 0\n";
                }
                for (const auto& e : errors) {
                    out << "v " << e.nx << " " << e.ny << " " << e.nz << " 0 1 0\n";
                }
                for (std::size_t i = 0; i < num_errors; i++) {
                    out << "l " << (i + 1) << " " << (num_errors + i + 1) << "\n";
                }
            } else {
                std::cerr << "Warning: cannot write " << combined_obj << std::endl;
            }
        }
    }

    if (errors.empty()) {
        std::cout << "[Validation] No mismatches found. Error OBJ files are empty." << std::endl;
    } else {
        std::cout << "[Validation] Wrote " << errors.size() << " mismatch pairs to OBJ files." << std::endl;
    }

    return true;
}

bool Inout::benchmark_mesh_inside_outside(
    const std::string& reconstructed_mesh_obj,
    const std::string& reference_mesh_obj,
    const std::string& seed_pair_csv,
    std::size_t num_tests,
    const std::string& csv_output
)
{
#ifndef USE_CGAL
    (void)reconstructed_mesh_obj;
    (void)reference_mesh_obj;
    (void)seed_pair_csv;
    (void)num_tests;
    (void)csv_output;
    std::cerr << "Error: --inout-benchmark requires USE_CGAL=ON." << std::endl;
    return false;
#else
    if (num_tests == 0) {
        std::cerr << "Error: num_tests must be greater than 0." << std::endl;
        return false;
    }

    std::ifstream seed_file(seed_pair_csv);
    if (!seed_file.is_open()) {
        std::cerr << "Error: cannot open seed pair csv: " << seed_pair_csv << std::endl;
        return false;
    }

    std::vector<SeedPoint> seed_points;
    SeedCsvLayout seed_layout = SeedCsvLayout::Unknown;
    bool saw_explicit_region = false;
    std::string seed_line;
    while (std::getline(seed_file, seed_line)) {
        if (!append_points_from_csv_line(seed_line, seed_layout, seed_points, saw_explicit_region)) {
            std::cerr << "Error: failed to parse seed csv line: " << seed_line << std::endl;
            return false;
        }
    }
    if (seed_points.empty()) {
        std::cerr << "Error: no valid seed points parsed from " << seed_pair_csv << std::endl;
        return false;
    }

    Eigen::MatrixXd reconstructed_vertices;
    Eigen::MatrixXi reconstructed_faces;
    if (!load_obj_as_triangles(reconstructed_mesh_obj, reconstructed_vertices, reconstructed_faces) ||
        reconstructed_vertices.rows() == 0 ||
        reconstructed_faces.rows() == 0 ||
        reconstructed_faces.cols() != 3) {
        std::cerr << "Error: cannot read triangular reconstructed mesh: "
                  << reconstructed_mesh_obj << std::endl;
        return false;
    }

    Eigen::MatrixXd reference_vertices;
    Eigen::MatrixXi reference_faces;
    if (!load_obj_as_triangles(reference_mesh_obj, reference_vertices, reference_faces) ||
        reference_vertices.rows() == 0 ||
        reference_faces.rows() == 0 ||
        reference_faces.cols() != 3) {
        std::cerr << "Error: cannot read triangular reference mesh: "
                  << reference_mesh_obj << std::endl;
        return false;
    }

    if (!saw_explicit_region) {
        if ((seed_points.size() % 2) != 0) {
            std::cerr << "Error: seed csv has no explicit inside/outside label and contains an odd number of seeds. "
                      << "For inside/outside accuracy, provide either a seed-pair CSV or a fifth region column."
                      << std::endl;
            return false;
        }
        std::cerr << "Warning: seed csv has no explicit region column; assuming alternating row pairs "
                  << "[inside, outside], matching columns inside_x inside_y inside_z outside_x outside_y outside_z."
                  << std::endl;
        for (std::size_t i = 0; i < seed_points.size(); ++i) {
            seed_points[i].region = (i % 2 == 0) ? 1 : 0;
        }
        saw_explicit_region = true;
    }

    // Move the original data mesh into the same position as the reconstructed
    // dedup mesh before the inside/outside test.  The seed coordinates are
    // produced together with the reconstructed model, so all inputs then receive
    // the reconstructed/dedup transform.
    BBox3d reference_bbox;
    reference_bbox.include_vertices(reference_vertices);
    BBox3d reconstructed_bbox;
    reconstructed_bbox.include_vertices(reconstructed_vertices);
    if (!reference_bbox.valid() || !reconstructed_bbox.valid()) {
        std::cerr << "Error: cannot compute valid mesh bbox for inout alignment." << std::endl;
        return false;
    }
    const Eigen::Vector3d reference_center(
        reference_bbox.center_x(),
        reference_bbox.center_y(),
        reference_bbox.center_z()
    );
    const Eigen::Vector3d reconstructed_center(
        reconstructed_bbox.center_x(),
        reconstructed_bbox.center_y(),
        reconstructed_bbox.center_z()
    );
    const Eigen::Vector3d reference_to_dedup_offset = reconstructed_center - reference_center;
    translate_vertices(reference_vertices, reference_to_dedup_offset);

    const SimilarityTransform3d reconstructed_transform = make_unit_bbox_transform(reconstructed_vertices);
    apply_transform_to_vertices(reference_vertices, reconstructed_transform);
    apply_transform_to_vertices(reconstructed_vertices, reconstructed_transform);
    apply_transform_to_seeds(seed_points, reconstructed_transform);

    std::cout << std::fixed << std::setprecision(9)
              << "[Align] original/reference bbox center=(" << reference_center.x() << ", "
              << reference_center.y() << ", " << reference_center.z() << ")"
              << ", dedup/reconstruction bbox center=(" << reconstructed_center.x() << ", "
              << reconstructed_center.y() << ", " << reconstructed_center.z() << ")"
              << ", original_translation=(" << reference_to_dedup_offset.x() << ", "
              << reference_to_dedup_offset.y() << ", " << reference_to_dedup_offset.z() << ")"
              << std::endl;
    std::cout << std::fixed << std::setprecision(9)
              << "[Normalize] aligned original/reconstruction/seed center=(" << reconstructed_transform.cx << ", "
              << reconstructed_transform.cy << ", " << reconstructed_transform.cz << ")"
              << ", scale=" << reconstructed_transform.inv_scale << std::endl;

    orient_faces_outward(reconstructed_vertices, reconstructed_faces);
    orient_faces_outward(reference_vertices, reference_faces);

    double min_x = std::numeric_limits<double>::infinity();
    double min_y = std::numeric_limits<double>::infinity();
    double min_z = std::numeric_limits<double>::infinity();
    double max_x = -std::numeric_limits<double>::infinity();
    double max_y = -std::numeric_limits<double>::infinity();
    double max_z = -std::numeric_limits<double>::infinity();

    const auto include_vertices_in_bbox = [&](const Eigen::MatrixXd& vertices) {
        for (Eigen::Index i = 0; i < vertices.rows(); ++i) {
            min_x = std::min(min_x, vertices(i, 0));
            min_y = std::min(min_y, vertices(i, 1));
            min_z = std::min(min_z, vertices(i, 2));
            max_x = std::max(max_x, vertices(i, 0));
            max_y = std::max(max_y, vertices(i, 1));
            max_z = std::max(max_z, vertices(i, 2));
        }
    };
    // Generate query points in the reconstruction/seed canonical space.  The
    // reconstructed GT is the consistency target for the seed Voronoi classifier;
    // using the original mesh bbox here can sample a different space and unfairly
    // dilute the reconstruction accuracy.
    include_vertices_in_bbox(reconstructed_vertices);
    for (const auto& p : seed_points) {
        min_x = std::min(min_x, p.x);
        min_y = std::min(min_y, p.y);
        min_z = std::min(min_z, p.z);
        max_x = std::max(max_x, p.x);
        max_y = std::max(max_y, p.y);
        max_z = std::max(max_z, p.z);
    }

    const double max_span = std::max({max_x - min_x, max_y - min_y, max_z - min_z, 1e-9});
    const double padding = 0.03 * max_span;
    min_x -= padding;
    min_y -= padding;
    min_z -= padding;
    max_x += padding;
    max_y += padding;
    max_z += padding;

    std::mt19937_64 gen(20260430);
    std::uniform_real_distribution<double> dist_x(min_x, max_x);
    std::uniform_real_distribution<double> dist_y(min_y, max_y);
    std::uniform_real_distribution<double> dist_z(min_z, max_z);

    Eigen::MatrixXd query_points(num_tests, 3);
    std::vector<double> query_x(num_tests);
    std::vector<double> query_y(num_tests);
    std::vector<double> query_z(num_tests);
    for (std::size_t i = 0; i < num_tests; ++i) {
        const double x = dist_x(gen);
        const double y = dist_y(gen);
        const double z = dist_z(gen);

        query_x[i] = x;
        query_y[i] = y;
        query_z[i] = z;

        query_points(static_cast<Eigen::Index>(i), 0) = x;
        query_points(static_cast<Eigen::Index>(i), 1) = y;
        query_points(static_cast<Eigen::Index>(i), 2) = z;
    }

    std::cout << "[InOut] reconstructed=" << reconstructed_mesh_obj
              << " (V=" << reconstructed_vertices.rows()
              << ", F=" << reconstructed_faces.rows() << ")" << std::endl;
    std::cout << "[InOut] reference=" << reference_mesh_obj
              << " (V=" << reference_vertices.rows()
              << ", F=" << reference_faces.rows() << ")" << std::endl;
    std::cout << "[InOut] seed_pairs=" << seed_pair_csv
              << " (seeds=" << seed_points.size() << ")" << std::endl;
    std::cout << "[InOut] seed-pair convention: inside=(cols 0..2), outside=(cols 3..5); "
              << "face columns 6..14 are ignored by nearest-seed query" << std::endl;
    std::cout << "[InOut] num_tests=" << num_tests << std::endl;

    const auto original_fast_winding_start = std::chrono::steady_clock::now();
    Eigen::VectorXd original_fast_winding_numbers;
    igl::fast_winding_number(reference_vertices, reference_faces, query_points, original_fast_winding_numbers);
    const auto original_fast_winding_end = std::chrono::steady_clock::now();
    const auto original_fast_winding_ns =
        std::chrono::duration_cast<std::chrono::nanoseconds>(
            original_fast_winding_end - original_fast_winding_start
        );

    const auto reconstructed_fast_winding_start = std::chrono::steady_clock::now();
    Eigen::VectorXd reconstructed_fast_winding_numbers;
    igl::fast_winding_number(
        reconstructed_vertices,
        reconstructed_faces,
        query_points,
        reconstructed_fast_winding_numbers
    );
    const auto reconstructed_fast_winding_end = std::chrono::steady_clock::now();
    const auto reconstructed_fast_winding_ns =
        std::chrono::duration_cast<std::chrono::nanoseconds>(
            reconstructed_fast_winding_end - reconstructed_fast_winding_start
        );

    const auto method_build_start = std::chrono::steady_clock::now();
    const FastSeedIndex seed_index(seed_points);
    const auto method_build_end = std::chrono::steady_clock::now();
    const auto method_build_ns =
        std::chrono::duration_cast<std::chrono::nanoseconds>(method_build_end - method_build_start);

    using CgalKernel = CGAL::Simple_cartesian<double>;
    using CgalPoint = CgalKernel::Point_3;
    using CgalMesh = CGAL::Surface_mesh<CgalPoint>;
    using CgalPrimitive = CGAL::AABB_face_graph_triangle_primitive<CgalMesh>;
    using CgalTraits = CGAL::AABB_traits_3<CgalKernel, CgalPrimitive>;
    using CgalTree = CGAL::AABB_tree<CgalTraits>;
    using CgalSide = CGAL::Side_of_triangle_mesh<CgalMesh, CgalKernel, CGAL::Default, CgalTree>;

    auto run_cgal_queries = [&](const Eigen::MatrixXd& vertices,
                                const Eigen::MatrixXi& faces,
                                std::vector<CGAL::Bounded_side>& results,
                                std::size_t& boundary_count,
                                std::size_t& skipped_faces,
                                std::chrono::nanoseconds& build_ns,
                                std::chrono::nanoseconds& query_ns) {
        CgalMesh mesh;
        std::vector<CgalMesh::Vertex_index> mesh_vertices;
        mesh_vertices.reserve(static_cast<std::size_t>(vertices.rows()));
        for (Eigen::Index i = 0; i < vertices.rows(); ++i) {
            mesh_vertices.push_back(mesh.add_vertex(CgalPoint(
                vertices(i, 0),
                vertices(i, 1),
                vertices(i, 2)
            )));
        }
        skipped_faces = 0;
        for (Eigen::Index i = 0; i < faces.rows(); ++i) {
            const CgalMesh::Face_index face = mesh.add_face(
                mesh_vertices[static_cast<std::size_t>(faces(i, 0))],
                mesh_vertices[static_cast<std::size_t>(faces(i, 1))],
                mesh_vertices[static_cast<std::size_t>(faces(i, 2))]
            );
            if (face == CgalMesh::null_face()) {
                skipped_faces++;
            }
        }

        const auto build_start = std::chrono::steady_clock::now();
        CgalTree tree(CGAL::faces(mesh).first, CGAL::faces(mesh).second, mesh);
        tree.build();
        CgalSide side(tree);
        const auto build_end = std::chrono::steady_clock::now();
        build_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(build_end - build_start);

        boundary_count = 0;
        const auto query_start = std::chrono::steady_clock::now();
        for (std::size_t i = 0; i < num_tests; ++i) {
            const Eigen::Index row = static_cast<Eigen::Index>(i);
            const CGAL::Bounded_side result = side(CgalPoint(
                query_points(row, 0), query_points(row, 1), query_points(row, 2)
            ));
            if (result == CGAL::ON_BOUNDARY) {
                boundary_count++;
            }
            results[i] = result;
        }
        const auto query_end = std::chrono::steady_clock::now();
        query_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(query_end - query_start);
    };

    std::vector<CGAL::Bounded_side> original_cgal_results(num_tests, CGAL::ON_UNBOUNDED_SIDE);
    std::vector<CGAL::Bounded_side> cgal_results(num_tests, CGAL::ON_UNBOUNDED_SIDE);
    std::size_t original_cgal_boundary_count = 0;
    std::size_t cgal_boundary_count = 0;
    std::size_t original_cgal_skipped_faces = 0;
    std::size_t cgal_skipped_faces = 0;
    std::chrono::nanoseconds original_cgal_build_ns(0);
    std::chrono::nanoseconds original_cgal_query_ns(0);
    std::chrono::nanoseconds cgal_build_ns(0);
    std::chrono::nanoseconds cgal_query_ns(0);

    run_cgal_queries(
        reference_vertices,
        reference_faces,
        original_cgal_results,
        original_cgal_boundary_count,
        original_cgal_skipped_faces,
        original_cgal_build_ns,
        original_cgal_query_ns
    );
    if (original_cgal_skipped_faces > 0) {
        std::cerr << "Warning: original CGAL skipped " << original_cgal_skipped_faces
                  << " non-manifold or duplicate faces while building Surface_mesh." << std::endl;
    }
    run_cgal_queries(
        reconstructed_vertices,
        reconstructed_faces,
        cgal_results,
        cgal_boundary_count,
        cgal_skipped_faces,
        cgal_build_ns,
        cgal_query_ns
    );
    if (cgal_skipped_faces > 0) {
        std::cerr << "Warning: reconstructed CGAL skipped " << cgal_skipped_faces
                  << " non-manifold or duplicate faces while building Surface_mesh." << std::endl;
    }

    CgalMesh cgal_mesh;
    std::vector<CgalMesh::Vertex_index> cgal_vertices;
    cgal_vertices.reserve(static_cast<std::size_t>(reconstructed_vertices.rows()));
    for (Eigen::Index i = 0; i < reconstructed_vertices.rows(); ++i) {
        cgal_vertices.push_back(cgal_mesh.add_vertex(CgalPoint(
            reconstructed_vertices(i, 0),
            reconstructed_vertices(i, 1),
            reconstructed_vertices(i, 2)
        )));
    }
    std::size_t seed_check_cgal_skipped_faces = 0;
    for (Eigen::Index i = 0; i < reconstructed_faces.rows(); ++i) {
        const CgalMesh::Face_index face = cgal_mesh.add_face(
            cgal_vertices[static_cast<std::size_t>(reconstructed_faces(i, 0))],
            cgal_vertices[static_cast<std::size_t>(reconstructed_faces(i, 1))],
            cgal_vertices[static_cast<std::size_t>(reconstructed_faces(i, 2))]
        );
        if (face == CgalMesh::null_face()) {
            seed_check_cgal_skipped_faces++;
        }
    }
    if (seed_check_cgal_skipped_faces > 0) {
        std::cerr << "Warning: seed-check CGAL skipped " << seed_check_cgal_skipped_faces
                  << " non-manifold or duplicate faces while building Surface_mesh." << std::endl;
    }

    const auto cgal_build_start = std::chrono::steady_clock::now();
    CgalTree cgal_tree(faces(cgal_mesh).first, faces(cgal_mesh).second, cgal_mesh);
    cgal_tree.build();
    CgalSide cgal_side(cgal_tree);
    const auto cgal_build_end = std::chrono::steady_clock::now();
    const auto seed_check_cgal_build_ns =
        std::chrono::duration_cast<std::chrono::nanoseconds>(cgal_build_end - cgal_build_start);
    (void)seed_check_cgal_build_ns;

    {
        // Sanity-check seed labels against reconstruction GT using CGAL.
        // If this is ~50%, it's usually not a KD-tree/query bug; it means the CSV seed
        // "inside/outside" semantics are inconsistent with the reconstructed mesh (often due
        // to normal/winding issues at seed-generation time or swapped columns in the CSV).
        std::size_t seed_effective = 0;
        std::size_t seed_boundary = 0;
        std::size_t seed_label_mismatches = 0;
        std::size_t seed_gt_inside = 0;
        std::size_t seed_gt_outside = 0;
        for (const auto& s : seed_points) {
            const CGAL::Bounded_side side = cgal_side(CgalPoint(s.x, s.y, s.z));
            if (side == CGAL::ON_BOUNDARY) {
                seed_boundary++;
                continue;
            }
            const bool gt_inside = (side == CGAL::ON_BOUNDED_SIDE);
            seed_effective++;
            if (gt_inside) seed_gt_inside++; else seed_gt_outside++;
            const bool label_inside = (s.region != 0);
            if (label_inside != gt_inside) {
                seed_label_mismatches++;
            }
        }
        const double seed_label_accuracy = (seed_effective > 0)
            ? (100.0 * (1.0 - static_cast<double>(seed_label_mismatches) / static_cast<double>(seed_effective)))
            : 100.0;
        std::cout << std::fixed << std::setprecision(6)
                  << "[Seed Label Check]"
                  << " seeds=" << seed_points.size()
                  << ", effective=" << seed_effective
                  << ", boundary=" << seed_boundary
                  << ", gt_inside=" << seed_gt_inside
                  << ", gt_outside=" << seed_gt_outside
                  << ", label_mismatches=" << seed_label_mismatches
                  << ", label_accuracy=" << seed_label_accuracy << "%"
                  << std::endl;
    }

    std::vector<unsigned char> method_inside(num_tests, 0);
    std::vector<std::size_t> method_nearest_seed(num_tests, 0);
    const auto method_query_start = std::chrono::steady_clock::now();
    for (std::size_t i = 0; i < num_tests; ++i) {
        const double qx = query_x[i];
        const double qy = query_y[i];
        const double qz = query_z[i];
        unsigned char nearest_region = 0;
        const std::size_t nearest_seed =
            seed_index.nearest_index_and_region(qx, qy, qz, nearest_region);
        if (nearest_seed >= seed_points.size()) {
            std::cerr << "Error: nearest-seed query failed at test " << i << std::endl;
            return false;
        }
        method_nearest_seed[i] = nearest_seed;
        method_inside[i] = nearest_region;
    }
    const auto method_query_end = std::chrono::steady_clock::now();
    const auto method_query_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(
        method_query_end - method_query_start
    );

    std::vector<unsigned char> original_embree_inside(num_tests, 0);
    std::vector<unsigned char> embree_inside(num_tests, 0);
    std::chrono::nanoseconds original_embree_build_ns(0);
    std::chrono::nanoseconds original_embree_query_ns(0);
    std::chrono::nanoseconds embree_build_ns(0);
    std::chrono::nanoseconds embree_query_ns(0);
    bool original_embree_available = false;
    bool embree_available = false;
#ifdef USE_EMBREE
    const auto run_embree_queries = [&](const Eigen::MatrixXd& vertices,
                                        const Eigen::MatrixXi& faces,
                                        std::vector<unsigned char>& inside,
                                        std::chrono::nanoseconds& build_ns,
                                        std::chrono::nanoseconds& query_ns,
                                        bool& available,
                                        const char* label) {
        try {
            const auto embree_build_start = std::chrono::steady_clock::now();
            EmbreeParityTester embree_tester(vertices, faces);
            const auto embree_build_end = std::chrono::steady_clock::now();
            build_ns =
                std::chrono::duration_cast<std::chrono::nanoseconds>(embree_build_end - embree_build_start);
            available = true;

            const auto embree_query_start = std::chrono::steady_clock::now();
            for (std::size_t i = 0; i < num_tests; ++i) {
                const Eigen::Index row = static_cast<Eigen::Index>(i);
                inside[i] = static_cast<unsigned char>(embree_tester.inside(
                    query_points(row, 0),
                    query_points(row, 1),
                    query_points(row, 2)
                ));
            }
            const auto embree_query_end = std::chrono::steady_clock::now();
            query_ns =
                std::chrono::duration_cast<std::chrono::nanoseconds>(embree_query_end - embree_query_start);
        } catch (const std::exception& exc) {
            std::cerr << "Warning: " << label
                      << " Embree inside/outside test unavailable: " << exc.what() << std::endl;
            available = false;
        }
    };

    run_embree_queries(
        reference_vertices,
        reference_faces,
        original_embree_inside,
        original_embree_build_ns,
        original_embree_query_ns,
        original_embree_available,
        "original"
    );
    run_embree_queries(
        reconstructed_vertices,
        reconstructed_faces,
        embree_inside,
        embree_build_ns,
        embree_query_ns,
        embree_available,
        "reconstructed"
    );
#else
    std::cerr << "Warning: Embree inside/outside test unavailable because USE_EMBREE is OFF." << std::endl;
#endif

    struct MethodStats {
        std::size_t effective = 0;
        std::size_t gt_boundary = 0;
        std::size_t mismatches = 0;
        double accuracy = 100.0;
    };

    struct MethodMismatch {
        double qx = 0.0;
        double qy = 0.0;
        double qz = 0.0;
        double sx = 0.0;
        double sy = 0.0;
        double sz = 0.0;
        std::size_t seed_index = 0;
        bool predicted_inside = false;
        bool seed_inside = false;
        bool gt_inside = false;
    };

    const auto write_method_mismatches_obj =
        [&](const std::filesystem::path& out_path, const std::vector<MethodMismatch>& mismatches) {
            std::error_code ec;
            std::filesystem::create_directories(out_path.parent_path(), ec);
            std::ofstream out(out_path);
            if (!out.is_open()) {
                std::cerr << "Warning: cannot write mismatch obj: " << out_path << std::endl;
                return;
            }
            out << std::fixed << std::setprecision(16);
            out << "# method mismatches (query point -> nearest seed)\n";
            out << "# mismatches=" << mismatches.size() << "\n";
            // Two vertices per mismatch: query then seed. We also output a line segment.
            for (std::size_t i = 0; i < mismatches.size(); ++i) {
                const MethodMismatch& m = mismatches[i];
                out << "# i=" << i
                    << " seed_index=" << m.seed_index
                    << " predicted_inside=" << (m.predicted_inside ? 1 : 0)
                    << " seed_inside=" << (m.seed_inside ? 1 : 0)
                    << " gt_inside=" << (m.gt_inside ? 1 : 0)
                    << "\n";
                out << "v " << m.qx << " " << m.qy << " " << m.qz << "\n";
                out << "v " << m.sx << " " << m.sy << " " << m.sz << "\n";
            }
            for (std::size_t i = 0; i < mismatches.size(); ++i) {
                const std::size_t v0 = 2 * i + 1;
                const std::size_t v1 = 2 * i + 2;
                out << "l " << v0 << " " << v1 << "\n";
            }
        };

    const auto accuracy_percent = [](std::size_t mismatches, std::size_t count) {
        if (count == 0) return 100.0;
        return 100.0 * (1.0 - static_cast<double>(mismatches) / static_cast<double>(count));
    };

    const auto compare_binary_to_cgal_gt = [&](const std::vector<unsigned char>& predicted_inside,
                                               const std::vector<CGAL::Bounded_side>& gt_results,
                                               std::vector<MethodMismatch>* method_mismatches_out) {
        MethodStats stats;
        for (std::size_t i = 0; i < num_tests; ++i) {
            if (gt_results[i] == CGAL::ON_BOUNDARY) {
                stats.gt_boundary++;
                continue;
            }

            const bool gt_inside = (gt_results[i] == CGAL::ON_BOUNDED_SIDE);
            stats.effective++;
            const bool predicted = predicted_inside[i] != 0;
            if (predicted != gt_inside) {
                stats.mismatches++;
                if (method_mismatches_out) {
                    const std::size_t seed_idx = method_nearest_seed[i];
                    const SeedPoint& seed = seed_points[seed_idx];
                    const Eigen::Index row = static_cast<Eigen::Index>(i);
                    method_mismatches_out->push_back(MethodMismatch{
                        query_points(row, 0), query_points(row, 1), query_points(row, 2),
                        seed.x, seed.y, seed.z,
                        seed_idx,
                        predicted,
                        (seed.region != 0),
                        gt_inside
                    });
                }
            }
        }
        stats.accuracy = accuracy_percent(stats.mismatches, stats.effective);
        return stats;
    };

    const auto compare_winding_to_cgal_gt = [&](const Eigen::VectorXd& winding_numbers,
                                                const std::vector<CGAL::Bounded_side>& gt_results) {
        std::vector<unsigned char> predicted_inside(num_tests, 0);
        for (std::size_t i = 0; i < num_tests; ++i) {
            predicted_inside[i] =
                static_cast<unsigned char>(std::abs(winding_numbers(static_cast<Eigen::Index>(i))) > 0.5);
        }
        return compare_binary_to_cgal_gt(predicted_inside, gt_results, nullptr);
    };

    const MethodStats original_ours_stats =
        compare_binary_to_cgal_gt(method_inside, original_cgal_results, nullptr);
    const MethodStats reconstructed_ours_stats =
        compare_binary_to_cgal_gt(method_inside, cgal_results, nullptr);
    const MethodStats original_fast_winding_stats =
        compare_winding_to_cgal_gt(original_fast_winding_numbers, original_cgal_results);
    const MethodStats reconstructed_fast_winding_stats =
        compare_winding_to_cgal_gt(reconstructed_fast_winding_numbers, cgal_results);
    const MethodStats original_embree_stats = original_embree_available
        ? compare_binary_to_cgal_gt(original_embree_inside, original_cgal_results, nullptr)
        : MethodStats{};
    const MethodStats embree_stats = embree_available
        ? compare_binary_to_cgal_gt(embree_inside, cgal_results, nullptr)
        : MethodStats{};

    std::vector<MethodMismatch> reconstructed_method_mismatches;
    reconstructed_method_mismatches.reserve(1024);
    (void)compare_binary_to_cgal_gt(method_inside, cgal_results, &reconstructed_method_mismatches);

    if (!reconstructed_method_mismatches.empty()) {
        const std::filesystem::path reconstructed_path(reconstructed_mesh_obj);
        const std::string stem = reconstructed_path.stem().string();
        const std::filesystem::path out_obj =
            std::filesystem::path("inout") / "results" / (stem + "_reconstructed_gt_method_mismatches.obj");
        write_method_mismatches_obj(out_obj, reconstructed_method_mismatches);
        std::cout << "[Reconstructed GT] wrote method mismatch OBJ: " << out_obj
                  << " (count=" << reconstructed_method_mismatches.size() << ")" << std::endl;
    } else {
        std::cout << "[Reconstructed GT] no method mismatches; no mismatch OBJ written." << std::endl;
    }

    const auto ns_to_ms = [](std::chrono::nanoseconds ns) {
        return std::chrono::duration<double, std::milli>(ns).count();
    };
    const auto ns_to_avg_us = [](std::chrono::nanoseconds ns, std::size_t count) {
        if (count == 0) return 0.0;
        return std::chrono::duration<double, std::micro>(ns).count() / static_cast<double>(count);
    };

    const double original_fast_winding_total_ms = ns_to_ms(original_fast_winding_ns);
    const double original_fast_winding_avg_us = ns_to_avg_us(original_fast_winding_ns, num_tests);
    const double reconstructed_fast_winding_total_ms = ns_to_ms(reconstructed_fast_winding_ns);
    const double reconstructed_fast_winding_avg_us =
        ns_to_avg_us(reconstructed_fast_winding_ns, num_tests);
    const double reconstructed_ours_build_ms = ns_to_ms(method_build_ns);
    const double reconstructed_ours_query_total_ms = ns_to_ms(method_query_ns);
    const double reconstructed_ours_query_avg_us = ns_to_avg_us(method_query_ns, num_tests);
    const double original_cgal_build_ms = ns_to_ms(original_cgal_build_ns);
    const double original_cgal_query_total_ms = ns_to_ms(original_cgal_query_ns);
    const double original_cgal_query_avg_us = ns_to_avg_us(original_cgal_query_ns, num_tests);
    const double reconstructed_cgal_build_ms = ns_to_ms(cgal_build_ns);
    const double reconstructed_cgal_query_total_ms = ns_to_ms(cgal_query_ns);
    const double reconstructed_cgal_query_avg_us = ns_to_avg_us(cgal_query_ns, num_tests);
    const double original_embree_build_ms = ns_to_ms(original_embree_build_ns);
    const double original_embree_query_total_ms = ns_to_ms(original_embree_query_ns);
    const double original_embree_query_avg_us = ns_to_avg_us(original_embree_query_ns, num_tests);
    const double reconstructed_embree_build_ms = ns_to_ms(embree_build_ns);
    const double reconstructed_embree_query_total_ms = ns_to_ms(embree_query_ns);
    const double reconstructed_embree_query_avg_us = ns_to_avg_us(embree_query_ns, num_tests);

    std::cout << std::fixed << std::setprecision(6)
              << "[InOut Summary]"
              << " total_tests=" << num_tests
              << ", gt=cgal_reconstructed_mesh"
              << ", original_fast_winding_total_ms=" << original_fast_winding_total_ms
              << ", original_fast_winding_avg_us=" << original_fast_winding_avg_us
              << ", reconstructed_fast_winding_total_ms=" << reconstructed_fast_winding_total_ms
              << ", reconstructed_fast_winding_avg_us=" << reconstructed_fast_winding_avg_us
              << ", reconstructed_ours_seed_tree_build_ms=" << reconstructed_ours_build_ms
              << ", reconstructed_ours_query_total_ms=" << reconstructed_ours_query_total_ms
              << ", reconstructed_ours_query_avg_us=" << reconstructed_ours_query_avg_us
              << ", original_cgal_build_ms=" << original_cgal_build_ms
              << ", original_cgal_query_total_ms=" << original_cgal_query_total_ms
              << ", original_cgal_query_avg_us=" << original_cgal_query_avg_us
              << ", reconstructed_cgal_build_ms=" << reconstructed_cgal_build_ms
              << ", reconstructed_cgal_query_total_ms=" << reconstructed_cgal_query_total_ms
              << ", reconstructed_cgal_query_avg_us=" << reconstructed_cgal_query_avg_us
              << ", original_cgal_boundary_queries=" << original_cgal_boundary_count
              << ", reconstructed_cgal_boundary_queries=" << cgal_boundary_count
              << ", original_embree_available=" << (original_embree_available ? 1 : 0)
              << ", original_embree_build_ms=" << original_embree_build_ms
              << ", original_embree_query_total_ms=" << original_embree_query_total_ms
              << ", original_embree_query_avg_us=" << original_embree_query_avg_us
              << ", reconstructed_embree_available=" << (embree_available ? 1 : 0)
              << ", reconstructed_embree_build_ms=" << reconstructed_embree_build_ms
              << ", reconstructed_embree_query_total_ms=" << reconstructed_embree_query_total_ms
              << ", reconstructed_embree_query_avg_us=" << reconstructed_embree_query_avg_us
              << std::endl;

    std::cout << std::fixed << std::setprecision(6)
              << "[CGAL GT Accuracy]"
              << " original_effective_tests=" << original_ours_stats.effective
              << ", original_gt_boundary_skipped=" << original_ours_stats.gt_boundary
              << ", original_ours_mismatches=" << original_ours_stats.mismatches
              << ", original_ours_accuracy=" << original_ours_stats.accuracy << "%"
              << ", reconstructed_effective_tests=" << reconstructed_ours_stats.effective
              << ", reconstructed_gt_boundary_skipped=" << reconstructed_ours_stats.gt_boundary
              << ", reconstructed_ours_mismatches=" << reconstructed_ours_stats.mismatches
              << ", reconstructed_ours_accuracy=" << reconstructed_ours_stats.accuracy << "%"
              << ", original_fast_winding_effective_tests=" << original_fast_winding_stats.effective
              << ", original_fast_winding_boundary_skipped=" << original_fast_winding_stats.gt_boundary
              << ", original_fast_winding_mismatches=" << original_fast_winding_stats.mismatches
              << ", original_fast_winding_accuracy=" << original_fast_winding_stats.accuracy << "%"
              << ", reconstructed_fast_winding_effective_tests=" << reconstructed_fast_winding_stats.effective
              << ", reconstructed_fast_winding_boundary_skipped=" << reconstructed_fast_winding_stats.gt_boundary
              << ", reconstructed_fast_winding_mismatches=" << reconstructed_fast_winding_stats.mismatches
              << ", reconstructed_fast_winding_accuracy=" << reconstructed_fast_winding_stats.accuracy << "%"
              << ", original_embree_mismatches=" << original_embree_stats.mismatches
              << ", original_embree_accuracy=" << original_embree_stats.accuracy << "%"
              << ", embree_mismatches=" << embree_stats.mismatches
              << ", embree_accuracy=" << embree_stats.accuracy << "%"
              << std::endl;

    if (!csv_output.empty()) {
        const bool exists = std::filesystem::exists(csv_output);
        std::ofstream out(csv_output, std::ios::app);
        if (!out.is_open()) {
            std::cerr << "Error: cannot write csv output: " << csv_output << std::endl;
            return false;
        }
        if (!exists) {
            out << "reconstructed,reference,seed_pair_csv,seeds,total_tests,"
                << "gt_method,"
                << "original_fast_winding_total_ms,original_fast_winding_avg_us,"
                << "reconstructed_fast_winding_total_ms,reconstructed_fast_winding_avg_us,"
                << "reconstructed_ours_seed_tree_build_ms,reconstructed_ours_query_total_ms,reconstructed_ours_query_avg_us,"
                << "original_cgal_build_ms,original_cgal_query_total_ms,original_cgal_query_avg_us,original_cgal_boundary_queries,"
                << "reconstructed_cgal_build_ms,reconstructed_cgal_query_total_ms,reconstructed_cgal_query_avg_us,reconstructed_cgal_boundary_queries,"
                << "original_embree_available,original_embree_build_ms,original_embree_query_total_ms,original_embree_query_avg_us,"
                << "reconstructed_embree_available,reconstructed_embree_build_ms,reconstructed_embree_query_total_ms,reconstructed_embree_query_avg_us,"
                << "original_effective_tests,original_gt_boundary_skipped,"
                << "original_ours_mismatches,original_ours_accuracy_percent,"
                << "reconstructed_effective_tests,reconstructed_gt_boundary_skipped,"
                << "reconstructed_ours_mismatches,reconstructed_ours_accuracy_percent,"
                << "original_fast_winding_effective_tests,original_fast_winding_boundary_skipped,"
                << "original_fast_winding_mismatches,original_fast_winding_accuracy_percent,"
                << "reconstructed_fast_winding_effective_tests,reconstructed_fast_winding_boundary_skipped,"
                << "reconstructed_fast_winding_mismatches,reconstructed_fast_winding_accuracy_percent,"
                << "original_embree_mismatches,original_embree_accuracy_percent,"
                << "embree_mismatches,embree_accuracy_percent\n";
        }
        out << std::fixed << std::setprecision(9)
            << reconstructed_mesh_obj << ','
            << reference_mesh_obj << ','
            << seed_pair_csv << ','
            << seed_points.size() << ','
            << num_tests << ','
            << "cgal_reconstructed_mesh" << ','
            << original_fast_winding_total_ms << ','
            << original_fast_winding_avg_us << ','
            << reconstructed_fast_winding_total_ms << ','
            << reconstructed_fast_winding_avg_us << ','
            << reconstructed_ours_build_ms << ','
            << reconstructed_ours_query_total_ms << ','
            << reconstructed_ours_query_avg_us << ','
            << original_cgal_build_ms << ','
            << original_cgal_query_total_ms << ','
            << original_cgal_query_avg_us << ','
            << original_cgal_boundary_count << ','
            << reconstructed_cgal_build_ms << ','
            << reconstructed_cgal_query_total_ms << ','
            << reconstructed_cgal_query_avg_us << ','
            << cgal_boundary_count << ','
            << (original_embree_available ? 1 : 0) << ','
            << original_embree_build_ms << ','
            << original_embree_query_total_ms << ','
            << original_embree_query_avg_us << ','
            << (embree_available ? 1 : 0) << ','
            << reconstructed_embree_build_ms << ','
            << reconstructed_embree_query_total_ms << ','
            << reconstructed_embree_query_avg_us << ','
            << original_ours_stats.effective << ','
            << original_ours_stats.gt_boundary << ','
            << original_ours_stats.mismatches << ','
            << original_ours_stats.accuracy << ','
            << reconstructed_ours_stats.effective << ','
            << reconstructed_ours_stats.gt_boundary << ','
            << reconstructed_ours_stats.mismatches << ','
            << reconstructed_ours_stats.accuracy << ','
            << original_fast_winding_stats.effective << ','
            << original_fast_winding_stats.gt_boundary << ','
            << original_fast_winding_stats.mismatches << ','
            << original_fast_winding_stats.accuracy << ','
            << reconstructed_fast_winding_stats.effective << ','
            << reconstructed_fast_winding_stats.gt_boundary << ','
            << reconstructed_fast_winding_stats.mismatches << ','
            << reconstructed_fast_winding_stats.accuracy << ','
            << original_embree_stats.mismatches << ','
            << original_embree_stats.accuracy << ','
            << embree_stats.mismatches << ','
            << embree_stats.accuracy << '\n';
    }

    (void)saw_explicit_region;
    return true;
#endif
}
