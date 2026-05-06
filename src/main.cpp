#include "Generator.h"
#include "Offset.h"
#include "Inout.h"
#include "MeshingTree.h"
#include <vector>
#include <iostream>
#include <filesystem>
#include <numeric>
#include <Eigen/Geometry>
#include <Eigen/Dense>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/bfs_orient.h> // 用于统一面片方向
#include <igl/fast_winding_number.h>
#include <igl/volume.h>     // 用于计算体积判断朝向

extern void optimizer(MeshingTree* seeds, MeshingTree* spheres, std::vector<int> face_flat);

double **points;
size_t **faces;
double *seeds_sizing, *seedes;
size_t *seeds_region_id;
std::vector<int> face_flat;

bool FixMeshNormals(const std::string& input_path, const std::string& output_path) {
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;

    // 1. 读取 OBJ
    if (!igl::readOBJ(input_path, V, F)) {
        std::cerr << "Error: Cannot read " << input_path << std::endl;
        return false;
    }

    // 2. 统一局部方向 (BFS Orient)
    // 这一步保证了所有面片的“正反”是连贯的，但可能整体朝里，也可能整体朝外
    Eigen::MatrixXi F_consistent;
    Eigen::VectorXi C;
    igl::bfs_orient(F, F_consistent, C);

    // 3. 计算带符号体积 (Signed Volume) 替代 igl::volume
    // 公式：Sum( (p1 x p2) . p3 ) / 6.0
    // 为了数值稳定性，最好先减去质心
    Eigen::Vector3d centroid = V.colwise().mean();
    
    double total_volume = 0.0;
    for (int i = 0; i < F_consistent.rows(); ++i) {
        // 获取三角形的三个顶点，并相对于质心平移
        Eigen::Vector3d v0 = V.row(F_consistent(i, 0)) - centroid.transpose();
        Eigen::Vector3d v1 = V.row(F_consistent(i, 1)) - centroid.transpose();
        Eigen::Vector3d v2 = V.row(F_consistent(i, 2)) - centroid.transpose();

        // 计算该面片与原点构成的四面体体积（带符号）
        total_volume += v0.dot(v1.cross(v2));
    }
    total_volume /= 6.0;

    std::cout << "当前模型体积: " << total_volume << std::endl;

    // 4. 判断方向并翻转
    if (total_volume < 0) {
        std::cout << "检测到法向向内 (Volume < 0)，正在翻转为向外..." << std::endl;
        // 翻转所有面的顶点顺序： (v0, v1, v2) -> (v2, v1, v0)
        F_consistent = F_consistent.rowwise().reverse().eval();
    } else {
        std::cout << "all points to outside" << std::endl;
    }

    // 5. 保存结果
    if (!igl::writeOBJ(output_path, V, F_consistent)) {
        std::cerr << "Error: Cannot write to " << output_path << std::endl;
        return false;
    } else {
        std::cout << "已保存修复后的模型至: " << output_path << std::endl;
    }

    return true;
}

bool ValidateSeedPairsByWinding(const std::string& surface_obj, MeshingTree* seeds) {
    if (!seeds) return false;

    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    if (!igl::readOBJ(surface_obj, V, F) || V.rows() == 0 || F.rows() == 0) {
        std::cerr << "[SeedSideCheck] Cannot read surface mesh: " << surface_obj << std::endl;
        return false;
    }

    const size_t num_seeds = static_cast<size_t>(seeds->get_num_tree_points());
    std::vector<size_t> active_indices;
    active_indices.reserve(num_seeds);
    std::vector<int> row_of_seed(num_seeds, -1);
    for (size_t i = 0; i < num_seeds; ++i) {
        if (!seeds->tree_point_is_active(i)) continue;
        row_of_seed[i] = static_cast<int>(active_indices.size());
        active_indices.push_back(i);
    }
    if (active_indices.empty()) {
        std::cerr << "[SeedSideCheck] No active seeds to validate." << std::endl;
        return true;
    }

    Eigen::MatrixXd Q(static_cast<Eigen::Index>(active_indices.size()), 3);
    for (size_t row = 0; row < active_indices.size(); ++row) {
        double p[4] = {0.0, 0.0, 0.0, 0.0};
        seeds->get_tree_point(active_indices[row], 4, p);
        Q(static_cast<Eigen::Index>(row), 0) = p[0];
        Q(static_cast<Eigen::Index>(row), 1) = p[1];
        Q(static_cast<Eigen::Index>(row), 2) = p[2];
    }

    Eigen::VectorXd winding;
    igl::fast_winding_number(V, F, Q, winding);

    constexpr double boundary_tol = 1e-6;
    std::vector<int> side(num_seeds, -1);
    for (size_t row = 0; row < active_indices.size(); ++row) {
        const double w = winding(static_cast<Eigen::Index>(row));
        if (std::abs(std::abs(w) - 0.5) <= boundary_tol) {
            side[active_indices[row]] = -1;
        } else {
            side[active_indices[row]] = (std::abs(w) > 0.5) ? 1 : 0;
        }
    }

    size_t checked_pairs = 0;
    size_t relabeled_pairs = 0;
    size_t deleted_pairs = 0;
    size_t boundary_pairs = 0;
    for (size_t i = 0; i < num_seeds; ++i) {
        if (!seeds->tree_point_is_active(i)) continue;
        size_t* attrib_i = seeds->get_tree_point_attrib(i);
        if (!attrib_i || attrib_i[0] < 6) continue;
        const size_t j = attrib_i[1];
        if (j <= i || j >= num_seeds || !seeds->tree_point_is_active(j)) continue;

        size_t* attrib_j = seeds->get_tree_point_attrib(j);
        if (!attrib_j || attrib_j[0] < 6 || attrib_j[1] != i) continue;

        checked_pairs++;
        const int side_i = side[i];
        const int side_j = side[j];
        if (side_i < 0 || side_j < 0) {
            seeds->lazy_delete_tree_point(i);
            seeds->lazy_delete_tree_point(j);
            boundary_pairs++;
            deleted_pairs++;
            continue;
        }
        if (side_i == side_j) {
            seeds->lazy_delete_tree_point(i);
            seeds->lazy_delete_tree_point(j);
            deleted_pairs++;
            continue;
        }

        const size_t old_i = attrib_i[5];
        const size_t old_j = attrib_j[5];
        attrib_i[5] = static_cast<size_t>(side_i);
        attrib_j[5] = static_cast<size_t>(side_j);
        if (old_i != attrib_i[5] || old_j != attrib_j[5]) {
            relabeled_pairs++;
        }
    }

    std::cout << "[SeedSideCheck] pairs=" << checked_pairs
              << " relabeled=" << relabeled_pairs
              << " deleted=" << deleted_pairs
              << " boundary=" << boundary_pairs
              << " active_seeds=" << active_indices.size()
              << std::endl;
    return true;
}

int RunSurfaceSpherePipeline(const char* remesh_obj, const char* spheres_csv) {
    Generator generator;
    MeshingTree *spheres = new MeshingTree();
    MeshingTree *upper_seeds = new MeshingTree();
    MeshingTree *lower_seeds = new MeshingTree();
    MeshingTree *seeds = new MeshingTree();

    double *local_seeds_sizing = nullptr;
    double *local_seedes = nullptr;
    size_t *local_seeds_region_id = nullptr;
    std::vector<int> local_face_flat;
    size_t num_faces1 = 0;

    const std::filesystem::path remesh_path(remesh_obj);
    const std::filesystem::path fixed_mesh_path =
        std::filesystem::current_path() / (remesh_path.stem().string() + "_Fixed.obj");
    const std::string fixed_mesh_path_string = fixed_mesh_path.string();

    if (!FixMeshNormals(remesh_path.string(), fixed_mesh_path_string)) {
        delete spheres;
        delete upper_seeds;
        delete lower_seeds;
        delete seeds;
        return 1;
    }

    generator.generate_spheres(spheres_csv, spheres);
    if (spheres->get_num_tree_points() == 0) {
        std::cerr << "[Error] Failed to load spheres or empty file: " << spheres_csv << std::endl;
        delete spheres;
        delete upper_seeds;
        delete lower_seeds;
        delete seeds;
        return 1;
    }

    generator.read_obj_faces(fixed_mesh_path_string.c_str(), local_face_flat, num_faces1);
    if (local_face_flat.empty()) {
        std::cerr << "[Error] Failed to read face indices or empty file: " << fixed_mesh_path_string << std::endl;
        delete spheres;
        delete upper_seeds;
        delete lower_seeds;
        delete seeds;
        return 1;
    }

    generator.generate_surface_seeds(spheres, upper_seeds, lower_seeds, num_faces1, local_face_flat);
    generator.color_surface_seeds(
        static_cast<int>(num_faces1),
        spheres,
        upper_seeds,
        lower_seeds,
        seeds,
        local_face_flat,
        local_seedes,
        local_seeds_region_id,
        local_seeds_sizing
    );
    if (!ValidateSeedPairsByWinding(fixed_mesh_path_string, seeds)) {
        delete spheres;
        delete upper_seeds;
        delete lower_seeds;
        delete seeds;
        return 1;
    }
    optimizer(seeds, spheres, local_face_flat);

    const size_t num_active_seeds = generator.collect_active_seeds(
        seeds,
        local_seedes,
        local_seeds_region_id,
        local_seeds_sizing
    );
    generator.generate_seed_csv("seeds.csv", 3, num_active_seeds, local_seedes, local_seeds_sizing, local_seeds_region_id);
    generator.generate_seed_csv_with_faces("seed_pairs_with_faces.csv", "seed_points_with_faces.csv", seeds, spheres);
    generator.generate_surface_mesh(seeds, spheres, "surface_mesh1.obj");

    delete[] local_seedes;
    delete[] local_seeds_region_id;
    delete[] local_seeds_sizing;
    delete spheres;
    delete upper_seeds;
    delete lower_seeds;
    delete seeds;
    return 0;
}

int main(int argc, char** argv)
{
    Generator generator;
    MeshingTree seed_storage;
    MeshingTree *seeds = &seed_storage;

    if (argc == 4 && std::string(argv[1]) == "--surface-spheres") {
        return RunSurfaceSpherePipeline(argv[2], argv[3]);
    }

    if (argc == 2 && std::string(argv[1]) == "--help") {
        std::cerr << "Usage:" << std::endl;
        std::cerr << "  " << argv[0] << std::endl;
        std::cerr << "  " << argv[0] << " paired_seeds.csv" << std::endl;
        std::cerr << "  " << argv[0] << " inner.txt outer.txt" << std::endl;
        std::cerr << "  " << argv[0] << " --surface-spheres remesh.obj spheres.csv" << std::endl;
        std::cerr << "  " << argv[0] << " --inout-benchmark reconstructed.obj reference.obj seed_pairs.csv [num_tests] [csv_output]" << std::endl;
        return 0;
    }

    if (argc >= 5 && std::string(argv[1]) == "--inout-benchmark") {
        const std::size_t num_tests = (argc >= 6)
            ? static_cast<std::size_t>(std::stoull(argv[5]))
            : 100000;
        const std::string csv_output = (argc >= 7) ? argv[6] : "";
        return Inout::benchmark_mesh_inside_outside(argv[2], argv[3], argv[4], num_tests, csv_output) ? 0 : 1;
    }

    if (argc == 2) {
        const char* paired_csv = argv[1];

        if (!generator.load_paired_seeds_from_csv(paired_csv, seeds)) {
            return 1;
        }

        const size_t num_active_seeds = generator.collect_active_seeds(seeds, seedes, seeds_region_id, seeds_sizing);
        generator.generate_seed_csv("seeds.csv", 3, num_active_seeds, seedes, seeds_sizing, seeds_region_id);
        generator.generate_seed_csv_with_faces("seed_pairs_with_faces.csv", "seed_points_with_faces.csv", seeds, nullptr);
        generator.generate_surface_mesh(seeds, nullptr, "surface_mesh1.obj");
        return 0;
    }

    if (argc == 3) {
        const char* inner_txt = argv[1];
        const char* outer_txt = argv[2];

        if (!generator.load_paired_seeds_from_txt(inner_txt, outer_txt, seeds)) {
            return 1;
        }

        const size_t num_active_seeds = generator.collect_active_seeds(seeds, seedes, seeds_region_id, seeds_sizing);
        generator.generate_seed_csv("seeds.csv", 3, num_active_seeds, seedes, seeds_sizing, seeds_region_id);
        generator.generate_seed_csv_with_faces("seed_pairs_with_faces.csv", "seed_points_with_faces.csv", seeds, nullptr);
        generator.generate_surface_mesh(seeds, nullptr, "surface_mesh1.obj");
        return 0;
    }

    if (argc != 1) {
        std::cerr << "Usage:" << std::endl;
        std::cerr << "  " << argv[0] << std::endl;
        std::cerr << "  " << argv[0] << " paired_seeds.csv" << std::endl;
        std::cerr << "  " << argv[0] << " inner.txt outer.txt" << std::endl;
        std::cerr << "  " << argv[0] << " --surface-spheres remesh.obj spheres.csv" << std::endl;
        std::cerr << "  " << argv[0] << " --inout-benchmark reconstructed.obj reference.obj seed_pairs.csv [num_tests] [csv_output]" << std::endl;
        return 1;
    }

    return RunSurfaceSpherePipeline(
        "../data/obj/QuadCover_8000_gear_Remesh_qem.obj",
        "../data/spheres/QuadCover_8000_gear_Spheres_qem.csv"
    );
}
