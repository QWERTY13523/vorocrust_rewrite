#include "Generator.h"
#include "MeshingTree.h"
#include <vector>
#include <iostream>
#include <map>
#include <Eigen/Dense>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/bfs_orient.h> // 用于统一面片方向
#include <igl/volume.h>     // 用于计算体积判断朝向
#include <numeric>
#include <Eigen/Geometry>

extern void optimizer(MeshingTree* seeds, MeshingTree* spheres, std::vector<int> face_flat);

double **points;
size_t **faces;
double *seeds_sizing, *seedes;
size_t *seeds_region_id;
std::vector<int> face_flat;

void FixMeshNormals(const std::string& input_path, const std::string& output_path) {
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;

    // 1. 读取 OBJ
    if (!igl::readOBJ(input_path, V, F)) {
        std::cerr << "Error: Cannot read " << input_path << std::endl;
        return;
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
    } else {
        std::cout << "已保存修复后的模型至: " << output_path << std::endl;
    }
}
int main()
{
    FixMeshNormals("../data/obj/Ours_6000_block_Remesh_density_modified.obj", "../data/obj/Ours_6000_block_Remesh_density_modified_Fixed.obj");
    Generator generator;
    MeshingTree *spheres = new MeshingTree();
    MeshingTree *upper_seeds = new MeshingTree();
    MeshingTree *lower_seeds = new MeshingTree();
    MeshingTree *seeds = new MeshingTree();
    size_t num_points, num_faces, num_faces1;
    generator.read_input_obj_file("../data/obj/block.obj",num_points,points,num_faces,faces);
    if (num_points == 0) {
        std::cerr << "[Error] Failed to read input points or empty file." << std::endl;
        return 1;
    }

    generator.generate_spheres("../data/spheres/Sphere_6000_82_block.csv", spheres);
    if (spheres->get_num_tree_points() == 0) {
        std::cerr << "[Error] Failed to load spheres or empty file." << std::endl;
        return 1;
    }

    generator.read_obj_faces("../data/obj/Ours_6000_block_Remesh_density_modified_Fixed.obj",face_flat, num_faces1);
    if (face_flat.empty()) {
        std::cerr << "[Error] Failed to read face indices or empty file." << std::endl;
        return 1;
    }
    generator.generate_surface_seeds(num_points, points, num_faces, faces,
        spheres, upper_seeds, lower_seeds, num_faces1, face_flat);
    generator.color_surface_seeds(num_faces1, spheres,upper_seeds, lower_seeds, seeds,
        face_flat, points, seedes, seeds_region_id, seeds_sizing);
    optimizer(seeds, spheres, face_flat);
    generator.generate_surface_mesh(seeds, spheres, "surface_mesh1.obj");
}