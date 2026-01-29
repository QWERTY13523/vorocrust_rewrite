#include "Geometry.h"

const double PI = 3.141592653589793;
Geometry::Geometry()
{

}

Geometry::~Geometry()
{

}

double Geometry::get_point_angle(double x, double y)
{
	#pragma region Get Point Angle:
	if (fabs(x) < 1E-10 && fabs(y) < 1E-10) return DBL_MAX; // this is a singularity

	double angle(0.0);
	if (fabs(x) > fabs(y)) angle = atan(fabs(y) / fabs(x));
	else                   angle = 0.5 * PI - atan(fabs(x) / fabs(y));

	if (x > -1E-10 && y > -1E-10) return angle;
	else if (x > -1E-10)          return 2 * PI - angle;
	else if (y > -1E-10)          return PI - angle;
	else                          return PI + angle;
	return 0.0;
	#pragma endregion
}

double Geometry::distance(size_t num_dim, double* x, double* y)
{
	#pragma region Distance between two points:
	double dst = 0.0;
	for (size_t idim = 0; idim < num_dim; idim++)
	{
		double dx = y[idim] - x[idim];
		dst += dx * dx;
	}
	dst = sqrt(dst);
	return dst;
	#pragma endregion
}

double Geometry::cos_angle(size_t num_dim, double* x, double* y, double* z)
{
	#pragma region Distance between two points:
	double hxy = distance(num_dim, x, y);
	double hyz = distance(num_dim, y, z);
	if (hxy < 1E-10 || hyz < 1E-10) return 1.0;
	double dot(0.0);
	for (size_t idim = 0; idim < num_dim; idim++)
	{
		double dxy = x[idim] - y[idim];
		double dzy = z[idim] - y[idim];
		dot += dxy * dzy;
	}
	return dot /= (hxy * hyz);
	#pragma endregion
}

bool Geometry::get_power_vertex(size_t num_dim, size_t num_points, double** centers, double* radii, double* pv)
{
    // 1. 基础检查
    if (num_dim == 0 || num_points < num_dim + 1) return false;

    // 2. 构建线性系统 Ax = b
    // 我们有 (num_points - 1) 个方程，num_dim 个变量
    // A 的每一行是两个球心的向量差: 2 * (c_i - c_0)
    Eigen::MatrixXd A(num_points - 1, num_dim);
    Eigen::VectorXd b(num_points - 1);

    // 预计算第0个球的参数
    Eigen::VectorXd c0(num_dim);
    for (size_t d = 0; d < num_dim; ++d) c0[d] = centers[0][d];
    double pow0 = c0.squaredNorm() - radii[0] * radii[0];

    for (size_t i = 1; i < num_points; ++i)
    {
        Eigen::VectorXd ci(num_dim);
        for (size_t d = 0; d < num_dim; ++d) ci[d] = centers[i][d];
        
        double pow_i = ci.squaredNorm() - radii[i] * radii[i];

        // 填充 A 的第 (i-1) 行: 2 * (c_i - c_0)
        A.row(i - 1) = 2.0 * (ci - c0);
        
        // 填充 b 的第 (i-1) 行: pow_i - pow0
        b(i - 1) = pow_i - pow0;
    }

    // 3. 使用列主元 Householder QR 分解求解
    // 这种方法对于超定方程组（Overdetermined）会自动给出最小二乘解
    // 对于病态矩阵（接近奇异），它比普通 LU 分解更稳定
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> solver(A);
    
    // 检查矩阵的秩，如果秩小于维度，说明球心共面/共线，无法确定唯一的根心
    if (solver.rank() < num_dim) {
        return false; // 无唯一解
    }

    Eigen::VectorXd x = solver.solve(b);

    // 4. (可选) 验证解的质量
    // 如果你希望在“无解”时严格返回 false，而不是返回近似解，
    // 可以检查残差 (Ax - b) 的模长。
    if ((A * x - b).norm() > 1e-5) {
        // 残差过大，说明这些球并没有公共的根心（Power Vertex）
        // 如果你的应用允许误差，可以去掉这个检查
        return false; 
    }

    // 5. 输出结果
    for (size_t d = 0; d < num_dim; ++d) pv[d] = x[d];

    return true;
}

int Geometry::get_normal_component(size_t num_dim, size_t num_basis, double** basis, double* vect, double &norm)
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

double Geometry::dot_product(size_t num_dim, double* v1, double* v2)
{
	#pragma region dot product:
	double dot = 0.0;
	for (size_t idim = 0; idim < num_dim; idim++)
	{
		dot += v1[idim] * v2[idim];
	}
	return dot;
	#pragma endregion
}

bool Geometry::normalize_vector(size_t num_dim, double* vec)
{
	#pragma region Normalize Vecotr:
	double norm(0.0);
	for (size_t idim = 0; idim < num_dim; idim++)
	{
		norm += vec[idim] * vec[idim];
	}
	norm = sqrt(norm);
	if (norm < 1E-10) return false;
	for (size_t idim = 0; idim < num_dim; idim++) vec[idim] /= norm;
	return true;
	#pragma endregion
}

bool Geometry::overlapping_spheres(size_t num_dim, double* sphere_i, double* sphere_j)
{
	#pragma region Overlapping spheres check:
	double h = distance(num_dim, sphere_i, sphere_j);
	double ri = sphere_i[num_dim];
	double rj = sphere_j[num_dim];
	if (h < ri + rj + 1E-10) return true;
	return false;
	#pragma endregion
}

double Geometry::get_3d_triangle_area(double** corners)
{
	#pragma region Get Area of a triangle:
	double* a = new double[3]; double* b = new double[3];
	for (size_t idim = 0; idim < 3; idim++)
	{
		a[idim] = corners[1][idim] - corners[0][idim];
		b[idim] = corners[2][idim] - corners[0][idim];
	}
	double nx = a[1] * b[2] - a[2] * b[1];
	double ny = a[2] * b[0] - a[0] * b[2];
	double nz = a[0] * b[1] - a[1] * b[0];
	double norm = sqrt(nx * nx + ny * ny + nz * nz);
	delete[] a; delete[] b;
	return 0.5 * norm;
	#pragma endregion
}

bool Geometry::get_3d_triangle_normal(double** corners, double* normal)
{
	#pragma region Get Face Area and Normal:
	double* a = new double[3]; double* b = new double[3];
	for (size_t idim = 0; idim < 3; idim++)
	{
		a[idim] = corners[1][idim] - corners[0][idim];
		b[idim] = corners[2][idim] - corners[0][idim];
	}
	if (!normalize_vector(3, a) || !normalize_vector(3, b))
	{
		for (size_t idim = 0; idim < 3; idim++) normal[idim] = 0.0;
		delete[] a; delete[] b;
		return false; // a degenerate face due to edge length
	}
	normal[0] = a[1] * b[2] - a[2] * b[1];
	normal[1] = a[2] * b[0] - a[0] * b[2];
	normal[2] = a[0] * b[1] - a[1] * b[0];
	double norm = normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2];
	norm = sqrt(norm);
	if (fabs(norm) < 1E-10)
	{
		delete[] a; delete[] b;
		return false; // a degenerate face due to angle
	}
	for (size_t idim = 0; idim < 3; idim++) normal[idim] /= norm;
	delete[] a; delete[] b;
	return true;
	#pragma endregion
}

int Geometry::project_to_3d_triangle(double* p, double** corners, double* q, double &proj_dist)
{
	#pragma region Project to 3d Triangle:
	double* n = new double[3];
	get_3d_triangle_normal(corners, n);
	double* a = new double[3]; double* b = new double[3];
	for (size_t idim = 0; idim < 3; idim++) a[idim] = corners[0][idim] - p[idim];
	double dot = dot_product(3, a, n);
	for (size_t idim = 0; idim < 3; idim++) q[idim] = p[idim] + dot * n[idim];

	bool in_triangle(true);
	double* m = new double[3];
	for (size_t i = 0; i < 3; i++)
	{
		size_t ip = i + 1; if (ip == 3) ip = 0;
		size_t im = 2; if (i > 0) im = i - 1;
		double* tmp = corners[im];
		corners[im] = q;
		if (!get_3d_triangle_normal(corners, m))
		{
			corners[im] = tmp;
			continue;
		}
		dot = dot_product(3, m, n);
		corners[im] = tmp;
		if (dot < -1E-10)
		{
			in_triangle = false; break;
		}
	}
	if (in_triangle)
	{
		proj_dist = distance(3, p, q);
		delete[] n; delete[] m; delete[] a; delete[] b;
		return 0;
	}

	// project to edges/corners
	double hclosest(DBL_MAX);
	for (size_t i = 0; i < 3; i++)
	{
		size_t ip = i + 1; if (ip == 3) ip = 0;
		size_t im = 2; if (i > 0) im = i - 1;
		for (size_t idim = 0; idim < 3; idim++)
		{
			a[idim] = p[idim] - corners[i][idim];
			b[idim] = corners[ip][idim] - corners[i][idim];
		}
		dot = dot_product(3, a, b);

		double b_norm(0.0);
		for (size_t idim = 0; idim < 3; idim++) b_norm += b[idim] * b[idim];

		double u = dot / b_norm;
		if (u < 0.0) u = 0.0;
		if (u > 1.0) u = 1.0;
		for (size_t idim = 0; idim < 3; idim++) a[idim] = corners[i][idim] + u * (corners[ip][idim] - corners[i][idim]);

		double h = distance(3, a, p);
		if (h < hclosest)
		{
			hclosest =  h;
			for (size_t idim = 0; idim < 3; idim++) q[idim] = a[idim];
		}
	}
	proj_dist = distance(3, p, q);
	delete[] n; delete[] m; delete[] a; delete[] b;
	return 0;
	#pragma endregion
}

int Geometry::project_to_3d_line_segment(double* p, double* xo, double* xn, double* q, double& proj_dist)
{
	#pragma region Project to 3d line segement:
	double* a = new double[3];
	double* b = new double[3];
	for (size_t idim = 0; idim < 3; idim++)
	{
		a[idim] = p[idim] - xo[idim];
		b[idim] = xn[idim] - xo[idim];
	}
	double dot = dot_product(3, a, b);

	double b_norm(0.0);
	for (size_t idim = 0; idim < 3; idim++) b_norm += b[idim] * b[idim];

	double u = dot / b_norm;
	if (u < 0.0) u = 0.0;
	if (u > 1.0) u = 1.0;
	for (size_t idim = 0; idim < 3; idim++) q[idim] = xo[idim] + u * (xn[idim] - xo[idim]);
	proj_dist = distance(3, p, q);

	delete[] a; delete[] b;
	return 0;
	#pragma endregion
}

int Geometry::project_to_3d_line(double* p, double* xo, double* edir, double* q, double &proj_dist)
{
	#pragma region Project to 3d line segement:
	double* a = new double[3];
	for (size_t idim = 0; idim < 3; idim++) a[idim] = p[idim] - xo[idim];

	double u = dot_product(3, a, edir);
	for (size_t idim = 0; idim < 3; idim++) q[idim] = xo[idim] + u * edir[idim];
	proj_dist = distance(3, p, q);

	delete[] a;
	return 0;
	#pragma endregion
}
