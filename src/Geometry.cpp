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

void Geometry::get_power_vertex(size_t num_dim, double* co, double ro, double* c1, double r1, double* pv)
{
    // 计算圆心距离
    double h_sq = 0.0;
    for (size_t idim = 0; idim < num_dim; idim++) {
        double dx = c1[idim] - co[idim];
        h_sq += dx * dx;
    }
    double h = sqrt(h_sq);

    // 【鲁棒性检查】如果球心重合或极度接近，防止除以零
    if (h < 1E-10) {
        // 球心重合时，中点没有良定义（整个空间都是等幂面，或者无解）
        // 这里简单返回 co 作为 fallback，或者你可以设为 NaN
        for (size_t idim = 0; idim < num_dim; idim++) pv[idim] = co[idim];
        return;
    }

    // 解析公式: x = (d^2 + r1^2 - r2^2) / 2d
    // ho 是从 co 出发，沿 (c1-co) 方向的距离
    double ho = (h_sq + ro * ro - r1 * r1) / (2.0 * h);

    // 线性插值计算坐标
    double ratio = ho / h; 
    for (size_t idim = 0; idim < num_dim; idim++) {
        pv[idim] = co[idim] + ratio * (c1[idim] - co[idim]);
    }
}

// ---------------------------------------------------------
// 2. 多球情况（主要针对3球）：求球心平面上的根心 (Radical Center)
// ---------------------------------------------------------
bool Geometry::get_power_vertex(size_t num_dim, size_t num_points, double** centers, double* radii, double* pv)
{
    // 目前 VoroCrust 逻辑主要调用 3 球情况，这里针对 3 球做优化
    // 使用局部坐标系投影法（Gram-Schmidt），数值最稳定
    if (num_points != 3) {
        // 如果未来有 4 球需求，需扩展此函数求解线性方程组 Ax=b
        return false; 
    }

    double* p0 = centers[0];
    double* p1 = centers[1];
    double* p2 = centers[2];
    double r0 = radii[0];
    double r1 = radii[1];
    double r2 = radii[2];

    // --- 构建局部坐标系 (Origin at p0) ---
    // X轴方向向量: p1 - p0
    double v10[3];
    double d10_sq = 0.0;
    for(int i=0; i<3; ++i) { v10[i] = p1[i] - p0[i]; d10_sq += v10[i]*v10[i]; }
    double d10 = sqrt(d10_sq);

    if (d10 < 1E-10) return false; // p0 和 p1 重合

    // 单位化 X 轴
    double axis_x[3];
    for(int i=0; i<3; ++i) axis_x[i] = v10[i] / d10;

    // 辅助向量: p2 - p0
    double v20[3];
    for(int i=0; i<3; ++i) v20[i] = p2[i] - p0[i];
    
    // 投影 v20 到 axis_x
    double proj = 0.0;
    for(int i=0; i<3; ++i) proj += v20[i] * axis_x[i];

    // Y轴方向向量: v20 - proj * axis_x (Gram-Schmidt 正交化)
    double axis_y[3];
    double d_y_sq = 0.0;
    for(int i=0; i<3; ++i) {
        axis_y[i] = v20[i] - proj * axis_x[i];
        d_y_sq += axis_y[i] * axis_y[i];
    }
    double d_y = sqrt(d_y_sq);

    // 【关键检查】三点共线判定
    // 如果 p2 到 p0p1 直线的距离极小，说明共线，无法构成平面，无唯一根心
    if (d_y < 1E-10) return false; 

    // 单位化 Y 轴
    for(int i=0; i<3; ++i) axis_y[i] /= d_y;

    // --- 转换到局部 2D 坐标 (x, y) ---
    // p0: (0, 0)
    // p1: (d10, 0)
    // p2: (proj, d_y) -> 记为 (x2, y2)
    double x2 = proj;
    double y2 = d_y;

    // --- 求解 2D 根心 (rx, ry) ---
    // 建立方程组：
    // 1. 对 p0, p1: 2 * d10 * rx = d10^2 + r0^2 - r1^2
    // 2. 对 p0, p2: 2 * x2 * rx + 2 * y2 * ry = x2^2 + y2^2 + r0^2 - r2^2

    double rx = (d10_sq + r0*r0 - r1*r1) / (2.0 * d10);
    
    double term2 = (x2*x2 + y2*y2 + r0*r0 - r2*r2);
    double ry = (term2 - 2.0 * x2 * rx) / (2.0 * y2);

    // --- 转换回世界坐标 ---
    // P = p0 + rx * axis_x + ry * axis_y
    for(size_t i=0; i<num_dim; ++i) {
        pv[i] = p0[i] + rx * axis_x[i] + ry * axis_y[i];
    }

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
