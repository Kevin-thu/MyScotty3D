
#include "shape.h"
#include "../geometry/util.h"

namespace Shapes {

Vec2 Sphere::uv(Vec3 dir) {
	float u = std::atan2(dir.z, dir.x) / (2.0f * PI_F);
	if (u < 0.0f) u += 1.0f;
	float v = std::acos(-1.0f * std::clamp(dir.y, -1.0f, 1.0f)) / PI_F;
	return Vec2{u, v};
}

BBox Sphere::bbox() const {
	BBox box;
	box.enclose(Vec3(-radius));
	box.enclose(Vec3(radius));
	return box;
}

PT::Trace Sphere::hit(Ray ray) const {
	//A3T2 - sphere hit

    // TODO (PathTracer): Task 2
    // Intersect this ray with a sphere of radius Sphere::radius centered at the origin.

    // If the ray intersects the sphere twice, ret should
    // represent the first intersection, but remember to respect
    // ray.dist_bounds! For example, if there are two intersections,
    // but only the _later_ one is within ray.dist_bounds, you should
    // return that one!

    PT::Trace ret;
    ret.origin = ray.point;
	float det = 4 * pow(dot(ray.point, ray.dir), 2) - 4 * ray.dir.norm_squared() * (ray.point.norm_squared() - radius * radius);
	if (det < 0) {
		ret.hit = false;
		return ret;
	}
	float t1 = (-2 * dot(ray.point, ray.dir) - sqrt(det)) / (2 * ray.dir.norm_squared());
	float t2 = (-2 * dot(ray.point, ray.dir) + sqrt(det)) / (2 * ray.dir.norm_squared());
	if (ray.dist_bounds.x <= t1) {
		if (ray.dist_bounds.y < t1) {
			ret.hit = false;
			return ret;
		}
		ret.distance = t1;
		ret.position = ray.at(t1);
		ret.normal = ret.position.unit(); // normal向量应该做归一化
		ret.uv = uv(ret.position);
	} else {
		if (ray.dist_bounds.y < t2 || ray.dist_bounds.x > t2) {
			ret.hit = false;
			return ret;
		}
		ret.distance = t2;
		ret.position = ray.at(t2);
		ret.normal = ret.position.unit();
		ret.uv = uv(ret.position);
	}
	ret.hit = true;
    return ret;
}

Vec3 Sphere::sample(RNG &rng, Vec3 from) const {
	die("Sampling sphere area lights is not implemented yet.");
}

float Sphere::pdf(Ray ray, Mat4 pdf_T, Mat4 pdf_iT) const {
	die("Sampling sphere area lights is not implemented yet.");
}

Indexed_Mesh Sphere::to_mesh() const {
	return Util::closed_sphere_mesh(radius, 2);
}

} // namespace Shapes

bool operator!=(const Shapes::Sphere& a, const Shapes::Sphere& b) {
	return a.radius != b.radius;
}

bool operator!=(const Shape& a, const Shape& b) {
	if (a.shape.index() != b.shape.index()) return false;
	return std::visit(
		[&](const auto& shape) {
			return shape != std::get<std::decay_t<decltype(shape)>>(b.shape);
		},
		a.shape);
}
