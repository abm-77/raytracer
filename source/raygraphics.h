#pragma once
#include "raymath.h"

typedef V3F Color;

class Material; // forward decalare
class Texture;
class SolidColorTexture;

class AABB {
public:
	AABB(){}
	AABB(const Point3& a, const Point3& b) : min(a), max(b) {}
	b32 hit (const Ray& ray, f32 t_min, f32 t_max) const;
public:
	Point3 min;
	Point3 max;
};

internal AABB surrounding_box (const AABB& box0, const AABB& box1) {
	Point3 small (Min (box0.min.x, box1.min.x), Min(box0.min.y, box1.min.y), Min(box0.min.z, box1.min.z));
	Point3 big (Max (box0.max.x, box1.max.x), Max(box0.max.y, box1.max.y), Max(box0.max.z, box1.max.z));
	return AABB(small, big);
}

typedef struct HitResult {
	Point3 p;
	V3F normal;
	f32 t;
	f32 u;
	f32 v;
	b32 front_face;
	Material* p_material;
} HitResult;
inline void set_front_face_normal (HitResult& hit, const Ray& ray, const V3F& outward_normal);

class Hittable {
public:
	virtual b32 hit (const Ray& ray, f32 t_min, f32 t_max, HitResult& result) const = 0;
	virtual b32 bounding_box(f32 time0, f32 time1, AABB& output_box) const = 0;
};

typedef b32 (*hittable_comparator)(const Hittable* a, const Hittable* b);

class Sphere : public Hittable {
public:
	Sphere() {}
	Sphere(Point3 center, f32 radius, Material* material) : center(center), radius(radius), p_material(material) {}
	virtual b32 hit (const Ray& ray, f32 t_min, f32 t_max, HitResult& result) const override; 
	virtual b32 bounding_box(f32 time0, f32 time1, AABB& output_box) const override;
private:
	static void calculate_uv(const Point3& p, f32& u, f32& v);
public:
	Point3 center;
	f32 radius;
	Material* p_material;
};

class MovingSphere : public Hittable {
public:
	MovingSphere() {}
	MovingSphere(Point3 cen0, Point3 cen1, f32 _time0, f32 _time1, f32 radius, Material* material) 
		: center0(cen0), center1(cen1), time0(_time0), time1(_time1), radius(radius), p_material(material) {}
	virtual b32 hit (const Ray& ray, f32 t_min, f32 t_max, HitResult& result) const override; 
	virtual b32 bounding_box(f32 _time0, f32 _time1, AABB& output_box) const override;
	Point3 get_center (f32 time) const;
public:
	Point3 center0, center1;
	f32 time0, time1;
	f32 radius;
	Material* p_material;
};

class XYRect : public Hittable {
public:
	XYRect(){}
	XYRect(f32 _x0, f32 _x1, f32 _y0, f32 _y1, f32 _k, Material* mat) : x0(_x0), x1(_x1), y0(_y0), y1(_y1), k(_k), p_material(mat) {}

	virtual b32 hit (const Ray& ray, f32 t_min, f32 t_max, HitResult& result) const override; 
	virtual b32 bounding_box(f32 _time0, f32 _time1, AABB& output_box) const override;
public:
	Material* p_material;
	f32 x0, x1, y0, y1, k;
};


enum {MAX_LIST_HITTABLES = 1024};
class HittableList : public Hittable {
public:
	HittableList() : count(0) {};
	HittableList(Hittable* object) { count = 0; append(object); }
	void clear ();
	void append (Hittable* object);
	virtual b32 hit (const Ray& ray, f32 t_min, f32 t_max, HitResult& result) const override; 
	virtual b32 bounding_box(f32 time0, f32 time1, AABB& output_box) const override;
	Point3 get_center (f32 time) const;
public:
	Hittable* objects [MAX_LIST_HITTABLES];
	u32 count;
};


class Camera {
public:
	Camera (Point3 lookfrom, Point3 lookat, V3F vup, f32 vfov, f32 aspect_ratio, f32 aperture, f32 focus_dist, f32 _time0, f32 _time1);
	Ray get_ray(f32 s, f32 t) const;
public:
	Point3 origin;
	Point3 lower_left;
	V3F horizontal;
	V3F vertical;
	V3F right, up, forward;
	f32 lens_radius;
	f32 time0, time1;
};


class Material {
public:
	virtual b32 scatter (const Ray& ray_in, const HitResult& result, Color& attenuation, Ray& scattered) const = 0;
	virtual Color emitted(f32 u, f32 v, const Point3& p) const { return Color(0, 0, 0); }
};

class LambertianMaterial : public Material {
public:
	LambertianMaterial (M_Arena* arena, const Color& color);
	LambertianMaterial (Texture* texture) : albedo(texture) {}

	virtual b32 scatter (const Ray& ray_in, const HitResult& result, Color& attenuation, Ray& scattered) const override;
public:
	Texture* albedo;
};

class MetalMaterial : public Material {
public:
	MetalMaterial(M_Arena* arena, const Color& color, f32 roughness); 
	MetalMaterial(Texture* texture, f32 roughness) : albedo(texture), roughness(roughness) {}
	virtual b32 scatter (const Ray& ray_in, const HitResult& result, Color& attenuation, Ray& scattered) const override;
public:
	Texture* albedo;
	f32 roughness; // makes metal less reflective
};

class DielectricMaterial : public Material {
public:
	DielectricMaterial(f32 ri) : refraction_index(ri) {}
	virtual b32 scatter (const Ray& ray_in, const HitResult& result, Color& attenuation, Ray& scattered) const override;
public:
	f32 refraction_index;
private:
	// using Schlick's approximation for reflectance
	static f64 reflectance (f64 cosine, f64 ref_idx) {
		f64 r0 = (1 - ref_idx) / (1 + ref_idx);
		r0 = r0*r0;
		return r0 + (1 - r0) * pow((1 - cosine), 5);
	}
};

class DiffuseLightMaterial : public Material {
public:
	DiffuseLightMaterial(Texture* text) : emit(text) {}
	DiffuseLightMaterial(M_Arena* arena, const Color& c);
	virtual b32 scatter (const Ray& ray_in, const HitResult& result, Color& attenuation, Ray& scattered) const override;
	virtual Color emitted(f32 u, f32 v, const Point3& p) const override;
public:
	Texture* emit;
};

class BVHNode : public Hittable {
public:
	BVHNode () {}
	BVHNode (HittableList& list, f32 time0, f32 time1, M_Arena* arena) { BVHNode(&list.objects[0], 0, list.count, time0, time1, arena); }
	BVHNode (Hittable** src_objects, i32 start, i32 end, f32 time0, f32 time1, M_Arena* arena);
	virtual b32 hit (const Ray& ray, f32 t_min, f32 t_max, HitResult& result) const override; 
	virtual b32 bounding_box(f32 time0, f32 time1, AABB& output_box) const override;
public:
	Hittable* left;
	Hittable* right;
	AABB box;
};

internal void write_color (const Color& pixel_color);
