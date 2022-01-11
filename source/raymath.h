#pragma once

#include "language_layer.h"

typedef union V3F {
	V3F (){}
	V3F (f32 f) : x(f), y(f), z(f) {}
	V3F (f32 x, f32 y, f32 z) : x(x), y(y), z(z) {};
	f32 length () const;
	f32 length_squared () const;
	b32 near_zero() const;
	V3F normalized() const;
	V3F operator- () const;
	V3F& operator+= (const V3F& other);
	V3F& operator*= (const f32 t);
	V3F& operator/= (const f32 t);
	struct {f32 x, y, z;};
	struct {f32 r, g, b;};
	f32 v[3];
} V3F;

typedef V3F Point3;

internal inline void print_vector (const V3F &v);
internal inline V3F operator+(const V3F &u, const V3F &v);
internal inline V3F operator-(const V3F &u, const V3F &v);
internal inline V3F operator*(const V3F &u, const V3F &v);
internal inline V3F operator*(f32 t, const V3F &v);
internal inline V3F operator*(const V3F &v, f32 t);
internal inline V3F operator/(const V3F& v, f32 t);
internal inline f32 dot(const V3F& u, const V3F &v);
internal inline V3F cross(const V3F& u, const V3F &v); 
internal V3F reflect (const V3F& v, const V3F& n);
internal V3F refract (const V3F& uv, const V3F& n, f64 etai_over_etat);

class Ray {
public:
	Ray () {}
	Ray (const Point3& o, const V3F& d, f32 tm) : origin(o), direction(d), time(tm) {}	
	Point3 eval (const f32 t) const;
public:
	Point3 origin;
	V3F direction;
	f32 time;
};

internal inline f64 deg_to_rad(f64 degrees);
internal inline f64 random_double();
internal inline i32 random_int();
internal inline i32 random_int(i32 min, i32 max);
internal inline f64 random_double(f64 min, f64 max);
internal inline V3F random_v3();
internal inline V3F random_v3(f64 min, f64 max);
internal V3F random_unit_v3();
internal V3F random_in_unit_disk();