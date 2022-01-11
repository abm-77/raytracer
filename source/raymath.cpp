#pragma region VECTOR
f32 V3F::length () const {
	return SquareRoot(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}
f32 V3F::length_squared () const {
	return v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
}

// return true if all dimensions are near zero
b32 V3F::near_zero() const {
	const double s = 1e-8;
	return AbsoluteValue(v[0]) < s && AbsoluteValue(v[1]) < s && AbsoluteValue(v[2]) < s;

}
V3F V3F::normalized() const {
	f32 l = 1 / this->length();
	return V3F(v[0] * l, v[1] * l, v[2] * l);
}
V3F V3F::operator- () const {
	return V3F(-v[0], -v[1], -v[2]);
}
V3F& V3F::operator+= (const V3F& other) {
	v[0] += other.v[0];
	v[1] += other.v[1];
	v[2] += other.v[2];
	return *this;
}
V3F& V3F::operator*= (const f32 t) {
	v[0] *= t; 
	v[1] *= t;
	v[2] *= t;
	return *this;
}
V3F& V3F::operator/= (const f32 t) {
	f32 inv = 1/t;
	v[0] *= inv;
	v[1] *= inv;
	v[2] *= inv;
	return *this;
}

internal inline void print_vector (const V3F &v) {
    printf("%f %f %f\n", v.v[0], v.v[1], v.v[2]);
}
internal inline V3F operator+(const V3F &u, const V3F &v) {
    return V3F (u.v[0] + v.v[0], u.v[1] + v.v[1], u.v[2] + v.v[2]);
}
internal inline V3F operator-(const V3F &u, const V3F &v) {
    return V3F (u.v[0] - v.v[0], u.v[1] - v.v[1], u.v[2] - v.v[2]);
}
internal inline V3F operator*(const V3F &u, const V3F &v) {
    return V3F (u.v[0] * v.v[0], u.v[1] * v.v[1], u.v[2] * v.v[2]);
}
internal inline V3F operator*(f32 t, const V3F &v) {
    return V3F (t*v.v[0], t*v.v[1], t*v.v[2]);
}
internal inline V3F operator*(const V3F &v, f32 t) {
    return t * v;
}
internal inline V3F operator/(const V3F& v, f32 t) {
    return (1/t) * v;
}
internal inline f32 dot(const V3F &u, const V3F &v) {
    return u.v[0] * v.v[0]
         + u.v[1] * v.v[1]
         + u.v[2] * v.v[2];
}
internal inline V3F cross(const V3F &u, const V3F &v) {
    return V3F (u.v[1] * v.v[2] - u.v[2] * v.v[1],
                u.v[2] * v.v[0] - u.v[0] * v.v[2],
                u.v[0] * v.v[1] - u.v[1] * v.v[0]);
}

internal V3F reflect (const V3F& v, const V3F& n) {
	return v - 2 * dot(v, n) * n;
}

internal V3F refract (const V3F& uv, const V3F& n, f64 etai_over_etat) {
	f64 cos_theta = Min(dot(-uv, n), 1.0);
	V3F r_out_perp = etai_over_etat * (uv + cos_theta*n);
	V3F r_out_parallel = -SquareRoot(AbsoluteValue(1 - r_out_perp.length_squared())) * n;
	return r_out_perp + r_out_parallel;
}

#pragma region RAY

Point3 Ray::eval (const f32 t) const {
	return origin + t * direction;
}

internal inline f64 deg_to_rad(f64 degrees) {
    return degrees * PI_F64 / 180.0;
}

// returns in range [0, 1)
internal inline f64 random_double () {
	return rand() / (RAND_MAX + 1.0);
}

internal inline i32 random_int() {
	return (int) (rand() / (RAND_MAX + 1.0));
}
internal inline i32 random_int(i32 min, i32 max) {
	return min + (max-min) * (rand() / (RAND_MAX + 1.0));
}

internal inline f64 random_double (f64 min, f64 max) {
	return min + (max-min) * (rand() / (RAND_MAX + 1.0));
}

internal inline V3F random_v3 () {
	return V3F (random_double(), random_double(), random_double());
}

internal inline V3F random_v3 (f64 min, f64 max) {
	return V3F (random_double(min, max), random_double(min, max), random_double(min, max));
}

internal V3F random_in_unit_sphere() {
	while (true) {
		V3F p = random_v3();
		if(p.length_squared() >= 1) continue;
		return p;
	}
}

internal V3F random_in_hemisphere(const V3F& normal) {
	V3F in_sphere = random_in_unit_sphere();
	if (dot(normal, in_sphere) > 0) // in the same hemisphere as normal
		return in_sphere;
	else
		return -in_sphere;
}

// makes a random unit vector (point on unit sphere) by normalizing a point in the unit sphere
internal V3F random_unit_v3() {
	return random_in_unit_sphere().normalized();
}

internal V3F random_in_unit_disk() {
	while (true) {
		Point3 p = V3F(random_double(-1,1), random_double(-1,1), 0);
		if (p.length_squared() >= 1) continue;
		return p;
	}
}