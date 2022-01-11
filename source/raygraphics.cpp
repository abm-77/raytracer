inline void set_front_face_normal (HitResult& hit, const Ray& ray, const V3F& outward_normal) {
	hit.front_face = dot(ray.direction, outward_normal) < 0;
	hit.normal = hit.front_face ? outward_normal : -outward_normal;
}

b32 Sphere::hit (const Ray& ray, f32 t_min, f32 t_max, HitResult& result) const {
	V3F oc = ray.origin - center;
	f32 a = ray.direction.length_squared();
	f32 half_b = dot(ray.direction, oc);
	f32 c = oc.length_squared() - radius*radius;
	f32 discriminant = (half_b*half_b) - (a*c);
	if (discriminant < 0) return false;

	f32 sqrtd = SquareRoot(discriminant);	

	// Find nearest root that lies in range
	f32 root = (-half_b - SquareRoot(discriminant)) / (a);
	if (root < t_min || root > t_max) {
		root = (-half_b + SquareRoot(discriminant)) / (a);
		if (root < t_min || root > t_max)  return false;
	}

	result.t = root;
	result.p = ray.eval(result.t);
	V3F outward_normal = (result.p - center) / radius;
	set_front_face_normal(result, ray, outward_normal);
	calculate_uv(outward_normal, result.u, result.v);
	result.p_material = p_material;

	return true;
}

b32 Sphere::bounding_box(f32 time0, f32 time1, AABB& output_box) const {
	V3F extents = V3F (radius, radius, radius);
	output_box.min = center - extents;
	output_box.max = center + extents;
	return true;
}

void Sphere::calculate_uv(const Point3& p, f32& u, f32& v) {
	f32 theta = acos(-p.y);
	f32 phi = atan2(-p.z, p.x) + PI_F32;

	u = phi / (2 * PI_F32);
	v = theta / PI_F32;
}

b32 MovingSphere::hit (const Ray& ray, f32 t_min, f32 t_max, HitResult& result) const {
	Point3 cen = get_center(ray.time);
	V3F oc = ray.origin - cen;
	f32 a = ray.direction.length_squared();
	f32 half_b = dot(ray.direction, oc);
	f32 c = oc.length_squared() - radius*radius;
	f32 discriminant = (half_b*half_b) - (a*c);
	if (discriminant < 0) return false;

	f32 sqrtd = SquareRoot(discriminant);	

	// Find nearest root that lies in range
	f32 root = (-half_b - SquareRoot(discriminant)) / (a);
	if (root < t_min || root > t_max) {
		root = (-half_b + SquareRoot(discriminant)) / (a);
		if (root < t_min || root > t_max)  return false;
	}

	result.t = root;
	result.p = ray.eval(result.t);
	set_front_face_normal(result, ray, (result.p - cen) / radius);
	result.p_material = p_material;

	return true;
}

b32 MovingSphere::bounding_box(f32 _time0, f32 _time1, AABB& output_box) const {
	V3F extents = V3F(radius, radius, radius);
	AABB box0 (get_center(_time0) - extents, get_center(_time0) + extents);
	AABB box1 (get_center(_time1) - extents, get_center(_time1) + extents);
	output_box = surrounding_box(box0, box1);
	return true;
}

Point3 MovingSphere::get_center (f32 time) const {
	return center0 + ((time - time0) / (time1 - time0)) * (center1 - center0);
}

b32 XYRect::hit (const Ray& ray, f32 t_min, f32 t_max, HitResult& result) const {
	f32 t = (k - ray.origin.z) / ray.direction.z;
	if (t < t_min || t > t_max) return false;

	f32 x = ray.origin.x + (t * ray.direction.x);
	f32 y = ray.origin.y + (t * ray.direction.y);
	if (x < x0 || x > x1 || y < y0 || y > y1) return false;

	result.u = (x - x0) / (x1 - x0);
	result.v = (y - y0) / (y1 - y0);
	result.t = t;
	V3F outward_normal = V3F(0,0,1.0);
	set_front_face_normal(result, ray, outward_normal);
	result.p_material = p_material;
	result.p = ray.eval(t);
	return true;
}

b32 XYRect::bounding_box(f32 _time0, f32 _time1, AABB& output_box) const {
	output_box = AABB(Point3(x0, y0, k-0.0001), Point3(x1,y1,k+0.0001));
	return true;
}

void HittableList::clear () {
	count = 0;
	memset(objects, 0, MAX_LIST_HITTABLES);
}

void HittableList::append (Hittable* object) {
	objects[count++] = object;
}

b32 HittableList::hit (const Ray& ray, f32 t_min, f32 t_max, HitResult& result) const {
	HitResult temp_result;
	b32 hit_anything = false;
	f32 closest_so_far = t_max;

	for (i32 i = 0; i < count; ++i) {
		if(objects[i]->hit(ray, t_min, closest_so_far, temp_result)) {
			hit_anything = true;
			closest_so_far = temp_result.t;
			result = temp_result;
		}
	}

	return hit_anything;
}

b32 HittableList::bounding_box(f32 time0, f32 time1, AABB& output_box) const {
	if (count == 0) return false;
	AABB temp_box;
	b32 first_box = true;
	for (i32 i = 0; i < count; ++i) {
		if (!objects[i]->bounding_box(time0, time1, temp_box)) return false;
		output_box = (first_box) ? temp_box : surrounding_box (output_box, temp_box);
		first_box = false;
	}
	return true;
}

b32 BVHNode::bounding_box(f32 time0, f32 time1, AABB& output_box) const {
	output_box = box;
	return true;
}

b32 BVHNode::hit (const Ray& ray, f32 t_min, f32 t_max, HitResult& result) const {
	if (!box.hit(ray, t_min, t_max)) return false;
	b32 hit_left = left->hit(ray, t_min, t_max, result);
	b32 hit_right = right->hit(ray, t_min, (hit_left) ? result.t : t_max, result);
	return hit_left || hit_right;
}


inline b32 box_compare (const Hittable* a, const Hittable* b, i32 axis) {
	AABB box_a;
	AABB box_b;

	if (!a->bounding_box(0, 0, box_a) || !b->bounding_box(0, 0, box_b))
		fprintf(stderr, "No bounding box in bvh_node constructor.\n");

	return box_a.min.v[axis] < box_b.min.v[axis];
}

b32 box_x_compare (const Hittable* a, const Hittable* b) {
	return box_compare(a, b, 0);
}
b32 box_y_compare (const Hittable* a, const Hittable* b) {
	return box_compare(a, b, 1);
}
b32 box_z_compare (const Hittable* a, const Hittable* b) {
	return box_compare(a, b, 2);
}

internal void hittable_list_sort (Hittable** list, i32 start, i32 end, hittable_comparator comparator) {
	Hittable* key;
  	i32 i, j;

    for (i = start; i < end; ++i) {
        key = list[i];
        j = i - 1;
 
        while (j >= 0 && comparator(key, list[j])) {
            list[j + 1] = list[j];
            j = j - 1;
        }
        list[j + 1] = key;
    }
}

BVHNode::BVHNode (Hittable** src_objects, i32 start, i32 end, f32 time0, f32 time1, M_Arena* arena) {
	Hittable** objects = src_objects;
	i32 axis = random_int(0,2);

	hittable_comparator comparator 	= 	(axis == 0) ? box_x_compare 
									: 	(axis == 1) ? box_y_compare 
													: box_z_compare;

	i32 object_span = end - start;

	if (object_span == 1) {
		left = right = objects[start];
	}
	else if (object_span == 2) {
		left = objects[start];
		right = objects[start+1];
	}
	else {
		hittable_list_sort(objects, start, end, comparator);
		i32 mid = start + object_span * 0.5f;
		left = M_ArenaPush<BVHNode>(arena, objects, start, mid, time0, time1, arena);
		right = M_ArenaPush<BVHNode>(arena, objects, mid, end, time0, time1, arena);
	}

	AABB box_left, box_right;
	if (!left->bounding_box(time0, time1, box_left) || !right->bounding_box(time0, time1, box_right))
		fprintf(stderr, "No bounding box in bvh_node constructor.\n");

	box = surrounding_box(box_left, box_right);
}

Camera::Camera (Point3 lookfrom, Point3 lookat, V3F vup, f32 vfov, f32 aspect_ratio, f32 aperture, f32 focus_dist, f32 _time0, f32 _time1) {
	f32 theta = deg_to_rad(vfov);
	f32 h = tan(theta / 2);
	f32 viewport_height = 2.0f * h;
	f32 viewport_width = aspect_ratio * viewport_height;

	forward = (lookfrom - lookat).normalized();
	right = cross(vup, forward).normalized();
	up = cross(forward, right);

	origin = lookfrom;
	horizontal = focus_dist * viewport_width * right;
	vertical = focus_dist * viewport_height * up;
	lower_left = origin - (horizontal*0.5f) - (vertical*0.5f) - (focus_dist*forward);

	lens_radius = aperture * 0.5f;

	time0 = _time0;
	time1 = _time1;
}

Ray Camera::get_ray(f32 s, f32 t) const {
	V3F rd = lens_radius * random_in_unit_disk();
	V3F offset = right*rd.x + up*rd.y;
	return Ray(origin + offset, lower_left + s*horizontal + t*vertical - origin - offset, random_double(time0, time1));
}



b32 LambertianMaterial::scatter (const Ray& ray_in, const HitResult& result, Color& attenuation, Ray& scattered) const {
	V3F scatter_direction = result.normal + random_unit_v3();

	// catch degerate scatter direction
	if (scatter_direction.near_zero()) 
		scatter_direction = result.normal;

	scattered = Ray (result.p, scatter_direction, ray_in.time);
	attenuation = albedo->value(result.u, result.v, result.p);	
	return true;
}

MetalMaterial::MetalMaterial(M_Arena* arena, const Color& color, f32 roughness) {
	this->albedo = M_ArenaPush<SolidColorTexture>(arena, color);
	this->roughness = roughness;
}

b32 MetalMaterial::scatter (const Ray& ray_in, const HitResult& result, Color& attenuation, Ray& scattered) const {
	V3F reflected = reflect(ray_in.direction.normalized(), result.normal);
	scattered = Ray(result.p, reflected + roughness * random_in_unit_sphere(), ray_in.time);
	attenuation = albedo->value(result.u, result.v, result.p);
	return dot(scattered.direction, result.normal) > 0;
}

b32 DielectricMaterial::scatter (const Ray& ray_in, const HitResult& result, Color& attenuation, Ray& scattered) const {
	attenuation = Color(1.0f);
	f32 refraction_ratio = (result.front_face) ? (1 / refraction_index) : refraction_index;

	V3F unit_direction = ray_in.direction.normalized();
	f32 cos_theta = Min(dot(-unit_direction, result.normal), 1.0f);
	f32 sin_theta = SquareRoot(1 - cos_theta*cos_theta);
	b32 cannot_refract = refraction_ratio * sin_theta > 1.0; // no solution to Snell's law, total internal refraction
	V3F direction = (cannot_refract || reflectance(cos_theta, refraction_ratio) > random_double()) ? 
		reflect(unit_direction, result.normal) : refract(unit_direction, result.normal, refraction_ratio);

	scattered = Ray(result.p, direction, ray_in.time);
	return true;
}

DiffuseLightMaterial::DiffuseLightMaterial(M_Arena* arena, const Color& c) {
	this->emit = M_ArenaPush<SolidColorTexture>(arena, c);
}

b32 DiffuseLightMaterial::scatter (const Ray& ray_in, const HitResult& result, Color& attenuation, Ray& scattered) const {
	return false;
}

Color DiffuseLightMaterial::emitted(f32 u, f32 v, const Point3& p) const {
	return emit->value(u, v, p);
}

b32 AABB::hit (const Ray& ray, f32 t_min, f32 t_max) const {
	for (i32 dim = 0; dim < 3; ++dim) {
		f32 inv_d = 1 / ray.direction.v[dim];
		f32 t0 = (min.v[dim] - ray.origin.v[dim]) * inv_d;
		f32 t1 = (max.v[dim] - ray.origin.v[dim]) * inv_d;

		if (inv_d <= 0) {
			f32 temp = t0;
			t0 = t1;
			t1 = temp;
		}

		t_min = (t0 > t_min) ? t0 : t_min;
		t_max = (t1 < t_max) ? t1 : t_max;

		if (t_max <= t_min) return false;
	}
	return true;
}


internal void write_color (const Color& pixel_color, i32 samples_per_pixel) {
	f32 r = pixel_color.r;
	f32 g = pixel_color.g;
	f32 b = pixel_color.b;

	// scale and gamma correct colors (gamma-2)
	f32 scale = 1.0f / samples_per_pixel; 
	r = SquareRoot(scale * r);
	g = SquareRoot(scale * g);
	b = SquareRoot(scale * b);

	printf("%d %d %d\n", i32(256 * Clamp(r, 0.0f, 0.999f)), i32(256 * Clamp(g, 0.0f, 0.999f)), i32(256 * Clamp(b, 0.0f, 0.999f)));
}