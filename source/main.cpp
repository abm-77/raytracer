#include <stdio.h>

#include "raytracer.inc"

M_Arena program_memory;

Color ray_color (const Ray& ray, const Color& background_color, const Hittable& world, i32 depth) {
    HitResult result;

    // If we've exceeded the ray bounce limit, no more light is gathered.
    if (depth <= 0)
        return Color(0,0,0);

    // If the ray hits nothing, return the background color.
    if (!world.hit(ray, 0.001, INF_F64, result))
        return background_color;

    Ray scattered;
    Color attenuation;
    Color emitted = result.p_material->emitted(result.u, result.v, result.p);

    if (!result.p_material->scatter(ray, result, attenuation, scattered))
        return emitted;

    return emitted + attenuation * ray_color(scattered, background_color, world, depth-1);
}

HittableList* earth() {
    ImageTexture* earth_texture = M_ArenaPush<ImageTexture>(&program_memory, "./earthmap.jpeg");
    LambertianMaterial* earth_surface = M_ArenaPush<LambertianMaterial>(&program_memory, earth_texture);
    Sphere* globe = M_ArenaPush<Sphere>(&program_memory, Point3(0,2,0), 2, earth_surface);
    return M_ArenaPush<HittableList>(&program_memory, globe);
}

HittableList* light_test() {
	SolidColorTexture* ground_texture = M_ArenaPush<SolidColorTexture>(&program_memory, Color(0.0, 1.0, 0.0));
	LambertianMaterial* ground = M_ArenaPush<LambertianMaterial>(&program_memory, ground_texture);

	SolidColorTexture* light_texture = M_ArenaPush<SolidColorTexture>(&program_memory, 4 * Color(1.0, 1.0, 1.0));
	DiffuseLightMaterial* light_mat =  M_ArenaPush<DiffuseLightMaterial>(&program_memory, light_texture);

	HittableList* world = M_ArenaPush<HittableList>(&program_memory);
	world->append(M_ArenaPush<Sphere>(&program_memory, Point3(0, -1000, 0), 1000, ground));
	world->append(M_ArenaPush<XYRect>(&program_memory, 3, 5, 1, 3, -2, light_mat));
	return world;
}

int main (void) {
	// Memory
	program_memory = M_ArenaMake();

	// Image 
	const f32 aspect_ratio = 16.0f / 9.0f;
	const i32 image_width = 400;
	const i32 image_height = image_width / aspect_ratio;
	const i32 samples_per_pixel = 100;
	const i32 max_depth = 50;

	// World
	HittableList* world = light_test();

	// Camera
	Point3 lookfrom(26,3,6);
	Point3 lookat(0,2,0);
	V3F world_up(0,1,0);
	f32 dist_to_focus = 10;
	f32 aperture = 0.0f;
	f32 vfov = 20.0;
	const Color background_color = Color(0.0, 0.0, 0.0);
	Camera camera (lookfrom, lookat, world_up, vfov, aspect_ratio, aperture, dist_to_focus, 0.0f, 1.0f);

	// Render
	printf("P3\n%d %d\n255\n", image_width, image_height);
	for (i32 y = image_height - 1; y >= 0; --y){
		fprintf(stderr, "\rScanlines Remaining: %d\n", y); fflush(stderr);
		for (i32 x = 0; x < image_width; ++x) {
			Color pixel_color (0, 0, 0);
			for (i32 s = 0; s < samples_per_pixel; ++s) {
				f64 u = (x + random_double()) / (image_width - 1);
				f64 v = (y + random_double()) / (image_height - 1);
				Ray r = camera.get_ray(u, v);
				pixel_color += ray_color(r, background_color, *world, max_depth);
			}
			write_color(pixel_color, samples_per_pixel);
		}
	}

	M_ArenaRelease(&program_memory);
	return 0;
}