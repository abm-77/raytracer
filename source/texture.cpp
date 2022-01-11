
Color SolidColorTexture::value(f32 u, f32 v, const Point3& p) const {
	return color_value;
}

CheckerTexture::CheckerTexture(M_Arena* arena, const Color& c1, const Color& c2) {
	this->even 	= M_ArenaPush<SolidColorTexture>(arena, c1);
	this->odd 	= M_ArenaPush<SolidColorTexture>(arena, c2);

}

Color CheckerTexture::value(f32 u, f32 v, const Point3& p) const {
	f32 sine = sin(10 * p.x) * sin(10 * p.y) * sin(10 * p.z);
	if (sine < 0)
		return odd->value(u, v, p);
	else
		return even->value(u, v, p);
}

ImageTexture::ImageTexture(const char* file_name) {
	i32 components_per_pixel = BYTES_PER_PIXEL;
	data = stbi_load(file_name, &width, &height, &components_per_pixel, components_per_pixel);

	if(!data) {
		printf("Could not load image: %s!\n", file_name);
		width = height = 0;
	}

	bytes_per_scanline = BYTES_PER_PIXEL * width;
}

ImageTexture::~ImageTexture() {
	free(data);
}

Color ImageTexture::value(f32 u, f32 v, const Point3& p) const {
	if (!data)
		return Color(1.0, 0.0, 0.0);
	
	u = Clamp(u, 0.0, 1.0);
	v = 1.0 - Clamp(v, 0.0, 1.0);

	i32 i = (i32) (u * width);
	i32 j = (i32) (v * height);

	if (i >= width) i = width - 1;
	if (j >= height) i = height - 1;

	f32 color_scale = 1.0 / 255.0;
	u8* pixel = data + (j * bytes_per_scanline) + (i * BYTES_PER_PIXEL);

	return Color(color_scale * pixel[0], color_scale * pixel[1], color_scale * pixel[2]);
}