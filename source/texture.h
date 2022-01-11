#pragma once

class Texture {
public:
	virtual Color value (f32 u, f32 v, const Point3& p) const = 0;
};

class SolidColorTexture : public Texture {
public:
	SolidColorTexture();
	SolidColorTexture(const Color& c) : color_value(c) {}
	SolidColorTexture(f32 r, f32 g, f32 b) : color_value(Color(r, g ,b)) {}

	virtual Color value(f32 u, f32 v, const Point3& p) const override;
private:
	Color color_value;
};

class CheckerTexture : public Texture {
public:
	CheckerTexture() {};
	CheckerTexture(M_Arena* arena, const Color& c1, const Color& c2);
	CheckerTexture(Texture* even, Texture* odd) : even(even), odd(odd) {}

	virtual Color value(f32 u, f32 v, const Point3& p) const override;
public:
	Texture* even;
	Texture* odd;
};

enum {BYTES_PER_PIXEL = 3};
class ImageTexture : public Texture {
public:
	ImageTexture() : data(NULL), width(0), height(0), bytes_per_scanline(0){}
	ImageTexture(const char* file_name);
	~ImageTexture();
	virtual Color value(f32 u, f32 v, const Point3& p) const override;
private:
	u8* data;
	i32 width, height;
	i32 bytes_per_scanline;

};