#include <stdlib.h>
#include <math.h>
#include "image.h"

float nn_interpolate(image im, float x, float y, int c)
{
    // TODO Fill in
	int rx = roundf(x);
	int ry = roundf(y);
	float pix_val = get_pixel(im, rx, ry, c);
    return pix_val;
}

image nn_resize(image im, int w, int h)
{
    // TODO Fill in (also fix that first line)

	float w_ratio = (float)im.w/w;
	float h_ratio = (float)im.h/h;	
	image resized_im = make_image(w, h, im.c);
	for(int i = 0; i < im.c; i++){
		for(int j = 0; j < h; j++){
			for(int k = 0; k < w; k++){
				float im_x = (-0.5) + w_ratio / 2.0 + (float)k * w_ratio;
				float im_y = (-0.5) + h_ratio / 2.0 + (float)j * h_ratio;
				im_x = fmin(im_x,im.w-1);
				im_y = fmin(im_y,im.h-1);
				set_pixel(resized_im, k, j , i, nn_interpolate(im, im_x, im_y, i));
			}
		}
	}
    return resized_im;
}

float bilinear_interpolate(image im, float x, float y, int c)
{
    // TODO
	int v1_x, v4_x, v1_y, v4_y;
	v1_x = floorf(x);
	v4_x = ceilf(x);
	v1_y = floorf(y);
	v4_y = ceilf(y);
	
	float d1, d2, d3, d4;
	d1 = x - v1_x;
	d2 = v4_x - x;
	d3 = y - v1_y;
	d4 = v4_y - y;
	
	float v1_val = get_pixel(im, v1_x, v1_y, c);
	float v2_val = get_pixel(im, v4_x, v1_y, c);
	float v3_val = get_pixel(im, v1_x, v4_y, c);
	float v4_val = get_pixel(im, v4_x, v4_y, c);
	
	float a1, a2, a3, a4;
	a1 = d2 * d4;
	a2 = d1 * d4;
	a3 = d2 * d3;
	a4 = d1 * d3;

	float p_val = a1 * v1_val + a2 * v2_val + a3 * v3_val + a4 * v4_val;

    return p_val;
}

image bilinear_resize(image im, int w, int h)
{
    // TODO
        float w_ratio = (float)im.w/(float)w;
	float h_ratio = (float)im.h/(float)h;	
	image bilinear_resized_im = make_image(w, h, im.c);
	// The actual co-ordinates of the starting points of edges of images begins at -0.5 + (width/height)ratio/2
	for(int i = 0; i < im.c; i++){
		for(int j = 0; j < h; j++){
			for(int k = 0; k < w; k++){
				float im_x = (-0.5) + w_ratio / 2.0 + (float)k * w_ratio;
				float im_y = (-0.5) + h_ratio / 2.0 + (float)j * h_ratio;             

				set_pixel(bilinear_resized_im, k, j , i, bilinear_interpolate(im, im_x, im_y, i));
			}
		}
	}
    return bilinear_resized_im;
}

