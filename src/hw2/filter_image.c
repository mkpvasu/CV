#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"
#define TWOPI 6.2831853

void l1_normalize(image im)
{
     float sum = 0.0;
    int i, j, k;

    // calculate sum
    for (i = 0; i < im.w; i++){
        for (j = 0; j < im.h; j++){
            for (k = 0; k < im.c; k++){
                sum += (get_pixel(im, i, j, k));
            }
        }
    }

    // to avoid divide by zero
    if (sum == 0) { return; }

    // divide each pixel value by the sum
    float val = 0.0;
    for (i = 0; i < im.w; i++){
        for (j = 0; j < im.h; j++){
            for (k = 0; k < im.c; k++){
                val = get_pixel(im, i, j, k);
                set_pixel(im, i, j, k, val / sum);
            }
        }
    }
}

image make_box_filter(int w)
{
    image filter = make_image(w, w, 1);
    for(int i = 0; i < w * w; i++){
        filter.data[i] = 1.0/(w*w);
    }
    return filter;
}

image convolve_image(image im, image filter, int preserve)
{
    assert(im.c == filter.c || filter.c == 1);
    int pad_w = (int)((filter.w -1)/ 2);
    int pad_h = (int)((filter.h -1)/ 2);
    image convolved = make_image(im.w,im.h,preserve ? im.c:1);

    if(im.c == filter.c){
    if(preserve == 1){
        for(int k = 0; k < im.c; k++){
            for(int j = 0 - pad_h; j < im.h - pad_h; j++){
                for(int i = 0 - pad_w; i < im.w - pad_w; i++){
                    float value = 0;

                    for(int fy = 0; fy < filter.h; fy++){
                        for(int fx = 0; fx < filter.w; fx++){
                            value += get_pixel(im, i+fx, j+fy, k) * get_pixel(filter, fx, fy, k);
                        }
                    }
                    set_pixel(convolved, i+pad_w, j+pad_h, k, value);
                }
            }
        }
    } 
    
    if (preserve != 1){
        for(int j = 0 - pad_h; j < im.h - pad_h; j++){
             for(int i = 0 - pad_w; i < im.w - pad_w; i++){
                float value = 0;
                
                for(int fy = 0; fy < filter.h; fy++){
                        for(int fx = 0; fx < filter.w; fx++){
			    for(int k = 0; k<im.c; k++){
                            value += get_pixel(im, i+fx, j+fy, k) * get_pixel(filter, fx, fy, k);
			    }
                        }
                    }
                set_pixel(convolved, i+pad_w, j+pad_h, 0, value);
            }
        }
    }
    }

    if(filter.c == 1){
    if(preserve == 1){
        for(int k = 0; k < im.c; k++){
            for(int j = 0 - pad_h; j < im.h - pad_h; j++){
                for(int i = 0 - pad_w; i < im.w - pad_w; i++){
                    float value = 0;

                    for(int fy = 0; fy < filter.h; fy++){
                        for(int fx = 0; fx < filter.w; fx++){
                            value += get_pixel(im, i+fx, j+fy, k) * get_pixel(filter, fx, fy, 0);
                        }
                    }
                    set_pixel(convolved, i+pad_w, j+pad_h, k, value);
                }
            }
        }
    } 
    
    if (preserve != 1){
        for(int j = 0 - pad_h; j < im.h - pad_h; j++){
             for(int i = 0 - pad_w; i < im.w - pad_w; i++){
                float value = 0;
                
                for(int fy = 0; fy < filter.h; fy++){
                        for(int fx = 0; fx < filter.w; fx++){
			    for(int k = 0; k<im.c; k++){
                            value += get_pixel(im, i+fx, j+fy, k) * get_pixel(filter, fx, fy, 0);
			    }
                        }
                    }
                set_pixel(convolved, i+pad_w, j+pad_h, 0, value);
            }
        }
    }
    }

return convolved;
}

image make_highpass_filter()
{
    image f = make_image(3,3,1);
    set_pixel(f,0,0,0,0.0);
    set_pixel(f,1,0,0,-1.0);
    set_pixel(f,2,0,0,0.0);
    set_pixel(f,0,1,0,-1.0);
    set_pixel(f,1,1,0,4.0);
    set_pixel(f,2,1,0,-1.0);
    set_pixel(f,0,2,0,0.0);
    set_pixel(f,1,2,0,-1.0);
    set_pixel(f,2,2,0,0.0);

    return f;
}

image make_sharpen_filter()
{
    image f = make_image(3,3,1);
    set_pixel(f,0,0,0,0.0);
    set_pixel(f,1,0,0,-1.0);
    set_pixel(f,2,0,0,0.0);
    set_pixel(f,0,1,0,-1.0);
    set_pixel(f,1,1,0,5.0);
    set_pixel(f,2,1,0,-1.0);
    set_pixel(f,0,2,0,0.0);
    set_pixel(f,1,2,0,-1.0);
    set_pixel(f,2,2,0,0.0);

    return f;
}

image make_emboss_filter()
{
    image f = make_image(3,3,1);
    set_pixel(f,0,0,0,-2.0);
    set_pixel(f,1,0,0,-1.0);
    set_pixel(f,2,0,0,0.0);
    set_pixel(f,0,1,0,-1.0);
    set_pixel(f,1,1,0,1.0);
    set_pixel(f,2,1,0,1.0);
    set_pixel(f,0,2,0,0.0);
    set_pixel(f,1,2,0,1.0);
    set_pixel(f,2,2,0,2.0);

    return f;
}

// Question 2.2.1: Which of these filters should we use preserve when we run our convolution and which ones should we not? Why?
// Answer: Box_filter should have preserve set to '1', because we are performing smoothing operation for the image.
//         Highpass_filter gives a '0' sum of weights, while sharpen_filter and emboss_filter gives a '1'
//         Highpass_filter has '0' sum of weights, the preserve should be set to '0'.
//         Sharpen_filter and Emboss_filter produces stylish images when preserve is set to '1' and if preserve is set to '0' it gives a normal image.


// Question 2.2.2: Do we have to do any post-processing for the above filters? Which ones and why?
// Answer: Highpass_filter, sharpen_filter or emboss_filter needs post-processing (clamp_image() or feature_normalize()) after using them to
//         ensure all the pixel values lies between 0 and 1. No post-processing is required for box_filter, because all the pixel
//         values are already between 0 and 1.

image make_gaussian_filter(float sigma)
{
    // TODO
    	int i, j, offset;
    	int w = 6 * sigma;
    	w = w > 0 ? ((w % 2 == 1) ? w : w + 1) : 1;
    	image gaussian_filter = make_image(w, w, 1);
    	offset = w / 2;
    	for (i = 0 - offset; i < w - offset; i++){
      	for (j = 0 - offset; j < w - offset; j++){
        float v = 1 / (TWOPI*sigma*sigma) * expf((- j*j - i*i) / (2*sigma*sigma));
        set_pixel(gaussian_filter, j+offset, i+offset, 0, v);
      	}
    	}
	l1_normalize(gaussian_filter);
    return gaussian_filter;

}

image add_image(image a, image b)
{
    // TODO
	image im_add = make_image(a.w,a.h,a.c);	
	assert(a.w == b.w && a.h == b.h && a.c == b.c);
		
		// Adding image a's info with image b's info to get a new image
		for(int i = 0; i<a.c*a.w*a.h; i++){
		im_add.data[i] = a.data[i]+b.data[i];
		}
    return im_add;
}

image sub_image(image a, image b)
{
    // TODO
	image im_sub = make_image(a.w,a.h,a.c);	
	assert(a.w == b.w && a.h == b.h && a.c == b.c);
		
		// Subtracting image a's info with image b's info to get a new image
		for(int i = 0; i<a.c*a.w*a.h; i++){
		im_sub.data[i] = a.data[i]-b.data[i];
		}
    return im_sub;
}

image make_gx_filter()
{
    // TODO
	image gx_filter = make_image(3,3,1);
	set_pixel(gx_filter, 0, 0, 0, -1.0);
	set_pixel(gx_filter, 1, 0, 0, 0.0);
	set_pixel(gx_filter, 2, 0, 0, 1.0);
	set_pixel(gx_filter, 0, 1, 0, -2.0);
	set_pixel(gx_filter, 1, 1, 0, 0.0);
	set_pixel(gx_filter, 2, 1, 0, 2.0);
	set_pixel(gx_filter, 0, 2, 0, -1.0);
	set_pixel(gx_filter, 1, 2, 0, 0.0);
	set_pixel(gx_filter, 2, 2, 0, 1.0);
    return gx_filter;
}

image make_gy_filter()
{
    // TODO
	image gy_filter = make_image(3,3,1);
	set_pixel(gy_filter, 0, 0, 0, -1.0);
	set_pixel(gy_filter, 1, 0, 0, -2.0);
	set_pixel(gy_filter, 2, 0, 0, -1.0);
	set_pixel(gy_filter, 0, 1, 0, 0.0);
	set_pixel(gy_filter, 1, 1, 0, 0.0);
	set_pixel(gy_filter, 2, 1, 0, 0.0);
	set_pixel(gy_filter, 0, 2, 0, 1.0);
	set_pixel(gy_filter, 1, 2, 0, 2.0);
	set_pixel(gy_filter, 2, 2, 0, 1.0);
    return gy_filter;
}

void feature_normalize(image im)
{
    // TODO

	// Calculating the min and max value in the image to find the range of values
	float min = 1.0;
	float max = 0.0;
	
	float range;
    for (int i = 0; i < im.c*im.h*im.w; ++i){
      if (im.data[i] > max) max = im.data[i];
      if (im.data[i] < min) min = im.data[i];
    }
    range = max - min;
    if (!range){
      for (int i = 0; i < im.c*im.h*im.w; ++i){
        im.data[i] = 0;
      }
    } else{
      for (int i = 0; i < im.c*im.h*im.w; ++i){
        im.data[i] = (im.data[i] - min) / range;
      }
    }
			
}

image *sobel_image(image im)
{
    // TODO

    image* images = calloc(2, sizeof(image));

    // calculate Gx, and Gy
    image gx_filter = make_gx_filter();
    image gy_filter = make_gy_filter();
    image gx = convolve_image(im, gx_filter, 0);
    image gy = convolve_image(im, gy_filter, 0);
    image mag = make_image(gx.w, gx.h, 1);
    image dir = make_image(gx.w, gx.h, 1);

    // calculate gradient and direction
     for (int i = 0; i < im.h*im.w; i++){
      mag.data[i] = sqrtf(gx.data[i]*gx.data[i] + gy.data[i]*gy.data[i]);
	dir.data[i] = atan2f(gy.data[i], gx.data[i]);
    }

    images[0] = mag;
    images[1] = dir;

return images;
}

image colorize_sobel(image im)
{
    // TODO
	image color_sobel = make_image(im.w,im.h,im.c);

	// Assigning the dir_sobel and mag_sobel images
	image *sobel_return = sobel_image(color_sobel);
    	feature_normalize(sobel_return[0]);
    	feature_normalize(sobel_return[1]);
    	memcpy(color_sobel.data, sobel_return[1].data, im.w*im.h*sizeof(float));
    	memcpy(color_sobel.data+im.h*im.w, sobel_return[0].data, im.w*im.h*sizeof(float));
    	memcpy(color_sobel.data+2*im.h*im.w, sobel_return[0].data, im.w*im.h*sizeof(float));
    	hsv_to_rgb(color_sobel);
	
	hsv_to_rgb(color_sobel);
	
    return color_sobel;
}

