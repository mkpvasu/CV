#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "image.h"


float get_pixel(image im, int x, int y, int c)
{
    int pix_value;
        if (c<0 || c>2 || x<0 || x>=im.w || y<0 || y>=im.h)
	{
        if (x<0){ x = 0; }
        if (x>=im.w){ x = im.w - 1; }
        if (y<0){ y = 0; }
        if (y>=im.h){ y = im.h - 1; }  
    	}

    	pix_value = (c*im.h*im.w) + (y*im.w) + x;

    return im.data[pix_value];
}

void set_pixel(image im, int x, int y, int c, float v)
{
        if (c<0 || c>2 || x<0 || x>=im.w || y<0 || y>=im.h){
	return;}
	else {
        im.data[(c*im.h*im.w) + (y*im.w) + x] = v;}
}

image copy_image(image im)
{
    int i;
    image copy = make_image(im.w, im.h, im.c);
    for (i = 0; i<(im.c*im.h*im.w); i++)
    {
        copy.data[i] = im.data[i];
    }
    return copy;
}

image rgb_to_grayscale(image im)
{
    assert(im.c == 3);
    image gray = make_image(im.w, im.h, 1);
    int i;
    for (i=0; i < im.w*im.h; i++)
    {
        gray.data[i] = 0.299 * im.data[i] + 0.587 * im.data [im.w*im.h + i] + 0.114 * im.data [2*im.w*im.h + i];
    }

    return gray;
}

void shift_image(image im, int c, float v)
{
    int i, k;
    for(i=0; i < im.w*im.h; i++)
    {
        k = (c*im.w*im.h) + i;
        im.data[k] = im.data[k] + v;
    }
}

void clamp_image(image im)
{
    int i;
    for (i = 0; i < im.c*im.h*im.w; ++i){
      im.data[i] = (im.data[i] > 1) ? 1 : ((im.data[i] < 0) ? 0 : im.data[i]);
    }
}


// These might be handy
float three_way_max(float a, float b, float c)
{
    return (a > b) ? ( (a > c) ? a : c) : ( (b > c) ? b : c) ;
}

float three_way_min(float a, float b, float c)
{
    return (a < b) ? ( (a < c) ? a : c) : ( (b < c) ? b : c) ;
}

void rgb_to_hsv(image im)
{
    // TODO Fill this in
    for(int i = 0; i < im.w; i ++) {
        for(int j = 0; j < im.h; j ++) {
            // get RGB values
            float Red = get_pixel(im, i, j, 0);
            float Gre = get_pixel(im, i, j, 1);
            float Blu = get_pixel(im, i, j, 2);

            // calculate the value of V
            float Val = three_way_max(Red, Gre, Blu);

            // calculate the value of S
            float m = three_way_min(Red, Gre, Blu);
            float C = Val - m;
            float Sat = 0;
            if(Val != 0) Sat = C / Val;

            // calculate the value of H'
            float H0 = 0;
            if(C != 0) {
                if(Val == Red) H0 = (Gre - Blu) / C;
                if(Val == Gre) H0 = (Blu - Red) / C + 2;
                if(Val == Blu) H0 = (Red - Gre) / C + 4;
            }

            // calculate the value of H
            float Hue;
            if(H0 < 0) Hue = H0 / 6 + 1;
            else Hue = H0 / 6.0;

            set_pixel(im, i, j, 0, Hue);
            set_pixel(im, i, j, 1, Sat);
            set_pixel(im, i, j, 2, Val);
        }
    }
}
void hsv_to_rgb(image im)
{
    for(int i = 0; i < im.w; i ++) {
        for(int j = 0; j < im.h; j ++) {
            float Hue = 360 * get_pixel(im, i, j, 0);
            float Sat = get_pixel(im, i, j, 1);
            float Val = get_pixel(im, i, j, 2);

            float C = Val * Sat;
            // X = C × (1 - |(H / 60°) mod 2 - 1|)
            float X = C * (1 - fabs(fmod(Hue / 60, 2) - 1));
            float m = Val - C;

            float Red = 0, Gre = 0, Blu = 0;

            if (Hue >= 0 && Hue < 60){
                Red = C;
                Gre = X;
                Blu = 0;
            } else if (Hue >= 60 && Hue < 120) {
                Red = X;
                Gre = C;
                Blu = 0;
            } else if (Hue >= 120 && Hue < 180) {
                Red = 0;
                Gre = C;
                Blu = X;
            } else if (Hue >= 180 && Hue < 240) {
                Red = 0;
                Gre = X;
                Blu = C;
            } else if (Hue >= 240 && Hue < 300) {
                Red = X;
                Gre = 0;
                Blu = C;
            } else if (Hue >= 300 && Hue < 360) {
                Red = C;
                Gre = 0;
                Blu = X;
            }

            set_pixel(im, i, j, 0, Red + m);
            set_pixel(im, i, j, 1, Gre + m);
            set_pixel(im, i, j, 2, Blu + m);
        }
    }
}

void scale_image(image im, int c, float v) {
    for(int i = 0; i < im.w; i ++) {
        for(int j = 0; j < im.h; j ++) {
            set_pixel(im, i, j, c, v * get_pixel(im, i, j, c));
        }
    }
}

void rgb_to_hcl(image im)
{
    double tm_xyz[3][3] = {{0.4124, 0.3576, 0.1805}, {0.2126, 0.7152, 0.0722}, {0.0193, 0.1192, 0.9505}};
    float var_U, var_V,var_Y,cie_l,cie_u,cie_v;

    // Converting RGB to sRGB
    // Channel 1 will store sR values
    // Channel 2 will store sG values
    // Channel 3 will store sB values
    for (int i = 0; i < im.c*im.w*im.h; i++)
    {
        im.data[i] = im.data[i]/255;

        if (im.data[i] <= 0.04045)
        {
            im.data[i] = im.data[i]/12.92;
        }
        else
        {
            im.data[i] = pow(((im.data[i] + 0.055)/1.055),2.5);
        }
    }

    // Converting sRGB to CIEXYZ
    // Channel 1 will store X values
    // Channel 2 will store Y values
    // Channel 3 will store Z values
    for (int i = 0; i < im.w*im.h; i++)
    {
        double a[3][1] = {{im.data[i]},{im.data[im.w*im.h + i]},{im.data[2*im.w*im.h + i]}};
        double result[3][1];
        for (int m = 0; m < 1; m++) {
            for (int j = 0; j < 3; m++) {
                for (int k = 0; k < 3; k++) {
                    result[m][j] += tm_xyz[m][k] * a[k][j];
                }
            }
        }
        im.data[i] = result[0][0];
        im.data[im.w*im.h + i] = result[1][0];
        im.data[2*im.w*im.h + i] = result[2][0];
    }

    // Converting CIEXYZ to CIELUV
    // Channel 1 will store Luminance values
    // Channel 2 will store U values
    // Channel 3 will store V values
    for (int i = 0; i < im.w*im.h; i++)
    {
        var_U = (4 * im.data[i]) / (im.data[i] + (15 * im.data[im.w*im.h + i]) + (3 * im.data[2*im.w*im.h + i]));
        var_V = (9 * im.data[im.w*im.h + i]) / (im.data[i] + (15 * im.data[im.w*im.h + i]) + (3 * im.data[2*im.w*im.h + i]));

        var_Y =  im.data[im.w*im.h + i]/100;

        if ( var_Y > 0.008856 ) 
        {
            cie_l = (116 * pow(var_Y,(1/3))) - 16;
        }
        else
        {
            cie_l = 903.2963 * var_Y;
        }

        cie_u = 13 * cie_l* (var_U - 0.2009);
        cie_v = 13 * cie_l* (var_V - 0.4610);

        im.data[i] = cie_l;
        im.data[im.w*im.h + i] = cie_u;
        im.data[2*im.w*im.h + i] = cie_v;
    }

    // Converting CIELUV to HCL
    // Channel 1 will store Hue values
    // Channel 2 will store Chroma values
    // Channel 3 will store Luminance values
    for (int i = 0; i < im.w*im.h; i++)
    {   
        float c1 = (pow(im.data[im.w*im.h + i],2) + pow(im.data[2*im.w*im.h + i],2));
        float temp = im.data[i];
        im.data[i] = atan2(im.data[2*im.w*im.h +i], im.data[im.w*im.h + i]);
        im.data[im.w*im.h + i] = sqrt(c1);
        im.data[2*im.w*im.h + i] = temp;
    }
}
