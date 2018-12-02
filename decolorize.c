// convert the input RGB colored image into grayscale
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include "imcore.h"
#include "prcore.h"

#define RANDOM_SAMPLE_COUNT   1000000

// compute the signed distance between the two color according to
// Color2Gray: Salience-Preserving Color Removal
float color_distance(struct color_t c1, struct color_t c2, float theta, float alpha)
{
    float L1, a1, b1;
    float L2, a2, b2;

    // Observer. = 2°, Illuminant = D65
    rgb2lab(c1, &L1, &a1, &b1);
    rgb2lab(c2, &L2, &a2, &b2);

#define crunch(x) (alpha * tanhf((x) / alpha))

    float lDiff = fabsf(L1 - L2);
    float crDiff = square(a1-a2) + square(b1-b2);

    // if the luminance difference is much higher than the chrominance, use the Luminance difference
    if(lDiff > crunch(crDiff))
    {
        return lDiff;
    }

    // find the angle between the dA and dB
    float diffSign = (a1 - a2) * cosf(theta) + (b1 - b2) * sinf(theta);

    if(diffSign > 0)
    {
        return crunch(crDiff);
    }
    else
    {
        return crunch(-crDiff);
    }
#undef crunch
}

// http://www.cescript.com/2015/10/imge-renksizlestirme-image.html
void decolorize(matrix_t *in, uint32_t sample_count, matrix_t *out)
{
    // for loop iterators
    uint32_t i, j;

    // allocate out before use it
    matrix_resize(out, height(in), width(in), 1);

    // TODO: create in_data pointer based on the input type
    uint8_t *in_data  = data(uint8_t, in);
    // create two color_t and distance arrays
    struct color_t *C1 = array_create(struct color_t, sample_count);
    struct color_t *C2 = array_create(struct color_t, sample_count);
    float *distance    = array_create(float,          sample_count);

    // select #sample_count random pixel difference
    for(i=0; i < sample_count; i++) {

        // select two random pixels inside image
        uint32_t x1 = random_int(0, width(in) - 1);
        uint32_t x2 = random_int(0, width(in) - 1);
        uint32_t y1 = random_int(0, height(in) - 1);
        uint32_t y2 = random_int(0, height(in) - 1);

        // compute the index of the pixel positions
        uint32_t idx1 = idx(in, y1, x1, 0);
        uint32_t idx2 = idx(in, y2, x2, 0);

        // get the RGB color of the input image
        C1[i] = RGB(in_data[idx1+2], in_data[idx1+1], in_data[idx1+0]);
        C2[i] = RGB(in_data[idx2+2], in_data[idx2+1], in_data[idx2+0]);

        // get the distance between the two colors
        distance[i] = color_distance(C1[i], C2[i], 3.14159f / 4.0f, 15);
    }

    // create a color weight table
    uint32_t table_length = 0;
    int32_t c1 = 0; int32_t c2 = 0;;
    float W[3][66] = {0.0f};

    for(c1 = 0; c1 < 11; c1++) {
        // sum of the weights must be less than 1.0
        for(c2 = 10 - c1; c2 >= 0; c2--) {
            // set the weights
            W[0][table_length] = c1 / 10.0f;
            W[1][table_length] = c2 / 10.0f;
            W[2][table_length] = (10 - c1 - c2) / 10.0f;
            // increase the table length
            table_length++;
        }
    }
    // find the best coefficients
    double minEg = INFINITY;
    float WR = 0.29f; float WG = 0.58f; float WB = 0.11f;

    for(j = 0; j < table_length; j++) {

        float wr = W[0][j];
        float wg = W[1][j];
        float wb = W[2][j];

        double Eg = 0;
        for(i = 0; i < sample_count; i++) {
            // find the projection
            float dX = wr*C1[i].red + wg*C1[i].green + wb*C1[i].blue;
            float dY = wr*C2[i].red + wg*C2[i].green + wb*C2[i].blue;
            // sum the error
            Eg += square( (dX - dY) - distance[i] );
        }
        // if the current error is the minimum save it
        if(Eg < minEg) {
            minEg = Eg;
            WR = wr;
            WG = wg;
            WB = wb;
        }
    }
    // print the found weights and errors
    printf("E[%3.2f %3.2f %3.2f]: %3.5f\n", WR, WG, WB, sqrt(minEg) / sample_count);

    // now convert the rgb image into grayscale using the weight minIdx
    uint8_t *out_data = data(uint8_t, out);
    //#pragma omp parallel for
    for(i = 0; i < width(in)*height(in); i++) {
        out_data[i] = (uint8_t) clamp( WB*in_data[3*i+0] + WG*in_data[3*i+1] + WR*in_data[3*i+2], 0,255);
    }

    // free up teh space
    array_free(C1);
    array_free(C2);
    array_free(distance);
}


int main(int argc, unsigned char *argv[]) {

    // read the test image
    unsigned char filename[256] = "..//data//impression_sunrise.bmp";
    if(argc > 1) {
        strncpy(filename, argv[1], 256);
    }
    // read the input image
    matrix_t *image = imread(filename);

    // create the gray images
    matrix_t *gray_image = matrix_create(uint8_t);
    matrix_t *decolorized_image = matrix_create(uint8_t);

    // convert the input into grascale
    rgb2gray(image, gray_image);
    // decolorize the image
    decolorize(image, RANDOM_SAMPLE_COUNT, decolorized_image);

    imwrite(gray_image, "gray_image.bmp");
    imwrite(decolorized_image, "decolorized_image.bmp");

    return 0;
}