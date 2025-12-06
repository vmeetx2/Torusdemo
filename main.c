#include <stdio.h>
#include <string.h>
#include <math.h>

#define WIDTH 80 // terminal columns
#define HEIGHT 24 // terminal rows
#define FRAMESIZE (WIDTH * HEIGHT) // ++++total number of “pixels”
#define THETA_STEP 0.05f
#define PHI_STEP 0.03f
#define R_MAJOR 2.0f
#define R_MINOR 1.0f

const char *shades = ".,-~:;=!*#$@"; // brightness rampy for ASCII shading

int clamp(int val, int min, int max) {
    if (val < min) return min;
    if (val > max) return max;
    return val;
}

// projs the first 2d stuff, opens for zbuf
void project_point(float theta, float phi, float A, float B,
                   int *x_out, int *y_out, float *D_out) {

    // compute basic trigonometry
    float cos_theta = cosf(theta); // position around tube
    float sin_theta = sinf(theta);
    float cos_phi = cosf(phi); // position around ring
    float sin_phi = sinf(phi);
    float cos_A = cosf(A); // rotation around X
    float sin_A = sinf(A);
    float cos_B = cosf(B); // rotation around Z
    float sin_B = sinf(B);

    // torus loc
    float circle_x = R_MAJOR + R_MINOR * cos_theta; // x of tube 
    float circle_y = R_MINOR * sin_theta;           // y of tube 

    // rot. matrix hell
    float x3d = circle_x * (cos_B * cos_phi + sin_A * sin_B * sin_phi) - circle_y * cos_A * sin_B;
    float y3d = circle_x * (sin_B * cos_phi - sin_A * cos_B * sin_phi) + circle_y * cos_A * cos_B;
    float z3d = cos_A * circle_x * sin_phi + circle_y * sin_A + 5.0f; // push forward so donut is in front of “camera”

    // depth factor = 1 / z to simulate perspective
    *D_out = 1.0f / z3d;

    // project to screen coordinates, scale and center
    *x_out = (int)(WIDTH / 2 + x3d * (*D_out) * 30);
    *y_out = (int)(HEIGHT / 2 + y3d * (*D_out) * 15);
}

int main() {
    float A = 0.0f; // rotation angle X
    float B = 0.0f; // rotation angle Z

    while (1) {
        float zbuf[FRAMESIZE]; // hset val
        char frame[FRAMESIZE];  // hset val

        memset(frame, ' ', FRAMESIZE); // clr
        memset(zbuf, 0, sizeof(zbuf)); // clr

        // per frame rendering post clearing in while, SAMPLES 
        for (float theta = 0; theta < 2 * M_PI; theta += THETA_STEP) {
            for (float phi = 0; phi < 2 * M_PI; phi += PHI_STEP) {
                int x, y;
                float D;    

                // compute 3D point rotated and projected
                project_point(theta, phi, A, B, &x, &y, &D);

                // linear index in frame buffer
                int idx = clamp(x, 0, WIDTH - 1) + WIDTH * clamp(y, 0, HEIGHT - 1); //array ko lia as -1 for l

                // depth check: closer points overwrite farther ones
                if (D > zbuf[idx]) {
                    zbuf[idx] = D;
                    // convert depth to brightness index (0-11)
                    int shade_idx = clamp((int)(D * 12), 0, 11);
                    frame[idx] = shades[shade_idx]; // store character depth bai
                }
            }
        }

        printf("\x1b[H"); // move cursor to top-left (ANSI escape)

        // print frame buffer row by row
        for (int i = 0; i < FRAMESIZE; i++) {
            putchar((i % WIDTH) ? frame[i] : '\n'); 
        }

        A += 0.04f; 
        B += 0.02f; 
    }

    return 0;
}
