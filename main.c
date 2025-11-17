#include <stdio.h>
#include <math.h>
#include <string.h>
#ifdef _WIN32
#include <windows.h>
#define sleep_ms(x) Sleep(x)
#else
#include <unistd.h>
#define sleep_ms(x) usleep((x) * 1000)
#endif
#define WIDTH 80
#define HEIGHT 40
#define R1 1.0
#define R2 0.5
#define K2 5.0
#define THETA_STEP 0.07
#define PHI_STEP 0.02
#define FRAME_DELAY 30
#define BRIGHTNESS ".,-~:;=!*#$@"

void render_frame(double A, double B) {
    char screen[HEIGHT][WIDTH + 1];
    double zbuffer[HEIGHT][WIDTH];

    for (int i = 0; i < HEIGHT; i++) {
        for (int j = 0; j < WIDTH; j++) {
            screen[i][j] = ' ';
            zbuffer[i][j] = 0.0;
        }
        screen[i][WIDTH] = '\0';
    }

    double K1 = WIDTH * K2 * 3 / (8 * (R1 + R2));

    for (double theta = 0; theta < 2 * M_PI; theta += THETA_STEP) {
        double costheta = cos(theta);
        double sintheta = sin(theta);
        for (double phi = 0; phi < 2 * M_PI; phi += PHI_STEP) {
            double cosphi = cos(phi);
            double sinphi = sin(phi);

            double circleX = R1 + R2 * costheta;
            double circleY = R2 * sintheta;

            double x = circleX * (cosphi * cos(B) + sinphi * sin(A) * sin(B)) - circleY * cos(A) * sin(B);
            double y = circleX * (cosphi * sin(B) - sinphi * sin(A) * cos(B)) + circleY * cos(A) * cos(B);
            double z = circleX * sinphi * cos(A) + circleY * sin(A);

            double ooz = 1.0 / (z + K2);

            int xp = (int)(WIDTH / 2 + K1 * ooz * x);
            int yp = (int)(HEIGHT / 2 - K1 * ooz * y);

            double nx = cosphi * costheta * cos(B) - costheta * sinphi * sin(A) * sin(B) - sintheta * cos(A) * sin(B);
            double ny = cosphi * costheta * sin(B) - costheta * sinphi * sin(A) * cos(B) + sintheta * cos(A) * cos(B);
            double nz = costheta * sinphi * cos(A) + sintheta * sin(A);

            double L = nx * 0.0 + ny * 1.0 + nz * -1.0;

            if (L > 0 && xp >= 0 && xp < WIDTH && yp >= 0 && yp < HEIGHT) {
                if (ooz > zbuffer[yp][xp]) {
                    zbuffer[yp][xp] = ooz;
                    int index = (int)(L * (strlen(BRIGHTNESS) - 1));
                    if (index < 0) index = 0;
                    if (index > (int)strlen(BRIGHTNESS) - 1) index = (int)strlen(BRIGHTNESS) - 1;
                    screen[yp][xp] = BRIGHTNESS[index];
                }
            }
        }
    }

    printf("\x1b[H");
    for (int i = 0; i < HEIGHT; i++) {
        printf("%s\n", screen[i]);
    }
}

int main() {
    double A = 0.0;
    double B = 0.0;

    printf("\x1b[2J");

    while (1) {
        render_frame(A, B);
        A += 0.04;
        B += 0.02;
        sleep_ms(FRAME_DELAY);
    }

    return 0;
}
