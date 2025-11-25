#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

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
#define THETA_STEP 0.02
#define PHI_STEP 0.01
#define FRAME_DELAY_MS 30
#define BRIGHTNESS_CHARS ".,-~:;=!*#$@"

#define LIGHT_X 0.0
#define LIGHT_Y 1.0
#define LIGHT_Z -1.0

void enable_vt100_on_windows(void) {
#ifdef _WIN32
    HANDLE hOut = GetStdHandle(STD_OUTPUT_HANDLE);
    DWORD dwMode = 0;
    GetConsoleMode(hOut, &dwMode);
    dwMode |= ENABLE_VIRTUAL_TERMINAL_PROCESSING;
    SetConsoleMode(hOut, dwMode);
#endif
}

static inline double clamp(double val, double low, double high) {
    return val < low ? low : (val > high ? high : val);
}

static inline int iround(double x) {
    return (int)(x >= 0 ? x + 0.5 : x - 0.5);
}

static char screen[HEIGHT][WIDTH + 1];
static double zbuffer[HEIGHT][WIDTH];

#define THETA_STEPS ((int)(2 * M_PI / THETA_STEP) + 1)
#define PHI_STEPS   ((int)(2 * M_PI / PHI_STEP) + 1)
static double cos_theta[THETA_STEPS], sin_theta[THETA_STEPS];
static double cos_phi[PHI_STEPS], sin_phi[PHI_STEPS];
static int theta_steps = 0, phi_steps = 0;

void precompute_trig_tables(void) {
    if (theta_steps > 0) return;

    theta_steps = 0;
    for (double t = 0; t < 2 * M_PI; t += THETA_STEP) {
        cos_theta[theta_steps] = cos(t);
        sin_theta[theta_steps] = sin(t);
        theta_steps++;
    }

    phi_steps = 0;
    for (double p = 0; p < 2 * M_PI; p += PHI_STEP) {
        cos_phi[phi_steps] = cos(p);
        sin_phi[phi_steps] = sin(p);
        phi_steps++;
    }
}

void clear_buffers(void) {
    for (int i = 0; i < HEIGHT; i++) {
        for (int j = 0; j < WIDTH; j++) {
            screen[i][j] = ' ';
            zbuffer[i][j] = 0.0;
        }
        screen[i][WIDTH] = '\0';
    }
}

void torus_point_normal(int ti, int pi, double *x, double *y, double *z,
                        double *nx, double *ny, double *nz) {
    double costheta = cos_theta[ti];
    double sintheta = sin_theta[ti];
    double cosphi = cos_phi[pi];
    double sinphi = sin_phi[pi];

    double circleX = R1 + R2 * costheta;
    double circleY = R2 * sintheta;

    *x = circleX * cosphi;
    *y = circleY;
    *z = circleX * sinphi;

    *nx = cosphi * costheta;
    *ny = sintheta;
    *nz = sinphi * costheta;
}

void rotate_and_project(double x, double y, double z, double A, double B,
                        double *xp_out, double *yp_out, double *ooz_out) {
    double x_rot = x * cos(B) + z * sin(B);
    double z_rot = -x * sin(B) + z * cos(B);

    double y_rot = y * cos(A) - z_rot * sin(A);
    double z_final = y * sin(A) + z_rot * cos(A);

    double ooz = 1.0 / (z_final + K2);
    double K1 = WIDTH * K2 * 3.0 / (8.0 * (R1 + R2));

    *xp_out = (double)WIDTH / 2.0 + K1 * ooz * x_rot;
    *yp_out = (double)HEIGHT / 2.0 - K1 * ooz * y_rot;
    *ooz_out = ooz;
}

void rotate_normal(double nx, double ny, double nz, double A, double B,
                   double *nx_out, double *ny_out, double *nz_out) {
    double nx1 = nx * cos(B) - nz * sin(B);
    double nz1 = nx * sin(B) + nz * cos(B);

    double ny1 = ny * cos(A) - nz1 * sin(A);
    double nz2 = ny * sin(A) + nz1 * cos(A);

    *nx_out = nx1;
    *ny_out = ny1;
    *nz_out = nz2;
}

double compute_lighting(double nx, double ny, double nz, double lx, double ly, double lz) {
    static int light_normalized = 0;
    static double lnx, lny, lnz;
    if (!light_normalized) {
        double len = sqrt(lx*lx + ly*ly + lz*lz);
        lnx = (len > 0) ? lx / len : 0;
        lny = (len > 0) ? ly / len : 0;
        lnz = (len > 0) ? lz / len : 0;
        light_normalized = 1;
    }

    double dot = nx * lnx + ny * lny + nz * lnz;
    return clamp(dot, 0.0, 1.0);
}

void render_frame(double A, double B) {
    clear_buffers();
    precompute_trig_tables();

    double light_dir[3] = {LIGHT_X, LIGHT_Y, LIGHT_Z};

    for (int ti = 0; ti < theta_steps; ti++) {
        for (int pi = 0; pi < phi_steps; pi++) {
            double x, y, z, nx, ny, nz;
            torus_point_normal(ti, pi, &x, &y, &z, &nx, &ny, &nz);

            double xp_f, yp_f, ooz;
            rotate_and_project(x, y, z, A, B, &xp_f, &yp_f, &ooz);

            double rx, ry, rz;
            rotate_normal(nx, ny, nz, A, B, &rx, &ry, &rz);

            double L = compute_lighting(rx, ry, rz, light_dir[0], light_dir[1], light_dir[2]);

            if (L <= 0) continue;
            if (ooz <= 0) continue;

            int xp = iround(xp_f);
            int yp = iround(yp_f);

            if (xp < 0 || xp >= WIDTH || yp < 0 || yp >= HEIGHT) continue;

            if (ooz > zbuffer[yp][xp]) {
                zbuffer[yp][xp] = ooz;

                const char* charset = BRIGHTNESS_CHARS;
                int nchars = (int)strlen(charset);
                int final_idx = (int)(L * (nchars - 1) + 0.5);
                final_idx = clamp(final_idx, 0, nchars - 1);

                screen[yp][xp] = charset[final_idx];
            }
        }
    }

    printf("\x1b[H");
    for (int i = 0; i < HEIGHT; i++) {
        printf("%s\n", screen[i]);
    }
}

int main(void) {
    enable_vt100_on_windows();
    printf("\x1b[2J");

    double A = 0.0, B = 0.0;
    clock_t last_time = clock();

    while (1) {
        clock_t now = clock();
        double delta = (double)(now - last_time) / CLOCKS_PER_SEC;
        last_time = now;

        A += 1.5 * delta;
        B += 0.75 * delta;

        render_frame(A, B);
        sleep_ms(FRAME_DELAY_MS);
    }

    return 0;
}

    return 0;
}
