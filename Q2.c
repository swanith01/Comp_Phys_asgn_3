#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>
#include <complex.h>

#define N 64
#define XMIN -50
#define XMAX 50

// Define the function f(x)
double f(double x) {
    if (x == 0) {
        return 1.0; // Handle singularity at x = 0
    }
    return sin(x) / x;
}

// Analytical Fourier Transform of sinc function
double analytical_fourier_transform(double k) {
    return sqrt(M_PI / 2);
}

int main() {
    double delta = (XMAX - XMIN) / (double)N;

    // Generate x values at regular intervals
    double x_values[N];
    for (int i = 0; i < N; i++) {
        x_values[i] = XMIN + i * delta;
    }

    // Calculate y values using the function f(x)
    double y_values[N];
    for (int i = 0; i < N; i++) {
        y_values[i] = f(x_values[i]);
    }

    // Calculate the DFT using FFTW
    fftw_complex *in, *out;
    fftw_plan plan;
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    for (int i = 0; i < N; i++) {
        in[i][0] = y_values[i]; // Real part
        in[i][1] = 0;           // Imaginary part
    }
    plan = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);

    // Print the DFT values
    printf("DFT Values:\n");
    for (int i = 0; i < N; i++) {
        printf("%f + %fi\n", out[i][0], out[i][1]);
    }

    // Calculate analytical Fourier transform
    printf("Analytical Fourier Transform of sinc function:\n");
    for (int i = 0; i < N; i++) {
        double k = 2 * M_PI * i / (delta * N);
        printf("%f\n", analytical_fourier_transform(k));
    }

    // Free memory
    fftw_free(in);
    fftw_free(out);

    return 0;
}
/*
DFT Values:
1.991343 + 0.000000i
-2.029988 + 0.000000i
1.990967 + -0.000000i
-2.030759 + -0.000000i
1.989763 + -0.000000i
-2.032463 + -0.000000i
1.987458 + 0.000000i
-2.035521 + 0.000000i
1.983415 + 0.000000i
-2.040912 + -0.000000i
1.976081 + -0.000000i
-2.051225 + -0.000000i
1.960819 + 0.000000i
-2.075688 + -0.000000i
1.915966 + -0.000000i
-2.183298 + -0.000000i
0.836066 + 0.000000i
0.173488 + -0.000000i
0.094857 + 0.000000i
0.065160 + -0.000000i
0.049854 + 0.000000i
0.040641 + 0.000000i
0.034564 + 0.000000i
0.030312 + -0.000000i
0.027220 + 0.000000i
0.024914 + 0.000000i
0.023173 + 0.000000i
0.021854 + -0.000000i
0.020866 + -0.000000i
0.020148 + -0.000000i
0.019660 + 0.000000i
0.019377 + -0.000000i
0.019284 + 0.000000i
0.019377 + 0.000000i
0.019660 + -0.000000i
0.020148 + 0.000000i
0.020866 + -0.000000i
0.021854 + 0.000000i
0.023173 + 0.000000i
0.024914 + 0.000000i
0.027220 + 0.000000i
0.030312 + 0.000000i
0.034564 + 0.000000i
0.040641 + 0.000000i
0.049854 + 0.000000i
0.065160 + 0.000000i
0.094857 + -0.000000i
0.173488 + 0.000000i
0.836066 + 0.000000i
-2.183298 + -0.000000i
1.915966 + -0.000000i
-2.075688 + 0.000000i
1.960819 + 0.000000i
-2.051225 + 0.000000i
1.976081 + 0.000000i
-2.040912 + 0.000000i
1.983415 + 0.000000i
-2.035521 + -0.000000i
1.987458 + -0.000000i
-2.032463 + 0.000000i
1.989763 + -0.000000i
-2.030759 + 0.000000i
1.990967 + 0.000000i
-2.029988 + 0.000000i
Analytical Fourier Transform of sinc function:
1.253314
1.253314
1.253314
1.253314
1.253314
1.253314
1.253314
1.253314
1.253314
1.253314
1.253314
1.253314
1.253314
1.253314
1.253314
1.253314
1.253314
1.253314
1.253314
1.253314
1.253314
1.253314
1.253314
1.253314
1.253314
1.253314
1.253314
1.253314
1.253314
1.253314
1.253314
1.253314
1.253314
1.253314
1.253314
1.253314
1.253314
1.253314
1.253314
1.253314
1.253314
1.253314
1.253314
1.253314
1.253314
1.253314
1.253314
1.253314
1.253314
1.253314
1.253314
1.253314
1.253314
1.253314
1.253314
1.253314
1.253314
1.253314
1.253314
1.253314
1.253314
1.253314
1.253314
1.253314

*/
