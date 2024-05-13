#include <stdio.h>
#include <math.h>
#include <gsl/gsl_fft_complex.h>

#define N 64 // Number of points in the signal
#define PI 3.14159265358979323846

double sinc(double x) {
    if (x == 0.0) {
        return 1.0;
    }
    return sin(PI * x) / (PI * x);
}

int main() {
    // Define the input signal (sinc function)
    double signal[N][2]; // 2 columns for real and imaginary parts

    // Fill the signal array with the sinc function values
    for (int i = 0; i < N; i++) {
        double x = (i - N / 2.0) / (double)N;
        signal[i][0] = sinc(x); // Real part
        signal[i][1] = 0.0; // Imaginary part
    }

    // Print the input signal
    printf("Input Signal:\n");
    for (int i = 0; i < N; i++) {
        printf("%.6f + %.6fi\n", signal[i][0], signal[i][1]);
    }

    // Initialize the workspace and wavetable for FFT
    gsl_fft_complex_wavetable *wavetable;
    gsl_fft_complex_workspace *workspace;
    wavetable = gsl_fft_complex_wavetable_alloc(N);
    workspace = gsl_fft_complex_workspace_alloc(N);

    // Perform forward FFT (transforming time domain signal to frequency domain)
    gsl_fft_complex_forward(signal[0], 1, N, wavetable, workspace);

    // Compute the normalization factor
    double normalization_factor = 1.0 / sqrt(2 * PI / N);

    // Print the normalized DFT values
    printf("Normalized DFT Values:\n");
    for (int i = 0; i < N; i++) {
        double real = signal[i][0] * normalization_factor;
        double imag = signal[i][1] * normalization_factor;
        printf("%.6f + %.6fi\n", real, imag);
    }

    // Print the analytical Fourier transform of sinc function
    printf("Analytical Fourier Transform of sinc function:\n");
    for (int i = 0; i < N; i++) {
        double k = (i - N / 2.0) * 2 * PI / N;
        if (fabs(k) <= PI) {
            printf("%.6f\n", normalization_factor);
        } else {
            printf("0.000000\n");
        }
    }

    // Free the allocated memory
    gsl_fft_complex_wavetable_free(wavetable);
    gsl_fft_complex_workspace_free(workspace);

    return 0;
}

