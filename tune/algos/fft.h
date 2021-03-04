#ifndef __FFT_H__
#define __FFT_H__
#include <inttypes.h>
#include <complex.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#define pi 3.14159265359

double max_fft_frequency(double _Complex* signal, size_t N);
#endif