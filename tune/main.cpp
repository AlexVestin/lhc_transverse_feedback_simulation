#include <iostream>
extern "C" {
    #include <fftw3.h>
}


#define N 512
int main() {
    fftw_complex *in, *out;
    fftw_plan p;
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p); 
    fftw_destroy_plan(p);  
    fftw_free(in); 
    fftw_free(out);

    std::cout << "hello" << std::endl;
    return 0;   
}