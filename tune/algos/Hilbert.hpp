#include <vector>
#include "Algorithm.hpp"
#include <iostream>
#include <math.h>
extern "C" {
    #include <fftw3.h>
}

// https://stackoverflow.com/questions/24518989/how-to-perform-1-dimensional-valid-convolution
template<typename T>
std::vector<T> conv(std::vector<T> const &f, std::vector<T> const &g) {
  int const nf = f.size();
  int const ng = g.size();
  int const n  = nf + ng - 1;
  std::vector<T> out(n, T());
  for(auto i(0); i < n; ++i) {
    int const jmn = (i >= ng - 1)? i - (ng - 1) : 0;
    int const jmx = (i <  nf - 1)? i            : nf - 1;
    for(auto j(jmn); j <= jmx; ++j) {
      out[i] += (f[j] * g[i - j]);
    }
  }
  return out; 
}

class Hilbert {
    public:
        Hilbert(int size) {     
              
            in = (double*)fftw_malloc(size * sizeof(double));
            out = (fftw_complex*)fftw_malloc(size * sizeof(fftw_complex)); 
            planForward = fftw_plan_dft_r2c_1d(size, in, out, FFTW_ESTIMATE);
            //planInverse = fftw_plan_dft_1d(size, in, out, FFTW_COMPLEX)
        };

        ~Hilbert() {
       
            // TODO fix segfautl when destroying plan
            //fftw_destroy_plan(p);

            fftw_free(in);
            fftw_free(out);
            
        }

        float performAnalysis(const std::vector<double>& data) {
            std::vector<double> I = conv(data, iKernel);
            std::vector<double> Q = conv(data, qKernel);
            double totalPhase = 0;

            const double PI = 3.14159265359;
            std::vector<double> angle(I.size());           
            for(int i = 0; i < I.size(); i++) {
                angle[i] = std::remainder(std::atan(I[i] / Q[i]), 2.0 * M_PI);
            }

            // Diff phase 
            std::vector<double> diff(I.size());
            for(int i = 0; i < I.size() - 1; i++) {
                diff[i] = (angle[i + 1] - angle[i]) / (2 * PI) * 11245.;
                totalPhase += diff[i]; 
            }
        };

    private:
        double* in;
        fftw_complex* out;
        fftw_plan planForward;
        //fftw_plan planInverse;
        std::vector<double> iKernel{ -0.0906,-0.0197,-0.5941,0,0.5941,0.0197,0.0906 };
        std::vector<double> qKernel{ 0,0,0,1,0,0,0 };
};

