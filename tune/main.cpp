#include <iostream>
#include "algos/Naff.hpp"
#include <unistd.h>
#include <cstring>
#include "HDFLib.h"
//#include "algos/Hilbert.hpp"

#include <math.h>


#ifdef SKIA
    #include "Plotting.hpp"
#endif

#ifdef SCIPLOT
    #include <sciplot/sciplot.hpp>
    using namespace sciplot;
#endif

extern "C" {
    #include <fftw3.h>
}

// Number of samples to use for analysis
#define NUM_PICKUPS 1
#define N 256*NUM_PICKUPS
#define PI 3.14159265359
// Sample data signal 

#define AMP 1.0
// #define FREQUENCY 82.13
#define AVG_NOISE_AMT 0.03
#define REVOLUTION_FREQUENCY 11245.f




std::vector<double> fillSampleData(std::vector<double>& data, float frequency) {
    float increment = frequency * 2. * PI;
    const std::vector<float> pickupOffsets{0, 0.001, 0.00013, 0.002 };
    for(int turn = 0; turn < N; turn++) {
        for(int j = 0; j < NUM_PICKUPS; j++) {
            float pos = (turn + pickupOffsets[j]) / REVOLUTION_FREQUENCY;
            data[turn] = AMP * std::sin(pos * increment);
        }
    }

    return data;
}

void fillMagnitude(std::vector<fftw_complex>& data, std::vector<double>& magnitude) {
    for(int i = 0; i < magnitude.size(); i++) {
        magnitude[i] = std::sqrt(data[i][0] * data[i][0] + data[i][1] * data[i][1]) / (double)N;  
    } 
}


int main() {
    Naff h{N};

    // TODO Aligned?
    const int outSize = N / 2 + 1; 
    std::vector<double> sampleData(N);
    std::vector<fftw_complex> out(outSize);
    std::vector<double> magnitude(outSize);

    fftw_plan p = fftw_plan_dft_r2c_1d(N, sampleData.data(), out.data(), FFTW_ESTIMATE);
     
    bool running = true;
    int counter = 0; 

    // std::string testFile = "tune_data/7343/match10/07343_64k_B1H_Q9_20181025_05h05m34s.h5";
    // HDFLib::HDFFile test = HDFLib::HDFFile(testFile);
    // test.open();
    // test.setTranspose(true);
    // std::cout<<"data: "<<test[255].get()[0]<<std::endl;
    
    #ifdef SKIA
        InitWindow();
        while(running) {
            float frequency = 82.145 + counter / 20.;
            fillSampleData(sampleData, frequency);
            h.performAnalysis(sampleData);
            fftw_execute(p);
            fillMagnitude(out, magnitude);
            running = DrawPoints(magnitude, frequency);
            //usleep(10000); // = 0.01 second.
            counter++;
        }
        CloseWindow();
    #endif



    #ifdef SCIPLOT
    Vec x = linspace(0.0, N, N);
    Plot plot;
    plot.palette("set2");
    plot.xrange(0, N);
    plot.yrange(-AMP, AMP);
    plot.drawPoints(x, sampleData).labelNone().pointSize(1).pointType(0);
    plot.show();
    #endif
    // Save the plot to a PDF file
    return 0;   
}