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
#define NUM_PICKUPS 4
#define N 1024*NUM_PICKUPS
#define PI 3.14159265359
// Sample data signal 

#define AMP 1.0
// #define FREQUENCY 82.13
#define AVG_NOISE_AMT 0.03
#define REVOLUTION_FREQUENCY 11245.f





class FFTContainer {
public:
    FFTContainer(int size): size { size } {
        int outSize = size / 2 + 1; 
        in          = std::vector<double>(size);
        window      = std::vector<double>(size);
        out         = std::vector<fftw_complex>(outSize);
        magnitude   = std::vector<double>(outSize);
        
        plan = fftw_plan_dft_r2c_1d(size, in.data(), out.data(), FFTW_ESTIMATE);
        
        for(int i = 0; i < in.size(); i++) {  
            window[i] = 0.5 * (1 - std::cos(2*PI*i / (in.size() - 1)));
        }
    }

    ~FFTContainer() {
        fftw_destroy_plan(plan);
    }

    void fillMagnitude() {
        float scale = 2.0 / (double)size;
        for(int i = 0; i < magnitude.size(); i++) {
            magnitude[i] = std::sqrt(out[i][0] * out[i][0] + out[i][1] * out[i][1]) * scale;  
        } 
    }

    void fillSampleData(float frequency) {
        float increment = frequency * 2. * PI;
        const std::vector<float> pickupOffsets{0, 0.001, 0.00013, 0.002 };

        double totalEnergy = 0.0;
        for(int turn = 0; turn < size; turn++) {
            for(int j = 0; j < NUM_PICKUPS; j++) {
                float pos = (turn + pickupOffsets[j]) / REVOLUTION_FREQUENCY;
                in[turn] = (AMP * std::sin(pos * increment));
            }
        }
    }
    std::vector<double>& analyse(float frequency) {
        fillSampleData(frequency);
        fftw_execute(plan);
        fillMagnitude();
        return magnitude;
    }

private: 
    int size;
    std::vector<double> magnitude;
    std::vector<double> in;
    std::vector<fftw_complex> out;
    std::vector<double> window;
    fftw_plan plan;
};


int main() {
    FFTContainer container1{N};
    FFTContainer container2{N*4};

     
    bool running = true;
    int counter = 0; 

    //std::string testFile = "tune_data/7343/match10/07343_64k_B1H_Q10_20181025_05h05m39s.h5";
    //HDFLib::HDFFile test = HDFLib::HDFFile(testFile);
    //test.open();
    //test.setTranspose(true);
    //std::cout<<"data: "<<test[255].get()[0]<<std::endl;
    
    #ifdef SKIA
        InitWindow();
        while(running) {
            float frequency = 82.145 + counter / 20.;
            auto data1 = container1.analyse(frequency);
            auto data2 = container2.analyse(frequency);
            ClearCanvas();
            DrawPoints(data1, frequency, N, REVOLUTION_FREQUENCY, { 1.0, 0.0, 1.0, 1.0 });
            DrawPoints(data2, frequency, N*4, REVOLUTION_FREQUENCY, { .1, 0.0, 1.0, 1.0 });
            running = FlushCanvas();
            usleep(1000000); // = 0.01 second.
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