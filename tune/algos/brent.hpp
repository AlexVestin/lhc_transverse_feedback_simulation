#ifndef __BRENT_H__
#define __BRENT_H__

#include <inttypes.h>
#include <math.h>
#include <complex>

#define cplx std::complex<double> 
#define cplxvec std::vector<cplx>


struct merit_args_cpp
{
    size_t N;
    cplxvec window;
    cplxvec signal;
};

typedef struct merit_args_cpp merit_args_cpp;



void hann_harm_window_cpp(cplxvec& window, const size_t N, const double n)
{
    double T1 = 0.;
    double T2 = N;
    double TM = (T2-T1)/2.;
    double PIST = M_PI/TM;

    int factorial_1 = 1;
    int factorial_2 = 1;

    for( size_t j = 1; j <= n; j++ )
        factorial_1 *= j;
    for( size_t j = 1; j <= 2*n; j++ )
        factorial_2 *= j;

    //double cn = pow(2., n)*(((1.*factorial_1)*factorial_1)/(1.*factorial_2));
    double cn = exp(n*log(2.))*(((1.*factorial_1)*factorial_1)/(1.*factorial_2));
    
    for( size_t i = 0; i < N; i++)
        //window[i] = cn*pow(1. + cos( (i-TM)*PIST ) ,n);
        window[i] = cn*exp(n*log(1. + cos( (i-TM)*PIST )));
    
    return;
}

std::complex<double> inner_product(const cplxvec& signal, double amplitude, double frequency,const cplxvec& window, size_t N)
{
    double omega = (2*M_PI)*frequency;
    std::complex<double> result = 0.;
    double theta = omega*N;
    for( size_t i = N; i--; )
    {
        theta -= omega;
        cplx val = cplx(cos(theta), -sin(theta));
        result += signal[i]*val*window[i];
    }
    return (amplitude*result) / (double)N;
}


double minus_magnitude_fourier_integral_v2(double frequency, const merit_args_cpp* S) {
    std::complex<double> amp = inner_product(S->signal, 1., frequency, S->window, S->N);
    return -(amp.real()*amp.real() + amp.imag()*amp.imag());
}

double brent_minimize_cpp(double (*f)(double,const merit_args_cpp*), double min, double max, const merit_args_cpp* S)
{

    const int max_iter = 10000;
    const double golden = 0.3819660;
    const double tolerance = 1.490116e-8; //ldexp(1.-25)

    double x, w, v, u; 
    double fu, fv, fw, fx;
    double mid;
    double delta1, delta2;  
    double tol0 = tolerance*0.25;
    double tol1, tol2;  // minimal relative movement in x


    x = w = v = max;

    // Merit function
    fw = fv = fx = (*f)(x, S);
    delta2 = delta1 = 0;


    for(int i = max_iter; i--;)
    {
        mid = 0.5 * (min + max);
        tol1 = tolerance * fabs(x) + tol0;
        tol2 = 2. * tol1;
        if( fabs(x - mid) <= (tol2 - 0.5 * (max - min)) )
            return x;

        if( fabs(delta2) > tol1 )
        {
            double r = (x - w) * (fx - fv);
            double q = (x - v) * (fx - fw);
            double p = (x - v) * q - (x - w) * r;
            q = 2.0 * (q - r);
            if( q > 0 )
                p = -p;
            q = fabs(q);
            double delta0 = delta2;
            delta2 = delta1;
            if( fabs(p) >= fabs((0.5 * q) * delta0) || p <= q * (min - x) || p >= q * (max - x) )
            {
                delta2 = (x >= mid) ? min - x : max - x;
                delta1 = golden * delta2;
            }
            else
            {   
                delta1 = p / q;
                u = x + delta1;
                if( (u - min) < tol2 || (max - u) < tol2)
                    delta1 = (mid - x) < 0 ? -fabs(tol1) : fabs(tol1);
            }
        }
        else
        {
            delta2 = (x >= mid) ? min - x : max - x;
            delta1 = golden * delta2;
        }
        u = (fabs(delta1) >= tol1) ? x + delta1 : x + ( (delta1 > 0) ? fabs(tol1) : -fabs(tol1) );
        fu = (*f)(u,S);
        if(fu <= fx)
        {
            if( u >= x )
                min = x;
            else
                max = x;
            v = w;
            w = x;
            x = u;
            fv = fw;
            fw = fx;
            fx = fu;
        }
        else
        {
            if( u < x )
                min = u;
            else
                max = u;
            if( fu <= fw || w == x )
            {
                v = w;
                w = u;
                fv = fw;
                fw = fu;
            }
            else if( fu <= fv || v == x || v == w)
            {
                v = u;
                fv = fu;
            }
        }
    }
    printf("WARNING: nafflib Brent minimization reached maximum number of iterations: %d.\n",max_iter);
    return x;
}


#endif