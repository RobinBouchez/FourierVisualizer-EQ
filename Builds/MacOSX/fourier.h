#pragma once
#include "frequency.h"
#include "iostream"
#include <math.h>
#include <vector>
#include <complex>

class Fourier
{
public:
    Fourier(std::vector<float> samples, int sampleRate);
    
    void DFT(std::vector<float>&);
    std::vector<Frequency*> spectrum;
    
private:
    std::vector<float> samples;
    int sampleRate;
    const float TAU = 2 * M_PI;
    

};
