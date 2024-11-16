#pragma once
#include "frequency.h"
#include "iostream"
#include <math.h>
#include <vector>

class Fourier
{
public:
    Fourier(std::vector<float> samples, int sampleRate);
    
    std::vector<Frequency*> DFT();
    
private:
    std::vector<float> samples;
    int sampleRate;
    const float TAU = 2 * M_PI;
};
