#pragma once
#include <vector>


class Signal
{
public:
    Signal(std::vector<float> samples, int sampleRate);
    
    
private:
    std::vector<float> Samples;
    int SampleRate;

    int NumSamples = (int)Samples.size();
    float Duration = NumSamples / (float)SampleRate;
};
