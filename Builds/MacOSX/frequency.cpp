#include "frequency.h"

Frequency::Frequency(float frequency, float amplitude, float offset)
:frequency(frequency), amplitude(amplitude), phase(offset)
{
    
}
    
float Frequency::getFrequency() const {
    return frequency;
}

float Frequency::getAmplitude() const {
    return amplitude;
}
