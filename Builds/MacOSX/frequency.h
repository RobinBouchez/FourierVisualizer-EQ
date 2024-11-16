
class Frequency {
public:
    Frequency(float frequency, float amplitude, float offset);
    
    float getFrequency() const;
    float getAmplitude() const;
    
private:
    float frequency;
    float amplitude;
    float phase;
};

