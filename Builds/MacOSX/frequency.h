
class Frequency {
public:
    Frequency(float frequency, float amplitude, float offset);
    
    float getFrequency() const;
    float getAmplitude() const;
    
    void setAmplitude(float);
    
private:
    float frequency;
    float amplitude;
    float phase;
};

