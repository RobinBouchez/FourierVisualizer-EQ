#include "fourier.h"

struct Vector2 {
    float x;
    float y;

    // Default constructor
    Vector2() : x(0), y(0) {}

    // Parameterized constructor
    Vector2(float x, float y) : x(x), y(y) {}

    // Vector addition
    Vector2 operator+(const Vector2& other) const {
        return {x + other.x, y + other.y};
    }

    // Vector subtraction
    Vector2 operator-(const Vector2& other) const {
        return {x - other.x, y - other.y};
    }

    // Scalar multiplication
    Vector2 operator*(float scalar) const {
        return {x * scalar, y * scalar};
    }

    // Scalar division
    Vector2 operator/(float scalar) const {
        if (scalar == 0) throw std::invalid_argument("Division by zero");
        return {x / scalar, y / scalar};
    }

    // Dot product
    float dot(const Vector2& other) const {
        return x * other.x + y * other.y;
    }

    // Magnitude (length)
    float magnitude() const {
        return std::sqrt(x * x + y * y);
    }

    // Normalize (make unit vector)
    Vector2 normalized() const {
        float mag = magnitude();
        if (mag == 0) throw std::invalid_argument("Cannot normalize a zero-length vector");
        return {x / mag, y / mag};
    }

    // Print vector
    void print() const {
        std::cout << "(" << x << ", " << y << ")" << std::endl;
    }
};



Fourier::Fourier(std::vector<float> samples, int sampleRate)
: samples(samples), sampleRate(sampleRate)
{
    
}

void Fourier::DFT(std::vector<float>& samples) {
    size_t n = samples.size();
    if (n <= 1) return;

    // Split even and odd elements
    std::vector<float> even(n / 2);
    std::vector<float> odd(n / 2);
    for (int i = 0; i < n / 2; ++i) {
        even[i] = samples[i * 2];
        odd[i] = samples[i * 2 + 1];
    }

    // Recursively perform FFT on both halves
    DFT(even);
    DFT(odd);

    // Combine
    for (int k = 0; k < n / 2; ++k) {
        std::complex<double> t = std::polar(1.0, -2 * M_PI * k / n).real() * odd[k];
        spectrum[k] = new Frequency(k * sampleRate / n, even[k] + t.real(), k);
        spectrum[k + n / 2] = new Frequency(k * sampleRate / n, even[k] - t.real(), k);
    }
}

//std::vector<Frequency*> Fourier::DFT()
//{
//    int numFrequencies = samples.size() / 2 + 1; // Add one since we want to start at 0 Hz
//    std::vector<Frequency*> spectrum(numFrequencies);
//    
//    float frequencyStep = sampleRate / (float)samples.size(); // Equivalent to 1 / duration
//    
//    for(int i = 0; i < numFrequencies; ++i)
//    {
//        Vector2 sampleSum;
//        sampleSum.x = 0;
//        sampleSum.y = 0;
//        for (int i = 0; i < samples.size(); i++)
//        {
//            float angle = i / (float)(samples.size()) * TAU * i;
//            Vector2 testPoint;
//            testPoint.x = cos(angle);
//            testPoint.y = sin(angle);
//            sampleSum.x += testPoint.x * samples[i];
//            sampleSum.y += testPoint.y * samples[i];
//        }
//        
//        Vector2 sampleCentre;
//        sampleCentre.x = sampleSum.x / samples.size();
//        sampleCentre.y = sampleSum.y / samples.size();
//        
//        bool is0Hz = i == 0;
//        // The last frequency is equal to samplerate/2 only if sample count is even
//        bool isNyquistFreq = i == spectrum.size() - 1 && samples.size() % 2 == 0;
//        float amplitudeScale = is0Hz || isNyquistFreq ? 1 : 2;
//        float amplitude = sampleCentre.magnitude() * amplitudeScale;
//        
//        float frequency = i * frequencyStep;
//        float phase = -atan2(sampleCentre.y, sampleCentre.x);
//        spectrum[i] = new Frequency(frequency, amplitude, phase);
//    }
//    
//    return spectrum;
//}
