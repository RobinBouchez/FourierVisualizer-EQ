#pragma once

#include <JuceHeader.h>

//==============================================================================
/*
*/
class SpectrogramComponent  : public juce::AudioAppComponent, private juce::Timer
{
public:
    SpectrogramComponent();
    ~SpectrogramComponent() override;

    void paint (juce::Graphics&) override;
    void prepareToPlay (int /*samplesPerBlockExpected*/, double /*newSampleRate*/) override;
    void resized() override;
    void timerCallback() override;
    void releaseResources() override;
    void getNextAudioBlock (const juce::AudioSourceChannelInfo&) override;

    enum
    {
        fftOrder  = 11,             // [1]
        fftSize   = 1 << fftOrder,  // [2]
        scopeSize = 512             // [3]
    };
 

private:
    juce::dsp::FFT forwardFFT;                          // [3]
    juce::dsp::WindowingFunction<float> window;     // [5]
    
    float fifo [fftSize];                           // [6]
    float fftData [2 * fftSize];                    // [7]
    int fifoIndex = 0;                              // [8]
    bool nextFFTBlockReady = false;                 // [9]
    float scopeData [scopeSize];
    
    void pushNextSampleIntoFifo (float sample) noexcept;
    void drawNextFrameOfSpectrum();
    void drawFrame (juce::Graphics& g);
    
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (SpectrogramComponent)
};


