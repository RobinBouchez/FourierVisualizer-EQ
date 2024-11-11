#include "SpectrogramComponent.h"

//==============================================================================
SpectrogramComponent::SpectrogramComponent()
: forwardFFT (fftOrder), window (fftSize, juce::dsp::WindowingFunction<float>::hann)
{
    setOpaque (true);
    setAudioChannels (2, 0);  // we want a couple of input channels but no outputs
    startTimerHz (60);
    setSize (800, 300);

}

SpectrogramComponent::~SpectrogramComponent()
{
    shutdownAudio();
}

void SpectrogramComponent::paint (juce::Graphics& g)
{
    g.setOpacity (1.0f);
    g.setColour (juce::Colours::white);
    drawFrame (g);
}

void SpectrogramComponent::drawFrame (juce::Graphics& g)
{
    for (int i = 1; i < scopeSize; ++i)
    {
        auto width  = getLocalBounds().getWidth();
        auto height = getLocalBounds().getHeight();
        

        g.drawLine ({ (float) juce::jmap (i - 1, 0, scopeSize - 1, 0, width),
                              juce::jmap (scopeData[i - 1], 0.0f, 1.0f, (float) height, 0.0f),
                      (float) juce::jmap (i,     0, scopeSize - 1, 0, width),
                              juce::jmap (scopeData[i],     0.0f, 1.0f, (float) height, 0.0f) });
    }
}

void SpectrogramComponent::resized()
{

}

void SpectrogramComponent::prepareToPlay (int samplesPerBlockExpected, double newSampleRate){
    
}

void SpectrogramComponent::releaseResources() {
    
}

void SpectrogramComponent::getNextAudioBlock (const juce::AudioSourceChannelInfo& bufferToFill)
{
    if (bufferToFill.buffer->getNumChannels() > 0)
    {
        auto* channelData = bufferToFill.buffer->getReadPointer (0, bufferToFill.startSample);
        
        for (auto i = 0; i < bufferToFill.numSamples; ++i) {
            pushNextSampleIntoFifo (channelData[i]);
        }
    }
}


void SpectrogramComponent::pushNextSampleIntoFifo (float sample) noexcept
{
    if (fifoIndex == fftSize)               // [11]
    {
        if (! nextFFTBlockReady)            // [12]
        {
            juce::zeromem (fftData, sizeof (fftData));
            memcpy (fftData, fifo, sizeof (fifo));
            nextFFTBlockReady = true;
        }

        fifoIndex = 0;
    }

    fifo[(size_t) fifoIndex++] = sample; // [9]
}

void SpectrogramComponent::drawNextFrameOfSpectrum()
{
    // first apply a windowing function to our data
    window.multiplyWithWindowingTable (fftData, fftSize);       // [1]

    // then render our FFT data..
    forwardFFT.performFrequencyOnlyForwardTransform (fftData);  // [2]

    auto mindB = -100.0f;
    auto maxdB =    0.0f;

    for (int i = 0; i < scopeSize; ++i)                         // [3]
    {
        auto skewedProportionX = 1.0f - std::exp (std::log (1.0f - (float) i / (float) scopeSize) * 0.2f);
        auto fftDataIndex = juce::jlimit (0, fftSize / 2, (int) (skewedProportionX * (float) fftSize * 0.5f));
        auto level = juce::jmap (juce::jlimit (mindB, maxdB, juce::Decibels::gainToDecibels (fftData[fftDataIndex])
                                                           - juce::Decibels::gainToDecibels ((float) fftSize)),
                                 mindB, maxdB, 0.0f, 1.0f);

        scopeData[i] = level;                                   // [4]
    }
}

void SpectrogramComponent::timerCallback()
{
    if (nextFFTBlockReady)
    {
        //std::cout << "callback\n";
        drawNextFrameOfSpectrum();
        nextFFTBlockReady = false;
        repaint();
    }
}
