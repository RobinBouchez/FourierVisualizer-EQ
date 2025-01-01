#include "MainComponent.h"
#include <fftw3.h>
#include <complex.h>
#include <stdio.h>


#define BUFFERSIZE 4096

int sampleBufferSize = 0;

float* sampleBuffer[BUFFERSIZE];
float* sampleArray[BUFFERSIZE];

double magnitudes[BUFFERSIZE];
Signal* signalOfSamples;
Fourier* fourier;
juce::AudioBuffer<float> preComputedAudioBuffer;
juce::File loadedFile;


class FFTProcessor {
private:
    int N; // Size of input
    fftw_complex *in, *out;
    fftw_plan forward_plan;
    fftw_plan inverse_plan;

public:
    FFTProcessor(int size) : N(size) {
        // Allocate memory for input and output
        in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
        out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
        
        // Create forward FFT plan
        forward_plan = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
        
        // Create inverse FFT plan
        inverse_plan = fftw_plan_dft_1d(N, out, in, FFTW_BACKWARD, FFTW_ESTIMATE);
    }

    ~FFTProcessor() {
        fftw_destroy_plan(forward_plan);
        fftw_destroy_plan(inverse_plan);
        fftw_free(in);
        fftw_free(out);
    }

    std::vector<std::complex<double>> computeForwardDFT(const std::vector<double>& input) {
        if (input.size() != static_cast<size_t>(N)) {
            throw std::runtime_error("Input size mismatch");
        }

        // Copy input data to FFTW input array
        for (int i = 0; i < N; i++) {
            in[i][0] = input[i];  // Real part
            in[i][1] = 0.0;       // Imaginary part
        }

        // Execute the forward plan
        fftw_execute(forward_plan);

        // Copy results to output vector
        std::vector<std::complex<double>> result(N);
        for (int i = 0; i < N; i++) {
            result[i] = std::complex<double>(out[i][0], out[i][1]);
        }
        

        return result;
    }

    std::vector<double> computeInverseDFT(const std::vector<std::complex<double>>& input) {
        if (input.size() != static_cast<size_t>(N)) {
            throw std::runtime_error("Input size mismatch");
        }

        // Copy input data to FFTW output array (since inverse FFT reads from out)
        for (int i = 0; i < N; i++) {
            out[i][0] = input[i].real();
            out[i][1] = input[i].imag();
        }

        // Execute the inverse plan
        fftw_execute(inverse_plan);

        // Copy results to output vector and normalize
        std::vector<double> result(N);
        for (int i = 0; i < N; i++) {
            // Normalize by dividing by N (FFTW does not normalize automatically)
            result[i] = in[i][0] / N;
        }

        return result;
    }

    // Utility method to perform forward FFT and then inverse FFT
    std::vector<double> performRoundTrip(const std::vector<double>& input) {
        auto fft_result = computeForwardDFT(input);
        return computeInverseDFT(fft_result);
    }
};

double rectangleFunction(double gain, double freq, double cutoffFreq, double offset) {
    return gain / (1 + std::pow((freq - offset) / cutoffFreq, 16));
}

double scaleValue(double value, double oldMin, double oldMax, double newMin, double newMax) {
    if (oldMin == oldMax) {
        throw std::invalid_argument("oldMin and oldMax must be different to avoid division by zero.");
    }

    // Normalize the value to a [0, 1] range within the oldMin-oldMax range
    double normalized = (value - oldMin) / (oldMax - oldMin);

    // Scale the normalized value to the newMin-newMax range
    return newMin + normalized * (newMax - newMin);
}

//==============================================================================
MainComponent::MainComponent()
: state(Stopped), openButton("Open"), playButton("Play"), stopButton("Stop"),
thumbnailCache (5),                            // [4]
thumbnail (BUFFERSIZE, formatManager, thumbnailCache)
{
    // Make sure you set the size of the component after
    // you add any child components.
    setSize (1600, 900);

    // Some platforms require permissions to open input channels so request that here
    if (juce::RuntimePermissions::isRequired (juce::RuntimePermissions::recordAudio)
        && ! juce::RuntimePermissions::isGranted (juce::RuntimePermissions::recordAudio))
    {
        juce::RuntimePermissions::request (juce::RuntimePermissions::recordAudio,
                                           [&] (bool granted) { setAudioChannels (granted ? 2 : 0, 2); });
    }
    else
    {
        // Specify the number of input and output channels that we want to open
        setAudioChannels (0, 2);
    }
    juce::AudioDeviceManager::AudioDeviceSetup currentAudioSetup;
    deviceManager.getAudioDeviceSetup (currentAudioSetup);
    currentAudioSetup.bufferSize = BUFFERSIZE;
    deviceManager.setAudioDeviceSetup (currentAudioSetup, true);
    
    addAndMakeVisible (&openButton);
    openButton.setButtonText ("Open...");
    openButton.onClick = [this] { openButtonClicked(); };
     
    addAndMakeVisible (&playButton);
    playButton.setButtonText ("Play");
    playButton.onClick = [this] { playButtonClicked(); };
    playButton.setColour (juce::TextButton::buttonColourId, juce::Colours::green);
    playButton.setEnabled (false);
     
    addAndMakeVisible (&stopButton);
    stopButton.setButtonText ("Stop");
    stopButton.onClick = [this] { stopButtonClicked(); };
    stopButton.setColour (juce::TextButton::buttonColourId, juce::Colours::red);
    stopButton.setEnabled (false);
    
    addAndMakeVisible(&lowPassButton);
    lowPassButton.onClick = [this] {
        isLowPassEnabled = !isLowPassEnabled;
        lowPassButton.setButtonText ((isLowPassEnabled) ? "enabled" : "disabled");
    };
    lowPassButton.setButtonText ((isLowPassEnabled) ? "enabled" : "disabled");
    lowPassButton.setColour (juce::TextButton::buttonColourId, juce::Colours::grey);
    
    addAndMakeVisible(&highPassButton);
    highPassButton.onClick = [this] {
        isHighPassEnabled = !isHighPassEnabled;
        highPassButton.setButtonText ((isHighPassEnabled) ? "enabled" : "disabled");
    };
    highPassButton.setButtonText ((isHighPassEnabled) ? "enabled" : "disabled");
    highPassButton.setColour (juce::TextButton::buttonColourId, juce::Colours::grey);
    
    addAndMakeVisible(&filterButton3);
    filterButton3.onClick = [this] {
        isFilter3Enabled = !isFilter3Enabled;
        filterButton3.setButtonText ((isFilter3Enabled) ? "enabled" : "disabled");
    };
    filterButton3.setButtonText ((isFilter3Enabled) ? "enabled" : "disabled");
    filterButton3.setColour (juce::TextButton::buttonColourId, juce::Colours::grey);
    
    addAndMakeVisible(&filterButton4);
    filterButton4.onClick = [this] {
        isFilter4Enabled = !isFilter4Enabled;
        filterButton4.setButtonText ((isFilter4Enabled) ? "enabled" : "disabled");
    };
    filterButton4.setButtonText ((isFilter4Enabled) ? "enabled" : "disabled");
    filterButton4.setColour (juce::TextButton::buttonColourId, juce::Colours::grey);
    
    addAndMakeVisible(&filterButton5);
    filterButton5.onClick = [this] {
        isFilter5Enabled = !isFilter5Enabled;
        filterButton5.setButtonText ((isFilter5Enabled) ? "enabled" : "disabled");
    };
    filterButton5.setButtonText ((isFilter5Enabled) ? "enabled" : "disabled");
    filterButton5.setColour (juce::TextButton::buttonColourId, juce::Colours::grey);
    
    addAndMakeVisible(&filterButton6);
    filterButton6.onClick = [this] {
        isFilter6Enabled = !isFilter6Enabled;
        filterButton6.setButtonText ((isFilter6Enabled) ? "enabled" : "disabled");
    };
    filterButton6.setButtonText ((isFilter6Enabled) ? "enabled" : "disabled");
    filterButton6.setColour (juce::TextButton::buttonColourId, juce::Colours::grey);
    
    addAndMakeVisible(&filterButton7);
    filterButton7.onClick = [this] {
        isFilter7Enabled = !isFilter7Enabled;
        filterButton7.setButtonText ((isFilter7Enabled) ? "enabled" : "disabled");
    };
    filterButton7.setButtonText ((isFilter7Enabled) ? "enabled" : "disabled");
    filterButton7.setColour (juce::TextButton::buttonColourId, juce::Colours::grey);
    
    addAndMakeVisible(&filterButton8);
    filterButton8.onClick = [this] {
        isFilter8Enabled = !isFilter8Enabled;
        filterButton8.setButtonText ((isFilter8Enabled) ? "enabled" : "disabled");
    };
    filterButton8.setButtonText ((isFilter8Enabled) ? "enabled" : "disabled");
    filterButton8.setColour (juce::TextButton::buttonColourId, juce::Colours::grey);
    
    addAndMakeVisible(&filterButton9);
    filterButton9.onClick = [this] {
        isFilter9Enabled = !isFilter9Enabled;
        filterButton9.setButtonText ((isFilter9Enabled) ? "enabled" : "disabled");
    };
    filterButton9.setButtonText ((isFilter9Enabled) ? "enabled" : "disabled");
    filterButton9.setColour (juce::TextButton::buttonColourId, juce::Colours::grey);
    
    addAndMakeVisible(&filterButton10);
    filterButton10.onClick = [this] {
        isFilter10Enabled = !isFilter10Enabled;
        filterButton10.setButtonText ((isFilter10Enabled) ? "enabled" : "disabled");
    };
    filterButton10.setButtonText ((isFilter10Enabled) ? "enabled" : "disabled");
    filterButton10.setColour (juce::TextButton::buttonColourId, juce::Colours::grey);
    
    addAndMakeVisible(&realtimeButton);
    realtimeButton.onClick = [this] {
        isRealtime = !isRealtime;
        realtimeButton.setButtonText ((isRealtime) ? "realtime enabled" : "realtime disabled");
    };
    realtimeButton.setButtonText ((isRealtime) ? "realtime enabled" : "realtime disabled");
    realtimeButton.setColour (juce::TextButton::buttonColourId, juce::Colours::green);
    
    
    addAndMakeVisible (lowPassFreqSlider);
    lowPassFreqSlider.setRange (minimumFrequency, maximumFrequency);
    lowPassFreqSlider.setTextValueSuffix (" Hz");
    lowPassFreqSlider.setSliderStyle(juce::Slider::Rotary);
    lowPassFreqSlider.setValue(maximumFrequency);
     
    addAndMakeVisible(lowPassFreqLabel);
    lowPassFreqLabel.setText ("Cutoff Freq", juce::dontSendNotification);
    lowPassFreqLabel.attachToComponent (&lowPassFreqSlider, true);
    lowPassFreqSlider.setTextBoxStyle(juce::Slider::TextBoxBelow, false, 160, lowPassFreqSlider.getTextBoxHeight());
    
    addAndMakeVisible (lowPassFilterSlider);
    lowPassFilterSlider.setRange(minFilterGain, maxFilterGain);
    lowPassFilterSlider.setTextValueSuffix (" dB");
    lowPassFilterSlider.setSliderStyle(juce::Slider::LinearVertical);
     
    addAndMakeVisible(lowPassFilterLabel);
    lowPassFilterLabel.setText ("Gain", juce::dontSendNotification);
    lowPassFilterLabel.attachToComponent (&lowPassFilterSlider, true);
    lowPassFilterSlider.setTextBoxStyle(juce::Slider::TextBoxBelow, false, 160, lowPassFilterSlider.getTextBoxHeight());
    
    addAndMakeVisible (frequencySlider3);
    frequencySlider3.setRange (50, 500.0);
    frequencySlider3.setTextValueSuffix (" Hz");
    frequencySlider3.setSliderStyle(juce::Slider::Rotary);
    frequencySlider3.setTextBoxStyle(juce::Slider::TextBoxBelow, false, 160, frequencySlider3.getTextBoxHeight());
    
    addAndMakeVisible (filterSlider3);
    filterSlider3.setRange (minFilterGain, maxFilterGain);
    filterSlider3.setTextValueSuffix (" Db");
    filterSlider3.setSliderStyle(juce::Slider::LinearVertical   );
    filterSlider3.setTextBoxStyle(juce::Slider::TextBoxBelow, false, 160, filterSlider3.getTextBoxHeight());
    
    addAndMakeVisible (frequencySlider4);
    frequencySlider4.setRange (100, 200.0);
    frequencySlider4.setTextValueSuffix (" Hz");
    frequencySlider4.setSliderStyle(juce::Slider::Rotary);
    frequencySlider4.setTextBoxStyle(juce::Slider::TextBoxBelow, false, 160, frequencySlider4.getTextBoxHeight());
    
    addAndMakeVisible (filterSlider4);
    filterSlider4.setRange (minFilterGain, maxFilterGain);
    filterSlider4.setTextValueSuffix (" Db");
    filterSlider4.setSliderStyle(juce::Slider::LinearVertical);
    filterSlider4.setTextBoxStyle(juce::Slider::TextBoxBelow, false, 160, filterSlider4.getTextBoxHeight());
    
    addAndMakeVisible (frequencySlider5);
    frequencySlider5.setRange (200, 500.0);
    frequencySlider5.setTextValueSuffix (" Hz");
    frequencySlider5.setSliderStyle(juce::Slider::Rotary);
    frequencyLabel4.attachToComponent (&frequencySlider5, true);
    frequencySlider5.setTextBoxStyle(juce::Slider::TextBoxBelow, false, 160, frequencySlider5.getTextBoxHeight());
    
    addAndMakeVisible (filterSlider5);
    filterSlider5.setRange (minFilterGain, maxFilterGain);
    filterSlider5.setTextValueSuffix (" Db");
    filterSlider5.setSliderStyle(juce::Slider::LinearVertical   );
    filterSlider5.setTextBoxStyle(juce::Slider::TextBoxBelow, false, 160, filterSlider5.getTextBoxHeight());
    
    addAndMakeVisible (frequencySlider6);
    frequencySlider6.setRange (500, 1000.0);
    frequencySlider6.setTextValueSuffix (" Hz");
    frequencySlider6.setSliderStyle(juce::Slider::Rotary);
    frequencySlider6.setTextBoxStyle(juce::Slider::TextBoxBelow, false, 160, frequencySlider6.getTextBoxHeight());
    
    addAndMakeVisible (filterSlider6);
    filterSlider6.setRange (minFilterGain, maxFilterGain);
    filterSlider6.setTextValueSuffix (" Db");
    filterSlider6.setSliderStyle(juce::Slider::LinearVertical   );
    filterSlider6.setTextBoxStyle(juce::Slider::TextBoxBelow, false, 160, filterSlider6.getTextBoxHeight());
    
    addAndMakeVisible (frequencySlider7);
    frequencySlider7.setRange (1000, 5000.0);
    frequencySlider7.setTextValueSuffix (" Hz");
    frequencySlider7.setSliderStyle(juce::Slider::Rotary);
    frequencySlider7.setTextBoxStyle(juce::Slider::TextBoxBelow, false, 160, frequencySlider7.getTextBoxHeight());
    
    addAndMakeVisible (filterSlider7);
    filterSlider7.setRange (minFilterGain, maxFilterGain);
    filterSlider7.setTextValueSuffix (" Db");
    filterSlider7.setSliderStyle(juce::Slider::LinearVertical   );
    filterSlider7.setTextBoxStyle(juce::Slider::TextBoxBelow, false, 160, filterSlider7.getTextBoxHeight());
    
    addAndMakeVisible (frequencySlider8);
    frequencySlider8.setRange (5000, 10000.0);
    frequencySlider8.setTextValueSuffix (" Hz");
    frequencySlider8.setSliderStyle(juce::Slider::Rotary);
    frequencySlider8.setTextBoxStyle(juce::Slider::TextBoxBelow, false, 160, frequencySlider8.getTextBoxHeight());
    
    addAndMakeVisible (filterSlider8);
    filterSlider8.setRange (minFilterGain, maxFilterGain);
    filterSlider8.setTextValueSuffix (" Db");
    filterSlider8.setSliderStyle(juce::Slider::LinearVertical   );
    filterSlider8.setTextBoxStyle(juce::Slider::TextBoxBelow, false, 160, filterSlider8.getTextBoxHeight());
    
    addAndMakeVisible (frequencySlider9);
    frequencySlider9.setRange (10000.0, 20000.0);
    frequencySlider9.setTextValueSuffix (" Hz");
    frequencySlider9.setSliderStyle(juce::Slider::Rotary);
    frequencySlider9.setTextBoxStyle(juce::Slider::TextBoxBelow, false, 160, frequencySlider9.getTextBoxHeight());
    
    addAndMakeVisible (filterSlider9);
    filterSlider9.setRange (minFilterGain, maxFilterGain);
    filterSlider9.setTextValueSuffix (" Db");
    filterSlider9.setSliderStyle(juce::Slider::LinearVertical   );
    filterSlider9.setTextBoxStyle(juce::Slider::TextBoxBelow, false, 160, filterSlider9.getTextBoxHeight());
    
    addAndMakeVisible (frequencySlider10);
    frequencySlider10.setRange (10000.0, 20000.0);
    frequencySlider10.setTextValueSuffix (" Hz");
    frequencySlider10.setSliderStyle(juce::Slider::Rotary);
    frequencySlider10.setTextBoxStyle(juce::Slider::TextBoxBelow, false, 160, frequencySlider10.getTextBoxHeight());
    
    addAndMakeVisible (filterSlider10);
    filterSlider10.setRange (minFilterGain, maxFilterGain);
    filterSlider10.setTextValueSuffix (" Db");
    filterSlider10.setSliderStyle(juce::Slider::LinearVertical   );
    filterSlider10.setTextBoxStyle(juce::Slider::TextBoxBelow, false, 160, filterSlider10.getTextBoxHeight());
    
    addAndMakeVisible (highPassFreqSlider);
    highPassFreqSlider.setRange (minimumFrequency, maximumFrequency);
    highPassFreqSlider.setTextValueSuffix (" Hz");
    highPassFreqSlider.setSliderStyle(juce::Slider::Rotary);
    highPassFreqSlider.setTextBoxStyle(juce::Slider::TextBoxBelow, false, 160, highPassFreqSlider.getTextBoxHeight());
    highPassFreqSlider.setValue(minimumFrequency);

    
    addAndMakeVisible (highPassFilterSlider);
    highPassFilterSlider.setRange (minFilterGain, maxFilterGain);
    highPassFilterSlider.setTextValueSuffix (" Db");
    highPassFilterSlider.setSliderStyle(juce::Slider::LinearVertical   );
    highPassFilterSlider.setTextBoxStyle(juce::Slider::TextBoxBelow, false, 160, highPassFilterSlider.getTextBoxHeight());
    highPassFilterSlider.setValue(0);
    
    formatManager.registerBasicFormats();
    transportSource.addChangeListener (this);
    thumbnail.addChangeListener (this);
    
    startTimer (40);
}

MainComponent::~MainComponent()
{
    // This shuts down the audio device and clears the audio source.
    shutdownAudio();
}

//==============================================================================
void MainComponent::prepareToPlay (int samplesPerBlockExpected, double sampleRate)
{
    transportSource.prepareToPlay(samplesPerBlockExpected, sampleRate);
    sampleRateVar = sampleRate;
}


void MainComponent::getNextAudioBlock (const juce::AudioSourceChannelInfo& bufferToFill)
{
    if (playSource.get() == nullptr)
    {
        bufferToFill.clearActiveBufferRegion();
        return;
    }

        channelInfo = juce::AudioSourceChannelInfo(bufferToFill);
        transportSource.getNextAudioBlock(channelInfo);
        
        if (sampleBufferSize >= BUFFERSIZE) {
            for (auto sample = 0; sample < sampleBufferSize; ++sample) {
                sampleArray[sample] = sampleBuffer[sample];
            }
            sampleBufferSize = 0;
        } else {
            auto channelData = bufferToFill.buffer->getReadPointer(0, bufferToFill.startSample);
            for (auto sample = 0; sample < bufferToFill.numSamples; ++sample) {
                sampleBuffer[sample + sampleBufferSize] = const_cast<float*>(&channelData[sample]);
            }
            sampleBufferSize += bufferToFill.numSamples;
        }
        
        if(sampleArray[0] != NULL) {
            std::vector<double> samples(BUFFERSIZE);
            for (int i = 0; i < BUFFERSIZE; ++i) {
                samples[i] = *sampleArray[i];
            }
            if(isRealtime) {
            FFTProcessor fft(BUFFERSIZE);
            std::vector<std::complex<double>> spectrum = fft.computeForwardDFT(samples);
            
            
            
            for (int n = 0; n < BUFFERSIZE; ++n) {
                float frequency = n * sampleRateVar/BUFFERSIZE;
                double gainLinear = juce::jmap((float)lowPassFilterSlider.getValue(), minFilterGain, maxFilterGain, 0.f, 2.f);
                double rectResult = rectangleFunction(gainLinear, frequency, lowPassFreqSlider.getValue() , 0);
                gainLinear = juce::jmap((float)highPassFilterSlider.getValue(), minFilterGain, maxFilterGain, 0.f, 2.f);
                double rectResultHighpass = rectangleFunction(gainLinear, highPassFreqSlider.getValue(), frequency, 0);
                magnitudes[n] = std::sqrt(std::pow(spectrum[n].real(), 2) + std::pow(spectrum[n].imag(), 2)) * rectResult * rectResultHighpass;
            }
            
            for (int i = 0; i < BUFFERSIZE; ++i) {
                float frequency = i * sampleRateVar/BUFFERSIZE;
                double gainLinear = juce::jmap((float)lowPassFilterSlider.getValue(), minFilterGain, maxFilterGain, 0.f, 2.f);
                double rectResult = rectangleFunction(gainLinear, frequency, lowPassFreqSlider.getValue(), 0);
                gainLinear = juce::jmap((float)highPassFilterSlider.getValue(), minFilterGain, maxFilterGain, 0.f, 2.f);
                double rectResultHighpass = rectangleFunction(gainLinear, highPassFreqSlider.getValue(), frequency, 0);
                spectrum[i] *= rectResult * rectResultHighpass;
            }
            
            // Inverse transform
            auto reconstructed = fft.computeInverseDFT(spectrum);
            
            
            for (auto sample = 0; sample < channelInfo.numSamples; ++sample) {
                channelInfo.buffer->setSample(leftChannel, sample, reconstructed[sample]);
                channelInfo.buffer->setSample(rightChannel, sample, reconstructed[sample]);
            }
            } else {
                for (int k = 0; k < BUFFERSIZE; ++k) {
                    std::complex<float> f = 0;
                    for (int n = 0; n < BUFFERSIZE; ++n) {
                        f += (std::complex<float>)samples[n] *
                        std::exp(std::complex<float>(0, -2 * M_PI * k * n / BUFFERSIZE));
                    }
                    magnitudes[k] = std::abs(f); // Store the magnitude
                }
            }
    }
}

void MainComponent::releaseResources() {
    transportSource.releaseResources();
}


// Helper function to scale magnitudes to 0-1 range
float normalizeMagnitude(float magnitude, float maxPossibleMagnitude) {
    // Apply scaling and optional dB conversion
    float scaledMag = magnitude / maxPossibleMagnitude;
    
    // Optional: Convert to dB scale (common for spectrum analyzers)
    float dbMag = 20 * std::log10(scaledMag + 1e-6f); // Add small value to avoid log(0)
    
    // Map dB range (-100dB to 0dB) to 0-1
    float normalizedDb = juce::jmap(dbMag, -100.f, 0.f, -18.f, 18.f);
    return normalizedDb;
}

juce::Line<float> MainComponent::createFilterLine(juce::Rectangle<int> bounds, float pos, float nextPos, float gain, float cutOffFreq, float nextCutOffFreq, float freq, float nextFreq, float offset) {
    juce::Line<float> line;
    
    double rectResult1 = rectangleFunction(juce::jmap(gain, minFilterGain, maxFilterGain, 0.f, 2.f), freq, cutOffFreq, offset);
    double rectResult2 = rectangleFunction(juce::jmap(gain, minFilterGain, maxFilterGain, 0.f, 2.f), nextFreq, nextCutOffFreq, offset);
    
    line.setStart(pos, bounds.getY() + bounds.getHeight() - juce::jmap((float)rectResult1, 0.f, 2.f, 0.f, (float)bounds.getHeight()));
    line.setEnd(nextPos, bounds.getY() + bounds.getHeight() - juce::jmap((float)rectResult2, 0.f, 2.f, 0.f, (float)bounds.getHeight()));
    
    return line;
}

//==============================================================================
void MainComponent::paint (juce::Graphics& g)
{
    // (Our component is opaque, so we must completely fill the background with a solid colour)
    g.fillAll(juce::Colours::black);
    
    juce::Rectangle<int> thumbnailBounds (windowBorder_x, getHeight() - thumbnailHeight - buttonHeight - offset - windowBorder_y, getWidth() - (2 * windowBorder_x), thumbnailHeight);
     
            if (thumbnail.getNumChannels() == 0)
                paintIfNoFileLoaded (g, thumbnailBounds);
            else
                paintIfFileLoaded (g, thumbnailBounds);
    
    juce::Rectangle<int> spectraBounds(windowBorder_x, windowBorder_y, getWidth() - windowBorder_x * 2, 300);
    
    g.setColour (juce::Colours::grey);
    g.drawLine(spectraBounds.getX(), spectraBounds.getY() + spectraBounds.getHeight(), spectraBounds.getX() + spectraBounds.getWidth(), spectraBounds.getY() + spectraBounds.getHeight(), 1);
    g.drawLine(spectraBounds.getX(), spectraBounds.getY(), spectraBounds.getX(), spectraBounds.getY() + spectraBounds.getHeight(), 1);
    
    
        juce::Path freqPath;
        juce::Path filterPath1;
        juce::Path filterPath2;
        juce::Path filterPath3;
        juce::Path filterPath4;
        juce::Path filterPath5;
        juce::Path filterPath6;
        juce::Path filterPath7;
        juce::Path filterPath8;
        juce::Path filterPath9;
        juce::Path filterPath10;
        
        float maxFreq = std::min(maximumFrequency, sampleRateVar / 2.0f);
        float minFreq = minimumFrequency;
        
        static const float labelDBs[] = {18, 12, 6, 0, -6, -12, -18};
            for (float labelDB: labelDBs) {
                juce::String label = juce::String(labelDB) + "dB";
                g.drawText(label, spectraBounds.getX() - 40, spectraBounds.getBottom() - 20  - (spectraBounds.getHeight() / 2)  + ((labelDB / -6) * spectraBounds.getHeight() / 6), 40, 20, juce::Justification::centred);
            }
        
        for (int i = 0; i < BUFFERSIZE; ++i) {
            float frequency = i * sampleRateVar/BUFFERSIZE;
            
            // Skip frequencies outside 20Hz-20kHz range
            if (frequency < minFreq || frequency > maxFreq) {
                continue;
            }
            
            float normalizedMagnitude = normalizeMagnitude(magnitudes[i], BUFFERSIZE / 2);
            
            // Scale frequency logarithmically for better visualization
            float logFreq = (isLogEnabled) ? std::log10(frequency) : frequency;
            float logMaxFreq = maxFreq;
            float logMinFreq = minFreq;
            
            // Normalize position to the visible frequency range
            float normalizedPos = (logFreq - logMinFreq) / (logMaxFreq - logMinFreq);
            float xPos = spectraBounds.getX() + normalizedPos * (getWidth() - 2 * spectraBounds.getX());
            float yPos = spectraBounds.getBottom() - juce::jmap(normalizedMagnitude, minFilterGain, maxFilterGain, 0.f, (float)(spectraBounds.getHeight()));//(normalizedMagnitude * spectraBounds.getHeight() /2);

            g.setColour(juce::Colours::white);
            if (yPos < spectraBounds.getBottom()) {
                g.drawLine(xPos, yPos, xPos, spectraBounds.getBottom());
            }

            float nextfrequency = (i + 1) * sampleRateVar/BUFFERSIZE;
            float nextlogFreq = (isLogEnabled) ? std::log10(nextfrequency) : nextfrequency;
            
            // Normalize position to the visible frequency range
            float nextnormalizedPos = (nextlogFreq - logMinFreq) / (logMaxFreq - logMinFreq);
            float nextxPos = spectraBounds.getX() + nextnormalizedPos * (getWidth() - 2 * spectraBounds.getX());

            if(isLowPassEnabled){
                filterPath1.addLineSegment(createFilterLine(spectraBounds,
                                                            xPos, nextxPos,
                                                            lowPassFilterSlider.getValue(),
                                                            lowPassFreqSlider.getValue(), lowPassFreqSlider.getValue(),
                                                            frequency, nextfrequency,
                                                            0), 1);
            }
            if(isHighPassEnabled) {
                filterPath2.addLineSegment(createFilterLine(spectraBounds,
                                                            xPos, nextxPos,
                                                            highPassFilterSlider.getValue(),
                                                            frequency, nextfrequency,
                                                            highPassFreqSlider.getValue(), highPassFreqSlider.getValue(),
                                                            0), 1);
            }
            if(isFilter3Enabled){
                filterPath3.addLineSegment(createFilterLine(spectraBounds,
                                                            xPos, nextxPos,
                                                            filterSlider3.getValue(),
                                                            (frequencySlider3.getMaximum() - frequencySlider3.getMinimum()) / 2,
                                                            (frequencySlider3.getMaximum() - frequencySlider3.getMinimum()) / 2,
                                                            frequency, nextfrequency,
                                                            frequencySlider3.getValue()), 1);
            }
            if(isFilter4Enabled){
                filterPath4.addLineSegment(createFilterLine(spectraBounds,
                                                            xPos, nextxPos,
                                                            filterSlider4.getValue(),
                                                            (frequencySlider4.getMaximum() - frequencySlider4.getMinimum()) / 2,
                                                            (frequencySlider4.getMaximum() - frequencySlider4.getMinimum()) / 2,
                                                            frequency, nextfrequency,
                                                            frequencySlider4.getValue()), 1);
            }
            if(isFilter5Enabled) {
                filterPath5.addLineSegment(createFilterLine(spectraBounds,
                                                            xPos, nextxPos,
                                                            filterSlider5.getValue(),
                                                            (frequencySlider5.getMaximum() - frequencySlider5.getMinimum()) / 2,
                                                            (frequencySlider5.getMaximum() - frequencySlider5.getMinimum()) / 2,
                                                            frequency, nextfrequency,
                                                            frequencySlider5.getValue()), 1);
            }
            if(isFilter6Enabled) {
                filterPath6.addLineSegment(createFilterLine(spectraBounds,
                                                            xPos, nextxPos,
                                                            filterSlider6.getValue(),
                                                            (frequencySlider6.getMaximum() - frequencySlider6.getMinimum()) / 2,
                                                            (frequencySlider6.getMaximum() - frequencySlider6.getMinimum()) / 2,
                                                            frequency, nextfrequency,
                                                            frequencySlider6.getValue()), 1);
            }
            if(isFilter7Enabled) {
                            filterPath7.addLineSegment(createFilterLine(spectraBounds,
                                                                        xPos, nextxPos,
                                                                        filterSlider7.getValue(),
                                                                        (frequencySlider7.getMaximum() - frequencySlider7.getMinimum()) / 2,
                                                                        (frequencySlider7.getMaximum() - frequencySlider7.getMinimum()) / 2,
                                                                        frequency, nextfrequency,
                                                                        frequencySlider7.getValue()), 1);
            }
            if(isFilter8Enabled) {
                filterPath8.addLineSegment(createFilterLine(spectraBounds,
                                                            xPos, nextxPos,
                                                            filterSlider8.getValue(),
                                                            (frequencySlider8.getMaximum() - frequencySlider8.getMinimum()) / 2,
                                                            (frequencySlider8.getMaximum() - frequencySlider8.getMinimum()) / 2,
                                                            frequency, nextfrequency,
                                                            frequencySlider8.getValue()), 1);
            }
            if(isFilter9Enabled) {
                filterPath9.addLineSegment(createFilterLine(spectraBounds,
                                                            xPos, nextxPos,
                                                            filterSlider9.getValue(),
                                                            (frequencySlider9.getMaximum() - frequencySlider9.getMinimum()) / 2,
                                                            (frequencySlider9.getMaximum() - frequencySlider9.getMinimum()) / 2,
                                                            frequency, nextfrequency,
                                                            frequencySlider9.getValue()), 1);
            }
            if(isFilter10Enabled) {
                filterPath10.addLineSegment(createFilterLine(spectraBounds,
                                                             xPos, nextxPos,
                                                             filterSlider10.getValue(),
                                                             frequencySlider10.getValue(), frequencySlider10.getValue(),
                                                             frequency, nextfrequency,
                                                             frequencySlider10.getValue()), 1);
            }


            // Draw frequency labels at logarithmically spaced intervals
            static const float labelFreqs[] = {20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000};
            for (float labelFreq : labelFreqs) {
                if (std::abs(frequency - labelFreq) < (sampleRateVar/BUFFERSIZE/2)) {
                    float displayFreq = labelFreq >= 1000 ? labelFreq/1000 : labelFreq;
                    juce::String label = juce::String(displayFreq) + (labelFreq >= 1000 ? "k" : "");
                    g.drawText(label, xPos - 20, spectraBounds.getBottom(), 40, 20, juce::Justification::left);
                }
            }
        }
        g.setColour (juce::Colours::orange);
        g.fillPath(filterPath1);
        g.setColour (juce::Colours::purple);
        g.fillPath(filterPath2);
        g.setColour (juce::Colours::red);
        g.fillPath(filterPath3);
        g.setColour (juce::Colours::blue);
        g.fillPath(filterPath4);
        g.setColour (juce::Colours::green);
        g.fillPath(filterPath5);
        g.setColour (juce::Colours::yellow);
        g.fillPath(filterPath6);
        g.setColour (juce::Colours::grey);
        g.fillPath(filterPath7);
        g.setColour (juce::Colours::cyan);
        g.fillPath(filterPath8);
        g.setColour (juce::Colours::magenta);
        g.fillPath(filterPath9);
        g.setColour (juce::Colours::aqua);
        g.fillPath(filterPath10);
        
    
}

void MainComponent::resized()
{
    openButton.setBounds(windowBorder_x,getHeight() - buttonHeight - windowBorder_y,buttonWidth, buttonHeight);
    playButton.setBounds(getWidth() / 2 - buttonWidth,getHeight() - buttonHeight - windowBorder_y,buttonWidth, buttonHeight);
    stopButton.setBounds(getWidth() / 2 + buttonWidth / 2,getHeight() - buttonHeight - windowBorder_y,buttonWidth, buttonHeight);
    
    auto sliderLeft = windowBorder_x + labelOffset;
    auto sliderOffset = (getWidth() - labelOffset - windowBorder_x * 2) / 10 ;
    auto sliderY = getHeight() - 320;
    
    realtimeButton.setBounds(getWidth() - windowBorder_x - 150, getHeight() - buttonHeight - windowBorder_y, 150, buttonHeight);
    
    lowPassButton.setBounds(sliderLeft, sliderY + frequencySliderHeight + 10, 100, buttonHeight);
    filterButton3.setBounds(sliderLeft + sliderOffset, sliderY + frequencySliderHeight + 10, 100, buttonHeight);
    filterButton4.setBounds(sliderLeft + (sliderOffset * 2), sliderY + frequencySliderHeight + 10, 100, buttonHeight);
    filterButton5.setBounds(sliderLeft + (sliderOffset * 3), sliderY + frequencySliderHeight + 10, 100, buttonHeight);
    filterButton6.setBounds(sliderLeft + (sliderOffset * 4), sliderY + frequencySliderHeight + 10, 100, buttonHeight);
    filterButton7.setBounds(sliderLeft + (sliderOffset * 5), sliderY + frequencySliderHeight + 10, 100, buttonHeight);
    filterButton8.setBounds(sliderLeft + (sliderOffset * 6), sliderY + frequencySliderHeight + 10, 100, buttonHeight);
    filterButton9.setBounds(sliderLeft + (sliderOffset * 7), sliderY + frequencySliderHeight + 10, 100, buttonHeight);
    filterButton10.setBounds(sliderLeft + (sliderOffset * 8), sliderY + frequencySliderHeight + 10, 100, buttonHeight);
    highPassButton.setBounds(sliderLeft + (sliderOffset * 9), sliderY + frequencySliderHeight + 10, 100, buttonHeight);
    
    lowPassFilterSlider.setBounds  (sliderLeft,                      sliderY - 210, filterSliderWidth, filterSliderHeight);
    filterSlider3.setBounds        (sliderLeft + sliderOffset,       sliderY - 210, filterSliderWidth, filterSliderHeight);
    filterSlider4.setBounds        (sliderLeft + (sliderOffset * 2), sliderY - 210, filterSliderWidth, filterSliderHeight);
    filterSlider5.setBounds        (sliderLeft + (sliderOffset * 3), sliderY - 210, filterSliderWidth, filterSliderHeight);
    filterSlider6.setBounds        (sliderLeft + (sliderOffset * 4), sliderY - 210, filterSliderWidth, filterSliderHeight);
    filterSlider7.setBounds        (sliderLeft + (sliderOffset * 5), sliderY - 210, filterSliderWidth, filterSliderHeight);
    filterSlider8.setBounds        (sliderLeft + (sliderOffset * 6), sliderY - 210, filterSliderWidth, filterSliderHeight);
    filterSlider9.setBounds        (sliderLeft + (sliderOffset * 7), sliderY - 210, filterSliderWidth, filterSliderHeight);
    filterSlider10.setBounds        (sliderLeft + (sliderOffset * 8), sliderY - 210, filterSliderWidth, filterSliderHeight);
    highPassFilterSlider.setBounds (sliderLeft + (sliderOffset * 9), sliderY - 210, filterSliderWidth, filterSliderHeight);
    
    lowPassFreqSlider.setBounds  (sliderLeft,                      sliderY, frequencySliderWidth, frequencySliderHeight);
    frequencySlider3.setBounds   (sliderLeft + sliderOffset,       sliderY, frequencySliderWidth, frequencySliderHeight);
    frequencySlider4.setBounds   (sliderLeft + (sliderOffset * 2), sliderY, frequencySliderWidth, frequencySliderHeight);
    frequencySlider5.setBounds   (sliderLeft + (sliderOffset * 3), sliderY, frequencySliderWidth, frequencySliderHeight);
    frequencySlider6.setBounds   (sliderLeft + (sliderOffset * 4), sliderY, frequencySliderWidth, frequencySliderHeight);
    frequencySlider7.setBounds   (sliderLeft + (sliderOffset * 5), sliderY, frequencySliderWidth, frequencySliderHeight);
    frequencySlider8.setBounds   (sliderLeft + (sliderOffset * 6), sliderY, frequencySliderWidth, frequencySliderHeight);
    frequencySlider9.setBounds   (sliderLeft + (sliderOffset * 7), sliderY, frequencySliderWidth, frequencySliderHeight);
    frequencySlider10.setBounds   (sliderLeft + (sliderOffset * 8), sliderY, frequencySliderWidth, frequencySliderHeight);
    highPassFreqSlider.setBounds (sliderLeft + (sliderOffset * 9), sliderY, frequencySliderWidth, frequencySliderHeight);
}

void MainComponent::readFile(juce::File file) {
    if (file != juce::File{} && file.existsAsFile())  {
        juce::AudioFormatReader* reader = formatManager.createReaderFor(file);
        
        if (reader != nullptr)
        {
            std::unique_ptr<juce::AudioFormatReaderSource> newSource (new juce::AudioFormatReaderSource(reader, true));
            
//
//                preComputedAudioBuffer.setSize(reader->numChannels, reader->lengthInSamples);
//                reader->read(&preComputedAudioBuffer, 0, reader->lengthInSamples, 0, true, true);
//
            
            transportSource.setSource(newSource.get(), 0, nullptr, reader->sampleRate);
            playButton.setEnabled(true);
            changeState(Stopped);
            thumbnail.setSource (new juce::FileInputSource (file));
            playSource.reset(newSource.release());
        }
    }
}

void MainComponent::openButtonClicked()
  {
    chooser = std::make_unique<juce::FileChooser> ("Please select the file you want to load...",
                                                   juce::File::getSpecialLocation (juce::File::userDesktopDirectory),
                                                "*.mp3;*.wav");
 
    
    auto folderChooserFlags = juce::FileBrowserComponent::openMode |  juce::FileBrowserComponent::canSelectFiles;
 
    chooser->launchAsync (folderChooserFlags, [this] (const juce::FileChooser& chsr)
    {
        juce::File file = chsr.getResult();
        loadedFile = file;
        readFile(file);
    });
  }




void MainComponent::changeListenerCallback (juce::ChangeBroadcaster* source)
{
    if (source == &transportSource)
    {
        if (!transportSource.isPlaying() && state == Playing){
            readFile(loadedFile);
        }
        
        if (transportSource.isPlaying())
            changeState (Playing);
        else
            if(state == Pausing) {
                changeState (Paused);
            } else if(state == Stopping) {
                changeState(Stopped);
            } else changeState(Stopping);
    }
    
    if (source == &thumbnail)       thumbnailChanged();
}

void MainComponent::changeState(TransportState newState)
 {
    if (state != newState)
    {
        state = newState;
        
        switch (state)
        {
            case Stopped:
                playButton.setButtonText ("Play");
                stopButton.setButtonText ("Stop");
                stopButton.setEnabled (false);
                transportSource.setPosition (0.0);
                break;
                
            case Starting:
                transportSource.start();
                break;
                
            case Playing:
                playButton.setButtonText ("Pause");
                stopButton.setButtonText ("Stop");
                stopButton.setEnabled (true);
                break;
                
            case Pausing:
                transportSource.stop();
                break;
                
            case Paused:
                playButton.setButtonText ("Resume");
                stopButton.setButtonText ("Return to Zero");
                break;
                
            case Stopping:
                transportSource.stop();
                break;
        }
    }
 }

void MainComponent::playButtonClicked()
{
    if ((state == Stopped) || (state == Paused))
        changeState (Starting);
    else if (state == Playing)
        changeState (Pausing);
}

void MainComponent::stopButtonClicked()
{
    if (state == Paused)
        changeState (Stopped);
    else
        changeState (Stopping);
}

void MainComponent::timerCallback()
 {
     repaint();
 }

void MainComponent::transportSourceChanged()
{
    changeState (transportSource.isPlaying() ? Playing : Stopped);
}

void MainComponent::thumbnailChanged()
{
    repaint();
}
void MainComponent::paintIfNoFileLoaded (juce::Graphics& g, const juce::Rectangle<int>& thumbnailBounds)
{
    g.setColour (juce::Colours::black);
    g.fillRect (thumbnailBounds);
    g.setColour (juce::Colours::white);
    g.drawFittedText ("No File Loaded", thumbnailBounds, juce::Justification::centred, 1);
}

void MainComponent::paintIfFileLoaded (juce::Graphics& g, const juce::Rectangle<int>& thumbnailBounds)
{
    g.setColour (juce::Colours::black);
    g.fillRect (thumbnailBounds);

    g.setColour (juce::Colours::white);

    auto audioLength = (float) thumbnail.getTotalLength();
    thumbnail.drawChannels (g, thumbnailBounds, 0.0, audioLength, 1.0f);

    g.setColour (juce::Colours::red);

    auto audioPosition = (float) transportSource.getCurrentPosition();
    auto drawPosition = (audioPosition / audioLength) * (float) thumbnailBounds.getWidth()
                        + (float) thumbnailBounds.getX();
    g.drawLine (drawPosition, (float) thumbnailBounds.getY(), drawPosition,
                (float) thumbnailBounds.getBottom(), 2.0f);
}
