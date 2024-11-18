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


// Function to perform in-place Cooley-Tukey FFT
void fftSplit(std::vector<std::complex<double>>& a) {
    int n = a.size();
    if (n <= 1) return;

    // Split even and odd elements
    std::vector<std::complex<double>> even(n / 2);
    std::vector<std::complex<double>> odd(n / 2);
    for (int i = 0; i < n / 2; ++i) {
        even[i] = a[i * 2];
        odd[i] = a[i * 2 + 1];
    }

    // Recursively perform FFT on both halves
    fftSplit(even);
    fftSplit(odd);

    // Combine
    for (int k = 0; k < n / 2; ++k) {
        std::complex<double> t = std::polar(1.0, -2 * M_PI * k / n) * odd[k];
        a[k] = even[k] + t;
        a[k + n / 2] = even[k] - t;
    }
}

void fft(std::vector<std::complex<double>>& x) {
    unsigned int n = x.size();
    if (n <= 1) return;
    
    
    std::vector<std::complex<double>> X(n);
    for (int k = 0; k < n; ++k) {
        X[k] = 0;
        for (int i = 0; i < n; ++i) {
            std::complex<double> exponent(0, -2 * M_PI * k * i / n);
            X[k] += std::exp(exponent) * x[i];
        }
        x[k] = X[k];
    }
}


void ifft(std::vector<std::complex<double>>& a) {
    int n = a.size();
    if (n <= 1) return;

    for (auto& x : a) {
        x = std::conj(x);
    }
    fftSplit(a);
    for (auto& x : a) {
        x = std::conj(x) / static_cast<double>(n);
    }
}

// Low-pass filter function using rectangular function
void lowPassFilter(std::vector<std::complex<double>>& frequencyDomainSignal, double cutoffFrequency, double gain, double sampleRate) {
    int n = frequencyDomainSignal.size();
    double nyquist = sampleRate / 2.0;

    // Calculate the cutoff bin index based on cutoff frequency
    int cutoffBin = static_cast<int>((cutoffFrequency / nyquist) * (n / 2));

    // Zero out frequencies beyond the cutoff frequency
    for (int i = cutoffBin; i < n - cutoffBin; ++i) {
        frequencyDomainSignal[i] = frequencyDomainSignal[i] * gain;
    }
}

// Low-pass filter function
void applyLowPassFilter(std::vector<std::complex<double>>& frequencyDomainSignal, double sampleRate, double cutoffFrequency, double slope_dB_per_octave) {
    size_t N = frequencyDomainSignal.size();
    if (N == 0 || sampleRate <= 0) {
        throw std::invalid_argument("Invalid input parameters");
    }

    // Compute frequency resolution
    double frequencyResolution = sampleRate / N;

    // Convert slope from dB per octave to attenuation factor per frequency step
    double slopeFactor = std::pow(10.0, -slope_dB_per_octave / 20.0);

    // Loop through the frequency domain signal
    for (size_t i = 0; i < N; ++i) {
        double frequency = i * frequencyResolution;
        
        // Compute attenuation based on the distance from the cutoff frequency
        if (frequency > cutoffFrequency) {
            double octaveDifference = std::log2(frequency / cutoffFrequency);
            double attenuationFactor = std::pow(slopeFactor, octaveDifference);

            frequencyDomainSignal[i] *= attenuationFactor;

            // Handle symmetry for real signals
            if (i > 0 && i < N / 2) {
                frequencyDomainSignal[N - i] *= attenuationFactor; // Mirror index
            }
        }
    }
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
    
    addAndMakeVisible (frequencySlider1);
    frequencySlider1.setRange (minimumFrequency, maximumFrequency);
    frequencySlider1.setTextValueSuffix (" Hz");
    frequencySlider1.setSliderStyle(juce::Slider::Rotary);
    frequencySlider1.setValue(maximumFrequency);
    
//            frequencySlider.onValueChange = [this] { durationSlider.setValue (1.0 / frequencySlider.getValue(), juce::dontSendNotification); };
     
    addAndMakeVisible(frequencyLabel1);
    frequencyLabel1.setText ("Frequency", juce::dontSendNotification);
    frequencyLabel1.attachToComponent (&frequencySlider1, true);
    frequencySlider1.setTextBoxStyle(juce::Slider::TextBoxBelow, false, 160, frequencySlider1.getTextBoxHeight());
    
    addAndMakeVisible (filterSlider1);
    filterSlider1.setRange (-24.0, 24.0);
    filterSlider1.setTextValueSuffix (" dB");
    filterSlider1.setSliderStyle(juce::Slider::LinearVertical);
//            frequencySlider.onValueChange = [this] { durationSlider.setValue (1.0 / frequencySlider.getValue(), juce::dontSendNotification); };
     
    addAndMakeVisible(filterLabel1);
    filterLabel1.setText ("Filter", juce::dontSendNotification);
    filterLabel1.attachToComponent (&filterSlider1, true);
    filterSlider1.setTextBoxStyle(juce::Slider::TextBoxBelow, false, 160, filterSlider1.getTextBoxHeight());
    
    addAndMakeVisible (frequencySlider2);
    frequencySlider2.setRange (50, 500.0);
    frequencySlider2.setTextValueSuffix (" Hz");
    frequencySlider2.setSliderStyle(juce::Slider::Rotary);
//            frequencySlider.onValueChange = [this] { durationSlider.setValue (1.0 / frequencySlider.getValue(), juce::dontSendNotification); };
    
    addAndMakeVisible(frequencyLabel2);
    frequencyLabel2.setText ("Frequency", juce::dontSendNotification);
    frequencyLabel2.attachToComponent (&frequencySlider2, true);
    frequencySlider2.setTextBoxStyle(juce::Slider::TextBoxBelow, false, 160, frequencySlider2.getTextBoxHeight());
    
    addAndMakeVisible (filterSlider2);
    filterSlider2.setRange (-12.0, 24.0);
    filterSlider2.setTextValueSuffix (" Db");
    filterSlider2.setSliderStyle(juce::Slider::LinearVertical   );
//            frequencySlider.onValueChange = [this] { durationSlider.setValue (1.0 / frequencySlider.getValue(), juce::dontSendNotification); };
     
    addAndMakeVisible(filterLabel2);
    filterLabel2.setText ("Filter", juce::dontSendNotification);
    filterLabel2.attachToComponent (&filterSlider2, true);
    filterSlider2.setTextBoxStyle(juce::Slider::TextBoxBelow, false, 160, filterSlider2.getTextBoxHeight());
    
    addAndMakeVisible (frequencySlider3);
    frequencySlider3.setRange (50, 500.0);
    frequencySlider3.setTextValueSuffix (" Hz");
    frequencySlider3.setSliderStyle(juce::Slider::Rotary);
//            frequencySlider.onValueChange = [this] { durationSlider.setValue (1.0 / frequencySlider.getValue(), juce::dontSendNotification); };
     
    addAndMakeVisible(frequencyLabel3);
    frequencyLabel3.setText ("Frequency", juce::dontSendNotification);
    frequencyLabel3.attachToComponent (&frequencySlider3, true);
    frequencySlider3.setTextBoxStyle(juce::Slider::TextBoxBelow, false, 160, frequencySlider3.getTextBoxHeight());
    
    addAndMakeVisible (frequencySlider4);
    frequencySlider4.setRange (50, 500.0);
    frequencySlider4.setTextValueSuffix (" Hz");
    frequencySlider4.setSliderStyle(juce::Slider::Rotary);
//            frequencySlider.onValueChange = [this] { durationSlider.setValue (1.0 / frequencySlider.getValue(), juce::dontSendNotification); };
     
    addAndMakeVisible(frequencyLabel4);
    frequencyLabel4.setText ("Frequency", juce::dontSendNotification);
    frequencyLabel4.attachToComponent (&frequencySlider4, true);
    frequencySlider4.setTextBoxStyle(juce::Slider::TextBoxBelow, false, 160, frequencySlider4.getTextBoxHeight());
    
    
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
        
        FFTProcessor fft(BUFFERSIZE);
        std::vector<std::complex<double>> spectrum = fft.computeForwardDFT(samples);
        
        applyLowPassFilter(spectrum, sampleRateVar, frequencySlider1.getValue(), filterSlider1.getValue() * -1);
        
        
        for (int i = 0; i < BUFFERSIZE; ++i) {
            magnitudes[i] = std::pow(std::abs(spectrum[i].real()), 2);
        }
        
        
        // Inverse transform
        auto reconstructed = fft.computeInverseDFT(spectrum);
        
    
       for (auto sample = 0; sample < channelInfo.numSamples; ++sample) {
           channelInfo.buffer->setSample(0, sample, reconstructed[sample]);
           channelInfo.buffer->setSample(1, sample, reconstructed[sample]);
        }
    }
}

void MainComponent::releaseResources()
{
    
}


// Helper function to scale magnitudes to 0-1 range
float normalizeMagnitude(float magnitude, float maxPossibleMagnitude) {
    // Apply scaling and optional dB conversion
    float scaledMag = magnitude / maxPossibleMagnitude;
    
    // Optional: Convert to dB scale (common for spectrum analyzers)
    float dbMag = 20 * std::log10(scaledMag + 1e-6f); // Add small value to avoid log(0)
    
    // Map dB range (-100dB to 0dB) to 0-1
    float normalizedDb = juce::jmap(dbMag, -100.0f, 0.0f, 0.0f, 1.0f);
    return juce::jlimit(0.0f, 1.0f, normalizedDb);
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
    
    juce::Rectangle<int> spectraBounds (windowBorder_x, windowBorder_y, getWidth() - 40, 300);
    
    g.setColour (juce::Colours::grey);
    g.drawLine(spectraBounds.getX(), spectraBounds.getY() + spectraBounds.getHeight(), spectraBounds.getX() + spectraBounds.getWidth(), spectraBounds.getY() + spectraBounds.getHeight(), 1);
    
    if(sampleArray[0] != NULL) {
        juce::Path freqPath;
        juce::Path filterPath1;
        
        float maxFreq = std::min(20000.0f, sampleRateVar / 2.0f);
        float minFreq = 20.0f;

        for (int i = 0; i < BUFFERSIZE / 2; ++i) {
            float frequency = i * sampleRateVar/BUFFERSIZE;
            
            // Skip frequencies outside 20Hz-20kHz range
            if (frequency < minFreq || frequency > maxFreq) {
                continue;
            }
            
            // Scale magnitude to 0-1 range
            float normalizedMagnitude = normalizeMagnitude(magnitudes[i], BUFFERSIZE);
            
            // Scale frequency logarithmically for better visualization
            float logFreq = std::log10(frequency);
            float logMaxFreq = std::log10(maxFreq);
            float logMinFreq = std::log10(minFreq);
            
            // Normalize position to the visible frequency range
            float normalizedPos = (logFreq - logMinFreq) / (logMaxFreq - logMinFreq);
            float xPos = spectraBounds.getX() + normalizedPos * (getWidth() - 2 * spectraBounds.getX());
            float yPos = spectraBounds.getBottom() - (normalizedMagnitude * spectraBounds.getHeight());

            g.setColour(juce::Colours::white);
            g.drawLine(xPos, yPos, xPos, spectraBounds.getBottom());

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
        juce::Line<float> filterLine1;
        filterLine1.setStart(spectraBounds.getX() + (frequencySlider1.getValue() / maxFreq)  * (getWidth() - 2 * spectraBounds.getX()),
                             scaleValue(filterSlider1.getValue() * -1, filterSlider1.getMinimum(), filterSlider1.getMaximum(), spectraBounds.getY(), spectraBounds.getY() + spectraBounds.getHeight()));
        filterLine1.setEnd(filterLine1.getStart().getX(), spectraBounds.getBottom());
        filterPath1.addLineSegment(filterLine1, 1);
        
        g.setColour (juce::Colours::orange);
        g.fillPath(filterPath1);
    }
}

void MainComponent::resized()
{
    openButton.setBounds(windowBorder_x,getHeight() - buttonHeight - windowBorder_y,buttonWidth, buttonHeight);
    playButton.setBounds(getWidth() / 2 - buttonWidth,getHeight() - buttonHeight - windowBorder_y,buttonWidth, buttonHeight);
    stopButton.setBounds(getWidth() / 2 + buttonWidth / 2,getHeight() - buttonHeight - windowBorder_y,buttonWidth, buttonHeight);
    
    auto sliderLeft = windowBorder_x + labelOffset;
    auto sliderOffset = frequencySliderWidth * 2;
    auto sliderY = getHeight() - 320;
    frequencySlider1.setBounds (sliderLeft, sliderY, frequencySliderWidth, frequencySliderHeight);
    filterSlider1.setBounds (sliderLeft, sliderY - 210, filterSliderWidth, filterSliderHeight);
    frequencySlider2.setBounds (sliderLeft + sliderOffset, sliderY, frequencySliderWidth, frequencySliderHeight);
    filterSlider2.setBounds (sliderLeft + sliderOffset, sliderY - 210, filterSliderWidth, filterSliderHeight);
    frequencySlider3.setBounds (sliderLeft + (sliderOffset * 2), sliderY, frequencySliderWidth, frequencySliderHeight);
    frequencySlider4.setBounds (sliderLeft + (sliderOffset * 3), sliderY, frequencySliderWidth, frequencySliderHeight);
    frequencySlider5.setBounds (sliderLeft + (sliderOffset * 4), sliderY, frequencySliderWidth, frequencySliderHeight);
    frequencySlider6.setBounds (sliderLeft + (sliderOffset * 5), sliderY, frequencySliderWidth, frequencySliderHeight);
    frequencySlider7.setBounds (sliderLeft + (sliderOffset * 6), sliderY, frequencySliderWidth, frequencySliderHeight);
    frequencySlider8.setBounds (sliderLeft + (sliderOffset * 7), sliderY, frequencySliderWidth, frequencySliderHeight);
    frequencySlider9.setBounds (sliderLeft + (sliderOffset * 8), sliderY, frequencySliderWidth, frequencySliderHeight);
    frequencySlider10.setBounds (sliderLeft + (sliderOffset * 9), sliderY, frequencySliderWidth, frequencySliderHeight);
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
        if (file != juce::File{} && file.existsAsFile())  {
            juce::AudioFormatReader* reader = formatManager.createReaderFor(file);
            
            if (reader != nullptr)
            {
                std::unique_ptr<juce::AudioFormatReaderSource> newSource (new juce::AudioFormatReaderSource(reader, true));
                
                transportSource.setSource(newSource.get(), 0, nullptr, reader->sampleRate);
                playButton.setEnabled(true);
                changeState(Stopped);
                thumbnail.setSource (new juce::FileInputSource (file));       
                playSource.reset(newSource.release());
            }
        }
    });
  }




void MainComponent::changeListenerCallback (juce::ChangeBroadcaster* source)
{
    if (source == &transportSource)
    {
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

void MainComponent::changeState (TransportState newState)
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
