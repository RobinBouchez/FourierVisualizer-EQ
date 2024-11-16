#include "MainComponent.h"

#define BUFFERSIZE 4096

int sampleBufferSize = 0;

float* sampleBuffer[BUFFERSIZE];
float* sampleArray[BUFFERSIZE];
double magnitudes[BUFFERSIZE];

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
void applyLowPassFilter(std::vector<std::complex<float>>& frequencyDomainSignal, double sampleRate, double cutoffFrequency, double slope_dB_per_octave) {
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
, signal(std::vector<float>(0),0)
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
//    auto channelData = bufferToFill.buffer->getReadPointer(0, bufferToFill.startSample);
//    for (auto sample = 0; sample < bufferToFill.numSamples; ++sample) {
//        sampleArray[sample] = const_cast<float*>(&channelData[sample]);
//    }
//    
    if(sampleArray[0] != NULL) {
        juce::dsp::FFT fftObj = juce::dsp::FFT(sqrt(BUFFERSIZE));
        std::vector<std::complex<float>> audioSignal(BUFFERSIZE);
        std::vector<float> samples(BUFFERSIZE);
        std::vector<std::complex<float>> FFTin(BUFFERSIZE);
        std::vector<std::complex<float>> FFTout(BUFFERSIZE);
        
        for (int i = 0; i < BUFFERSIZE; ++i) {
            samples[i] = *sampleArray[i];
            //audioSignal[i] = std::complex<double>(*sampleArray[i], 0);
            FFTin[i] = std::complex<float>(*sampleArray[i], 0);
            
            fftObj.perform(&FFTin[i], &FFTout[i], false);
        }
        Signal signal = Signal(samples, bufferToFill.numSamples);
        Fourier fourier = Fourier();
        
        // Perform FFT
        //fftSplit(audioSignal);
        
        //applyLowPassFilter(FFTout, sampleRateVar, frequencySlider1.getValue(), filterSlider1.getValue() * -1);

        
        for (int i = 0; i < BUFFERSIZE; ++i) {
            magnitudes[i] = std::pow(std::abs(FFTout[i]), 2);
        }
        
        for (int i = 0; i < BUFFERSIZE; ++i) {
            fftObj.perform(&FFTout[i], &FFTin[i], true);
        }
        
        
        
        
        
//        ifft(audioSignal);
//
//        
       for (auto sample = 0; sample < channelInfo.numSamples; ++sample) {
           channelInfo.buffer->setSample(0, sample, FFTin[sample].real());
           channelInfo.buffer->setSample(1, sample, FFTin[sample].real());
        }
    }
}

void MainComponent::releaseResources()
{
    
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
        float maxFreq = sampleRateVar / 2.0f;
        
        for (int i = 0; i < BUFFERSIZE; ++i) {
            float frequency = i * sampleRateVar/BUFFERSIZE;
            
            // Scale frequency logarithmically for better visualization of lower frequencies
            //float logFreq = std::log10(frequency + 1);  // Add 1 to avoid log(0)
            //float logMaxFreq = std::log10(maxFreq + 1);
            float xPos = spectraBounds.getX() + i;//(frequency / maxFreq) * (getWidth() - 2 * spectraBounds.getX());
            float yPos = spectraBounds.getBottom() - magnitudes[i] * 300; //scaleValue(magnitudes[i], 0, 1, spectraBounds.getY(), spectraBounds.getY() - spectraBounds.getHeight());

            g.setColour (juce::Colours::white);
            g.drawLine(xPos, yPos, xPos, spectraBounds.getBottom());
 //           g.drawText(std::to_string(i), xPos,spectraBounds.getBottom(), 80, 20, juce::Justification::left);
            if (i % 50 == 0){
                float mhz = frequency / 1000;
                std::string rounded = std::to_string((int) std::round(mhz));
                g.drawText(rounded, xPos,spectraBounds.getBottom(), 80, 20, juce::Justification::left);
            }
//            juce::Line<float> line;
//            line.setStart(xPos, yPos);
//            line.setEnd(xPos, spectraBounds.getBottom());
//            freqPath.addLineSegment(line, 1);
            
        }
        
        juce::Line<float> filterLine1;
        filterLine1.setStart(spectraBounds.getX() + (frequencySlider1.getValue() / maxFreq)  * (getWidth() - 2 * spectraBounds.getX()),
                             scaleValue(filterSlider1.getValue() * -1, filterSlider1.getMinimum(), filterSlider1.getMaximum(), spectraBounds.getY(), spectraBounds.getY() + spectraBounds.getHeight()));
        filterLine1.setEnd(filterLine1.getStart().getX(), spectraBounds.getBottom());
        filterPath1.addLineSegment(filterLine1, 1);
        
        g.setColour (juce::Colours::orange);
        g.fillPath(filterPath1);
//        g.setColour (juce::Colours::white);
//        g.fillPath(freqPath);
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
