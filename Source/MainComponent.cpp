#include "MainComponent.h"

int sampleBufferSize = 0;

const float* sampleBuffer[512];

#include <cmath>
#include <complex>
#include <vector>

// Function to perform in-place Cooley-Tukey FFT
void fft(std::vector<std::complex<double>>& a) {
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
    fft(even);
    fft(odd);

    // Combine
    for (int k = 0; k < n / 2; ++k) {
        std::complex<double> t = std::polar(1.0, -2 * M_PI * k / n) * odd[k];
        a[k] = even[k] + t;
        a[k + n / 2] = even[k] - t;
    }
}


//==============================================================================
MainComponent::MainComponent()
: state(Stopped), openButton("Open"), playButton("Play"), stopButton("Stop"),
thumbnailCache (5),                            // [4]
        thumbnail (512, formatManager, thumbnailCache)
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
    frequencySlider1.setRange (0, 100.0);
    frequencySlider1.setTextValueSuffix (" Hz");
    frequencySlider1.setSliderStyle(juce::Slider::Rotary);
//            frequencySlider.onValueChange = [this] { durationSlider.setValue (1.0 / frequencySlider.getValue(), juce::dontSendNotification); };
     
    addAndMakeVisible(frequencyLabel1);
    frequencyLabel1.setText ("Frequency", juce::dontSendNotification);
    frequencyLabel1.attachToComponent (&frequencySlider1, true);
    frequencySlider1.setTextBoxStyle(juce::Slider::TextBoxBelow, false, 160, frequencySlider1.getTextBoxHeight());
    
    addAndMakeVisible (filterSlider1);
    filterSlider1.setRange (-12.0, 24.0);
    filterSlider1.setTextValueSuffix (" Db");
    filterSlider1.setSliderStyle(juce::Slider::LinearBarVertical);
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
    filterSlider2.setSliderStyle(juce::Slider::LinearBarVertical);
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
    
    sampleBufferSize = bufferToFill.numSamples;
    auto* channelData = bufferToFill.buffer->getReadPointer(0, bufferToFill.startSample);

    for (auto sample = 0; sample < sampleBufferSize; ++sample) {
        sampleBuffer[sample] = &channelData[sample];
    }

    
    transportSource.getNextAudioBlock (bufferToFill);
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
    
    
    g.setColour (juce::Colours::white);
    
    
    // Example signal of size 512 (can be filled with actual data)
    std::vector<std::complex<double>> signal(sampleBufferSize);
    for (int i = 0; i < sampleBufferSize; ++i) {
        signal[i] = std::complex<double>(*sampleBuffer[i], 0);  // Example: Sine wave
    }

    // Perform FFT
    fft(signal);

    for (int i = 0; i < sampleBufferSize / 2; ++i) {
        double magnitude = std::abs(signal[i]);
        float frequency = i * sampleRateVar/sampleBufferSize;
        
        // Scale frequency logarithmically for better visualization of lower frequencies
        float maxFreq = sampleRateVar / 2.0f;
        float logFreq = std::log10(frequency + 1);  // Add 1 to avoid log(0)
        float logMaxFreq = std::log10(maxFreq + 1);
        float xPos = spectraBounds.getX() + (logFreq / logMaxFreq) * 1500;
        
        float yPos = spectraBounds.getBottom() - (magnitude * 2);
        // Draw the point
        g.fillEllipse(xPos, yPos, 1, 1);
    }
}

void MainComponent::resized()
{
    openButton.setBounds(windowBorder_x,getHeight() - buttonHeight - windowBorder_y,buttonWidth, buttonHeight);
    playButton.setBounds(getWidth() / 2 - buttonWidth,getHeight() - buttonHeight - windowBorder_y,buttonWidth, buttonHeight);
    stopButton.setBounds(getWidth() / 2 + buttonWidth / 2,getHeight() - buttonHeight - windowBorder_y,buttonWidth, buttonHeight);
    
    auto sliderLeft = windowBorder_x + labelOffset;
    auto sliderOffset = frequencySliderWidth * 2;
    auto sliderY = getHeight() - 340;
    frequencySlider1.setBounds (sliderLeft, sliderY, frequencySliderWidth, frequencySliderHeight);
    filterSlider1.setBounds (sliderLeft, sliderY - 220, filterSliderWidth, filterSliderHeight);
    frequencySlider2.setBounds (sliderLeft + sliderOffset, sliderY, frequencySliderWidth, frequencySliderHeight);
    filterSlider2.setBounds (sliderLeft + sliderOffset, sliderY - 220, filterSliderWidth, filterSliderHeight);
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
                                                "*.wav","*.mp3");
 
    
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
