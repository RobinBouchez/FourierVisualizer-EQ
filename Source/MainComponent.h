#pragma once

#include <JuceHeader.h>
#include "audioSignal.h"
#include "fourier.h"

class MainComponent  : public juce::AudioAppComponent, public juce::ChangeListener, private juce::Timer
{
public:
    //==============================================================================
    MainComponent();
    ~MainComponent() override;

    //==============================================================================
    void prepareToPlay (int samplesPerBlockExpected, double sampleRate) override;
    void getNextAudioBlock (const juce::AudioSourceChannelInfo& bufferToFill) override;
    void releaseResources() override;
    void changeListenerCallback (juce::ChangeBroadcaster* source) override;

    //==============================================================================
    void paint (juce::Graphics& g) override;
    void resized() override;
    void timerCallback() override;
    //void processBlock(juce::AudioBuffer< float > &buffer, juce::MidiBuffer &midiMessages) override;

private:
    enum TransportState
       {
           Stopped,
           Starting,
           Playing,
           Paused,
           Pausing,
           Stopping
       };
    
    void openButtonClicked();

    void changeState (TransportState newState);
 
    void playButtonClicked();
    void stopButtonClicked();
    
    void transportSourceChanged();
    void thumbnailChanged();
    void paintIfNoFileLoaded (juce::Graphics& g, const juce::Rectangle<int>& thumbnailBounds);
    void paintIfFileLoaded (juce::Graphics& g, const juce::Rectangle<int>& thumbnailBounds);
    juce::Line<float> createFilterLine(juce::Rectangle<int> bounds, float pos, float nextPos, float gain, float cutOffFreq, float nextCutOffFreq, float freq, float nextFreq, float offset);
    
    void readFile(juce::File file);
    
    juce::AudioFormatManager formatManager;
    std::unique_ptr<juce::AudioFormatReaderSource> playSource;
    juce::AudioTransportSource transportSource;
    TransportState state;
    
    std::unique_ptr<juce::FileChooser> chooser;
    
    juce::TextButton openButton;
    juce::TextButton playButton;
    juce::TextButton stopButton;
    
    juce::AudioThumbnailCache thumbnailCache;
    juce::AudioThumbnail thumbnail;
    
    juce::TextButton logButton;
    bool isLogEnabled = false;
    
    juce::TextButton lowPassButton;
    bool isLowPassEnabled = false;
    
    juce::TextButton highPassButton;
    bool isHighPassEnabled = false;
    
    juce::TextButton filterButton3;
    bool isFilter3Enabled = false;
    
    juce::TextButton filterButton4;
    bool isFilter4Enabled = false;
    
    juce::TextButton filterButton5;
    bool isFilter5Enabled = false;
    
    juce::TextButton filterButton6;
    bool isFilter6Enabled = false;
    
    juce::TextButton filterButton7;
    bool isFilter7Enabled = false;
    
    juce::TextButton filterButton8;
    bool isFilter8Enabled = false;
    
    juce::TextButton filterButton9;
    bool isFilter9Enabled = false;
    
    juce::TextButton filterButton10;
    bool isFilter10Enabled = false;
    
    juce::TextButton realtimeButton;
    bool isRealtime = true;
    
    
    juce::Slider lowPassFreqSlider;
    juce::Label  lowPassFreqLabel;
    juce::Slider lowPassFilterSlider;
    juce::Label  lowPassFilterLabel;
    
    juce::Slider frequencySlider3;
    juce::Label  frequencyLabel2;
    juce::Slider filterSlider3;
    juce::Label  filterLabel2;
    
    juce::Slider frequencySlider4;
    juce::Label  frequencyLabel3;
    juce::Slider filterSlider4;
    juce::Label  filterLabel3;
    
    juce::Slider frequencySlider5;
    juce::Label  frequencyLabel4;
    juce::Slider filterSlider5;
    juce::Label  filterLabel4;
    
    juce::Slider frequencySlider6;
    juce::Label  frequencyLabel5;
    juce::Slider filterSlider6;
    juce::Label  filterLabel5;
    
    juce::Slider frequencySlider7;
    juce::Label  frequencyLabel6;
    juce::Slider filterSlider7;
    juce::Label  filterLabel6;
    
    juce::Slider frequencySlider8;
    juce::Label  frequencyLabel7;
    juce::Slider filterSlider8;
    juce::Label  filterLabel7;
    
    juce::Slider frequencySlider9;
    juce::Label  frequencyLabel8;
    juce::Slider filterSlider9;
    juce::Label  filterLabel8;
    
    juce::Slider frequencySlider10;
    juce::Label  frequencyLabel9;
    juce::Slider filterSlider10;
    juce::Label  filterLabel9;
    
    juce::Slider highPassFreqSlider;
    juce::Label  highPassFreqLabel;
    juce::Slider highPassFilterSlider;
    juce::Label  highPassFilterLabel;
    
    juce::AudioSourceChannelInfo channelInfo;
    
    float sampleRateVar;
    
    const float windowBorder_x = 50.f;
    const float windowBorder_y = 20.f;
    
    const float frequencySliderWidth = 120;
    const float frequencySliderHeight = 80;
    const float filterSliderWidth = 100;
    const float filterSliderHeight = 200;
    
    const float thumbnailWidth = 150;
    const float thumbnailHeight = 130;
    
    const float buttonWidth = 100;
    const float buttonHeight = 25;
    const float offset = 20;
    const float labelOffset = 100;
    
    const float minimumFrequency = 20;
    const float maximumFrequency = 20000;
    
    const float minFilterGain = -18.f;
    const float maxFilterGain = 18.f;
    
    const int leftChannel = 0;
    const int rightChannel = 1;
    
    


    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (MainComponent)
};
