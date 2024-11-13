#pragma once

#include <JuceHeader.h>
#include "SpectrogramComponent.h"

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
    
    juce::Slider frequencySlider1;
    juce::Label  frequencyLabel1;
    juce::Slider filterSlider1;
    juce::Label  filterLabel1;
    juce::Slider frequencySlider2;
    juce::Label  frequencyLabel2;
    juce::Slider filterSlider2;
    juce::Label  filterLabel2;
    juce::Slider frequencySlider3;
    juce::Label  frequencyLabel3;
    juce::Slider filterSlider3;
    juce::Label  filterLabel3;
    juce::Slider frequencySlider4;
    juce::Label  frequencyLabel4;
    juce::Slider frequencySlider5;
    juce::Label  frequencyLabel5;
    juce::Slider frequencySlider6;
    juce::Label  frequencyLabel6;
    juce::Slider frequencySlider7;
    juce::Label  frequencyLabel7;
    juce::Slider frequencySlider8;
    juce::Label  frequencyLabel8;
    juce::Slider frequencySlider9;
    juce::Label  frequencyLabel9;
    juce::Slider frequencySlider10;
    juce::Label  frequencyLabel10;
    
    float sampleRateVar;
    
    const float windowBorder_x = 20.f;
    const float windowBorder_y = 20.f;
    
    const float frequencySliderWidth = 100;
    const float frequencySliderHeight = 100;
    const float filterSliderWidth = 100;
    const float filterSliderHeight = 200;
    
    const float thumbnailWidth = 150;
    const float thumbnailHeight = 150;
    
    const float buttonWidth = 100;
    const float buttonHeight = 25;
    const float offset = 20;
    const float labelOffset = 100;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (MainComponent)
};
