#pragma once

#include <JuceHeader.h>
#include "SpectrogramComponent.h"

//==============================================================================
/*
    This component lives inside our window, and this is where you should put all
    your controls and content.
*/
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
    
    enum
    {
        fftOrder  = 11,             // [1]
        fftSize   = 1 << fftOrder,  // [2]
        scopeSize = 512             // [3]
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
    

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (MainComponent)
};
