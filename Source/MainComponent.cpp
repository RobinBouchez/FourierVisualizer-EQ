#include "MainComponent.h"

//==============================================================================
MainComponent::MainComponent()
: state(Stopped), openButton("Open"), playButton("Play"), stopButton("Stop"),
thumbnailCache (5),                            // [4]
        thumbnail (512, formatManager, thumbnailCache)
{
    // Make sure you set the size of the component after
    // you add any child components.
    setSize (1600, 800);

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
}

void MainComponent::getNextAudioBlock (const juce::AudioSourceChannelInfo& bufferToFill)
{
    if (playSource.get() == nullptr)
    {
        bufferToFill.clearActiveBufferRegion();
        return;
    }
    transportSource.getNextAudioBlock (bufferToFill);
}

void MainComponent::releaseResources()
{
    // This will be called when the audio device stops, or when it is being
    // restarted due to a setting change.

    // For more details, see the help for AudioProcessor::releaseResources()
}

//==============================================================================
void MainComponent::paint (juce::Graphics& g)
{
    // (Our component is opaque, so we must completely fill the background with a solid colour)
    g.fillAll(juce::Colours::black);
    
    juce::Rectangle<int> thumbnailBounds (20, 550, getWidth() - 40, 180);
     
            if (thumbnail.getNumChannels() == 0)
                paintIfNoFileLoaded (g, thumbnailBounds);
            else
                paintIfFileLoaded (g, thumbnailBounds);
    
    juce::Rectangle<int> spectraBounds (20, 20, getWidth() - 40, 300);
    
    g.setColour (juce::Colours::grey);
    g.drawLine(spectraBounds.getX(), spectraBounds.getY() + spectraBounds.getHeight(), spectraBounds.getX() + spectraBounds.getWidth(), spectraBounds.getY() + spectraBounds.getHeight(), 1);
    
    
//    juce::Slider s;
//    
//    
//    getLookAndFeel().drawLinearSlider(g,
//    20, 60,
//                                     200,200,
//                                     .5f, 0.0f, 1.0f,
//                                      juce::Slider::SliderStyle::LinearVertical,
//                                     s);
    
    
}

void MainComponent::resized()
{
    int borderOffset = 20;
    int width = 100;
    int height = 25;
    openButton.setBounds(borderOffset,getHeight() - height - borderOffset,width, height);
    playButton.setBounds(getWidth() / 2 - width,getHeight() - height - borderOffset,width, height);
    stopButton.setBounds(getWidth() / 2 + width / 2,getHeight() - height - borderOffset,width, height);
}


void MainComponent::openButtonClicked()
  {
    chooser = std::make_unique<juce::FileChooser> ("Please select the file you want to load...",
                                                   juce::File::getSpecialLocation (juce::File::userDesktopDirectory),
                                                "*.mp3");
 
    
    //juce::FileChooser c ("", )
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

    auto audioLength = (float) thumbnail.getTotalLength();                               // [12]
    thumbnail.drawChannels (g, thumbnailBounds, 0.0, audioLength, 1.0f);

    g.setColour (juce::Colours::red);

    auto audioPosition = (float) transportSource.getCurrentPosition();
    auto drawPosition = (audioPosition / audioLength) * (float) thumbnailBounds.getWidth()
                        + (float) thumbnailBounds.getX();                                // [13]
    g.drawLine (drawPosition, (float) thumbnailBounds.getY(), drawPosition,
                (float) thumbnailBounds.getBottom(), 2.0f);                              // [14]
}
