#pragma once

#include <JuceHeader.h>

class equalizer  : public juce::Component
{
public:
    equalizer();
    ~equalizer() override;

    void paint (juce::Graphics&) override;
    void resized() override;

private:
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (equalizer)
};
