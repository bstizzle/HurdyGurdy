#pragma once

#include <JuceHeader.h>
#include "BowedString.h"

//==============================================================================
/*
    This component lives inside our window, and this is where you should put all
    your controls and content.
*/
class MainComponent  : public juce::AudioAppComponent, public juce::Timer, public juce::Slider::Listener, public juce::KeyListener
{
public:
    //==============================================================================
    MainComponent();
    ~MainComponent() override;

    //==============================================================================
    void prepareToPlay (int samplesPerBlockExpected, double sampleRate) override;
    void getNextAudioBlock (const juce::AudioSourceChannelInfo& bufferToFill) override;
    void releaseResources() override;

    //==============================================================================
    void paint (juce::Graphics& g) override;
    void resized() override;

    //double limit(double val); // limiter for your ears

    void timerCallback() override;

    void sliderValueChanged(juce::Slider* slider) override
    {
        if (slider == &melodySlider)
        {
            //melodyString->tune(melodySlider.getValue());
            melodyString->refreshParameters(parameters, melodySlider.getValue());
        }
    }

    bool keyPressed(const juce::KeyPress& key, juce::Component* originatingComponent) override {
        if (key == juce::KeyPress('c'))
        {
            melodyString->refreshParameters(parameters, 233.08);
            return true; // Indicates the key was handled
        }
        else if (key == juce::KeyPress('a')) 
        {
            melodyString->refreshParameters(parameters, 246.94);
            return true; // Indicates the key was handled
        }
        else if (key == juce::KeyPress('w'))
        {
            melodyString->refreshParameters(parameters, 261.63);
            return true; // Indicates the key was handled
        }
        else if (key == juce::KeyPress('s'))
        {
            melodyString->refreshParameters(parameters, 277.18);
            return true; // Indicates the key was handled
        }
        else if (key == juce::KeyPress('e'))
        {
            melodyString->refreshParameters(parameters, 293.66);
            return true; // Indicates the key was handled
        }
        else if (key == juce::KeyPress('d'))
        {
            melodyString->refreshParameters(parameters, 311.13);
            return true; // Indicates the key was handled
        }
        else if (key == juce::KeyPress('f'))
        {
            melodyString->refreshParameters(parameters, 329.63);
            return true; // Indicates the key was handled
        }
        else if (key == juce::KeyPress('t'))
        {
            melodyString->refreshParameters(parameters, 349.23);
            return true; // Indicates the key was handled
        }
        else if (key == juce::KeyPress('g'))
        {
            melodyString->refreshParameters(parameters, 369.99);
            return true; // Indicates the key was handled
        }
        else if (key == juce::KeyPress('y'))
        {
            melodyString->refreshParameters(parameters, 392);
            return true; // Indicates the key was handled
        }
        else if (key == juce::KeyPress('h'))
        {
            melodyString->refreshParameters(parameters, 415.30);
            return true; // Indicates the key was handled
        }

        return false;
    }


private:
    //==============================================================================
    // Your private member variables go here...
    //std::unique_ptr<StiffString> myStiffString;
    //std::unique_ptr<Wheel> myWheel;
    //std::unique_ptr<BowedString> bowedString;
    std::unique_ptr<BowedString> droneString1;
    std::unique_ptr<BowedString> droneString2;
    std::unique_ptr<BowedString> melodyString;

    //// Set the paramters ///
    juce::NamedValueSet parameters;

    juce::Slider melodySlider;
    juce::Label melodyLabel;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (MainComponent)
};
