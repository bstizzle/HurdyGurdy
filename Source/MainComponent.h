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
            //melodyString1->tune(melodySlider.getValue());
            melodyString1->refreshParameters(parameters, melodySlider.getValue());
        }
    }

    bool keyPressed(const juce::KeyPress& key, juce::Component* originatingComponent) override {
        if (key == juce::KeyPress('a'))
        {
            melodyString1->refreshParameters(parameters, rootFreq * pow(1.05946, 12));
            melodyString2->refreshParameters(parameters, rootFreq * pow(1.05946, 24));
            return true; // Indicates the key was handled
        }
        else if (key == juce::KeyPress('w')) 
        {
            melodyString1->refreshParameters(parameters, rootFreq * pow(1.05946, 13));
            melodyString2->refreshParameters(parameters, rootFreq * pow(1.05946, 25));
            return true; // Indicates the key was handled
        }
        else if (key == juce::KeyPress('s'))
        {
            melodyString1->refreshParameters(parameters, rootFreq * pow(1.05946, 14));
            melodyString2->refreshParameters(parameters, rootFreq * pow(1.05946, 26));
            return true; // Indicates the key was handled
        }
        else if (key == juce::KeyPress('e'))
        {
            melodyString1->refreshParameters(parameters, rootFreq * pow(1.05946, 15));
            melodyString2->refreshParameters(parameters, rootFreq * pow(1.05946, 27));
            return true; // Indicates the key was handled
        }
        else if (key == juce::KeyPress('d'))
        {
            melodyString1->refreshParameters(parameters, rootFreq * pow(1.05946, 16));
            melodyString2->refreshParameters(parameters, rootFreq * pow(1.05946, 28));
            return true; // Indicates the key was handled
        }
        else if (key == juce::KeyPress('r'))
        {
            melodyString1->refreshParameters(parameters, rootFreq * pow(1.05946, 17));
            melodyString2->refreshParameters(parameters, rootFreq * pow(1.05946, 29));
            return true; // Indicates the key was handled
        }
        else if (key == juce::KeyPress('f'))
        {
            melodyString1->refreshParameters(parameters, rootFreq * pow(1.05946, 18));
            melodyString2->refreshParameters(parameters, rootFreq * pow(1.05946, 30));
            return true; // Indicates the key was handled
        }
        else if (key == juce::KeyPress('t'))
        {
            melodyString1->refreshParameters(parameters, rootFreq * pow(1.05946, 19));
            melodyString2->refreshParameters(parameters, rootFreq * pow(1.05946, 31));
            return true; // Indicates the key was handled
        }
        else if (key == juce::KeyPress('g'))
        {
            melodyString1->refreshParameters(parameters, rootFreq * pow(1.05946, 20));
            melodyString2->refreshParameters(parameters, rootFreq * pow(1.05946, 32));
            return true; // Indicates the key was handled
        }
        else if (key == juce::KeyPress('y'))
        {
            melodyString1->refreshParameters(parameters, rootFreq * pow(1.05946, 21));
            melodyString2->refreshParameters(parameters, rootFreq * pow(1.05946, 33));
            return true; // Indicates the key was handled
        }
        else if (key == juce::KeyPress('h'))
        {
            melodyString1->refreshParameters(parameters, rootFreq * pow(1.05946, 22));
            melodyString2->refreshParameters(parameters, rootFreq * pow(1.05946, 34));
            return true; // Indicates the key was handled
        }
        else if (key == juce::KeyPress('u'))
        {
            melodyString1->refreshParameters(parameters, rootFreq * pow(1.05946, 23));
            melodyString2->refreshParameters(parameters, rootFreq * pow(1.05946, 35));
            return true; // Indicates the key was handled
        }
        else if (key == juce::KeyPress('j'))
        {
            melodyString1->refreshParameters(parameters, rootFreq * pow(1.05946, 24));
            melodyString2->refreshParameters(parameters, rootFreq * pow(1.05946, 36));
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
    std::unique_ptr<BowedString> melodyString1;
    std::unique_ptr<BowedString> melodyString2;

    //// Set the paramters ///
    juce::NamedValueSet parameters;

    juce::Slider melodySlider;
    juce::Label melodyLabel;

    double rootFreq = 98.00;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (MainComponent)
};
