#pragma once

#include <JuceHeader.h>
#include "BowedString.h"

//==============================================================================
/*
    This component lives inside our window, and this is where you should put all
    your controls and content.
*/
class MainComponent  : public juce::AudioAppComponent, public juce::Timer, public juce::KeyListener
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

    bool keyPressed(const juce::KeyPress& key, juce::Component* originatingComponent) override {
        if (key == juce::KeyPress('a'))
        {
            //melodyString1->refreshParameters(parameters, rootFreq * pow(1.05946, 12));
            //melodyString2->refreshParameters(parameters, rootFreq * pow(1.05946, 24));
            melodyString1->refreshParameters(parameters, lowMelodyFreq);
            melodyString2->refreshParameters(parameters, highMelodyFreq);
            return true; // Indicates the key was handled
        }
        else if (key == juce::KeyPress('w')) 
        {
            //melodyString1->refreshParameters(parameters, rootFreq * pow(1.05946, 13));
            //melodyString2->refreshParameters(parameters, rootFreq * pow(1.05946, 25));
            melodyString1->refreshParameters(parameters, lowMelodyFreq * pow(1.05946, 1) - pow(1.15, 1));
            melodyString2->refreshParameters(parameters, highMelodyFreq * pow(1.05946, 1) - pow(1.3, 1));
            return true; // Indicates the key was handled
        }
        else if (key == juce::KeyPress('s'))
        {
            //melodyString1->refreshParameters(parameters, rootFreq * pow(1.05946, 14));
            //melodyString2->refreshParameters(parameters, rootFreq * pow(1.05946, 26));
            melodyString1->refreshParameters(parameters, lowMelodyFreq * pow(1.05946, 2) - pow(1.155, 2));
            melodyString2->refreshParameters(parameters, highMelodyFreq * pow(1.05946, 2) - pow(1.32, 2));
            return true; // Indicates the key was handled
        }
        else if (key == juce::KeyPress('e'))
        {
            //melodyString1->refreshParameters(parameters, rootFreq * pow(1.05946, 15));
            //melodyString2->refreshParameters(parameters, rootFreq * pow(1.05946, 27));
            melodyString1->refreshParameters(parameters, lowMelodyFreq * pow(1.05946, 3) - pow(1.16, 3));
            melodyString2->refreshParameters(parameters, highMelodyFreq * pow(1.05946, 3) - pow(1.335, 3));
            return true; // Indicates the key was handled
        }
        else if (key == juce::KeyPress('d'))
        {
            //melodyString1->refreshParameters(parameters, rootFreq * pow(1.05946, 16));
            //melodyString2->refreshParameters(parameters, rootFreq * pow(1.05946, 28));
            melodyString1->refreshParameters(parameters, lowMelodyFreq * pow(1.05946, 4) - pow(1.165, 4));
            melodyString2->refreshParameters(parameters, highMelodyFreq * pow(1.05946, 4) - pow(1.345, 4));
            return true; // Indicates the key was handled
        }
        else if (key == juce::KeyPress('r'))
        {
            //melodyString1->refreshParameters(parameters, rootFreq * pow(1.05946, 17));
            //melodyString2->refreshParameters(parameters, rootFreq * pow(1.05946, 29));
            melodyString1->refreshParameters(parameters, lowMelodyFreq * pow(1.05946, 5) - pow(1.17, 5));
            melodyString2->refreshParameters(parameters, highMelodyFreq * pow(1.05946, 5) - pow(1.355, 5));
            return true; // Indicates the key was handled
        }
        else if (key == juce::KeyPress('f'))
        {
            //melodyString1->refreshParameters(parameters, rootFreq * pow(1.05946, 18));
            //melodyString2->refreshParameters(parameters, rootFreq * pow(1.05946, 30));
            melodyString1->refreshParameters(parameters, lowMelodyFreq * pow(1.05946, 6) - pow(1.175, 6));
            melodyString2->refreshParameters(parameters, highMelodyFreq * pow(1.05946, 6) - pow(1.365, 6));
            return true; // Indicates the key was handled
        }
        else if (key == juce::KeyPress('t'))
        {
            //melodyString1->refreshParameters(parameters, rootFreq * pow(1.05946, 19));
            //melodyString2->refreshParameters(parameters, rootFreq * pow(1.05946, 31));
            melodyString1->refreshParameters(parameters, lowMelodyFreq * pow(1.05946, 7) - pow(1.18, 7));
            melodyString2->refreshParameters(parameters, highMelodyFreq * pow(1.05946, 7) - pow(1.375, 7));
            return true; // Indicates the key was handled
        }
        else if (key == juce::KeyPress('g'))
        {
            //melodyString1->refreshParameters(parameters, rootFreq * pow(1.05946, 20));
            //melodyString2->refreshParameters(parameters, rootFreq * pow(1.05946, 32));
            melodyString1->refreshParameters(parameters, lowMelodyFreq * pow(1.05946, 8) - pow(1.185, 8));
            melodyString2->refreshParameters(parameters, highMelodyFreq * pow(1.05946, 8) - pow(1.385, 8));
            return true; // Indicates the key was handled
        }
        else if (key == juce::KeyPress('y'))
        {
            //melodyString1->refreshParameters(parameters, rootFreq * pow(1.05946, 21));
            //melodyString2->refreshParameters(parameters, rootFreq * pow(1.05946, 33));
            melodyString1->refreshParameters(parameters, lowMelodyFreq * pow(1.05946, 9) - pow(1.19, 9));
            melodyString2->refreshParameters(parameters, highMelodyFreq * pow(1.05946, 9) - pow(1.4, 9));
            return true; // Indicates the key was handled
        }
        else if (key == juce::KeyPress('h'))
        {
            //melodyString1->refreshParameters(parameters, rootFreq * pow(1.05946, 22));
            //melodyString2->refreshParameters(parameters, rootFreq * pow(1.05946, 34));
            melodyString1->refreshParameters(parameters, lowMelodyFreq * pow(1.05946, 10) - pow(1.195, 10));
            melodyString2->refreshParameters(parameters, highMelodyFreq * pow(1.05946, 10) - pow(1.38, 10));
            return true; // Indicates the key was handled
        }
        else if (key == juce::KeyPress('u'))
        {
            //melodyString1->refreshParameters(parameters, rootFreq * pow(1.05946, 23));
            //melodyString2->refreshParameters(parameters, rootFreq * pow(1.05946, 35));
            melodyString1->refreshParameters(parameters, lowMelodyFreq * pow(1.05946, 11) - pow(1.2, 11));
            melodyString2->refreshParameters(parameters, highMelodyFreq * pow(1.05946, 11) - pow(1.38, 11));
            return true; // Indicates the key was handled
        }
        else if (key == juce::KeyPress('j'))
        {
            //melodyString1->refreshParameters(parameters, rootFreq * pow(1.05946, 24));
            //melodyString2->refreshParameters(parameters, rootFreq * pow(1.05946, 36));
            melodyString1->refreshParameters(parameters, lowMelodyFreq * pow(1.05946, 12) - pow(1.205, 12));
            melodyString2->refreshParameters(parameters, highMelodyFreq * pow(1.05946, 12) - pow(1.362, 12));
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

    juce::Label melodyLabel;

    //double rootFreq = 98.00;
    double rootFreq = 96.9;

    //double lowDroneFreq = 70.9;

    //double lowMelodyFreq = 196.00;
    double lowMelodyFreq = 190.71;

    //double highMelodyFreq = 392.00;
    double highMelodyFreq = 372.1;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (MainComponent)
};
