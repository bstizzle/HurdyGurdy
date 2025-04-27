#include "MainComponent.h"
#include "Globals.h"
#include <JuceHeader.h>

//==============================================================================
MainComponent::MainComponent()
{
    // Make sure you set the size of the component after
    // you add any child components.
    setSize (800, 600);

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
}

MainComponent::~MainComponent()
{
    stopTimer();

    // This shuts down the audio device and clears the audio source.
    shutdownAudio();
}

//==============================================================================
void MainComponent::prepareToPlay (int samplesPerBlockExpected, double sampleRate)
{
    //// Set the paramters ///
    juce::NamedValueSet parameters;

    // parameters you'll use to initialise more than one other parameter should be defined here
    double r = 0.0005;

    // string parameters
    parameters.set("L", 1);
    parameters.set("rho", 7850);
    parameters.set("A", r * r * 3.1415);
    parameters.set("T", 1000);
    parameters.set("E", 2e11);
    parameters.set("I", r * r * r * r * 3.1415 * 0.25);
    parameters.set("sigma0", 1);
    parameters.set("sigma1", 0.005);
    // bow parameters

    parameters.set("xB", 0.125);
    parameters.set("fB", 1);
    parameters.set("vB", 0.2);
    parameters.set("a", 100);

    //// Initialise an instance of the SimpleString class ////
    //myWheel = std::make_unique<Wheel>(parameters);
    //myStiffString = std::make_unique<StiffString>(parameters);
    bowedString = std::make_unique<BowedString>(parameters, 1/sampleRate);

    //addAndMakeVisible(myStiffString.get()); // add the string to the application
    //addAndMakeVisible(myWheel.get()); // add the string to the application
    addAndMakeVisible(bowedString.get());

    // Call resized again as our components need a sample rate before they can get initialised.
    resized();

    startTimerHz(15); // start the timer (15 Hz is a nice tradeoff between CPU usage and update speed)
}

void MainComponent::getNextAudioBlock (const juce::AudioSourceChannelInfo& bufferToFill)
{
    // Your audio-processing code goes here!

    // For more details, see the help for AudioProcessor::getNextAudioBlock()

    // Right now we are not producing any data, in which case we need to clear the buffer
    // (to prevent the output of random noise)
    bufferToFill.clearActiveBufferRegion();

    int numChannels = bufferToFill.buffer->getNumChannels();

    // Get pointers to output locations
    float* const channelData1 = bufferToFill.buffer->getWritePointer(0, bufferToFill.startSample);
    float* const channelData2 = numChannels > 1 ? bufferToFill.buffer->getWritePointer(1, bufferToFill.startSample) : nullptr;

    float output = 0.0;

    std::vector<float* const*> curChannel{ &channelData1, &channelData2 };

    // only do control stuff out of the buffer (at least work with flags so that control doesn't interfere with the scheme calculation)
    //if (myStiffString->shouldExcite())
        //myStiffString->excitePluck();
        //myStiffString->exciteBow();
    if (bowedString->shouldExcite())
        bowedString->exciteBow();
        //bowedString->excitePluck();

    for (int i = 0; i < bufferToFill.numSamples; ++i)
    {
        //myStiffString->calculate();
        //myStiffString->updateStates();
        bowedString->calculateBow();
        //bowedString->calculatePluck();
        bowedString->updateStates();

        //output = myStiffString->getOutput(0.8); // get output at 0.8L of the string
        output = bowedString->getOutput(0.3); // get output at 0.8L of the string

        for (int channel = 0; channel < numChannels; ++channel)
            curChannel[channel][0][i] = Globals::limit(output);
    }
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

}

void MainComponent::resized()
{
    // This is called when the MainContentComponent is resized.
    // If you add any child components, this is where you should
    // update their positions.
    //if (myStiffString != nullptr)
      //  myStiffString->setBounds(getLocalBounds());
    if (bowedString != nullptr)
        bowedString->setBounds(getLocalBounds());
}
/*
double MainComponent::limit(double val)
{
    if (val < -1)
    {
        val = -1;
        return val;
    }
    else if (val > 1)
    {
        val = 1;
        return val;
    }
    return val;
}
*/
void MainComponent::timerCallback()
{
    repaint(); // update the graphics X times a second
}