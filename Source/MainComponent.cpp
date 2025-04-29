#include "MainComponent.h"
#include "Globals.h"
#include <JuceHeader.h>

//==============================================================================
MainComponent::MainComponent()
{
    // Make sure you set the size of the component after
    // you add any child components.
    setSize (800, 600);
    setWantsKeyboardFocus(true);
    addKeyListener(this);

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
    // parameters you'll use to initialise more than one other parameter should be defined here
    double r = 0.0005;

    // string parameters
    //parameters.set("L", 1.0);
    parameters.set("rho", 7850.0);
    parameters.set("A", r * r * juce::double_Pi);
    parameters.set("T", 2000.0);
    parameters.set("E", 2e11); // maybe increase E for more inharmonies
    parameters.set("I", r * r * r * r * juce::double_Pi * 0.25);
    parameters.set("sigma0", 1.0);
    parameters.set("sigma1", 0.005); // steel to nylon relationship
    // bow parameters

    parameters.set("xB", 0.13); // closer to edge for more irregular
    parameters.set("fB", 1.0);
    parameters.set("vB", 0.2);
    parameters.set("a", 100.0);

    //// Initialise an instance of the SimpleString class ////
    //myWheel = std::make_unique<Wheel>(parameters);
    //myStiffString = std::make_unique<StiffString>(parameters);
    //bowedString = std::make_unique<BowedString>(parameters, 1/sampleRate);

    droneString1 = std::make_unique<BowedString>(parameters, 1 / sampleRate, 116.54); //Bb - root
    droneString2 = std::make_unique<BowedString>(parameters, 1 / sampleRate, 174.61); //F - perfect 5th
    melodyString = std::make_unique<BowedString>(parameters, 1 / sampleRate, 466.16); //Bb4 - initialize melody string one octave up
    melodyString->refreshParameters(parameters, 233.08); // before play, tune down the melody string to the roob Bb3

    //addAndMakeVisible(myStiffString.get()); // add the string to the application
    //addAndMakeVisible(myWheel.get()); // add the string to the application
    //addAndMakeVisible(bowedString.get());
    addAndMakeVisible(droneString1.get());
    addAndMakeVisible(droneString2.get());
    addAndMakeVisible(melodyString.get());


    addAndMakeVisible(melodySlider);
    melodySlider.setRange(50.0, 1000.0); // [1]
    melodySlider.addListener(this); // [3]

    addAndMakeVisible(melodyLabel);
    melodyLabel.setText("melody", juce::dontSendNotification);
    melodyLabel.attachToComponent(&melodySlider, true); // [4]


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
    //if (bowedString->shouldExcite())
        //bowedString->exciteBow();
        //bowedString->excitePluck();

    for (int i = 0; i < bufferToFill.numSamples; ++i)
    {
        //myStiffString->calculate();
        //myStiffString->updateStates();
        //if (bowedString->shouldExcite())
        //{
        //    bowedString->exciteBow();
        //}
        if (droneString1->shouldExcite())
        {
            droneString1->exciteBow();
            droneString2->exciteBow();
            melodyString->exciteBow();
        }
        if (droneString2->shouldExcite())
        {
            droneString2->exciteBow();
        }
        if (melodyString->shouldExcite())
        {
            melodyString->exciteBow();
        }

        //bowedString->calculateBow();
        //bowedString->updateStates();
        droneString1->calculateBow();
        droneString1->updateStates();
        droneString2->calculateBow();
        droneString2->updateStates();
        melodyString->calculateBow();
        melodyString->updateStates();

        //output = myStiffString->getOutput(0.8); // get output at 0.8L of the string
        //output = bowedString->getOutput(0.08); // get output at 0.8L of the string
        // move  outpout location closer to edge for more irregularities

        double outputLoc = 0.08;
        output = (droneString1->getOutput(outputLoc) + droneString2->getOutput(outputLoc)
            + melodyString->getOutput(outputLoc));

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
    auto sliderLeft = 120;
    melodySlider.setBounds(sliderLeft, 20, getWidth() - sliderLeft - 10, 20);

    auto bounds = getLocalBounds();
    int numStrings = 3;
    int stringHeight = bounds.getHeight() / numStrings;

    if (droneString1 != nullptr)
        droneString1->setBounds(bounds.removeFromTop(stringHeight));
    if (droneString2 != nullptr)
        droneString2->setBounds(bounds.removeFromTop(stringHeight));
    if (melodyString != nullptr)
        melodyString->setBounds(bounds.removeFromTop(stringHeight));
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