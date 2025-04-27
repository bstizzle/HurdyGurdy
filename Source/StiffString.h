/*
  ==============================================================================

    StiffString.h
    Created: 9 Apr 2025 2:31:13pm
    Author:  where

  ==============================================================================
*/

#pragma once
#include <JuceHeader.h>
#include "Wheel.h"
#include "Globals.h"

class StiffString : public juce::Component
{
public:
    StiffString(juce::NamedValueSet parameters);
    ~StiffString();

    long calcCounter = 0;
    void calculate();
    void updateStates();
    void excitePluck();
    void exciteBow();

    void mouseDown(const juce::MouseEvent& e) override;

    bool shouldExcite() { return excitationFlag; };

    double getOutput(double Lratio)
    {
        return u[1][static_cast<int> (round(N * Lratio))];
    }

    void paint(juce::Graphics&) override;
    juce::Path visualiseState(juce::Graphics& g, double visualScaling);

    void setConnectionDivisionTerm(double cDT) { connectionDivisionTerm = cDT; };
    double getConnectionDivisionTerm() { return connectionDivisionTerm; };

private:
    //model params
    double L, r, rho, A, T, E, I, c, kappa, sigma0, sigma1, lambda, lambdaSq, h, k, fs, muSq, a, fB, FB;
    int N;
    //scheme vars
    double Adiv, B0, B1, B2, C0, C1, S0, S1, Bss, phi, Jl0;

    // An (N+1) x 3 'matrix' containing the state of the system at all time-steps
    std::vector<std::vector<double>> uStates;

    // vector of pointers that point to state vectors
    std::vector<double*> u;

    // flag to tell MainComponent whether to excite the scheme or not
    bool excitationFlag = true;

    // initialise location of excitation
    double excitationLoc = 0.5;

    bool clamped = true;

    double connectionDivisionTerm = -1;

    // bow parameters
    double vRel;
    std::unique_ptr<Wheel> myWheel;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(StiffString)
};