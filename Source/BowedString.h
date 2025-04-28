/*
  ==============================================================================

    BowedString.h
    Created: 25 Apr 2025 7:39:31pm
    Author:  where

  ==============================================================================
*/

#pragma once
#include <JuceHeader.h>
#include "Globals.h"

class BowedString : public juce::Component
{
public:
    BowedString(juce::NamedValueSet parameters, double K);
    ~BowedString();

    long calcCounter = 0;
    void calculateBow();
    void calculatePluck();
    void updateStates();
    void excitePluck();
    void exciteBow();

    void mouseDown(const juce::MouseEvent& e) override;
    void mouseUp(const juce::MouseEvent& e) override;

    bool shouldExcite() { return excitationFlag; };

    double getOutput(double Lratio)
    {
        int index = round(N * Lratio);
        return u[1][index]*100;
    }

    void paint(juce::Graphics&) override;
    juce::Path visualiseState(juce::Graphics& g, double visualScaling);

    void setConnectionDivisionTerm(double cDT) { connectionDivisionTerm = cDT; };
    double getConnectionDivisionTerm() { return connectionDivisionTerm; };

    double Jl0(int l);

private:
    //BOW VARIABLES
        //model params
        double L, rho, A, T, E, I, c, kappa, sigma0, sigma1, lambda, lambdaSq, h, k, muSq, cSq, kappaSq;
        int N;
        //scheme vars
        double Adiv, B0, B1, B2, C0, C1, S0, S1, Bss, phi;

        // An (N+1) x 3 'matrix' containing the state of the system at all time-steps
        std::vector<std::vector<double>> uStates;

        // vector of pointers that point to state vectors
        std::vector<double*> u;

        // flag to tell MainComponent whether to excite the scheme or not
        bool excitationFlag = false;

        // initialise location of excitation
        double excitationLoc = 0.5;

        bool clamped = true;

        double connectionDivisionTerm = -1;

    //BOW VARIABLES
        // bow variables
        double xB, vB, fB, a, FB, BM;

        // NR variables
        double vRel, vPrev, tol, bL, eps, xNext, b, uI, uIPrev, uI1, uI2, uIM1, uIM2, uIPrev1, uIPrevM1;
        int i;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(BowedString)
};