/*
  ==============================================================================

    Wheel.h
    Created: 11 Apr 2025 2:39:52pm
    Author:  where

  ==============================================================================
*/

#pragma once
#include <JuceHeader.h>
#include "Globals.h"

class Wheel : public juce::Component
{
public:
    Wheel(juce::NamedValueSet parameters);
    ~Wheel();

    double calculate(std::vector<double*>& u);

private:
    // bow variables
    double xB, vB, fB, a, FB, k;

    // NR variables
    double vRel, vPrev, tol, bL, eps, xNext, b, g, gD, uI, uIPrev, uI1, uI2, uIM1, uIM2, uIPrev1, uIPrevM1;
    int maxIterations, i;

    // string variables still needed in the NR solve
    double rho, A, sigma0, sigma1, c, kappa, T, I, E, L, h, cSq, kappaSq, connectionDivisionTerm;
    int N;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(Wheel)
};