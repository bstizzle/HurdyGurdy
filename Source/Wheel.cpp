/*
  ==============================================================================

    Wheel.cpp
    Created: 11 Apr 2025 2:40:04pm
    Author:  where

  ==============================================================================
*/

#include <JuceHeader.h>
#include "Wheel.h"

Wheel::Wheel(juce::NamedValueSet parameters)
{
    k = *parameters.getVarPointer("k");

    // bow parameters
    xB = *parameters.getVarPointer("xB");
    fB = *parameters.getVarPointer("fB");
    vB = *parameters.getVarPointer("vB");
    a = *parameters.getVarPointer("a");

    // parameters from string
    L = *parameters.getVarPointer("L");
    rho = *parameters.getVarPointer("rho");
    A = *parameters.getVarPointer("A");
    T = *parameters.getVarPointer("T");
    E = *parameters.getVarPointer("E");
    I = *parameters.getVarPointer("I");
    sigma0 = *parameters.getVarPointer("sigma0");
    sigma1 = *parameters.getVarPointer("sigma1");

    FB = (fB / (rho * A));
    c = sqrt(T / (rho * A));
    cSq = pow(c, 2);
    kappa = sqrt((E * I) / (rho * A));
    kappaSq = pow(kappa, 2);

    double stabilityTerm = c * c * k * k + 4.0 * sigma1 * k;

    //h = sqrt(((pow(c,2) * pow(k,2)) + 4 * sigma1 * k + sqrt(pow((pow(c,2) * pow(k,2) + 4 * sigma1 * k),2) + (16 * pow(kappa,2) * pow(k,2))) / 2));
    h = sqrt(0.5 * (stabilityTerm + sqrt((stabilityTerm * stabilityTerm) + 16.0 * pow(kappa, 2) * k * k)));
    N = floor(L / h);
    h = L / N;

    // NR parameters
    vRel = 0; // starting point
    vPrev = 0;
    tol = 1e-4; // threshold
    maxIterations = 100;

    connectionDivisionTerm = k * k / (rho * A * h * (1.0 + sigma0 * k));
}

Wheel::~Wheel()
{
}

double Wheel::calculate(std::vector<double*>& u)
{
    bL = round(xB / h);

    uI = Globals::interpolation(u[1], bL);
    uIPrev = Globals::interpolation(u[2], bL);
    uI1 = Globals::interpolation(u[1], bL + 1);
    uI2 = Globals::interpolation(u[1], bL + 2);
    uIM1 = Globals::interpolation(u[1], bL - 1);
    uIM2 = Globals::interpolation(u[1], bL - 2);
    uIPrev1 = Globals::interpolation(u[2], bL + 1);
    uIPrevM1 = Globals::interpolation(u[2], bL - 1);

    b = (-2 / pow(k, 2)) * (uI - uIPrev) - (cSq / pow(h, 2) * (uI1 - 2 * uI + uIM1))
        + (kappaSq / pow(h, 4) * (uI2 - 4 * uI1 + 6 * uI - 4 * uIM1 + uIM2))
        - (2 * sigma1 / (k * pow(h, 2)) * (uI1 - 2 * uI + uIM1 - uIPrev1 + 2 * uIPrev - uIPrevM1))
        + (2 / k + 2 * sigma0) * vB;

    // g = (2 / k + 2 * sigma0) * vRel + FB * (1 / h) * sqrt(2 * a) * vRel * exp(-a * pow(vRel, 2) + 0.5) + b;
    // gD = 2 / k + 2 * sigma0 + (1 / h) * FB * sqrt(2 * a) * (1 - 2 * a * pow(vRel, 2)) * exp(-a * pow(vRel, 2) + 0.5);

    // loop until a maximum number of iterations
    i = 0;
    eps = 1;
    while (eps > tol&& i < 100)
    {
        vRel = vPrev - ((2 / k + 2 * sigma0) * vPrev + FB * (1 / h) * sqrt(2 * a) * vPrev * exp(-a * pow(vPrev, 2) + 0.5) + b)
                         / (2 / k + 2 * sigma0 + (1 / h) * FB * sqrt(2 * a) * (1 - 2 * a * pow(vPrev, 2)) * exp(-a * pow(vPrev, 2) + 0.5));

        eps = abs(vRel - vPrev);
        vRel = vPrev;
        i = i + 1;
    }

    double excitation = connectionDivisionTerm * FB * sqrt(2.0*a) * vRel * exp(-a * vRel * vRel) * exp(0.5);
    Globals::extrapolation(u[0], bL, -excitation);

    return vRel;
}