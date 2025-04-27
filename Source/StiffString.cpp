/*
  ==============================================================================

    StiffString.cpp
    Created: 9 Apr 2025 2:31:13pm
    Author:  where

  ==============================================================================
*/

#include <JuceHeader.h>
#include "StiffString.h"
#include <iostream>
#include <cmath>

StiffString::StiffString(juce::NamedValueSet parameters)
{
    myWheel = std::make_unique<Wheel>(parameters);

    fs = *parameters.getVarPointer("fs");
    k = *parameters.getVarPointer("k");
    a = *parameters.getVarPointer("a");
    fB = *parameters.getVarPointer("fB");

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
    kappa = sqrt((E * I) / (rho * A));

    double stabilityTerm = c * c * k * k + 4.0 * sigma1 * k;

    //h = sqrt(((pow(c,2) * pow(k,2)) + 4 * sigma1 * k + sqrt(pow((pow(c,2) * pow(k,2) + 4 * sigma1 * k),2) + (16 * pow(kappa,2) * pow(k,2))) / 2));
    h = sqrt(0.5 * (stabilityTerm + sqrt((stabilityTerm * stabilityTerm) + 16.0 * pow(kappa,2) * k * k)));
    N = floor(L / h);
    h = L / N;

    // Initialise vectors
    uStates = std::vector<std::vector<double>>(3,
        std::vector<double>(N + 1, 0));

    /*  Make u pointers point to the first index of the state vectors.
        To use u (and obtain a vector from the state vectors) use indices like u[n][l] where,
             - n = 0 is u^{n+1},
             - n = 1 is u^n, and
             - n = 2 is u^{n-1}.
        Also see calculateScheme()
     */

     // Initialise pointer vector
    u.resize(3, nullptr);

    // Make set memory addresses to first index of the state vectors.
    for (int i = 0; i < 3; ++i)
        u[i] = &uStates[i][0];

    lambda = c * k / h;
    lambdaSq = pow(lambda, 2);
    muSq = ((kappa * kappa) * (k * k)) / (h * h * h * h);

    S0 = sigma0 * k;
    S1 = (2.0 * sigma1 * k) / (h * h);

    // Scheme coefficients
    B0 = 2.0 - 2.0 * lambdaSq - 6.0 * muSq - 2.0 * S1; // u_l^n
    B1 = lambdaSq + 4.0 * muSq + S1;                   // u_{l+-1}^n
    B2 = -muSq;                                        // u_{l+-2}^n
    C0 = -1.0 + S0 + 2.0 * S1;                         // u_l^{n-1}
    C1 = -S1;                                          // u_{l+-1}^{n-1}

    Bss = 2.0 - 2.0 * lambdaSq - 5.0 * muSq - 2.0 * S1;

    Adiv = 1.0 / (1.0 + S0);                           // u_l^{n+1}

    // Divide by u_l^{n+1} term
    B0 *= Adiv;
    B1 *= Adiv;
    B2 *= Adiv;
    C0 *= Adiv;
    C1 *= Adiv;
    Bss *= Adiv;

    phi = (pow(k, 2) * (1/h) * FB * (sqrt(2 * a) * vRel * exp(-a * pow(vRel, 2) + 0.5)));
    phi *= Adiv;

    setConnectionDivisionTerm(k * k / (rho * A * h * (1.0 + sigma0 * k)));
}

StiffString::~StiffString()
{

}

void StiffString::calculate()
{
    for (int l = 2; l < N - 1; ++l) // clamped boundaries
        u[0][l] = B0 * u[1][l] + B1 * (u[1][l + 1] + u[1][l - 1]) + B2 * (u[1][l + 2] + u[1][l - 2])
            + C0 * u[2][l] + C1 * (u[2][l + 1] + u[2][l - 1]);
    /*
    for (int l = 2; 1 < N - 1; ++l)
    {
        //u[0][l] = (2 * u[1][l] - u[2][l] + lambdaSq * (u[1][l + 1] - 2 * u[1][l] + u[1][l-1])
        //    - (pow(k,2)*pow(kappa,2)/pow(h,4))*(u[1][l+2] - 4*u[1][l+1] + 6*u[1][l] - 4*u[1][l-1] + u[1][l-2]) + sigma0*k*u[2][l]
        //    + (2 * sigma1 * k / pow(h, 2)) * (u[1][l + 1] - 2 * u[1][l] + u[1][l - 1] - u[2][l + 1] + 2 * u[2][l] - u[2][l - 1])
        //    - (pow(k, 2) * (1 / h) * FB * (sqrt(2 * a) * vRel * exp(-a * pow(vRel, 2) + 0.5)))
        //    ) / (1 + sigma0 * k);
        u[0][l] = B0 * u[1][l] + B1 * (u[1][l + 1] + u[1][l - 1]) + B2 * (u[1][l + 2] + u[1][l - 2])
            + C0 * u[2][l] + C1 * (u[2][l + 1] + u[2][l - 1]) + phi;
    }
    */
    ++calcCounter;
}

void StiffString::updateStates()
{
    // Do a pointer-switch. MUCH quicker than copying two entire state vectors every time-step.
    double* uTmp = u[2];
    u[2] = u[1];
    u[1] = u[0];
    u[0] = uTmp;
}

void StiffString::excitePluck()
{
    //// Arbitrary excitation function (raised cosine) ////

    // width (in grid points) of the excitation
    double width = 10;

    // make sure we're not going out of bounds at the left boundary
    int start = std::max(floor((N + 1) * excitationLoc) - floor(width * 0.5), 1.0);

    for (int l = 0; l < width; ++l)
    {
        // make sure we're not going out of bounds at the right boundary (this does 'cut off' the raised cosine)
        if (l + start > (clamped ? N - 2 : N - 1))
            break;

        u[1][l + start] += 0.5 * (1 - cos(2.0 * juce::double_Pi * l / (width - 1.0)));
        u[2][l + start] += 0.5 * (1 - cos(2.0 * juce::double_Pi * l / (width - 1.0)));
    }
    // Disable the excitation flag to only excite once
    excitationFlag = false;
}

void StiffString::mouseDown(const juce::MouseEvent& e)
{
    // Get the excitation location as a ratio between the x-location of the mouse-click and the width of the app
    excitationLoc = e.x / static_cast<double> (getWidth());

    // Activate the excitation flag to be used by the MainComponent to excite the string
    excitationFlag = true;
}

void StiffString::paint(juce::Graphics& g)
{
    // clear the background
    g.fillAll(getLookAndFeel().findColour(juce::ResizableWindow::backgroundColourId));

    // choose your favourite colour
    g.setColour(juce::Colours::cyan);

    // draw the state
    g.strokePath(visualiseState(g, 100), juce::PathStrokeType(2.0f));

}

juce::Path StiffString::visualiseState(juce::Graphics& g, double visualScaling)
{
    // String-boundaries are in the vertical middle of the component
    double stringBoundaries = getHeight() / 2.0;

    // initialise path
    juce::Path stringPath;

    // start path
    stringPath.startNewSubPath(0, -u[1][0] * visualScaling + stringBoundaries);

    double spacing = getWidth() / static_cast<double>(N);
    double x = spacing;

    for (int l = 1; l <= N; l++) // if you don't save the boundaries use l < N
    {
        // Needs to be -u, because a positive u would visually go down
        float newY = -u[1][l] * visualScaling + stringBoundaries;

        // if we get NAN values, make sure that we don't get an exception
        if (isnan(newY))
            newY = 0;

        stringPath.lineTo(x, newY);
        x += spacing;
    }
    // if you don't save the boundaries, and add a stringPath.lineTo (x, getWidth()) here to end the statedrawing

    return stringPath;
}

void StiffString::exciteBow()
{
    vRel = myWheel->calculate(u);
}