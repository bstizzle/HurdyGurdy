/*
  ==============================================================================

    Globals.h
    Created: 15 Apr 2025 11:34:24am
    Author:  where

  ==============================================================================
*/

#pragma once
#include <JuceHeader.h>
#include "Wheel.h"

namespace Globals
{
    static double limit(double val)
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

    static double interpolation(double* uVec, int bp)
    {
        return uVec[bp];
    }

    static void extrapolation(double* uVec, int bp, double val)
    {
        uVec[bp] = uVec[bp] + val;
    }
}