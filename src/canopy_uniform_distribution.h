#ifndef CANOPY_UNIFORM_H
#define CANOPY_UNIFORM_H

#include <math.h>
#include <cstring>
#include <plstream.h>
#include "canopy.h"

/*
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Class that stores canopy density, etc. information for a uniform distribution.
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*/

class canopy_uniform_distribution : public canopy
{

public:
    canopy_uniform_distribution();
    canopy_uniform_distribution(double crownRatio_);
    canopy_uniform_distribution(canopy_uniform_distribution &rhs);
    canopy_uniform_distribution &operator=(const canopy_uniform_distribution &rhs);
    ~canopy_uniform_distribution();

    //inputs
    double crownRatio;  //fraction of canopy height with foliage [0 1]

protected:
    void compute_haz();
};

#endif // CANOPY_UNIFORM_H
