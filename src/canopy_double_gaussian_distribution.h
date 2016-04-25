#ifndef CANOPY_DOUBLE_GAUSS_H
#define CANOPY_DOUBLE_GAUSS_H

#include <math.h>
#include <cstring>
#include <plstream.h>
#include "canopy.h"

/*
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Class that stores canopy density, etc. information for a double Gaussian distribution.
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*/

class canopy_double_gaussian_distribution : public canopy
{

public:
    canopy_double_gaussian_distribution();
    canopy_double_gaussian_distribution(double heightMaxFoliageDist_, double standardDevFoliageUpper_, double standardDevFoliageLower_);
    canopy_double_gaussian_distribution(canopy_double_gaussian_distribution &rhs);
    canopy_double_gaussian_distribution &operator=(const canopy_double_gaussian_distribution &rhs);
    ~canopy_double_gaussian_distribution();

    //inputs
    double heightMaxFoliage;    //height of maximum foliage distribution for the normal distribution, normalized from 0 to 1
    double standardDevFoliageUpper;  //standard deviation of foliage distribution for the upper part
    double standardDevFoliageLower;  //standard deviation of foliage distribution for the lower part

protected:
    void compute_haz();
};

#endif // CANOPY_DOUBLE_GAUSS_H
