

#include <iostream>
#include <stdio.h>
#include "canopyFlow.h"

int main() {

    static const double pi = acos(-1.0);

    double canopyHeight = 10.0;
    double canopyCover = 0.5;
    double crownRatio = 0.7;
    double dragCoef = 0.2;

    double LAI = canopyHeight * canopyCover * 10.6955 / (6 * pi);


    canopyFlow wind;
    //-------Uniform Distribution-----------------------
    wind.C = new canopy_uniform_distribution(crownRatio);
    wind.C->leafAreaIndex = LAI;
    wind.C->canopyHeight = canopyHeight;                        //canopy height (m)
    wind.C->dragCoefAth = dragCoef;

    wind.C->z0g = 0.0075;                      //ground roughness length (m), model is somewhat sensitive to this
    wind.C->numNodes = 10001;                 //number of cells to use for numerical integration
    wind.C->initialize();
    wind.computeWind();

    //std::cout << wind.get_windAdjustmentFactorShelteredMidFlame(6.096, 1.0) << std::endl;
    //std::cout << wind.get_windAdjustmentFactorShelteredIntegral(6.096, 2.0) << std::endl;
    //std::cout << wind.get_windAdjustmentFactorUnshelteredIntegral(6.096, 22.19) << std::endl;

    double inputHeight = wind.C->canopyHeight + 6.096;
    double midFlameHeight = 2.5;
    double inputSpeed = 10.0;
    double lowLAI = 0.001;
    double highLAI = 10.0;
    int profileType = 1;    //  0 => sheltered;  1 => unsheltered;

    wind.plotDimensionalWind(inputSpeed, inputHeight);
//    wind.plotWAFvsCdLAI(inputHeight, midFlameHeight, lowLAI, highLAI, profileType);
//    wind.plotz0ohvsCdLAI(inputHeight, lowLAI, highLAI);
//    wind.plotdohvsCdLAI(inputHeight, lowLAI, highLAI);
//    wind.plotz0ohvsone_doh(inputHeight, lowLAI, highLAI);
//    wind.plotz0ohvsdoh(inputHeight, lowLAI, highLAI);

    std::cout << "Done!" << std::endl;
    return 0;
}
