

#include <iostream>
#include <stdio.h>
#include "canopyFlow.h"
#include "windAdjustmentFactor.h"

int main() {

    static const double pi = acos(-1.0);
    double albiniWAF, massmanWAF;

    double canopyHeight = 10.0;
    double canopyCover = 0.5;
    double crownRatio = 0.7;
    double dragCoef = 0.2;
    double groundRoughness = 0.0075;
    double fuelBedDepth = 0.5;
    double inputHeight = canopyHeight + 6.096;
    double midFlameHeight = fuelBedDepth*1.5;   //This assumes the flame height is equal to the fuel bed depth height (above the top of the fuel bed)

    double LAI = canopyHeight * canopyCover * 10.6955 / (6.0 * pi);


    canopyFlow wind;
    //-------Uniform Distribution-----------------------
    wind.C = new canopy_uniform_distribution(crownRatio);
    wind.C->leafAreaIndex = LAI;
    wind.C->canopyHeight = canopyHeight;                        //canopy height (m)
    wind.C->dragCoefAth = dragCoef;

    wind.C->z0g = groundRoughness;                      //ground roughness length (m), model is somewhat sensitive to this
    wind.C->numNodes = 10001;                 //number of cells to use for numerical integration
    wind.C->initialize();
    wind.computeWind();
    massmanWAF = wind.get_windAdjustmentFactorShelteredMidFlame(inputHeight, midFlameHeight);

    WindAjustmentFactor albini;
    albiniWAF = albini.calculateWindAdjustmentFactor(canopyCover, canopyHeight, crownRatio, fuelBedDepth);

    std::cout << "Massman WAF \t= \t" << massmanWAF << std::endl;
    std::cout << "Albini WAF \t= \t" << albiniWAF << std::endl;





    double inputSpeed = 10.0;
    double lowLAI = 0.001;
    double highLAI = 10.0;
    int profileType = 1;    //  0 => sheltered;  1 => unsheltered;

//    wind.plotDimensionalWind(inputSpeed, inputHeight);
//    wind.plotWAFvsCdLAI(inputHeight, midFlameHeight, lowLAI, highLAI, profileType);
//    wind.plotz0ohvsCdLAI(inputHeight, lowLAI, highLAI);
//    wind.plotdohvsCdLAI(inputHeight, lowLAI, highLAI);
//    wind.plotz0ohvsone_doh(inputHeight, lowLAI, highLAI);
//    wind.plotz0ohvsdoh(inputHeight, lowLAI, highLAI);

    std::cout << "Done!" << std::endl;
    return 0;
}
