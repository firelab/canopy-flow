

#include <iostream>
#include "canopyFlow.h"

int main() {

    canopyFlow wind;

    //Set inputs
    wind.C.leafAreaIndex = 1.0;
    wind.C.canopyHeight = 10.0;                        //canopy height (m)
    wind.C.z0g = 0.025;                      //ground roughness length (m)
    wind.C.dragCoefAth = 0.2;
    wind.C.heightMaxFoliageDist = 0.5;     //height of maximum foliage distribution for the normal distribution (m)
    wind.C.standardDevFoliageDist = 0.25;   //standard deviation of foliage distribution for the normal distribution (m)
    wind.C.numNodes = 10001;                 //number of cells to use for numerical integration

    wind.C.initialize();
    wind.computeWind();

    //wind.C.plot();
    wind.plotWind(15.0, wind.C.canopyHeight + 5.0);

    std::cout << "Done!" << std::endl;
    return 0;
}
