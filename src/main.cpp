

#include <iostream>
#include "canopyFlow.h"

int main() {

    canopyFlow wind;

    //-------Normal Distribution------------------------
    //double heightMaxFoliageDist = 0.5;
    //double standardDevFoliageDist = 0.25;
    //wind.C = new canopy_normal_distribution(heightMaxFoliageDist, standardDevFoliageDist);

    //-------Triangle Distribution---------------------
    double A1 = 0.0;       //density of top
    double Ax = 3.0;          //density at max point
    double Ab = 0.1;       //density of bottom (trunk space)
    double zmax = 0.7;      //height to Ax (0 < zmax < 1)
    double zbot = 0.3;      //height to bottom of triangular part (0 < zbot < 1; zbot < zmax)
    wind.C = new canopy_triangle_distribution(A1, Ax, Ab, zmax, zbot);

    //Set inputs

    wind.C->leafAreaIndex = 1.0;
    wind.C->canopyHeight = 10.0;                        //canopy height (m)
    wind.C->z0g = 0.025;                      //ground roughness length (m)
    wind.C->dragCoefAth = 0.2;
    wind.C->numNodes = 10001;                 //number of cells to use for numerical integration

    wind.C->initialize();
    wind.computeWind();

    //wind.C->plot();
    wind.plotWind(15.0, wind.C->canopyHeight + 5.0);

    std::cout << "Done!" << std::endl;
    return 0;
}
