#ifndef CANOPY_H
#define CANOPY_H

#include <math.h>
#include <cstring>
#include <plstream.h>

/*
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Class that stores canopy density, etc. information for a normal distribution.
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*/

class canopy
{

public:
    canopy();
    canopy(canopy &rhs);
    canopy &operator=(const canopy &rhs);
    ~canopy();

    void initialize();
    void plot_haz();

    //inputs
    double leafAreaIndex;
    double h;                       //canopy height (m)
    double z0g;                     //ground roughness length (m)
    double dragCoefAth;
    double heightMaxFoliageDist;    //height of maximum foliage distribution for the normal distribution (m)
    double standardDevFoliageDist;  //standard deviation of foliage distribution for the normal distribution (m)
    int numNodes;                   //number of cells to use for numerical integration

    double* cumulativeLeafDragArea; //cumulative leaf drag area (m^2/m^2)
    double* haz;                    //nondimensional leaf area density
    double* hacpz;                  //nondimensional drag area density
    double* zetaz;                  //normalized mapped vertical coordinate
    double zetah;                //this is hacpn at the top node

protected:
    double get_dragCoef(double zOverh);
    double get_shelterFactor(double hazn);
    void compute_haz();
    void compute_foliage_drag_area_index();

    double cellsize;                //cellsize for integration, this is computed, NOT INPUT
    double totalDragAreaIndex;
};


#endif // CANOPY_H
