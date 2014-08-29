
#ifndef CANOPYFLOW_H_
#define CANOPYFLOW_H_

#include <iostream>
#include <cmath>
//#include "boost/math/special_functions.hpp"
#include <plstream.h>
#include "canopy.h"

/*
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Class that computes wind flow in a canopy.
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*/

class canopyFlow
{

public:
    canopyFlow();
    canopyFlow(canopyFlow &rhs);
    canopyFlow &operator=(const canopyFlow &rhs);
    ~canopyFlow();

    void plot();
    void plotWind(double inputSpeed, double inputHeight);

    canopy C;

    void computeWind(); //Precomputes some necessary wind stuff, run before calling get_windspeed
    double get_windspeed(double inputSpeed, double inputHeight, double desiredHeight);    //Computes windspeed for given input speed, input height, and desiredHeight (speed and height units must be same as inputs)

    double doh;     //non-dimensional canopy displacement height
    double z0oh;    //non-dimensional canopy roughness length

protected:
    double usuh;
    double* uzc;
    double* uzcs;

    //static const double K = 0.4;    //von Karman constant
    //static const double c1 = 0.38;
    //static const double c3 = 15.0;
    //static const double STRSC = 0.5;    //adjustment factor to get the stress profile closer to observations of canopy Reynolds stresses
    //static const double rough = 1.07;   //correction factor for roughness sublayer

    double K;    //von Karman constant
    double c1;
    double c3;
    double STRSC;    //adjustment factor to get the stress profile closer to observations of canopy Reynolds stresses
    double rough;   //correction factor for roughness sublayer
    double logRough;    //precompute log of roughness to speed computations
};

#endif /* CANOPYFLOW_H_ */
