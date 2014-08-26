
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

    canopy C;

    void computeWind();

protected:
    double H;   //height above canopy top.
    double UH;  //wind speed at H
};

#endif /* CANOPYFLOW_H_ */
