
#ifndef CANOPYFLOW_H_
#define CANOPYFLOW_H_

#include <iostream>
#include <cmath>
#include "boost/math/special_functions.hpp"
#include "eigen3/Eigen/Dense"

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

    void computeWind();

protected:
    double FAI; //Foliage Area Index
    double z0h; //ratio of ground surface roughness length (z_0) to canopy height (h)
    double Uph; //ratio mean background pressure gradient to wind speed at canopy height (U_h) - This is a dimensionless parameter: Wang's model
    double H;   //height above canopy top.
    double h;   //canopy height
    double UH;  //wind speed at H

    /*NOTE:  The canopy is assumed to be comprised of several layers, each
     * with a fraction of the total of FAI - The foliage in each of these
     * layers is assumed constant and each layer is bounded by an upper and
     * a lower node - The physical height of each node is expressed as a
     * fraction (zhfrac) of the total canopy height (h)
     *
     * zhfrac.rows() = folfrac.rows() (always)
     * folfrac(end) = trunk space fraction of FAI
     * zhfrac(1) = 1 (always) and zhfrac(end) = z0h (always)
     * zhfrac(1) and zhfrac(end) are the top and bottom canopy nodes
     */

    Eigen::VectorXd folfrac;    //vertical array of foliage distribution (fraction of FAI)
    Eigen::VectorXd zhfrac;     //vertical array of height of canopy nodes (fraction of h)

    double delh;
    Eigen::VectorXd z01;

};



#endif /* CANOPYFLOW_H_ */
