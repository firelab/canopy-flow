
#include "canopy_asymmetric_gaussian_distribution.h"

canopy_asymmetric_gaussian_distribution::canopy_asymmetric_gaussian_distribution() : canopy()
{
    distributionType = double_gaussian;
    heightMaxFoliage = 0.6;
    standardDevFoliageUpper = 0.31;
    standardDevFoliageLower = 0.31;
}

canopy_asymmetric_gaussian_distribution::canopy_asymmetric_gaussian_distribution(double heightMaxFoliageDist_, double standardDevFoliageUpper_, double standardDevFoliageLower_) : canopy()
{
    distributionType = double_gaussian;
    heightMaxFoliage = heightMaxFoliageDist_;
    standardDevFoliageUpper = standardDevFoliageUpper_;
    standardDevFoliageLower = standardDevFoliageLower_;
}

canopy_asymmetric_gaussian_distribution::canopy_asymmetric_gaussian_distribution(canopy_asymmetric_gaussian_distribution &rhs) : canopy(rhs)
{
    distributionType = rhs.distributionType;
    heightMaxFoliage = rhs.heightMaxFoliage;
    standardDevFoliageUpper = rhs.standardDevFoliageUpper;
    standardDevFoliageLower = rhs.standardDevFoliageLower;
}

canopy_asymmetric_gaussian_distribution &canopy_asymmetric_gaussian_distribution::operator=(const canopy_asymmetric_gaussian_distribution &rhs)
{
    if(&rhs != this)
    {
        canopy::operator=(rhs);
        heightMaxFoliage = rhs.heightMaxFoliage;
        standardDevFoliageUpper = rhs.standardDevFoliageUpper;
        standardDevFoliageLower = rhs.standardDevFoliageLower;
    }
    return *this;
}

canopy_asymmetric_gaussian_distribution::~canopy_asymmetric_gaussian_distribution()
{

}

void canopy_asymmetric_gaussian_distribution::compute_haz()
{
    double norm;
    double integHazn = 0.0;

    //---------------FIX THIS------------------!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //double intermediate = cellsize*2.0 / 3.0;
    double intermediate = 2.0 / 3.0;
    //---------------FIX THIS------------------!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    for(int i=0; i<numNodes; i++)   //integrate using extended simpson's rule
    {
        if(i*cellsize < heightMaxFoliage)
            norm = (heightMaxFoliage - i*cellsize) / standardDevFoliageLower;
        else
            norm = (i*cellsize - heightMaxFoliage) / standardDevFoliageUpper;

        haz[i] = exp(-norm * norm); //temporarily store this here

        if(i%2 == 0)    //if even numbers
            integHazn += haz[i];
        else            //if odd numbers
            integHazn += 2.0 * haz[i];
    }

    integHazn = integHazn - 0.5 * (haz[0] + haz[numNodes-1]);
    integHazn *= intermediate;

    for(int i=0; i<numNodes; i++)
        haz[i] = leafAreaIndex * haz[i] /integHazn;
}


