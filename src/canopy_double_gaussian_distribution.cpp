
#include "canopy_double_gaussian_distribution.h"

canopy_double_gaussian_distribution::canopy_double_gaussian_distribution() : canopy()
{
    distributionType = normal_distribution;
    heightMaxFoliageDist = 0.5;
    standardDevFoliageUpper = 0.25;
    standardDevFoliageLower = 0.25;
}

canopy_double_gaussian_distribution::canopy_double_gaussian_distribution(double heightMaxFoliageDist_, double standardDevFoliageUpper_, double standardDevFoliageLower_)
{
    distributionType = normal_distribution;
    heightMaxFoliageDist = heightMaxFoliageDist_;
    standardDevFoliageUpper = standardDevFoliageUpper_;
    standardDevFoliageLower = standardDevFoliageLower_;
}

canopy_double_gaussian_distribution::canopy_double_gaussian_distribution(canopy_double_gaussian_distribution &rhs) : canopy(rhs)
{
    distributionType = rhs.distributionType;
    heightMaxFoliageDist = rhs.heightMaxFoliageDist;
    standardDevFoliageUpper = rhs.standardDevFoliageUpper;
    standardDevFoliageLower = rhs.standardDevFoliageLower;
}

canopy_double_gaussian_distribution &canopy_double_gaussian_distribution::operator=(const canopy_double_gaussian_distribution &rhs)
{
    if(&rhs != this)
    {
        canopy::operator=(rhs);
        heightMaxFoliageDist = rhs.heightMaxFoliageDist;
        standardDevFoliageUpper = rhs.standardDevFoliageUpper;
        standardDevFoliageLower = rhs.standardDevFoliageLower;
    }
    return *this;
}

canopy_double_gaussian_distribution::~canopy_double_gaussian_distribution()
{

}

void canopy_double_gaussian_distribution::compute_haz()
{
    double norm;
    double hazn;
    double integHazn = 0.0;

    //---------------FIX THIS------------------!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //double intermediate = cellsize*2.0 / 3.0;
    double intermediate = 2.0 / 3.0;
    //---------------FIX THIS------------------!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    for(int i=0; i<numNodes; i++)   //integrate using extended simpson's rule
    {
        if(i*cellsize < heightMaxFoliageDist)
            norm = (heightMaxFoliageDist - i*cellsize) / standardDevFoliageLower;
        else
            norm = (i*cellsize - heightMaxFoliageDist) / standardDevFoliageUpper;

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


