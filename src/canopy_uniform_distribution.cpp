
#include "canopy_uniform_distribution.h"

canopy_uniform_distribution::canopy_uniform_distribution() : canopy()
{
    distributionType = uniform;
    crownRatio = 1.0;
}

canopy_uniform_distribution::canopy_uniform_distribution(double crownRatio_) : canopy()
{
    distributionType = uniform;
    crownRatio = crownRatio_;
}

canopy_uniform_distribution::canopy_uniform_distribution(canopy_uniform_distribution &rhs) : canopy(rhs)
{
    distributionType = rhs.distributionType;
    crownRatio = rhs.crownRatio;
}

canopy_uniform_distribution &canopy_uniform_distribution::operator=(const canopy_uniform_distribution &rhs)
{
    if(&rhs != this)
    {
        canopy::operator=(rhs);
        crownRatio = rhs.crownRatio;
    }
    return *this;
}

canopy_uniform_distribution::~canopy_uniform_distribution()
{

}

void canopy_uniform_distribution::compute_haz()
{
    double integHazn = 0.0;

    //---------------FIX THIS------------------!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //double intermediate = cellsize*2.0 / 3.0;
    double intermediate = 2.0 / 3.0;
    //---------------FIX THIS------------------!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    for(int i=0; i<numNodes; i++)   //integrate using extended simpson's rule
    {
        if(i*cellsize < (1.0-crownRatio))
            haz[i] = 0.0; //temporarily store this here
        else
            haz[i] = 1.0; //temporarily store this here

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


