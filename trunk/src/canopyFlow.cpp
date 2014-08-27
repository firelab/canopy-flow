#include "canopyFlow.h"

canopyFlow::canopyFlow()
{
    H = 20.0;
    UH = 5.0;
}

canopyFlow::canopyFlow(canopyFlow &rhs)
{
    H = rhs.H;
    UH = rhs.UH;
}

canopyFlow &canopyFlow::operator=(const canopyFlow &rhs)
{
    if(&rhs != this)
    {
        H = rhs.H;
        UH = rhs.UH;
    }
    return *this;
}

canopyFlow::~canopyFlow()
{

}

void canopyFlow::computeWind()
{
    //C.initialize();
}

