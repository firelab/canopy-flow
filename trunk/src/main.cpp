

#include <iostream>
#include "canopyFlow.h"

int main() {

    canopyFlow wind;

    wind.computeWind();

    wind.C.plot_haz();

    std::cout << "Done!" << std::endl;
    return 0;
}
