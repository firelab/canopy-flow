/*
canopyFlowSwig.h
Description: function definitions that we want exposed in python. Subject to change on how this works
Main file at the botton used for testing functions before using SWIG
Command for running the main
g++ -g canopyFlowSwig.cpp ../src/canopyFlow.cpp ../src/canopy_uniform_distribution.cpp ../src/canopy_asymmetric_gaussian_distribution.cpp ../src/canopy_normal_distribution.cpp ../src/canopy_triangle_distribution.cpp ../src/massman_distribution.cpp ../src/canopy.cpp ../src/measured_distribution.cpp -o canopyFlowSwigExecutable
*/


#include "canopyFlowSwig.h"

void uniformDistribution(double crownRatio, double leafAreaIndex, double canopyHeight, double dragCoefAth, double z0g, double numNodes, double inputSpeed, double inputHeight, std::vector<double>& windSpeed, std::vector<double>& z)
{
    canopyFlow wind;

    // Setup for Uniform Distribution
    wind.C = new canopy_uniform_distribution(crownRatio);
    wind.C->leafAreaIndex = leafAreaIndex;
    wind.C->canopyHeight = canopyHeight;    // canopy height (m)
    wind.C->dragCoefAth = dragCoefAth;

    wind.C->z0g = z0g;                      // ground roughness length (m), model is somewhat sensitive to this
    wind.C->numNodes = numNodes;            // number of cells to use for numerical integration

    // Call the wind profile data function
    wind.computeWind();

    wind.calculateDimensionalWind(inputSpeed, inputHeight, z, windSpeed);

}

void asymmetricGaussianDistribution(double heightMaxFoliageDist, double standardDevFoliageUpper, double standardDevFoliageLower, double leafAreaIndex, double canopyHeight, double dragCoefAth, double z0g, double numNodes, double inputSpeed, double inputHeight, std::vector<double>& windSpeed, std::vector<double>& z)
{
    canopyFlow wind;

    // Setup for Asymmetric Gaussian Distribution
    wind.C = new canopy_asymmetric_gaussian_distribution(heightMaxFoliageDist, standardDevFoliageUpper, standardDevFoliageLower);
    wind.C->leafAreaIndex = leafAreaIndex;
    wind.C->canopyHeight = canopyHeight;    // canopy height (m)
    wind.C->dragCoefAth = dragCoefAth;

    wind.C->z0g = z0g;                      // ground roughness length (m), model is somewhat sensitive to this
    wind.C->numNodes = numNodes;            // number of cells to use for numerical integration

    // Call the wind profile data function
    wind.computeWind();
    
    wind.calculateDimensionalWind(inputSpeed, inputHeight, z, windSpeed);\
}

void normalDistribution(double heightMaxFoliageDist, double standardDevFoliageDist, double leafAreaIndex, double canopyHeight, double dragCoefAth, double z0g, double numNodes, double inputSpeed, double inputHeight, std::vector<double>& windSpeed, std::vector<double>& z)
{
    canopyFlow wind;

    // Setup for Normal Distribution
    wind.C = new canopy_normal_distribution(heightMaxFoliageDist, standardDevFoliageDist);
    wind.C->leafAreaIndex = leafAreaIndex;
    wind.C->canopyHeight = canopyHeight;    // canopy height (m)
    wind.C->dragCoefAth = dragCoefAth;

    wind.C->z0g = z0g;                      // ground roughness length (m), model is somewhat sensitive to this
    wind.C->numNodes = numNodes;            // number of cells to use for numerical integration
    
    // Call the wind profile data function
    wind.computeWind();
    
    wind.calculateDimensionalWind(inputSpeed, inputHeight, z, windSpeed);

}

void triangleDistribution(double A1, double Ax, double Ab, double zmax, double zbot, double leafAreaIndex, double canopyHeight, double dragCoefAth, double z0g, double numNodes, double inputSpeed, double inputHeight, std::vector<double>& windSpeed, std::vector<double>& z)
{
    canopyFlow wind;

    // Setup for Triangle Distribution
    wind.C = new canopy_triangle_distribution(A1, Ax, Ab, zmax, zbot);
    wind.C->leafAreaIndex = leafAreaIndex;
    wind.C->canopyHeight = canopyHeight;    // canopy height (m)
    wind.C->dragCoefAth = dragCoefAth;

    wind.C->z0g = z0g;                      // ground roughness length (m), model is somewhat sensitive to this
    wind.C->numNodes = numNodes;            // number of cells to use for numerical integration

    // Call the wind profile data function
    wind.computeWind();

    wind.calculateDimensionalWind(inputSpeed, inputHeight, z, windSpeed);

}

void massmanDistribution(double A1, double A2, double A3, double zmax, double leafAreaIndex, double canopyHeight, double dragCoefAth, double z0g, double numNodes, double inputSpeed, double inputHeight, std::vector<double>& windSpeed, std::vector<double>& z)
{
    canopyFlow wind;

    // Setup for Massman Distribution
    wind.C = new massman_distribution(A1, A2, A3, zmax);
    wind.C->leafAreaIndex = leafAreaIndex;
    wind.C->canopyHeight = canopyHeight;    // canopy height (m)
    wind.C->dragCoefAth = dragCoefAth;

    wind.C->z0g = z0g;                      // ground roughness length (m), model is somewhat sensitive to this
    wind.C->numNodes = numNodes;            // number of cells to use for numerical integration

    // Call the wind profile data function
    wind.computeWind();

    wind.calculateDimensionalWind(inputSpeed, inputHeight, z, windSpeed);

}
/*
int main() {
    std::vector<double> windSpeed;
    std::vector<double> z;

    // Test Uniform Distribution
    std::cout << "Testing Uniform Distribution..." << std::endl;
    uniformDistribution(0.7, 1.0, 8.0, 0.2, 0.0075, 10001, 10.0, 8.0 + 6.096, windSpeed, z);
    for (size_t i = 0; i < windSpeed.size(); ++i) {
        std::cout << "z[" << i << "] = " << z[i] << ", windSpeed[" << i << "] = " << windSpeed[i] << std::endl;
    }

    // Clear vectors before next test
    windSpeed.clear();
    z.clear();

    // Test Asymmetric Gaussian Distribution
    std::cout << "Testing Asymmetric Gaussian Distribution..." << std::endl;
    asymmetricGaussianDistribution(0.36, 0.6, 0.2, 1.0, 8.0, 0.2, 0.0075, 10001, 10.0, 8.0 + 6.096, windSpeed, z);
    for (size_t i = 0; i < windSpeed.size(); ++i) {
        std::cout << "z[" << i << "] = " << z[i] << ", windSpeed[" << i << "] = " << windSpeed[i] << std::endl;
    }

    // Clear vectors before next test
    windSpeed.clear();
    z.clear();

    // Test Normal Distribution
    std::cout << "Testing Normal Distribution..." << std::endl;
    normalDistribution(0.5, 0.3, 1.0, 3.0, 0.2, 0.0075, 10001, 10.0, 3.0 + 6.096, windSpeed, z);
    for (size_t i = 0; i < windSpeed.size(); ++i) {
        std::cout << "z[" << i << "] = " << z[i] << ", windSpeed[" << i << "] = " << windSpeed[i] << std::endl;
    }

    // Clear vectors before next test
    windSpeed.clear();
    z.clear();

    // Test Triangle Distribution
    std::cout << "Testing Triangle Distribution..." << std::endl;
    triangleDistribution(0.32, 1.0, 0.02, 0.36, 0.12, 3.28, 10.0, 0.2, 0.0075, 10001, 10.0, 10.0 + 6.096, windSpeed, z);
    for (size_t i = 0; i < windSpeed.size(); ++i) {
        std::cout << "z[" << i << "] = " << z[i] << ", windSpeed[" << i << "] = " << windSpeed[i] << std::endl;
    }

    // Clear vectors before next test
    windSpeed.clear();
    z.clear();

    // Test Massman Distribution
    std::cout << "Testing Massman Distribution..." << std::endl;
    massmanDistribution(1.10, 2.0, 1.0, 0.7, 1.0, 10.0, 0.2, 0.0075, 10001, 10.0, 10.0 + 6.096, windSpeed, z);
    for (size_t i = 0; i < windSpeed.size(); ++i) {
        std::cout << "z[" << i << "] = " << z[i] << ", windSpeed[" << i << "] = " << windSpeed[i] << std::endl;
    }

    std::cout << "All tests completed successfully!" << std::endl;
    return 0;
}
*/