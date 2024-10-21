#ifndef CANOPYFLOWSWIG_H
#define CANOPYFLOWSWIG_H

#include <vector>
#include "../src/canopyFlow.cpp"

void uniformDistribution(double crownRatio, double leafAreaIndex, double canopyHeight, double dragCoefAth, double z0g, double numNodes, double inputSpeed, double inputHeight, std::vector<double>& windSpeed, std::vector<double>& z);

void asymmetricGaussianDistribution(double heightMaxFoliageDist, double standardDevFoliageUpper, double standardDevFoliageLower, double leafAreaIndex, double canopyHeight, double dragCoefAth, double z0g, double numNodes, double inputSpeed, double inputHeight, std::vector<double>& windSpeed, std::vector<double>& z);

void normalDistribution(double heightMaxFoliageDist, double standardDevFoliageDist, double leafAreaIndex, double canopyHeight, double dragCoefAth, double z0g, double numNodes, double inputSpeed, double inputHeight, std::vector<double>& windSpeed, std::vector<double>& z);

void triangleDistribution(double A1, double Ax, double Ab, double zmax, double zbot, double leafAreaIndex, double canopyHeight, double dragCoefAth, double z0g, double numNodes, double inputSpeed, double inputHeight, std::vector<double>& windSpeed, std::vector<double>& z);

void massmanDistribution(double A1, double A2, double A3, double zmax, double leafAreaIndex, double canopyHeight, double dragCoefAth, double z0g, double numNodes, double inputSpeed, double inputHeight, std::vector<double>& windSpeed, std::vector<double>& z);


#endif