

#include <iostream>
#include <stdio.h>
#include "canopyFlow.h"
#include <cstdlib>
#include <vector>
#include <fstream>

using namespace std;

/**
 * @brief calc_singlePoiunt
 *
 * Use "Compute wind"
 * to calculate a single data point that is passed through via a subprocess pipe
 * or terminal
 * @param initialWindSpeed
 * @param initialHeight
 * @param canopy_height
 * @param z0g
 * @param leafAreaIndex
 * @param dragCoeff
 * @param crownRatio
 * @param outputWindHeight
 * @return
 */
double calc_singlePoint(double initialWindSpeed, double initialHeight,
                         double canopy_height, double z0g, double leafAreaIndex,
                         double dragCoeff,double crownRatio,double outputWindHeight)
{
    canopyFlow wind;
    wind.C = new canopy_uniform_distribution(crownRatio);
    wind.C->z0g = z0g;
    wind.C->leafAreaIndex = leafAreaIndex;
    wind.C->canopyHeight = canopy_height;
    wind.C->dragCoefAth = dragCoeff;
    wind.C->numNodes = 100001;                 //number of cells to use for numerical integration
    wind.computeWind();
    double outSpeed = wind.get_windspeed(initialWindSpeed,initialHeight,outputWindHeight);
    return outSpeed;
}
/**
 * @brief calc_pointArray
 * Calculate an array of points passed in as a *.csv
 *output is written to csv as well
 * @param initialWindSpeed
 * @param initialHeight
 * @param canopy_height
 * @param z0g
 * @param leafAreaIndex
 * @param dragCoeff
 * @param crownRatio
 * @param inFile
 * @param outputFile
 */
void calc_pointArray(double initialWindSpeed, double initialHeight,
                    double canopy_height, double z0g, double leafAreaIndex,
                    double dragCoeff,double crownRatio, std::string inFile,std::string outputFile)
{
  vector<double> vResult;
  vector<double> inputHgts;

  // std::string cfg_path=FindDataPath(arrayPath); //Find the csv file with all the config options
  std::ifstream dFile(inFile.c_str());
  std::string line;

  if(dFile.is_open())
  {
      while(getline(dFile,line))
      {
//          cout<<line<<endl;
          if(line=="CanopyFlow output:")
          {
              cout<<"No New Input Detected!"<<endl;
              return;
          }
          double dLine = atof(line.c_str());
          inputHgts.push_back(dLine);
      }
      dFile.close();
  }
  else
  {
      cout<<"cannot open file.."<<endl;
  }
  canopyFlow wind;
  wind.C = new canopy_uniform_distribution(crownRatio);
  wind.C->z0g = z0g;
  wind.C->leafAreaIndex = leafAreaIndex;
  wind.C->canopyHeight = canopy_height;
  wind.C->dragCoefAth = dragCoeff;
  wind.C->numNodes = 100001;                 //number of cells to use for numerical integration
  wind.computeWind();
  std::ofstream outFile;
  outFile.open(outputFile.c_str());
  outFile<<"CanopyFlow output:"<<endl;
  for(unsigned long i=0;i<inputHgts.size();i++)
  {
    double outSpeed = wind.get_windspeed(initialWindSpeed,initialHeight,inputHgts[i]);
    outFile<<outSpeed<<endl;
  }
  outFile.close();
  cout<<"-:Done:-"<<endl;
}
/**
 * @brief main
 * Parse command line input and
 * direct to function
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char* argv[]) {

    std::string eval = std::string(argv[1]);
    if(eval=="vec")
    {
        double initialWindSpeed = atof(argv[2]);
        double initialHeight = atof(argv[3]);
        // double outWindHeight = atof(argv[4]);
        double canopy_height = atof(argv[4]);
        double z0_ground = atof(argv[5]);
        double PAI = atof(argv[6]);
        double drag_coeff = atof(argv[7]);
        double crownRatio = atof(argv[8]);
        std::string inputFile = std::string(argv[9]);
        std::string outputFile = std::string(argv[10]);
        calc_pointArray(initialWindSpeed,
                          initialHeight,
                          canopy_height,
                          z0_ground,
                          PAI,
                          drag_coeff,
                          crownRatio,
                          inputFile,
                          outputFile);
    }
    if(eval=="single")
    {
        double initialWindSpeed = atof(argv[2]);
        double initialHeight = atof(argv[3]);
        double canopy_height = atof(argv[4]);
        double z0_ground = atof(argv[5]);
        double PAI = atof(argv[6]);
        double drag_coeff = atof(argv[7]);
        double crownRatio = atof(argv[8]);
        double outWindHeight = atof(argv[9]);
        double output = calc_singlePoint(initialWindSpeed,
                                          initialHeight,
                                          canopy_height,
                                          z0_ground,
                                          PAI,
                                          drag_coeff,
                                          crownRatio,
                                          outWindHeight);
        cout<<"WindOutL-:"<<output<<":-"<<endl;
    }

    bool oldWay = false;

    if(oldWay==true){
        canopyFlow wind;
        double initialWindSpeed = atof(argv[1]);
        double initialHeight = atof(argv[2]);
        double outWindHeight = atof(argv[3]);
        double canopy_height = atof(argv[4]);
        double z0_ground = atof(argv[5]);
        double PAI = atof(argv[6]);
        double drag_coeff = atof(argv[7]);
        double cRat = atof(argv[8]);
        wind.C->z0g = z0_ground;
        wind.C->leafAreaIndex = PAI;
        wind.C->canopyHeight = canopy_height;                        //canopy height (m)
        wind.C->dragCoefAth = drag_coeff;
        double crownRatio = cRat;
        wind.C = new canopy_uniform_distribution(crownRatio);
        wind.C->numNodes = 100001;                 //number of cells to use for numerical integration
        wind.computeWind();
        double outSpeed = 0;
        outSpeed = wind.get_windspeed(initialWindSpeed,initialHeight,outWindHeight);
        cout<<"init Height: "<<initialHeight<<endl;
        cout<<"Canopy: "<<canopy_height<<endl;
        cout<<"Out Wind Height: "<<outWindHeight<<endl;
        cout<<"windIn: "<<initialWindSpeed<<endl;
        cout<<"WindOut-:"<<outSpeed<<":-"<<endl;
    }
    return 0;
}
