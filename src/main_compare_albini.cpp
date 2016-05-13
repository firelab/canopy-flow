

#include <iostream>
#include <stdio.h>
#include "canopyFlow.h"
#include "windAdjustmentFactor.h"

static const double pi = acos(-1.0);

int main() {

    //----Settings for flames below canopy top case----
    double canopyHeight = 10.0;                 //canopy height (m)
    //double canopyCover = 0.7;                   //canopy cover [0..1]
    double crownRatio = 0.5;                    //crown ratio [0..1]
    double dragCoef = 0.2;                      //drag coefficient (Albini used 1.0, Massman likes 0.2)
    double fuelBedDepth = 0.5 / 3.28084;                  //fuel bed depth (m)
    double inputHeight = canopyHeight + 6.096;  //input wind height (m)
    double betaSigma = 10.6955;                 //10.6955 is the value albini uses for beta*sigma, but in m^-1
    double midFlameHeight = fuelBedDepth;   //This assumes the flame height is equal to twice the fuel bed depth height (measured from the ground)
    //double groundRoughness = 0.0075;            //ground roughness (m)
    double groundRoughness = fuelBedDepth * 0.13;            //ground roughness (m), using Albini method

    int plotNodes = 1000;
    double* canopyCoverArray = new double[plotNodes];
    double* WAFarrayAlbiniBelow = new double[plotNodes];
    double* WAFarrayMassmanBelow = new double[plotNodes];
    double* WAFarrayAlbiniAbove = new double[plotNodes];
    double* WAFarrayMassmanAbove = new double[plotNodes];

    WindAjustmentFactor albini;
    canopyFlow massman;

    //-------Uniform Distribution, Below Canopy Flames-----------------------
    massman.C = new canopy_uniform_distribution(crownRatio);
    massman.C->canopyHeight = canopyHeight;
    massman.C->dragCoefAth = dragCoef;
    massman.C->z0g = groundRoughness;
    massman.C->numNodes = 10001;
    massman.computeWind();

    double canopyCoverStepSize = 1.0 / (plotNodes-1);

    //build plot arrays
    for(int i=0; i<plotNodes; i++)
    {
        canopyCoverArray[i] = i * canopyCoverStepSize;
        massman.C->leafAreaIndex = canopyHeight * crownRatio * canopyCoverArray[i] * betaSigma / (6.0 * pi);
        massman.computeWind();
        WAFarrayMassmanBelow[i] = massman.get_windAdjustmentFactorShelteredMidFlame(inputHeight, midFlameHeight);
        WAFarrayAlbiniBelow[i] = albini.calculateWindAdjustmentFactor(canopyCoverArray[i], canopyHeight*3.28084, crownRatio, fuelBedDepth*3.28084);  //arguments converted to feet here
    }

    //-------Uniform Distribution, Above Canopy Flames-----------------------
    //----Settings for flames above canopy top case----
    canopyHeight = 3.0;                 //canopy height (m)
    inputHeight = canopyHeight + 6.096;  //input wind height (m)

    crownRatio = 1.0;   //use 1.0 for crown ratio?
    massman.C = new canopy_uniform_distribution(crownRatio);
    massman.C->canopyHeight = canopyHeight;
    massman.C->dragCoefAth = dragCoef;
    massman.C->z0g = 0.0025;   //not sure what to use here...  fairly smooth roughness under canopy...?
    massman.C->numNodes = 10001;
    massman.computeWind();

    //build more plot arrays
    for(int i=0; i<plotNodes; i++)
    {
        massman.C->leafAreaIndex = canopyHeight * crownRatio * canopyCoverArray[i] * betaSigma / (6.0 * pi);  //10.6955 is the value albini uses for beta*sigma, but in m^-1
        massman.computeWind();
        WAFarrayMassmanAbove[i] = massman.get_windAdjustmentFactorUnshelteredIntegral(inputHeight, massman.C->canopyHeight*2.0); //Assume flame height is double the canopy height
        WAFarrayAlbiniAbove[i] = albini.calculateWindAdjustmentFactor(0.0, 0.0, crownRatio, massman.C->canopyHeight*3.28084);  //arguments converted to feet here
    }


    PLFLT xmin = canopyCoverArray[0], ymin = 0.0, xmax = canopyCoverArray[plotNodes-1], ymax = 1.0;

    PLINT just=0, axis=0;
    plstream *pls;

    // plplot initialization

    pls = new plstream();  // declare plplot object

    plsdev("svg"); //cairo uses the same color scheme as on screen - black on
    //red on black default
    plsfnam("massman_vs_albini_WAF.svg");// sets the names of the output file

    plscolbg(255,255,255);  //change background color

    pls->init();           // start plplot object
    plscol0(1, 0, 0, 0);        //change the color pallet0 color #1 to be black (0,0,0)
    plscol0(2, 255, 0, 0);      //change the color pallet0 color #2 to be red (255,0,0)
    plscol0(3, 255, 255, 255);  //change the color pallet0 color #3 to be white (255,255,255)
    plcol0(1);              //now change our font color to be the color #1

    pls->adv( 0 );
    pls->vpor( 0.15, 0.85, 0.1, 0.8 );
    pls->wind( xmin, xmax, ymin, ymax);
    pls->box( "bcnst", 0.0, 0, "bcnst", 0.0, 0 );

    pls->line(plotNodes, (PLFLT*) canopyCoverArray,(PLFLT*) WAFarrayAlbiniBelow);   //plot Albini
    plcol0(2);              //now change our font color to be the color #2
    pls->line(plotNodes, (PLFLT*) canopyCoverArray,(PLFLT*) WAFarrayMassmanBelow);   //plot Massman
    plcol0(1);              //now change our font color to be the color #1

    pls->lsty(2);
    pls->line(plotNodes, (PLFLT*) canopyCoverArray,(PLFLT*) WAFarrayAlbiniAbove);   //plot Albini
    plcol0(2);              //now change our font color to be the color #2
    pls->line(plotNodes, (PLFLT*) canopyCoverArray,(PLFLT*) WAFarrayMassmanAbove);   //plot Massman
    plcol0(1);              //now change our font color to be the color #1
    pls->lsty(1);

    pls->schr(0, 1.6);  //change font size
    pls->mtex( "t", 4.0, 0.5, 0.5, "Albini vs. Massman" );
    pls->schr(0, 1.0);  //change font size
    pls->mtex( "b", 3.0, 0.5, 0.5, "Canopy Cover" );
    pls->mtex( "l", 3.0, 0.5, 0.5, "WAF" );

    // Draw a legend
    PLINT        nlegend = 4;
    const char   *text[4], *symbols[4];
    PLINT        opt_array[4];
    PLINT        text_colors[4];
    PLINT        line_colors[4];
    PLINT        line_styles[4];
    PLFLT        line_widths[4];
    PLINT        symbol_numbers[4], symbol_colors[4];
    PLFLT        symbol_scales[4];
    PLFLT        legend_width, legend_height;

    // First legend entry.
    opt_array[0]   = PL_LEGEND_LINE;
    text_colors[0] = 1;
    text[0]        = "Albini, flames under canopy";
    line_colors[0] = 1;
    line_styles[0] = 1;
    line_widths[0] = 1.;
    symbols[0] = "";

    // Second legend entry.
    opt_array[1]      = PL_LEGEND_LINE;
    text_colors[1]    = 2;
    text[1]           = "Massman,flames under canopy";
    line_colors[1]    = 2;
    line_styles[1]    = 1;
    line_widths[1]    = 1.;
    //symbol_colors[1]  = 3;
    //symbol_scales[1]  = 1.;
    //symbol_numbers[1] = 4;
    symbols[1]        = "";
    // from the above opt_arrays we can completely ignore everything
    // to do with boxes.

    // First legend entry.
    opt_array[2]   = PL_LEGEND_LINE;
    text_colors[2] = 1;
    text[2]        = "Albini, flames above canopy";
    line_colors[2] = 1;
    line_styles[2] = 2;
    line_widths[2] = 1.;
    symbols[2] = "";

    // Second legend entry.
    opt_array[3]      = PL_LEGEND_LINE;
    text_colors[3]    = 2;
    text[3]           = "Massman, flames above canopy";
    line_colors[3]    = 2;
    line_styles[3]    = 2;
    line_widths[3]    = 1.;
    //symbol_colors[3]  = 3;
    //symbol_scales[3]  = 1.;
    //symbol_numbers[3] = 4;
    symbols[3]        = "";
    // from the above opt_arrays we can completely ignore everything
    // to do with boxes.

    plscol0a( 15, 32, 32, 32, 0.70 );

    pls->legend(&legend_width, &legend_height,
                PL_LEGEND_BACKGROUND | PL_LEGEND_BOUNDING_BOX, 0,
                0.03, 0.03, 0.1, 3,
                1, 1, 1, 1,
                nlegend, opt_array,
                -3.0, 1.0, 2.0,
                1., text_colors, (const char **) text,
                NULL, NULL, NULL, NULL,
                line_colors, line_styles, line_widths,
                symbol_colors, symbol_scales, symbol_numbers, (const char **) symbols );

    delete pls; // close plot

    delete canopyCoverArray;
    canopyCoverArray = NULL;
    delete WAFarrayAlbiniBelow;
    WAFarrayAlbiniBelow = NULL;
    delete WAFarrayMassmanBelow;
    WAFarrayMassmanBelow = NULL;
    delete WAFarrayAlbiniAbove;
    WAFarrayAlbiniAbove = NULL;
    delete WAFarrayMassmanAbove;
    WAFarrayMassmanAbove = NULL;


    std::cout << "Done!" << std::endl;
    return 0;

}
