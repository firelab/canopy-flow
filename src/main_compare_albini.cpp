

#include <iostream>
#include <stdio.h>
#include "canopyFlow.h"
#include "windAdjustmentFactor.h"

void plotWAFvsCanopyCover(double canopyHeight, double crownRatio, double dragCoef, double groundRoughness, double fuelBedDepth, double inputHeight, double midFlameHeight);

static const double pi = acos(-1.0);

int main() {

    double canopyHeight = 20.0;                 //canopy height (m)
    //double canopyCover = 0.7;                   //canopy cover [0..1]
    double crownRatio = 0.5;                    //crown ratio [0..1]
    double dragCoef = 0.2;                      //drag coefficient (Albini used 1.0, Massman likes 0.2)
    double fuelBedDepth = 0.5 / 3.28084;                  //fuel bed depth (m)
    double inputHeight = canopyHeight + 6.096;  //input wind height (m)
    double midFlameHeight = fuelBedDepth*1.5;   //This assumes the flame height is equal to the fuel bed depth height (above the top of the fuel bed)
    //double groundRoughness = 0.0075;            //ground roughness (m)
    double groundRoughness = fuelBedDepth * 0.13;            //ground roughness (m), using Albini method

    plotWAFvsCanopyCover(canopyHeight, crownRatio, dragCoef, groundRoughness, fuelBedDepth, inputHeight, midFlameHeight);

    std::cout << "Done!" << std::endl;
    return 0;
}

void plotWAFvsCanopyCover(double canopyHeight, double crownRatio, double dragCoef, double groundRoughness, double fuelBedDepth, double inputHeight, double midFlameHeight)
{

    WindAjustmentFactor albini;
    canopyFlow wind;
    //-------Uniform Distribution-----------------------
    wind.C = new canopy_uniform_distribution(crownRatio);
    wind.C->canopyHeight = canopyHeight;
    wind.C->dragCoefAth = dragCoef;
    wind.C->z0g = groundRoughness;
    wind.C->numNodes = 10001;
    wind.C->initialize();
    wind.computeWind();

    int plotNodes = 1000;
    double* canopyCoverArray = new double[plotNodes];
    double* WAFarrayAlbini = new double[plotNodes];
    double* WAFarrayMassman = new double[plotNodes];

    double canopyCoverStepSize = 1.0 / (plotNodes-1);

    for(int i=0; i<plotNodes; i++)
    {
        canopyCoverArray[i] = i * canopyCoverStepSize;
        wind.C->leafAreaIndex = canopyHeight * canopyCoverArray[i] * 10.6955 / (6.0 * pi);  //10.6955 is the value albini uses for beta*sigma, but in m^-1
        wind.C->initialize();
        wind.computeWind();
        WAFarrayMassman[i] = wind.get_windAdjustmentFactorShelteredMidFlame(inputHeight, midFlameHeight);
        WAFarrayAlbini[i] = albini.calculateWindAdjustmentFactor(canopyCoverArray[i], canopyHeight*3.28084, crownRatio, fuelBedDepth*3.28084);  //arguments converted to feet here
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

    pls->line(plotNodes, (PLFLT*) canopyCoverArray,(PLFLT*) WAFarrayAlbini);   //plot Albini
    plcol0(2);              //now change our font color to be the color #2
    pls->line(plotNodes, (PLFLT*) canopyCoverArray,(PLFLT*) WAFarrayMassman);   //plot Massman
    plcol0(1);              //now change our font color to be the color #1

    pls->schr(0, 1.6);  //change font size
    pls->mtex( "t", 4.0, 0.5, 0.5, "Albini vs. Massman" );
    pls->schr(0, 1.0);  //change font size
    pls->mtex( "b", 3.0, 0.5, 0.5, "Canopy Cover" );
    pls->mtex( "l", 3.0, 0.5, 0.5, "WAF" );

    // Draw a legend
    PLINT        nlegend = 2;
    const char   *text[2], *symbols[2];
    PLINT        opt_array[2];
    PLINT        text_colors[2];
    PLINT        line_colors[2];
    PLINT        line_styles[2];
    PLFLT        line_widths[2];
    PLINT        symbol_numbers[2], symbol_colors[2];
    PLFLT        symbol_scales[2];
    PLFLT        legend_width, legend_height;

    // First legend entry.
    opt_array[0]   = PL_LEGEND_LINE;
    text_colors[0] = 1;
    text[0]        = "Albini";
    line_colors[0] = 1;
    line_styles[0] = 1;
    line_widths[0] = 1.;
    symbols[0] = "";

    // Second legend entry.
    opt_array[1]      = PL_LEGEND_LINE;
    text_colors[1]    = 2;
    text[1]           = "Massman";
    line_colors[1]    = 2;
    line_styles[1]    = 1;
    line_widths[1]    = 1.;
    //symbol_colors[1]  = 3;
    //symbol_scales[1]  = 1.;
    //symbol_numbers[1] = 4;
    symbols[1]        = "";
    // from the above opt_arrays we can completely ignore everything
    // to do with boxes.

    plscol0a( 15, 32, 32, 32, 0.70 );

    pls->legend(&legend_width, &legend_height,
                PL_LEGEND_BACKGROUND | PL_LEGEND_BOUNDING_BOX, 0,
                0.03, 0.03, 0.1, 3,
                1, 1, 1, 1,
                nlegend, opt_array,
                1.0, 1.0, 2.0,
                1., text_colors, (const char **) text,
                NULL, NULL, NULL, NULL,
                line_colors, line_styles, line_widths,
                symbol_colors, symbol_scales, symbol_numbers, (const char **) symbols );

    delete pls; // close plot

    delete canopyCoverArray;
    canopyCoverArray = NULL;
    delete WAFarrayAlbini;
    WAFarrayAlbini = NULL;
    delete WAFarrayMassman;
    WAFarrayMassman = NULL;
}
