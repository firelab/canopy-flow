

#include <iostream>
#include <stdio.h>
#include "canopyFlow.h"
#include "behaveRun.h"
#include "fuelModelSet.h"
#include "windAdjustmentFactor.h"


static const double pi = acos(-1.0);

int main() {

    FuelModelSet fuelModelSet;
    BehaveRun behave(fuelModelSet);

    //----Settings for fire spread----------------------
    int fuelModelNumber = 2;
    double moistureOneHour = 5.0;
    double moistureTenHour = 6.0;
    double moistureHundredHour = 7.0;
    double moistureLiveHerbaceous = 100.0;
    double moistureLiveWoody = 100.0;
    double windSpeed = 5.0;
    double windDirection = 0;
    double slope = 0.0;
    double aspect = 0.0;
    double directionOfMaxSpread = 0;
    double directionOfInterest = 0;

    //----Settings for flames below canopy top case----
    double canopyHeight = 10.0;                 //canopy height (m)
    double canopyCover = 0.7;                   //canopy cover [0..1]
    double crownRatio = 0.5;                    //crown ratio [0..1]
    double dragCoef = 0.2;                      //drag coefficient (Albini used 1.0, Massman likes 0.2)
    double fuelBedDepth = fuelModelSet.getFuelbedDepth(fuelModelNumber) / 3.28084;                  //fuel bed depth (m)
    double inputHeight = canopyHeight + 6.096;  //input wind height (m)
    double betaSigma = 10.6955;                 //10.6955 is the value albini uses for beta*sigma, but in m^-1
    double midFlameHeight = fuelBedDepth;   //This assumes the flame height is equal to twice the fuel bed depth height (measured from the ground)
    //double groundRoughness = 0.0075;            //ground roughness (m)
    double groundRoughness = fuelBedDepth * 0.13;            //ground roughness (m), using Albini method

    WindAjustmentFactor albini;
    canopyFlow massman;

    int plotNodes = 1000;
    double* canopyCoverArray = new double[plotNodes];
    double* WAFarrayAlbiniBelow = new double[plotNodes];
    double* WAFarrayMassmanBelowUniZ0gHfuel = new double[plotNodes];
    double* WAFarrayMassmanBelowUniZ0g01 = new double[plotNodes];
    double* WAFarrayMassmanBelowGaussZ0gHfuel = new double[plotNodes];
    double* WAFarrayMassmanBelowGaussZ0g01 = new double[plotNodes];
    double* WAFarrayAlbiniAbove = new double[plotNodes];
    double* WAFarrayMassmanAbove = new double[plotNodes];

    double* albiniBelowSpread = new double[plotNodes];
    double* massmanBelowSpreadUniZ0gHfuel = new double[plotNodes];
    double* massmanBelowSpreadUniZ0g01 = new double[plotNodes];
    double* massmanBelowSpreadGaussZ0gHfuel = new double[plotNodes];
    double* massmanBelowSpreadGaussZ0g01 = new double[plotNodes];


    double* albiniAboveSpread = new double[plotNodes];
    double* massmanAboveSpread = new double[plotNodes];


    double* flameHeightBelow = new double[plotNodes];
    double* flameHeightAbove = new double[plotNodes];
    double* WAFBelowIntegral = new double[plotNodes];
    double* WAFBelowMidFlame = new double[plotNodes];
    double* WAFAboveIntegral = new double[plotNodes];


    //-------Below Canopy Flames-----------------------
    massman.C = new canopy_uniform_distribution(crownRatio);
    massman.C->canopyHeight = canopyHeight;
    massman.C->dragCoefAth = dragCoef;
    massman.C->z0g = groundRoughness;
    massman.C->numNodes = 10001;
    massman.computeWind();

    double canopyCoverStepSize = 1.0 / (plotNodes-1);
    double flameHeightStepSizeAbove = (inputHeight-canopyHeight) / (plotNodes);
    double flameHeightStepSizeBelow = canopyHeight / (plotNodes-1);

    //compute Albini below canopy
    for(int i=0; i<plotNodes; i++)
    {
        canopyCoverArray[i] = i * canopyCoverStepSize;
        WAFarrayAlbiniBelow[i] = albini.calculateWindAdjustmentFactor(canopyCoverArray[i], canopyHeight*3.28084, crownRatio, fuelBedDepth*3.28084);  //arguments converted to feet here

        behave.updateSurfaceInputs(fuelModelNumber, moistureOneHour, moistureTenHour, moistureHundredHour, moistureLiveHerbaceous, moistureLiveWoody, WindHeightInputMode::DIRECT_MIDFLAME, WAFarrayAlbiniBelow[i]*windSpeed, windDirection, slope, aspect, canopyCover, canopyHeight, crownRatio);
        albiniBelowSpread[i] = behave.calculateSurfaceFireForwardSpreadRate(directionOfInterest);
    }

    //compute Massman below canopy, uniform canopy profile, z0g=0.13*Hfuel
    for(int i=0; i<plotNodes; i++)
    {
        //build plot arrays for Massman vs Albini WAF plot below canopy
        canopyCoverArray[i] = i * canopyCoverStepSize;
        massman.C->leafAreaIndex = canopyHeight * crownRatio * canopyCoverArray[i] * betaSigma / (6.0 * pi);
        massman.computeWind();
        WAFarrayMassmanBelowUniZ0gHfuel[i] = massman.get_windAdjustmentFactorShelteredMidFlame(inputHeight, midFlameHeight);

        behave.updateSurfaceInputs(fuelModelNumber, moistureOneHour, moistureTenHour, moistureHundredHour, moistureLiveHerbaceous, moistureLiveWoody, WindHeightInputMode::DIRECT_MIDFLAME, WAFarrayMassmanBelowUniZ0gHfuel[i]*windSpeed, windDirection, slope, aspect, canopyCover, canopyHeight, crownRatio);
        massmanBelowSpreadUniZ0gHfuel[i] = behave.calculateSurfaceFireForwardSpreadRate(directionOfInterest);

        //reset canopy cover
        massman.C->leafAreaIndex = canopyHeight * crownRatio * canopyCover * betaSigma / (6.0 * pi);
        massman.computeWind();
        //build plot arrays for Massman flame height plot
        flameHeightAbove[i] = flameHeightStepSizeAbove + canopyHeight + i * flameHeightStepSizeAbove;   //Don't include the canopy height, start above that
        flameHeightBelow[i] = i * flameHeightStepSizeBelow;
        WAFBelowIntegral[i] = massman.get_windAdjustmentFactorShelteredIntegral(inputHeight, flameHeightBelow[i]);
        WAFBelowMidFlame[i] = massman.get_windAdjustmentFactorShelteredMidFlame(inputHeight, flameHeightBelow[i]/2.0);
        WAFAboveIntegral[i] = massman.get_windAdjustmentFactorUnshelteredIntegral(inputHeight, flameHeightAbove[i]);
    }

    //compute Massman below canopy, uniform canopy profile, z0g=0.01m
    for(int i=0; i<plotNodes; i++)
    {
        //build plot arrays for Massman vs Albini WAF plot below canopy
        canopyCoverArray[i] = i * canopyCoverStepSize;
        massman.C->leafAreaIndex = canopyHeight * crownRatio * canopyCoverArray[i] * betaSigma / (6.0 * pi);
        massman.C->z0g = 0.01;
        massman.computeWind();
        WAFarrayMassmanBelowUniZ0g01[i] = massman.get_windAdjustmentFactorShelteredMidFlame(inputHeight, midFlameHeight);

        behave.updateSurfaceInputs(fuelModelNumber, moistureOneHour, moistureTenHour, moistureHundredHour, moistureLiveHerbaceous, moistureLiveWoody, WindHeightInputMode::DIRECT_MIDFLAME, WAFarrayMassmanBelowUniZ0g01[i]*windSpeed, windDirection, slope, aspect, canopyCover, canopyHeight, crownRatio);
        massmanBelowSpreadUniZ0g01[i] = behave.calculateSurfaceFireForwardSpreadRate(directionOfInterest);
    }

    //change to asymmetric Gaussian canopy profile
    delete massman.C;
    massman.C = new canopy_asymmetric_gaussian_distribution(0.6, 0.2, 0.3);
    massman.C->canopyHeight = canopyHeight;
    massman.C->dragCoefAth = dragCoef;
    massman.C->z0g = groundRoughness;
    massman.C->numNodes = 10001;
    massman.computeWind();

    //compute Massman below canopy, asymmetric Gaussian canopy profile, z0g=0.13*Hfuel
    for(int i=0; i<plotNodes; i++)
    {
        //build plot arrays for Massman vs Albini WAF plot below canopy
        canopyCoverArray[i] = i * canopyCoverStepSize;
        massman.C->leafAreaIndex = canopyHeight * crownRatio * canopyCoverArray[i] * betaSigma / (6.0 * pi);
        massman.computeWind();
        WAFarrayMassmanBelowGaussZ0gHfuel[i] = massman.get_windAdjustmentFactorShelteredMidFlame(inputHeight, midFlameHeight);

        behave.updateSurfaceInputs(fuelModelNumber, moistureOneHour, moistureTenHour, moistureHundredHour, moistureLiveHerbaceous, moistureLiveWoody, WindHeightInputMode::DIRECT_MIDFLAME, WAFarrayMassmanBelowGaussZ0gHfuel[i]*windSpeed, windDirection, slope, aspect, canopyCover, canopyHeight, crownRatio);
        massmanBelowSpreadGaussZ0gHfuel[i] = behave.calculateSurfaceFireForwardSpreadRate(directionOfInterest);
    }

    //compute Massman below canopy, asymmetric Gaussian canopy profile, z0g=0.01m
    for(int i=0; i<plotNodes; i++)
    {
        //build plot arrays for Massman vs Albini WAF plot below canopy
        canopyCoverArray[i] = i * canopyCoverStepSize;
        massman.C->leafAreaIndex = canopyHeight * crownRatio * canopyCoverArray[i] * betaSigma / (6.0 * pi);
        massman.C->z0g = 0.01;
        massman.computeWind();
        WAFarrayMassmanBelowGaussZ0g01[i] = massman.get_windAdjustmentFactorShelteredMidFlame(inputHeight, midFlameHeight);

        behave.updateSurfaceInputs(fuelModelNumber, moistureOneHour, moistureTenHour, moistureHundredHour, moistureLiveHerbaceous, moistureLiveWoody, WindHeightInputMode::DIRECT_MIDFLAME, WAFarrayMassmanBelowGaussZ0g01[i]*windSpeed, windDirection, slope, aspect, canopyCover, canopyHeight, crownRatio);
        massmanBelowSpreadGaussZ0g01[i] = behave.calculateSurfaceFireForwardSpreadRate(directionOfInterest);
    }

    //-------Above Canopy Flames-----------------------
    fuelModelNumber = 4;
    canopyHeight = fuelModelSet.getFuelbedDepth(fuelModelNumber) / 3.28084;                 //canopy height (m)
    inputHeight = canopyHeight + 6.096;  //input wind height (m)

    crownRatio = 1.0;   //use 1.0 for crown ratio?
    delete massman.C;
    massman.C = new canopy_uniform_distribution(crownRatio);
    massman.C->canopyHeight = canopyHeight;
    massman.C->dragCoefAth = dragCoef;
    massman.C->z0g = 0.01;   //not sure what to use here...  fairly smooth roughness under canopy...?
    massman.C->numNodes = 10001;
    massman.computeWind();


    for(int i=0; i<plotNodes; i++)
    {
        //build plot arrays for Massman vs Albini WAF plot above canopy
        massman.C->leafAreaIndex = canopyHeight * crownRatio * canopyCoverArray[i] * betaSigma / (6.0 * pi);  //10.6955 is the value albini uses for beta*sigma, but in m^-1
        massman.computeWind();
        WAFarrayMassmanAbove[i] = massman.get_windAdjustmentFactorUnshelteredIntegral(inputHeight, massman.C->canopyHeight*2.0); //Assume flame height is double the canopy height
        WAFarrayAlbiniAbove[i] = albini.calculateWindAdjustmentFactor(0.0, 0.0, crownRatio, massman.C->canopyHeight*3.28084);  //arguments converted to feet here

        behave.updateSurfaceInputs(fuelModelNumber, moistureOneHour, moistureTenHour, moistureHundredHour, moistureLiveHerbaceous, moistureLiveWoody, WindHeightInputMode::DIRECT_MIDFLAME, WAFarrayAlbiniAbove[i]*windSpeed, windDirection, slope, aspect, canopyCover, canopyHeight, crownRatio);
        albiniAboveSpread[i] = behave.calculateSurfaceFireForwardSpreadRate(directionOfInterest);

        behave.updateSurfaceInputs(fuelModelNumber, moistureOneHour, moistureTenHour, moistureHundredHour, moistureLiveHerbaceous, moistureLiveWoody, WindHeightInputMode::DIRECT_MIDFLAME, WAFarrayMassmanAbove[i]*windSpeed, windDirection, slope, aspect, canopyCover, canopyHeight, crownRatio);
        massmanAboveSpread[i] = behave.calculateSurfaceFireForwardSpreadRate(directionOfInterest);
    }

    //----Build Albini vs Massman WAF plot-----------------

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
    plscol0(1, 0, 0, 0);        //change the color pallet0 color #1 to be black
    plscol0(2, 255, 0, 0);      //change the color pallet0 color #2 to be red
    plscol0(3, 255, 255, 255);  //change the color pallet0 color #3 to be white
    plscol0(4, 0, 0, 255);  //change the color pallet0 color #4 to be blue
    plscol0(5, 0, 255, 0);  //change the color pallet0 color #5 to be green
    plscol0(6, 255, 165, 0);  //change the color pallet0 color #6 to be orange
    plcol0(1);              //now change our font color to be the color #1

    pls->adv( 0 );
    pls->vpor( 0.15, 0.85, 0.1, 0.8 );
    pls->wind( xmin, xmax, ymin, ymax);
    pls->box( "bcnst", 0.0, 0, "bcnst", 0.0, 0 );

    pls->line(plotNodes, (PLFLT*) canopyCoverArray,(PLFLT*) WAFarrayAlbiniBelow);   //plot Albini
    plcol0(2);              //now change our font color to be the color #2
    pls->line(plotNodes, (PLFLT*) canopyCoverArray,(PLFLT*) WAFarrayMassmanBelowUniZ0gHfuel);   //plot Massman Uniform, z0g=0.13*Hfuel
    plcol0(4);              //now change our font color to be the color #4
    pls->line(plotNodes, (PLFLT*) canopyCoverArray,(PLFLT*) WAFarrayMassmanBelowUniZ0g01);   //plot Massman Uniform, z0g=0.01
    plcol0(5);              //now change our font color to be the color #5
    pls->line(plotNodes, (PLFLT*) canopyCoverArray,(PLFLT*) WAFarrayMassmanBelowGaussZ0gHfuel);   //plot Massman Asymmetric Gauss, z0g=0.13*Hfuel
    plcol0(6);              //now change our font color to be the color #6
    pls->line(plotNodes, (PLFLT*) canopyCoverArray,(PLFLT*) WAFarrayMassmanBelowGaussZ0g01);   //plot Massman Asymmetric Gauss, z0g=0.01
    plcol0(1);              //now change our font color to be the color #1

    pls->lsty(2);
    pls->line(plotNodes, (PLFLT*) canopyCoverArray,(PLFLT*) WAFarrayAlbiniAbove);   //plot Albini
    plcol0(2);              //now change our font color to be the color #2
    pls->line(plotNodes, (PLFLT*) canopyCoverArray,(PLFLT*) WAFarrayMassmanAbove);   //plot Massman
    plcol0(1);              //now change our font color to be the color #1
    pls->lsty(1);

    pls->schr(0, 1.6);  //change font size
    pls->mtex( "t", 4.0, 0.5, 0.5, "Albini vs. Massman WAF" );
    pls->schr(0, 1.0);  //change font size
    pls->mtex( "b", 3.0, 0.5, 0.5, "Canopy Cover" );
    pls->mtex( "l", 3.0, 0.5, 0.5, "WAF" );

    // Draw a legend
    PLINT        nlegend = 7;
    const char   *text[7], *symbols[7];
    PLINT        opt_array[7];
    PLINT        text_colors[7];
    PLINT        line_colors[7];
    PLINT        line_styles[7];
    PLFLT        line_widths[7];
    PLINT        symbol_numbers[7], symbol_colors[7];
    PLFLT        symbol_scales[7];
    PLFLT        legend_width, legend_height;

    // First legend entry.
    opt_array[0]   = PL_LEGEND_LINE;
    text_colors[0] = 1;
    text[0]        = "Albini, sub-canopy flames";
    line_colors[0] = 1;
    line_styles[0] = 1;
    line_widths[0] = 1.;
    symbols[0] = "";

    // Second legend entry.
    opt_array[1]      = PL_LEGEND_LINE;
    text_colors[1]    = 2;
    text[1]           = "Massman, sub-canopy flames, uniform profile, z0g=Hfuel";
    line_colors[1]    = 2;
    line_styles[1]    = 1;
    line_widths[1]    = 1.;
    //symbol_colors[1]  = 3;
    //symbol_scales[1]  = 1.;
    //symbol_numbers[1] = 4;
    symbols[1]        = "";
    // from the above opt_arrays we can completely ignore everything
    // to do with boxes.

    // Third legend entry.
    opt_array[2]      = PL_LEGEND_LINE;
    text_colors[2]    = 4;
    text[2]           = "Massman, sub-canopy flames, uniform profile, z0g=0.01";
    line_colors[2]    = 4;
    line_styles[2]    = 1;
    line_widths[2]    = 1.;
    //symbol_colors[2]  = 3;
    //symbol_scales[2]  = 1.;
    //symbol_numbers[2] = 4;
    symbols[2]        = "";
    // from the above opt_arrays we can completely ignore everything
    // to do with boxes.

    // Fourth legend entry.
    opt_array[3]      = PL_LEGEND_LINE;
    text_colors[3]    = 5;
    text[3]           = "Massman, sub-canopy flames, asymmetric Gaussian profile, z0g=Hfuel";
    line_colors[3]    = 5;
    line_styles[3]    = 1;
    line_widths[3]    = 1.;
    //symbol_colors[3]  = 3;
    //symbol_scales[3]  = 1.;
    //symbol_numbers[3] = 4;
    symbols[3]        = "";
    // from the above opt_arrays we can completely ignore everything
    // to do with boxes.

    // Fifth legend entry.
    opt_array[4]      = PL_LEGEND_LINE;
    text_colors[4]    = 6;
    text[4]           = "Massman, sub-canopy flames, asymmetric Gaussian profile, z0g=0.01";
    line_colors[4]    = 6;
    line_styles[4]    = 1;
    line_widths[4]    = 1.;
    //symbol_colors[4]  = 3;
    //symbol_scales[4]  = 1.;
    //symbol_numbers[4] = 4;
    symbols[4]        = "";
    // from the above opt_arrays we can completely ignore everything
    // to do with boxes.

    // Sixth legend entry.
    opt_array[5]   = PL_LEGEND_LINE;
    text_colors[5] = 1;
    text[5]        = "Albini, above canopy flames";
    line_colors[5] = 1;
    line_styles[5] = 2;
    line_widths[5] = 1.;
    symbols[5] = "";

    // Seventh legend entry.
    opt_array[6]      = PL_LEGEND_LINE;
    text_colors[6]    = 2;
    text[6]           = "Massman, above canopy flames";
    line_colors[6]    = 2;
    line_styles[6]    = 2;
    line_widths[6]    = 1.;
    //symbol_colors[6]  = 3;
    //symbol_scales[6]  = 1.;
    //symbol_numbers[6] = 4;
    symbols[6]        = "";
    // from the above opt_arrays we can completely ignore everything
    // to do with boxes.

    plscol0a( 15, 32, 32, 32, 0.70 );

    pls->legend(&legend_width, &legend_height,
                PL_LEGEND_BACKGROUND | PL_LEGEND_BOUNDING_BOX, 0,
                0.03, 0.03, 0.1, 3,
                1, 1, 1, 1,
                nlegend, opt_array,
                -3.0, 0.6, 1.0,
                1., text_colors, (const char **) text,
                NULL, NULL, NULL, NULL,
                line_colors, line_styles, line_widths,
                symbol_colors, symbol_scales, symbol_numbers, (const char **) symbols );

    delete pls; // close plot


    //----Build Albini vs Massman spread rate plot-------------------

    xmin = canopyCoverArray[0], ymin = 0.0, xmax = canopyCoverArray[plotNodes-1], ymax = std::max(albiniAboveSpread[0], massmanAboveSpread[1]);

    just=0, axis=0;
    //plstream *pls;

    // plplot initialization

    pls = new plstream();  // declare plplot object

    plsdev("svg"); //cairo uses the same color scheme as on screen - black on
    //red on black default
    plsfnam("massman_vs_albini_spread_rate.svg");// sets the names of the output file

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

    pls->line(plotNodes, (PLFLT*) canopyCoverArray,(PLFLT*) albiniBelowSpread);   //plot Albini
    plcol0(2);              //now change our font color to be the color #2
    pls->line(plotNodes, (PLFLT*) canopyCoverArray,(PLFLT*) massmanBelowSpreadUniZ0gHfuel);   //plot Massman
    plcol0(1);              //now change our font color to be the color #1

    pls->lsty(2);
    pls->line(plotNodes, (PLFLT*) canopyCoverArray,(PLFLT*) albiniAboveSpread);   //plot Albini
    plcol0(2);              //now change our font color to be the color #2
    pls->line(plotNodes, (PLFLT*) canopyCoverArray,(PLFLT*) massmanAboveSpread);   //plot Massman
    plcol0(1);              //now change our font color to be the color #1
    pls->lsty(1);

    pls->schr(0, 1.6);  //change font size
    pls->mtex( "t", 4.0, 0.5, 0.5, "Albini vs. Massman Spread Rate" );
    pls->schr(0, 1.0);  //change font size
    pls->mtex( "b", 3.0, 0.5, 0.5, "Canopy Cover" );
    pls->mtex( "l", 3.0, 0.5, 0.5, "Spread Rate (ch/h)" );

    // Draw a legend
//    PLINT        nlegend = 4;
//    const char   *text[4], *symbols[4];
//    PLINT        opt_array[4];
//    PLINT        text_colors[4];
//    PLINT        line_colors[4];
//    PLINT        line_styles[4];
//    PLFLT        line_widths[4];
//    PLINT        symbol_numbers[4], symbol_colors[4];
//    PLFLT        symbol_scales[4];
//    PLFLT        legend_width, legend_height;

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

    // Third legend entry.
    opt_array[2]   = PL_LEGEND_LINE;
    text_colors[2] = 1;
    text[2]        = "Albini, flames above canopy";
    line_colors[2] = 1;
    line_styles[2] = 2;
    line_widths[2] = 1.;
    symbols[2] = "";

    // Fourth legend entry.
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





    //----Build Massman flame height comparison plot--------------------------------

    xmin = flameHeightBelow[0], ymin = 0.0, xmax = flameHeightAbove[plotNodes-1], ymax = 1.0;

    just=0, axis=0;
    //plstream *pls;

    // plplot initialization

    pls = new plstream();  // declare plplot object

    plsdev("svg"); //cairo uses the same color scheme as on screen - black on
    //red on black default
    plsfnam("massman_vs_flameHeight.svg");// sets the names of the output file

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

    pls->line(plotNodes, (PLFLT*) flameHeightBelow,(PLFLT*) WAFBelowIntegral);   //plot below integral
    pls->line(plotNodes, (PLFLT*) flameHeightAbove,(PLFLT*) WAFAboveIntegral);   //plot above integral
    plcol0(2);              //now change our font color to be the color #2
    pls->line(plotNodes, (PLFLT*) flameHeightBelow,(PLFLT*) WAFBelowMidFlame);   //plot below midflame
    plcol0(1);              //now change our font color to be the color #1

    pls->schr(0, 1.6);  //change font size
    pls->mtex( "t", 4.0, 0.5, 0.5, "Massman flame height comparison" );
    pls->schr(0, 1.0);  //change font size
    pls->mtex( "b", 3.0, 0.5, 0.5, "Flame Height (m)" );
    pls->mtex( "l", 3.0, 0.5, 0.5, "WAF" );

    // Draw a legend
   nlegend = 2;
//    const char   *text[4], *symbols[4];
//    PLINT        opt_array[4];
//    PLINT        text_colors[4];
//    PLINT        line_colors[4];
//    PLINT        line_styles[4];
//    PLFLT        line_widths[4];
//    PLINT        symbol_numbers[4], symbol_colors[4];
//    PLFLT        symbol_scales[4];
//    PLFLT        legend_width, legend_height;

    // First legend entry.
    opt_array[0]   = PL_LEGEND_LINE;
    text_colors[0] = 1;
    text[0]        = "Integral method";
    line_colors[0] = 1;
    line_styles[0] = 1;
    line_widths[0] = 1.;
    symbols[0] = "";

    // Second legend entry.
    opt_array[1]      = PL_LEGEND_LINE;
    text_colors[1]    = 2;
    text[1]           = "Midflame method";
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
                0.44, 0.03, 0.1, 3,
                1, 1, 1, 1,
                nlegend, opt_array,
                2.0, 1.0, 2.0,
                1., text_colors, (const char **) text,
                NULL, NULL, NULL, NULL,
                line_colors, line_styles, line_widths,
                symbol_colors, symbol_scales, symbol_numbers, (const char **) symbols );

    delete pls; // close plot










    delete canopyCoverArray;
    canopyCoverArray = NULL;
    delete WAFarrayAlbiniBelow;
    WAFarrayAlbiniBelow = NULL;
    delete WAFarrayMassmanBelowUniZ0gHfuel;
    WAFarrayMassmanBelowUniZ0gHfuel = NULL;
    delete WAFarrayMassmanBelowUniZ0g01;
    WAFarrayMassmanBelowUniZ0g01 = NULL;
    delete WAFarrayMassmanBelowGaussZ0gHfuel;
    WAFarrayMassmanBelowGaussZ0gHfuel = NULL;
    delete WAFarrayMassmanBelowGaussZ0g01;
    WAFarrayMassmanBelowGaussZ0g01 = NULL;
    delete WAFarrayAlbiniAbove;
    WAFarrayAlbiniAbove = NULL;
    delete WAFarrayMassmanAbove;
    WAFarrayMassmanAbove = NULL;

    delete albiniBelowSpread;
    albiniBelowSpread = NULL;
    delete massmanBelowSpreadUniZ0gHfuel;
    massmanBelowSpreadUniZ0gHfuel = NULL;
    delete massmanBelowSpreadUniZ0g01;
    massmanBelowSpreadUniZ0g01 = NULL;
    delete massmanBelowSpreadGaussZ0gHfuel;
    massmanBelowSpreadGaussZ0gHfuel = NULL;
    delete massmanBelowSpreadGaussZ0g01;
    massmanBelowSpreadGaussZ0g01 = NULL;
    delete albiniAboveSpread;
    albiniAboveSpread = NULL;
    delete massmanAboveSpread;
    massmanAboveSpread = NULL;

    delete flameHeightBelow;
    flameHeightBelow = NULL;
    delete flameHeightAbove;
    flameHeightAbove = NULL;
    delete WAFBelowIntegral;
    WAFBelowIntegral = NULL;
    delete WAFBelowMidFlame;
    WAFBelowMidFlame = NULL;
    delete WAFAboveIntegral;
    WAFAboveIntegral = NULL;

    std::cout << "Done!" << std::endl;
    return 0;

}
