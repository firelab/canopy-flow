#include "canopyFlow.h"

canopyFlow::canopyFlow()
{
    C = NULL;
    usuh = 0.0;
    uzc = NULL;
    uzcs = NULL;
    doh = 0.0;
    z0oh = 0.0;
    K = 0.4;
    c1 = 0.38;
    c3 = 15.0;
    STRSC = 0.5;
    rough = 1.07;
    logRough = log(rough);
}

canopyFlow::canopyFlow(canopyFlow &rhs)
{
    if(rhs.C)
        make_canopy(rhs.C);
    else
        C = NULL;
    usuh = rhs.usuh;
    if(rhs.uzc)
    {
        uzc = new double[C->numNodes];
        memcpy(uzc, rhs.uzc, sizeof(double) * C->numNodes);
    }else{
        uzc = NULL;
    }
    if(rhs.uzcs)
    {
        uzcs = new double[C->numNodes];
        memcpy(uzcs, rhs.uzcs, sizeof(double) * C->numNodes);
    }else{
        uzcs = NULL;
    }
    doh = rhs.doh;
    z0oh = rhs.z0oh;
    K = rhs.K;
    c1 = rhs.c1;
    c3 = rhs.c3;
    STRSC = rhs.STRSC;
    rough = rhs.rough;
    logRough = rhs.logRough;
}

canopyFlow &canopyFlow::operator=(const canopyFlow &rhs)
{
    if(&rhs != this)
    {
        if(rhs.C)
            make_canopy(rhs.C);
        else
            C = NULL;
        usuh = rhs.usuh;
        if(rhs.uzc)
        {
            delete[] uzc;
            uzc = new double[C->numNodes];
            memcpy(uzc, rhs.uzc, sizeof(double) * C->numNodes);
        }else{
            uzc = NULL;
        }
        if(rhs.uzcs)
        {
            delete[] uzcs;
            uzcs = new double[C->numNodes];
            memcpy(uzcs, rhs.uzcs, sizeof(double) * C->numNodes);
        }else{
            uzcs = NULL;
        }
        doh = rhs.doh;
        z0oh = rhs.z0oh;
        K = rhs.K;
        c1 = rhs.c1;
        c3 = rhs.c3;
        STRSC = rhs.STRSC;
        rough = rhs.rough;
        logRough = rhs.logRough;
    }
    return *this;
}

canopyFlow::~canopyFlow()
{
    delete C;
    C = NULL;
    delete[] uzc;
    uzc = NULL;
    delete[] uzcs;
    uzcs = NULL;
}

void canopyFlow::plot()
{
    double* y = new double[C->numNodes];
    for(int i=0; i<C->numNodes; i++)
    {
        y[i] = i*C->cellsize;
    }

    double* x = new double[C->numNodes];
    for(int i=0; i<C->numNodes; i++)
    {
        x[i] = uzcs[i];
    }

    PLFLT xmin =0, ymin=0, xmax=1, ymax=1;

    PLINT just=0, axis=0;
    plstream *pls;

    // plplot initialization

    pls = new plstream();  // declare plplot object

    //plsdev("wxwidgets"); // sets the plot device to WX Widget which
    // allows for viewing and saving the plot to a file
    // Note that saving postscript from within widgets is buggy.
    // other useful values in place of wxwidgets:
    // xwin - X-window display to screen
    // ps - postscript file
    // psc - color postscript file
    // Or just comment out line to get a list of choices
    //plsdev("psc");
    plsdev("pdf"); //cairo uses the same color scheme as on screen - black on
    //red on black default
    plsfnam("uzcs.pdf");// sets the names of the output file

    // Parse and process command line arguments.
    //pls->parseopts( &argc, argv, PL_PARSE_FULL ); // device and other options
    // can be set from the command-line:
    // -dev devname        sets the output device to "devname"
    // -o output_file      sets the output file name to output_file                // -h                  gives a list of all possible options


    pls->init();           // start plplot object
    pls->env(xmin, xmax, ymin, ymax, just, axis );
    //Setup window size
    // - just=0 sets axis so they scale indepedently
    // - axis=0 draw axis box, ticks, and numeric labels
    //   see "man plenv" for details
    pls->lab( "uzcs", "z/h", "uzcs plot");

    // Plot the data points - (num_points, x, y. plot_symbol)
    //  - plot_symbol=9 sets a circle with a dot in the
    // middle for the plot symbol - see "man plpoin"
    //pls->poin( numNodes, (PLFLT*) x,(PLFLT*) y, 9 );
    pls->line( C->numNodes, (PLFLT*) x,(PLFLT*) y);

    delete pls; // close plot

    delete x;
    x = NULL;
    delete y;
    y = NULL;
}

void canopyFlow::plotWind(double inputSpeed, double inputHeight)
{
    int cells = 10000;
    double cellsize = inputHeight/cells;

    double* z = new double[cells];
    for(int i=0; i<cells; i++)
        z[i] = i*cellsize;

    double* windSpeed = new double[cells];
    for(int i=0; i<cells; i++)
        windSpeed[i] = get_windspeed(inputSpeed, inputHeight, i*cellsize);

    double* zCanopy = new double[C->numNodes];   //set up canopy density z axis
    for(int i=0; i<C->numNodes; i++)
        zCanopy[i] = i * C->cellsize * C->canopyHeight;

    PLFLT xmin =0, ymin=0, xmax=windSpeed[cells-1], ymax=z[cells-1];

    double maxHaz = C->haz[0];
    for(int i=0; i<C->numNodes; i++)
    {
        if(C->haz[i]>maxHaz)
            maxHaz = C->haz[i];
    }

    double maxWind = windSpeed[0];
    for(int i=0; i<cells; i++)
    {
        if(windSpeed[i]>maxWind)
            maxWind = windSpeed[i];
    }

    PLINT just=0, axis=0;
    plstream *pls;

    // plplot initialization

    pls = new plstream();  // declare plplot object

    //plsdev("wxwidgets"); // sets the plot device to WX Widget which
    // allows for viewing and saving the plot to a file
    // Note that saving postscript from within widgets is buggy.
    // other useful values in place of wxwidgets:
    // xwin - X-window display to screen
    // ps - postscript file
    // psc - color postscript file
    // Or just comment out line to get a list of choices
    //plsdev("psc");
    plsdev("pdf"); //cairo uses the same color scheme as on screen - black on
    //red on black default
    plsfnam("windProfile.pdf");// sets the names of the output file

    // Parse and process command line arguments.
    //pls->parseopts( &argc, argv, PL_PARSE_FULL ); // device and other options
    // can be set from the command-line:
    // -dev devname        sets the output device to "devname"
    // -o output_file      sets the output file name to output_file                // -h                  gives a list of all possible options

    plscolbg(255,255,255);  //change background color

    pls->init();           // start plplot object
    plscol0(1, 0, 0, 0);    //first change the color pallet0 first color(#1) to be black (0,0,0)
    plcol0(1);              //now change our font color to be the color #1

    //pls->env(xmin, xmax, ymin, ymax, just, axis );
    pls->adv( 0 );
    pls->vpor( 0.15, 0.85, 0.1, 0.8 );
    //pls->vsta();
    pls->wind( xmin, xmax, ymin, ymax);
    pls->box( "bnst", 0.0, 0, "bcnst", 0.0, 0 );
    //Setup window size
    // - just=0 sets axis so they scale indepedently
    // - axis=0 draw axis box, ticks, and numeric labels
    //   see "man plenv" for details

    //pls->lab( "wind speed", "height", "wind speed plot");

    // Plot the data points - (num_points, x, y. plot_symbol)
    //  - plot_symbol=9 sets a circle with a dot in the
    // middle for the plot symbol - see "man plpoin"
    //pls->poin( numNodes, (PLFLT*) x,(PLFLT*) y, 9 );

    //plot fill area for canopy
    double* xFill = new double[4];
    double* yFill = new double[4];
    xFill[0] = 0;
    yFill[0] = 0;
    xFill[1] = 0;
    yFill[1] = C->canopyHeight;
    xFill[2] = maxWind;
    yFill[2] = C->canopyHeight;
    xFill[3] = maxWind;
    yFill[3] = 0;
    pls->scol0a( 0, 0, 255, 0, 0.5);
    pls->col0( 0 );
    //plscol0(1, 0, 255, 0);    //first change the color pallet0 first color(#1)
    //plcol0(1);              //now change our font color to be the color #1
    pls->fill(4, xFill, yFill);

    plscol0(1, 0, 0, 0);    //first change the color pallet0 first color(#1)
    plcol0(1);              //now change our font color to be the color #1
    pls->line(cells, (PLFLT*) windSpeed,(PLFLT*) z);   //plot wind speed
    pls->schr(0, 1.6);  //change font size
    pls->mtex( "t", 4.0, 0.5, 0.5, "Canopy Wind Flow" );
    pls->schr(0, 1.0);  //change font size
    pls->mtex( "b", 3.0, 0.5, 0.5, "wind speed" );
    pls->mtex( "l", 3.0, 0.5, 0.5, "height" );
    //pls->wind(0.0, max, 0.0, zCanopy[C->numNodes-1]);  //reset window to haz coordinates
    plscol0(1, 14, 143, 14);    //first change the color pallet0 first color(#1) to be black (0,0,0)
    plcol0(1);              //now change our font color to be the color #1
    pls->wind(0.0, maxHaz, ymin, ymax);  //reset window to haz coordinates
    pls->box("cmstv", 0.0, 0, "", 0.0, 0);
    pls->line(C->numNodes, (PLFLT*) C->haz,(PLFLT*) zCanopy);   //plot haz
    pls->mtex( "t", 3.0, 0.5, 0.5, "Non-dimenstional leaf area density" );

    delete pls; // close plot

    delete windSpeed;
    windSpeed = NULL;
    delete z;
    z = NULL;
}

void canopyFlow::make_canopy(canopy::eCanopyType t)
{
    if(t == canopy::normal_distribution)
        C = new canopy_normal_distribution;
    else if(t == canopy::triangle)
        C = new canopy_normal_distribution;
    else if(t == canopy::Massman)
        C = new canopy_normal_distribution;
    else if(t == canopy::stepwise)
        C = new canopy_normal_distribution;
    else
        throw std::logic_error("ERROR: Cannot determine canopy type in canopyFlow::make_canopy().\n");
}

void canopyFlow::make_canopy(canopy* X)
{
    delete C;
    if(X->distributionType == canopy::normal_distribution)
        C = new canopy_normal_distribution(*((canopy_normal_distribution*)X));
    else if(X->distributionType == canopy::triangle)
        C = new canopy_normal_distribution(*((canopy_normal_distribution*)X));
    else if(X->distributionType == canopy::Massman)
        C = new canopy_normal_distribution(*((canopy_normal_distribution*)X));
    else if(X->distributionType == canopy::stepwise)
        C = new canopy_normal_distribution(*((canopy_normal_distribution*)X));
    else
        throw std::logic_error("ERROR: Cannot determine canopy type in canopyFlow::make_canopy().\n");
}

void canopyFlow::computeWind()
{
    //C->initialize();
    usuh = c1 - (c1 + K / log(C->z0gh)) * exp(-c3 *C->zetah);
    double Csurf = 2.0 * usuh * usuh;
    double nexp = C->zetah / Csurf;

    double* nzet = new double[C->numNodes];
    uzc =  new double[C->numNodes];
    uzcs =  new double[C->numNodes];

    double max_uzc = 0.0;
    for(int i=0; i<C->numNodes; i++)
    {
        nzet[i] = nexp * C->zetaz[i];
        if(i*C->cellsize <= C->z0gh)
            uzc[i] = 0.0;
        else
            uzc[i] = log(i*C->cellsize / C->z0gh) * cosh(nzet[i]) / cosh(nexp);
        if(uzc[i] > max_uzc)
            max_uzc = uzc[i];
    }

    delete nzet;

    for(int i=0; i<C->numNodes; i++)
        uzc[i] /= max_uzc;

    uzcs[0] = 0.0;
    for(int i=1;i<C->numNodes;i++) //calculate Reynolds stress using cumulative trapazoid
        uzcs[i] = uzcs[i-1] + (C->hacpz[i-1]*uzc[i-1]*uzc[i-1] + C->hacpz[i]*uzc[i]*uzc[i]) * 0.5;

    for(int i=1;i<C->numNodes;i++)
        uzcs[i] /= uzcs[C->numNodes-1];

    //compute lowest inflection point
    int begin = (int) 2.0 * (C->z0gh/C->cellsize);
    double difU;
    double difUOld = uzc[begin] - 2.0*uzc[begin-1] + uzc[begin-2];   //start by computing the value of the second derivative at node = begin-1
    int end;
    for(int i=begin; i<C->numNodes-1; i++)
    {
        difU = uzc[i+1] - 2.0*uzc[i] + uzc[i-1];
        if((difU*difUOld) < 0.0)   //if we've switched sign, break
        {
            end = i;
            break;
        }
        difUOld = difU;
    }

    //now set the stress below the lowest inflection point to the value at the inflection point
    for(int i=0; i<end; i++)
        uzcs[i] = uzcs[end-1];

     double inter1;

    for(int i=0; i<C->numNodes; i++)   //integrate using extened Simpson's rule
    {
        if(i%2 == 0)    //if even numbers
            inter1 += uzcs[i];
        else            //if odd numbers
            inter1 += 2.0 * uzcs[i];
    }

    inter1 = inter1 - 0.5 * (uzcs[0] + uzcs[C->numNodes-1]);
    inter1 *= C->cellsize * 2.0/3.0;
    doh = 1.0 - inter1;
    z0oh = rough * inter1 * exp(-K/usuh);
    printf("doh = %lf\nz0oh = %lf\n", doh, z0oh);
}

double canopyFlow::get_windspeed(double inputSpeed, double inputHeight, double desiredHeight)
{
    double uCanopyHeight = K * inputSpeed / (usuh * (log((inputHeight/C->canopyHeight - doh) / z0oh) + logRough));
    //double uCanopyHeight = inputSpeed * log((1.0 - doh) / z0oh) / log((inputHeight/C->canopyHeight - doh) / z0oh);
    if(desiredHeight <= C->z0g)  //if below ground roughness height
    {
        return 0.0;
    }else if(desiredHeight < C->canopyHeight)  //below canopy
    {
        return uCanopyHeight * uzc[(int)(desiredHeight/(C->canopyHeight*C->cellsize))];
    }else if(desiredHeight == C->canopyHeight)   //exactly at canopy height
    {
        return uCanopyHeight;
    }else   //above canopy height
    {
        //return uCanopyHeight * usuh * log((desiredHeight/C->canopyHeight - doh) / (z0oh + C->z0g/C->canopyHeight) / K;
        return uCanopyHeight * usuh * (log((desiredHeight/C->canopyHeight - doh) / z0oh) + logRough) / K;
    }
}


