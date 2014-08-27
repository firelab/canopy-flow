#include "canopyFlow.h"

canopyFlow::canopyFlow()
{
    H = 20.0;
    UH = 5.0;
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
}

canopyFlow::canopyFlow(canopyFlow &rhs)
{
    H = rhs.H;
    UH = rhs.UH;
    usuh = rhs.usuh;
    if(rhs.uzc)
        memcpy(uzc, rhs.uzc, sizeof(double) * C.numNodes);
    if(rhs.uzcs)
        memcpy(uzcs, rhs.uzcs, sizeof(double) * C.numNodes);
    doh = rhs.doh;
    z0oh = rhs.z0oh;
}

canopyFlow &canopyFlow::operator=(const canopyFlow &rhs)
{
    if(&rhs != this)
    {
        H = rhs.H;
        UH = rhs.UH;
        usuh = rhs.usuh;
        if(rhs.uzc)
            memcpy(uzc, rhs.uzc, sizeof(double) * C.numNodes);
        if(rhs.uzcs)
            memcpy(uzcs, rhs.uzcs, sizeof(double) * C.numNodes);
        doh = rhs.doh;
        z0oh = rhs.z0oh;
    }
    return *this;
}

canopyFlow::~canopyFlow()
{
    delete[] uzc;
    uzc = NULL;
    delete[] uzcs;
    uzcs = NULL;
}

void canopyFlow::plot()
{
    double* y = new double[C.numNodes];
    for(int i=0; i<C.numNodes; i++)
    {
        y[i] = i*C.cellsize;
    }

    double* x = new double[C.numNodes];
    for(int i=0; i<C.numNodes; i++)
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
    pls->line( C.numNodes, (PLFLT*) x,(PLFLT*) y);

    delete pls; // close plot

    delete x;
    x = NULL;
    delete y;
    y = NULL;
}

void canopyFlow::computeWind()
{
    //C.initialize();
    c2 = c1 + K / log(C.z0gh);
    usuh = c1 - c2 * exp(-c3 *C.zetah);
    double Csurf = 2.0 * usuh * usuh;
    double nexp = C.zetah / Csurf;

    double* nzet = new double[C.numNodes];
    uzc =  new double[C.numNodes];
    uzcs =  new double[C.numNodes];

    double max_uzc = 0.0;
    for(int i=0; i<C.numNodes; i++)
    {
        nzet[i] = nexp * C.zetaz[i];
        if(i*C.cellsize <= C.z0gh)
            uzc[i] = 0.0;
        else
            uzc[i] = log(i*C.cellsize / C.z0gh) * cosh(nzet[i]) / cosh(nexp);
        if(uzc[i] > max_uzc)
            max_uzc = uzc[i];
    }

    delete nzet;

    for(int i=0; i<C.numNodes; i++)
        uzc[i] /= max_uzc;

    uzcs[0] = 0.0;
    for(int i=1;i<C.numNodes;i++) //calculate Reynolds stress using cumulative trapazoid
        uzcs[i] = uzcs[i-1] + (C.hacpz[i-1]*uzc[i-1]*uzc[i-1] + C.hacpz[i]*uzc[i]*uzc[i]) * 0.5;

    for(int i=1;i<C.numNodes;i++)
        uzcs[i] /= uzcs[C.numNodes-1];

     double inter1;

    for(int i=0; i<C.numNodes; i++)   //integrate using extened simpson's rule
    {
        if(i%2 == 0)    //if even numbers
            inter1 += uzcs[i];
        else            //if odd numbers
            inter1 += 2.0 * uzcs[i];
    }

    inter1 = inter1 - 0.5 * (uzcs[0] + uzcs[C.numNodes-1]);
    inter1 *= C.cellsize * 2.0/3.0;
    doh = 1.0 - inter1;
    z0oh = rough * inter1 * exp(-K/usuh);
    printf("doh = %lf\nz0oh = %lf\n", doh, z0oh);
}


