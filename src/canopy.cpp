#include "canopy.h"

canopy::canopy()
{
    leafAreaIndex = 1.0;
    canopyHeight = 10.0;
    z0g = 0.025;
    dragCoefAth = 0.2;
    heightMaxFoliageDist = 0.5;
    standardDevFoliageDist = 0.25;
    numNodes = 101;   //MUST BE ODD NUMBER
    cellsize = 0.0;    //cellsize here is normalized from 0 to 1
    cumulativeLeafDragArea = NULL;
    haz = NULL;
    hacpz = NULL;
    zetaz = NULL;
    zetah = 0.0;
    z0gh = 0.0;
}

canopy::canopy(canopy &rhs)
{
    leafAreaIndex = rhs.leafAreaIndex;
    canopyHeight = rhs.canopyHeight;
    z0g = rhs.z0g;
    dragCoefAth = rhs.dragCoefAth;
    heightMaxFoliageDist = rhs.heightMaxFoliageDist;
    standardDevFoliageDist = rhs.standardDevFoliageDist;
    numNodes = rhs.numNodes;
    cellsize = rhs.cellsize;
    totalDragAreaIndex = rhs.totalDragAreaIndex;
    if(rhs.cumulativeLeafDragArea)
        memcpy(cumulativeLeafDragArea, rhs.cumulativeLeafDragArea, sizeof(double) * numNodes);
    if(rhs.haz)
        memcpy(haz, rhs.haz, sizeof(double) * numNodes);
    if(rhs.hacpz)
        memcpy(hacpz, rhs.hacpz, sizeof(double) * numNodes);
    if(rhs.zetaz)
        memcpy(zetaz, rhs.zetaz, sizeof(double) * numNodes);
    zetah = rhs.zetah;
    z0gh = rhs.z0gh;
}

canopy &canopy::operator=(const canopy &rhs)
{
    if(&rhs != this)
    {
        leafAreaIndex = rhs.leafAreaIndex;
        canopyHeight = rhs.canopyHeight;
        z0g = rhs.z0g;
        dragCoefAth = rhs.dragCoefAth;
        heightMaxFoliageDist = rhs.heightMaxFoliageDist;
        standardDevFoliageDist = rhs.standardDevFoliageDist;
        numNodes = rhs.numNodes;
        cellsize = rhs.cellsize;
        totalDragAreaIndex = rhs.totalDragAreaIndex;
        if(rhs.cumulativeLeafDragArea)
            memcpy(cumulativeLeafDragArea, rhs.cumulativeLeafDragArea, sizeof(double) * numNodes);
        if(rhs.haz)
            memcpy(haz, rhs.haz, sizeof(double) * numNodes);
        if(rhs.hacpz)
            memcpy(hacpz, rhs.hacpz, sizeof(double) * numNodes);
        if(rhs.zetaz)
            memcpy(zetaz, rhs.zetaz, sizeof(double) * numNodes);
        z0gh = rhs.z0gh;
    }
    return *this;
}

canopy::~canopy()
{
    delete[] cumulativeLeafDragArea;
    cumulativeLeafDragArea = NULL;
    delete[] haz;
    haz = NULL;
    delete[] hacpz;
    hacpz = NULL;
    delete[] zetaz;
    zetaz = NULL;
}

void canopy::initialize()
{
    cumulativeLeafDragArea = new double[numNodes];
    haz = new double[numNodes];
    hacpz = new double[numNodes];
    zetaz = new double[numNodes];
    z0gh = z0g/canopyHeight;
    cellsize = 1.0 / (numNodes - 1);    //cellsize here is normalized from 0 to 1

    compute_haz();
    compute_foliage_drag_area_index();
}

void canopy::plot()
{
    double max = haz[0];
    for(int i=0; i<numNodes; i++)
    {
        if(haz[i]>max)
            max = haz[i];
    }
    printf("max = %lf\n", max);

//    PLFLT xmin =0, ymin=0, xmax=1, ymax=1,
//            x[6]= {0.0, 1.0, 2.0, 3.0, 4.0, 5.0},
//            y[6] = {0., 1.0, 4.0, 9.1, 15.5, 25.3};

    double* y = new double[numNodes];
    for(int i=0; i<numNodes; i++)
    {
        y[i] = i*cellsize;
    }

    double* x = new double[numNodes];
    for(int i=0; i<numNodes; i++)
    {




        //---------------FIX THIS------------------!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        //x[i] = haz[i]/max;
        x[i] = haz[i];
        //---------------FIX THIS------------------!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





        if(i%100 == 0)
            printf("%lf\t%lf\n",x[i], y[i]);
    }

    PLFLT xmin =0, ymin=0, xmax=0.001, ymax=1;

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
    plsfnam("haz.pdf");// sets the names of the output file

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
    pls->lab( "haz/max", "z/h", "haz plot");

    // Plot the data points - (num_points, x, y. plot_symbol)
    //  - plot_symbol=9 sets a circle with a dot in the
    // middle for the plot symbol - see "man plpoin"
    //pls->poin( numNodes, (PLFLT*) x,(PLFLT*) y, 9 );
    pls->line( numNodes, (PLFLT*) x,(PLFLT*) y);

    delete pls; // close plot

    delete x;
    x = NULL;
    delete y;
    y = NULL;
}

void canopy::compute_haz()
{
    double norm;
    double hazn;
    double integHazn = 0.0;





    //---------------FIX THIS------------------!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //double intermediate = cellsize*2.0 / 3.0;
    double intermediate = 2.0 / 3.0;
    //---------------FIX THIS------------------!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!









    for(int i=0; i<numNodes; i++)   //integrate using extened simpson's rule
    {
        norm = (i*cellsize - heightMaxFoliageDist) / standardDevFoliageDist;
        haz[i] = exp(-norm * norm); //temporarily store this here
        if(i%2 == 0)    //if even numbers
            integHazn += haz[i];
        else            //if odd numbers
            integHazn += 2.0 * haz[i];
    }

    integHazn = integHazn - 0.5 * (haz[0] + haz[numNodes-1]);
    integHazn *= intermediate;

    for(int i=0; i<numNodes; i++)
        haz[i] = leafAreaIndex * haz[i] /integHazn;
}

double canopy::get_dragCoef(double zOverh)
{
    double d1 = 0.0;    //This is drag parameter #1.  It's value might be -1 <= d1 <= +1
    return dragCoefAth * exp(-d1 * (1.0 - zOverh));
}
double canopy::get_shelterFactor(double hazi)
 {
     double p1 = 0.0;   //This is the momentum shelter factor parameter.  Likely values are <= 1/3 * numcells
     return 1.0 + p1 * hazi;
 }

void canopy::compute_foliage_drag_area_index()
{

    hacpz[0] = get_dragCoef(0) * haz[0] / get_shelterFactor(0);
    zetaz[0] = 0.0;
    for(int i=1;i<numNodes;i++)
    {
        hacpz[i] = get_dragCoef(i*cellsize) * haz[i] / get_shelterFactor(i*cellsize);
        zetaz[i] = zetaz[i-1] + (hacpz[i-1] + hacpz[i]) * 0.5;
    }
    zetah = zetaz[numNodes-1];
    for(int i=1;i<numNodes;i++)
        zetaz[i] /= zetah;
}



