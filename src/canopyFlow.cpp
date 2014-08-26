#include "canopyFlow.h"

canopyFlow::canopyFlow()
{
    H = 20.0;
    UH = 5.0;
}

canopyFlow::canopyFlow(canopyFlow &rhs)
{
    H = rhs.H;
    UH = rhs.UH;
}

canopyFlow &canopyFlow::operator=(const canopyFlow &rhs)
{
    if(&rhs != this)
    {
        H = rhs.H;
        UH = rhs.UH;
    }
    return *this;
}

canopyFlow::~canopyFlow()
{

}

void canopyFlow::computeWind()
{

    C.initialize();

/////////////////////////////////////////////////////////

    PLFLT xmin =0, ymin=0, xmax=5, ymax=30,
        x[6]= {0.0, 1.0, 2.0, 3.0, 4.0, 5.0},
        y[6] = {0., 1.0, 4.0, 9.1, 15.5, 25.3};
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
      plsdev("svg"); //cairo uses the same color scheme as on screen - black on
                            //red on black default
      plsfnam("test2.svg");// sets the names of the output file

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
      pls->lab( "(x)", "(y)", "PlPlot example title");

      // Plot the data points - (num_points, x, y. plot_symbol)
      //  - plot_symbol=9 sets a circle with a dot in the
      // middle for the plot symbol - see "man plpoin"
      pls->poin( 6, x, y, 9 );

      delete pls; // close plot









}

