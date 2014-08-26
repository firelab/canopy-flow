#include "canopyFlow.h"

canopyFlow::canopyFlow()
{
    FAI = 1.0;
    z0h = 0.0025;
    Uph = 0.0;
    H = 20.0;
    h = 10.0;
    UH = 5.0;

    folfrac.resize(10);
    folfrac << 0.01,0.02,0.03,0.04,0.05,0.20,0.30,0.25,0.10,0.02;
    zhfrac.resize(11);
    zhfrac << 1.0,0.98,0.96,0.94,0.92,0.90,0.75,0.45,0.40,0.35,z0h;

    delh = 0.0001;
    z01.resize(((1.0-z0h)/delh)+2);
    //need to check here that there is no remainder...
    z01(0) = 1.0;
    for(int i=1; i<z01.rows(); i++)
    {
        z01(i) = z01(i-1) - delh;
    }

    //std::cout << z01;
}

canopyFlow::canopyFlow(canopyFlow &rhs)
{
    FAI = rhs.FAI;
    z0h = rhs.z0h;
    Uph = rhs.Uph;
    H = rhs.H;
    h = rhs.h;
    UH = rhs.UH;

    folfrac = rhs.folfrac;
    zhfrac = rhs.zhfrac;

    delh = rhs.delh;
    z01 = rhs.z01;

}

canopyFlow &canopyFlow::operator=(const canopyFlow &rhs)
{
    if(&rhs != this)
    {
        FAI = rhs.FAI;
        z0h = rhs.z0h;
        Uph = rhs.Uph;
        H = rhs.H;
        h = rhs.h;
        UH = rhs.UH;

        folfrac = rhs.folfrac;
        zhfrac = rhs.zhfrac;

        delh = rhs.delh;
        z01 = rhs.z01;
    }
    return *this;
}

canopyFlow::~canopyFlow()
{

}

void canopyFlow::computeWind()
{
    //Calculate the Bessel function arrays at the nodes
    Eigen::VectorXd Aj(folfrac.rows());
    Aj = (4.52 + 0.62*FAI*folfrac.array())*folfrac.array()*FAI;
    //for(int i=0; i<Aj.rows(); i++)
    //    Aj(i) = (4.52 + 0.62*FAI*folfrac(i))*folfrac(i)*FAI;

    Eigen::VectorXd aru(folfrac.rows());
    Eigen::VectorXd arl(folfrac.rows());

    //for(int i=0; i<aru.rows(); i++)
    //    aru(i) = 2.0 * std::sqrt(Aj(i) * zhfrac(i));
    aru = 2.0 * std::sqrt(Aj.array() * zhfrac.segment(0,zhfrac.rows()-1).array());

    //for(int i=0; i<arl.rows(); i++)
    //    arl(i) = 2.0 * std::sqrt((double)(Aj(i) * zhfrac(i+1)));
    arl = 2.0 * std::sqrt(Aj.array() * zhfrac.segment(1,zhfrac.rows()-1).array());

    Eigen::VectorXd I0u(aru.rows());
    Eigen::VectorXd I0l(arl.rows());
    Eigen::VectorXd K0u(aru.rows());
    Eigen::VectorXd K0l(arl.rows());
    Eigen::VectorXd I1u(aru.rows());
    Eigen::VectorXd I1l(arl.rows());
    Eigen::VectorXd K1u(aru.rows());
    Eigen::VectorXd K1l(arl.rows());

    for(int i=0; i<aru.rows(); i++)
        I0u(i) = boost::math::cyl_bessel_i(0.0, (double)aru(i));

    for(int i=0; i<arl.rows(); i++)
        I0l(i) = boost::math::cyl_bessel_i(0.0, (double)arl(i));

    for(int i=0; i<aru.rows(); i++)
        K0u(i) = boost::math::cyl_bessel_k(0.0, (double)aru(i));

    for(int i=0; i<arl.rows(); i++)
        K0l(i) = boost::math::cyl_bessel_k(0.0, (double)arl(i));

    for(int i=0; i<aru.rows(); i++)
        I1u(i) = boost::math::cyl_bessel_i(1.0, (double)aru(i));

    for(int i=0; i<arl.rows(); i++)
        I1l(i) = boost::math::cyl_bessel_i(1.0, (double)arl(i));

    for(int i=0; i<aru.rows(); i++)
        K1u(i) = boost::math::cyl_bessel_k(1.0, (double)aru(i));

    for(int i=0; i<arl.rows(); i++)
        K1l(i) = boost::math::cyl_bessel_k(1.0, (double)arl(i));

    //Initialize arrays for the vertical wind profile

    Eigen::VectorXd xar(2*zhfrac.rows()-4+2);
    xar.setZero(xar.rows());
    xar(0) = 1.0-Uph;
    xar(xar.rows()-1) = -Uph;

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Aar(2*zhfrac.rows()-2, 2*zhfrac.rows()-2);
    Aar.setZero(Aar.rows(), Aar.cols());

    //Populate array (Aar) for the vertical wind profile
    Aar(0,0) = I0u(0);
    Aar(0,1) = K0u(1);
    double ql, qu;
    for(int i=1; i<zhfrac.rows()-1; i++)
    {
        ql = std::sqrt((double)Aj(i-1));
        qu = std::sqrt((double)Aj(i));

        Aar(2*i-1,2*i-2) = I0l(i-1);
        Aar(2*i-1,2*i-1) = K0l(i-1);
        Aar(2*i-1,2*i) = -I0u(i);
        Aar(2*i-1,2*i+1) = -K0u(i);

        Aar(2*i,2*i-2) = ql*I1l(i-1);
        Aar(2*i,2*i-1) = -ql*K1l(i-1);
        Aar(2*i,2*i) = -qu*I1u(i);
        Aar(2*i,2*i+1) = qu*K1u(i);
    }

    Aar(Aar.rows()-1,Aar.cols()-2) = I0l(I0l.rows()-1);
    Aar(Aar.rows()-1,Aar.cols()-1) = K0l(K0l.rows()-1);

    //std::cout << Aar;

    Eigen::VectorXd CD(xar.rows());

    CD = Aar.colPivHouseholderQr().solve(xar);

    Eigen::VectorXi Ncan(zhfrac.rows()-1);

    double temp1 = zhfrac(0);   //have to make this temporary to evaluate expression below for some reason...  Eigen problem...

    Ncan = ((1.0/delh)*(temp1-zhfrac.segment(1,zhfrac.rows()-1).array()) + 0.5).cast<int>();
    //std::cout << Ncan;

    Eigen::VectorXd ZA(Ncan(zhfrac.rows()-2)+2);
    Eigen::VectorXd UzUh(Ncan(zhfrac.rows()-2)+2);
    Eigen::VectorXd dUun(Ncan(zhfrac.rows()-2)+2);

    temp1 = Aj(0);
    ZA.segment(0,static_cast<int> (Ncan(0)+1)) = 2.0*(temp1*z01.segment(0,static_cast<int>(Ncan(0))+1)).cwiseSqrt();
    for(int m=0; m<=Ncan(0); m++)
        UzUh(m) = CD(0)*boost::math::cyl_bessel_i(0.0,(double)ZA(m)) + CD(1)*boost::math::cyl_bessel_k(0.0,(double)ZA(m)) + Uph;
    for(int m=0; m<=Ncan(0); m++)
        dUun(m) = ZA(m)*(CD(0)*boost::math::cyl_bessel_i(1.0,(double)ZA(m))-CD(1)*boost::math::cyl_bessel_k(1.0,(double)ZA(m))) / z01(m);

    for(int n=1; n<zhfrac.rows()-2; n++)
    {
        for(int m=Ncan(n-1)+1; m<Ncan(n); m++)
        {
            ZA(m) = 2.0*sqrt((double)(Aj(n)*z01(m)));
            UzUh(m) = CD(2*n)*boost::math::cyl_bessel_i(0.0,(double)ZA(m)) + CD(2*n+1)*boost::math::cyl_bessel_k(0.0,(double)ZA(m)) + Uph;
            dUun(m) = ZA(m)*(CD(2*n)*boost::math::cyl_bessel_i(1.0,(double)ZA(m))-CD(2*n+1)*boost::math::cyl_bessel_k(1.0,(double)ZA(m))) / z01(m);
        }
    }

    temp1 = Aj(folfrac.rows()-1);
    int start = Ncan(zhfrac.rows()-3)+1;
    int sizeVec = Ncan(zhfrac.rows()-2)  -  Ncan(zhfrac.rows()-3);
    ZA.segment(start,sizeVec) = 2.0*(temp1*z01.segment(start,sizeVec)).cwiseSqrt();
    //for(int m=Ncan(zhfrac.rows()-3)+1; m<=Ncan(zhfrac.rows()-2); m++)
    for(int m=start; m<=sizeVec; m++)
        UzUh(m) = CD(CD.rows()-2)*boost::math::cyl_bessel_i(0.0,(double)ZA(m)) + CD(CD.rows()-1)*boost::math::cyl_bessel_k(0.0,(double)ZA(m)) + Uph;
    //for(int m=Ncan(zhfrac.rows()-3)+1; m<=Ncan(zhfrac.rows()-2); m++)
    for(int m=start; m<=sizeVec; m++)
        dUun(m) = ZA(m)*(CD(CD.rows()-2)*boost::math::cyl_bessel_i(1.0,(double)ZA(m))-CD(CD.rows()-1)*boost::math::cyl_bessel_k(1.0,(double)ZA(m))) / z01(m);

    Eigen::VectorXd dUdz(dUun.rows());
    dUdz = dUun.array()/dUun(0);

    double vkar = 0.4;  //Von Karman constant
    double b1 = 0.35 + vkar/log(z0h);
    double beta = 0.35 - b1*exp(-4.0*FAI);

    Eigen::VectorXd Ust2(dUdz.rows());
    Ust2 = z01.array() * dUdz.array();

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

