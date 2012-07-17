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

    std::cout << Aar;

    Eigen::VectorXd CD(xar.rows());

    CD = Aar.colPivHouseholderQr().solve(xar);

    Eigen::VectorXi Ncan(zhfrac.rows()-1);
    Ncan.array() = (int)((1.0/delh)*(zhfrac(0)-zhfrac.segment(1,zhfrac.rows()-1).array()) + 0.5);

    std::cout << Ncan;

}

