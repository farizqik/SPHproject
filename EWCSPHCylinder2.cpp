#include<iostream>
#include<cmath>
#include <string>
#include <iomanip>
#include <fstream>
#include <chrono>
#include <filesystem>
#include <sstream>
#include <vector>
#include <limits>


using namespace std;

// --------------------------------------------
// Timestep
// --------------------------------------------

const double dtMax = 0.001;
const double dtMin = 1.0e-10;
const double CFL = 0.3;
const double Totaltime = 10.0;

// ------------------------------------------------------------
// Geometry
// ------------------------------------------------------------

const double waterheight = 2.0;


// ------------------------------------------------------------
// Some Constants
// ------------------------------------------------------------
const double PI = 3.14159265358979323846;
const double g = 9.81;
const double rho0 = 1000.0;

const double velcoefX = 0.1;
const double velcoefY = 0.0;

const double inletVelocity = 0.1;  // m/s

const double c0 = 10.0*sqrt(g*(waterheight));
//const double c0 = 20.0*inletVelocity;

const double gammaEOS = 7.0;
const double B = c0*c0*rho0/gammaEOS;

const double viscosity = 1.0e-4;
const double tvis = 1.0;
const double Dcylinder = 0.1;
const double r0 = 0.5*Dcylinder;
const double R = 10*Dcylinder;
const double dr0 = sqrt(2.0*viscosity*tvis);
//const double dr0 = 0.1*Dcylinder;
//double drhodtexact = -100;
const double alphaAV = 0.01;
const double deltadifussion = 0.1;

const int Nr = 30;
const int Nbuffer = 8;      // open-boundary buffer rings
int NrTotal = Nr + Nbuffer;

double Reynolds =
    inletVelocity*Dcylinder/viscosity;



string kernel; 
string type; 




// -----------------------------------------------------------------------------------------------------------------------------
// Kernel functions
// -----------------------------------------------------------------------------------------------------------------------------

struct KernelResult {
    double Weight;
    double dWeightX;
    double dWeightY;
};


// -----------------------------------------------------------------------------------------------------------------------------
// Gaussian kernel
// -----------------------------------------------------------------------------------------------------------------------------
KernelResult gaussian(double q, double h, double dirX, double dirY)
{
    KernelResult result;
    double alpha = 1.0 / (PI * h * h);

    result.Weight = alpha*exp(-q*q);
    result.dWeightX = -2.0*alpha*q*exp(-q*q)/h * dirX;
    result.dWeightY = -2.0*alpha*q*exp(-q*q)/h * dirY;

    return result;
}


// -----------------------------------------------------------------------------------------------------------------------------
// Cubic kernel
// -----------------------------------------------------------------------------------------------------------------------------
KernelResult cubicSpline(double q, double h, double dirX, double dirY)
{
    double alpha = 10.0 / (7.0 * PI * h * h);
    KernelResult result;
    if (q<=1.0)
    {
        result.Weight = alpha*(1.0-1.5*pow(q,2)+0.75*pow(q,3));
        result.dWeightX = alpha*(-3.0*q+2.25*pow(q,2))/h * dirX;
        result.dWeightY = alpha*(-3.0*q+2.25*pow(q,2))/h * dirY;
    }
    else if (q<=2.0)
    {
        result.Weight = alpha*0.25*pow(2.0-q,3);
        result.dWeightX = -0.75*alpha*pow(2.0-q,2)/h * dirX;
        result.dWeightY = -0.75*alpha*pow(2.0-q,2)/h * dirY;
    }
    else
    {
        result.Weight = 0.0;
        result.dWeightX = 0.0;
        result.dWeightY = 0.0;
    }

    return result;
}

// -----------------------------------------------------------------------------------------------------------------------------
// Wendland kernel
// -----------------------------------------------------------------------------------------------------------------------------
KernelResult Wendland(double q, double h, double dirX, double dirY)
{
    double alpha = 7.0 / (4.0 * PI * h * h);
    KernelResult result;
    if (0.0<=q && q<=2.0)
    {
        result.Weight = alpha*pow(1.0-0.5*q,4)*(2.0*q+1.0);
        result.dWeightX = alpha*(4*pow(1.0-0.5*q,3)*(-0.5)*(2.0*q+1.0)+pow(1.0-0.5*q,4)*2.0)/h * dirX;
        result.dWeightY = alpha*(4*pow(1.0-0.5*q,3)*(-0.5)*(2.0*q+1.0)+pow(1.0-0.5*q,4)*2.0)/h * dirY;
    }
    else
    {
        result.Weight = 0.0;
        result.dWeightX = 0.0;
        result.dWeightY = 0.0;
    }

    return result;
}


// --------------------------------------------------
// Robust 3x3 linear solver with partial pivoting.
// Returns false when A is singular / nearly singular.
// --------------------------------------------------
bool solveLinear(double A[][3], double b[], double X[], int n)
{
    double maxA = 0.0;
    for (int row = 0; row < n; row++)
    {
        for (int col = 0; col < n; col++)
        {
            maxA = max(maxA, abs(A[row][col]));
        }
    }

    if (maxA == 0.0)
    {
        return false;
    }

    const double pivotTolerance = 1e-10 * maxA;

    for (int k = 0; k < n; k++)
    {
        // Find the largest pivot in this column.
        int pivotRow = k;
        double pivotValue = abs(A[k][k]);

        for (int row = k + 1; row < n; row++)
        {
            if (abs(A[row][k]) > pivotValue)
            {
                pivotValue = abs(A[row][k]);
                pivotRow = row;
            }
        }

        if (pivotValue < pivotTolerance)
        {
            return false;
        }

        if (pivotRow != k)
        {
            for (int col = k; col < n; col++)
            {
                swap(A[k][col], A[pivotRow][col]);
            }
            swap(b[k], b[pivotRow]);
        }

        // Elimination below the pivot.
        for (int row = k + 1; row < n; row++)
        {
            double factor = A[row][k] / A[k][k];

            for (int col = k; col < n; col++)
            {
                A[row][col] -= factor * A[k][col];
            }
            b[row] -= factor * b[k];
        }
    }

    // Back substitution.
    for (int row = n - 1; row >= 0; row--)
    {
        if (abs(A[row][row]) < pivotTolerance)
        {
            return false;
        }

        double sum = b[row];
        for (int col = row + 1; col < n; col++)
        {
            sum -= A[row][col] * X[col];
        }

        X[row] = sum / A[row][row];

        if (!isfinite(X[row]))
        {
            return false;
        }
    }

    return true;
}


// --------------------------------------------------
// Calculate adaptive timestep
// 
// --------------------------------------------------

double calculateDt(int Nboundary,
    int Nparticles, int NfluidEnd,
    const vector<double>& x,
    const vector<double>& y,
    const vector<double>& u,
    const vector<double>& v,
    const vector<double>& rho,
    const vector<double>& pressure,
    const vector<double>& dudt,
    const vector<double>& dvdt,
    const vector<double>& h,
    double c0,
    const string& kernel)
    {
        double dtForce = numeric_limits<double>::infinity(); //initialise dtForce to positive infinity (will be replaced if there are any fluid particles)
        double dtCFL = numeric_limits<double>::infinity(); //initialise dtCFL to positive infinity (will be replaced if there are any fluid particles)

        for (int i = Nboundary; i < NfluidEnd; ++i)
        {
            // Compute dtForce based on acceleration magnitude
            double ax = dudt[i];

            double ay = dvdt[i];

            double accelerationMagnitude =
                sqrt(ax*ax + ay*ay);

            if (accelerationMagnitude > 1e-14)
            {
                double particleDTForce = sqrt(h[i] / accelerationMagnitude);
                dtForce = min(dtForce, particleDTForce);
            }

            // Compute dtCFL based on velocity magnitude
            double maxRelativeSpeed = 0.0;

            for (int j = 0; j < Nparticles; ++j)
            {
                double sr = x[i] - x[j];
                double sz = y[i] - y[j];

                double s2 = sr*sr + sz*sz;

                if (s2 < 1e-14)
                {
                    continue;
                }

                if ((kernel == "cubic" || kernel == "wendland") && s2 > 4.0*h[i]*h[i])
                {
                    continue;
                }

                double du = u[i] - u[j];

                double dv = v[i] - v[j];

                double VabRab = du*sr + dv*sz;

                double RelativeSpeed = abs(h[i]*VabRab/(s2 + 0.01*h[i]*h[i]));

                maxRelativeSpeed = max(maxRelativeSpeed, RelativeSpeed);

            }

            double particleDTCFL = h[i] / (maxRelativeSpeed + c0);

            dtCFL = min(dtCFL, particleDTCFL);


        }
        //cout << "Acceleration-based dt = "<< dtForce << endl;
        //cout << "Velocity-based dt = "<< dtCFL << endl;

        double adaptivedt = min(dtForce, dtCFL);
        

        return adaptivedt;
        
    }


// --------------------------------------------------
// Compute the interaction between particle i and j
// --------------------------------------------------
void accumulateInteraction(
    int i,
    int j,
    const vector<double>& xs,
    const vector<double>& ys,
    const vector<double>& us,
    const vector<double>& vs,
    const vector<double>& rhos,
    const vector<double>& ps,
    double Vj,
    double mj,
    double h,
    const string& kernel,
    double& drhoAcc,
    double& duAcc,
    double& dvAcc)
{
    double rx = xs[i] - xs[j];
    double ry = ys[i] - ys[j];
    double r2 = rx*rx + ry*ry;

    if ((kernel == "cubic" || kernel == "wendland") &&
        r2 > 4.0*h*h)
    {
        return;
    }

    if (r2 < 1e-14)
    {
        return;
    }

    double r = sqrt(r2);
    double q = r/h;
    double dirX = rx/r;
    double dirY = ry/r;

    KernelResult result;

    if (kernel == "gaussian")
        result = gaussian(q, h, dirX, dirY);
    else if (kernel == "cubic")
        result = cubicSpline(q, h, dirX, dirY);
    else if (kernel == "wendland")
        result = Wendland(q, h, dirX, dirY);
    else
        return;

    double du = us[i] - us[j];
    double dv = vs[i] - vs[j];

    double vijrij = du*rx + dv*ry;
    double Piij = 0.0;

    if (vijrij < 0.0)
    {
        double muij =
            h*vijrij/(r2 + 0.01*h*h);

        double rhoij =
            0.5*(rhos[i] + rhos[j]);

        Piij =
            -alphaAV*c0*muij/rhoij;
    }

    double velocityDotGradient =
    us[i]*result.dWeightX
    + vs[i]*result.dWeightY;

    double relativeVelocityDotGradient =
        du*result.dWeightX
        + dv*result.dWeightY;

    double diffusion =
        2.0*deltadifussion*h*c0*Vj
        *(rhos[i]-rhos[j])
        *(rx*result.dWeightX + ry*result.dWeightY)
        /(r2 + 0.01*h*h);

    // Eulerian continuity
    drhoAcc +=
        -Vj*(rhos[j]-rhos[i])*velocityDotGradient
        + rhos[i]*Vj*relativeVelocityDotGradient
        + diffusion;


    // Symmetric pressure + artificial viscosity
    double pressureTerm =
        ps[i]/(rhos[i]*rhos[i])
        + ps[j]/(rhos[j]*rhos[j])
        + Piij;

    // Eulerian momentum
    duAcc -=
        Vj*(us[j]-us[i])*velocityDotGradient
        + Vj*rhos[j]*pressureTerm*result.dWeightX;

    dvAcc -=
        Vj*(vs[j]-vs[i])*velocityDotGradient
        + Vj*rhos[j]*pressureTerm*result.dWeightY;

// --------------------------------------------------
// Physical viscosity
// nu * Laplacian(u)
// --------------------------------------------------

    double rijDotGradW =
        rx*result.dWeightX
        + ry*result.dWeightY;

    double viscousX =
        2.0*viscosity*Vj
        *(us[i]-us[j])
        *rijDotGradW
        /(r2 + 0.01*h*h);

    double viscousY =
        2.0*viscosity*Vj
        *(vs[i]-vs[j])
        *rijDotGradW
        /(r2 + 0.01*h*h);

    duAcc += viscousX;
    dvAcc += viscousY;
}


// ------------------------------------------------------------
// solve Q in radial function
// ------------------------------------------------------------

double target = (R-r0)/dr0;

double radialSum(double qA, int N)
{
    if (abs(qA - 1.0) < 1.0e-12)
    {
        return N;
    }

    return (pow(qA, N) - 1.0) / (qA - 1.0);
}



// --------------------------------------------------
// Stationary free-slip cylinder: u_n = 0 and du_t/dn ~ 0.
// The boundary point and the first fluid point share the same angle.


void applyCylinderFreeSlip(
    int i, int NTheta,
    const vector<double>& theta,
    vector<double>& u,
    vector<double>& v)
{
    const int j = i + NTheta;
    const double tx = -sin(theta[i]);
    const double ty =  cos(theta[i]);
    const double ut = u[j]*tx + v[j]*ty;
    u[i] = ut*tx;
    v[i] = ut*ty;
}


// --------------------------------------------------
// Compute mDBC interpolation for a ghost particle
// --------------------------------------------------

void interpolateMDBC(
    int i,
    int Nboundary,
    int NfluidEnd,
    int NTheta,

    const vector<double>& xState,
    const vector<double>& yState,
    const vector<double>& rhoState,

    const vector<double>& xGhost,
    const vector<double>& yGhost,

    const vector<double>& Vi,
    const vector<double>& hp,

    const string& kernel,

    double& rhoGhostOut,
    double& drhoGhostXOut,
    double& drhoGhostYOut,
    double& rhoBoundaryOut)
{
    const int matrixSize = 3;

    double A[matrixSize][matrixSize] = {0.0};
    double b[matrixSize] = {0.0};
    double X[matrixSize] = {0.0};

    int neighborCount = 0;

    double shepardNumerator = 0.0;
    double shepardDenominator = 0.0;

    // h representative of the fluid near this boundary particle
    double h = hp[i + NTheta];

    for (int j = Nboundary;
         j < NfluidEnd;
         ++j)
    {
        double rx =
            xGhost[i] - xState[j];

        double ry =
            yGhost[i] - yState[j];

        double r2 =
            rx*rx + ry*ry;

        if ((kernel == "cubic" ||
             kernel == "wendland") &&
            r2 > 4.0*h*h)
        {
            continue;
        }

        double distance = sqrt(r2);
        double q = distance/h;

        // Optional Gaussian truncation
        if (kernel == "gaussian" && q > 3.0)
        {
            continue;
        }

        double dirX = 0.0;
        double dirY = 0.0;

        if (distance > 1.0e-14)
        {
            dirX = rx/distance;
            dirY = ry/distance;
        }

        KernelResult result;

        if (kernel == "gaussian")
        {
            result =
                gaussian(q,h,dirX,dirY);
        }
        else if (kernel == "cubic")
        {
            result =
                cubicSpline(q,h,dirX,dirY);
        }
        else if (kernel == "wendland")
        {
            result =
                Wendland(q,h,dirX,dirY);
        }
        else
        {
            throw invalid_argument(
                "Kernel must be gaussian, cubic, or wendland.");
        }

        if (abs(result.Weight) < 1.0e-14)
        {
            continue;
        }

        if (!isfinite(rhoState[j]) ||
            rhoState[j] <= 0.0)
        {
            continue;
        }

        neighborCount++;

        // Fixed Eulerian particle/control volume
        double Vj = Vi[j];

        // Shepard interpolation
        shepardNumerator +=
            rhoState[j]*
            result.Weight*
            Vj;

        shepardDenominator +=
            result.Weight*
            Vj;

        // x_j - x_g and y_j - y_g
        double dx =
            -rx;

        double dy =
            -ry;

        double dA[matrixSize][matrixSize] =
        {
            {
                result.Weight*Vj,
                result.Weight*Vj*dx,
                result.Weight*Vj*dy
            },
            {
                result.dWeightX*Vj,
                result.dWeightX*Vj*dx,
                result.dWeightX*Vj*dy
            },
            {
                result.dWeightY*Vj,
                result.dWeightY*Vj*dx,
                result.dWeightY*Vj*dy
            }
        };

        double db[matrixSize] =
        {
            rhoState[j]*
            result.Weight*Vj,

            rhoState[j]*
            result.dWeightX*Vj,

            rhoState[j]*
            result.dWeightY*Vj
        };

        for (int row = 0;
             row < matrixSize;
             ++row)
        {
            for (int column = 0;
                 column < matrixSize;
                 ++column)
            {
                A[row][column] +=
                    dA[row][column];
            }

            b[row] +=
                db[row];
        }
    }

    const int minMdbcNeighbors = 4;
    const double dryTolerance = 1.0e-12;
    const double supportTolerance = 0.4;

    bool hasFluid =
        shepardDenominator >
        dryTolerance;

    bool hasGoodSupport =
        shepardDenominator >
        supportTolerance;

    bool solved = false;

    if (hasGoodSupport &&
        neighborCount >= minMdbcNeighbors)
    {
        solved =
            solveLinear(
                A,b,X,matrixSize);
    }

    if (!hasFluid)
    {
        rhoGhostOut = rho0;

        drhoGhostXOut = 0.0;
        drhoGhostYOut = 0.0;
    }
    else if (solved)
    {
        rhoGhostOut = X[0];

        drhoGhostXOut = X[1];
        drhoGhostYOut = X[2];
    }
    else
    {
        rhoGhostOut =
            shepardNumerator/
            shepardDenominator;

        drhoGhostXOut = 0.0;
        drhoGhostYOut = 0.0;
    }

    rhoBoundaryOut =
        rhoGhostOut
        + drhoGhostXOut*
          (xState[i] - xGhost[i])
        + drhoGhostYOut*
          (yState[i] - yGhost[i]);

    // only prevent nonphysical density
    rhoBoundaryOut =
        max(rhoBoundaryOut,
            1.0e-6*rho0);
}
// ------------------------------------------------------------
// Riemann Invariants
// ------------------------------------------------------------

void applyRiemannBufferBoundary(
    vector<double>& u,
    vector<double>& v,
    vector<double>& rho,
    const vector<double>& theta,
    int NTheta,
    int Nr,
    int NrTotal)
{
    for (int ir = Nr + 1; ir <= NrTotal; ++ir)
    {
        for (int k = 0; k < NTheta; ++k)
        {
            // Buffer particle
            int i = ir*NTheta + k;

            // Last physical-fluid particle
            // at the same theta
            int j = Nr*NTheta + k;

            double nx = cos(theta[i]);
            double ny = sin(theta[i]);

            double tx = -ny;
            double ty =  nx;

            // Interior physical state
            double unInterior =
                u[j]*nx + v[j]*ny;

            double utInterior =
                u[j]*tx + v[j]*ty;

            double cInterior =
                c0*pow(
                    rho[j]/rho0,
                    0.5*(gammaEOS-1.0));

            // Leaving computational domain
            double Jplus =
                unInterior
                + 2.0*cInterior/
                  (gammaEOS-1.0);

            // Far field
            double unFar =
                inletVelocity*nx;

            double utFar =
                inletVelocity*tx;

            double Jminus =
                unFar
                - 2.0*c0/
                  (gammaEOS-1.0);

            // Reconstruct
            double unBuffer =
                0.5*(Jplus + Jminus);

            double cBuffer =
                0.25*(gammaEOS-1.0)
                *(Jplus-Jminus);

            cBuffer =
                max(cBuffer,1.0e-8);

            rho[i] =
                rho0*pow(
                    cBuffer/c0,
                    2.0/(gammaEOS-1.0));

            double utBuffer;

            if (unFar < 0.0)
            {
                // inflow
                utBuffer = utFar;
            }
            else
            {
                // outflow
                utBuffer = utInterior;
            }

            u[i] =
                unBuffer*nx
                + utBuffer*tx;

            v[i] =
                unBuffer*ny
                + utBuffer*ty;
        }
    }
}

// --------------------------------------------------
// Validate cylinder pressure-force integration
// using analytical inviscid potential flow
// --------------------------------------------------

void validatePotentialFlowForce(
    const vector<double>& theta,
    int Nboundary,
    double r0,
    double dTheta,
    double rhoInf,
    double UInf,
    double Dcylinder,
    const string& foldername)
{
    double dynamicPressure =
        0.5*rhoInf*UInf*UInf;

    double forceReference =
        dynamicPressure*Dcylinder;

    double Fx = 0.0;
    double Fy = 0.0;

    string filename =
        foldername
        + "/PotentialFlowValidation.csv";

    ofstream file(filename);

    file
        << "ID,theta,CpExact,pExact,"
        << "nx,ny,ds,dFx,dFy"
        << endl;

    for (int i = 0;
         i < Nboundary;
         ++i)
    {
        double th =
            theta[i];

        // Outward normal from cylinder into fluid
        double nx =
            cos(th);

        double ny =
            sin(th);

        // Analytical potential-flow pressure coefficient
        double CpExact =
            1.0
            - 4.0*sin(th)*sin(th);

        // Gauge pressure: p_inf = 0
        double pExact =
            dynamicPressure*CpExact;

        // Arc represented by this boundary particle
        double ds =
            r0*dTheta;

        // Pressure force on cylinder
        double dFx =
            -pExact*nx*ds;

        double dFy =
            -pExact*ny*ds;

        Fx += dFx;
        Fy += dFy;

        file
            << i << ","
            << th << ","
            << CpExact << ","
            << pExact << ","
            << nx << ","
            << ny << ","
            << ds << ","
            << dFx << ","
            << dFy
            << endl;
    }

    double CD =
        Fx/forceReference;

    double CL =
        Fy/forceReference;


    file << endl;

    file
        << "Fx" << ","
        << Fx << endl;

    file
        << "Fy" << ","
        << Fy << endl;

    file
        << "CD" << ","
        << CD << endl;

    file
        << "CL" << ","
        << CL << endl;

    file.close();

    cout << endl;
    cout << "====================================" << endl;
    cout << "Potential-flow force validation" << endl;
    cout << "====================================" << endl;

    cout << "Fx = "
         << Fx << endl;

    cout << "Fy = "
         << Fy << endl;

    cout << "CD = "
         << CD << endl;

    cout << "CL = "
         << CL << endl;

    cout << "Exact CD = 0" << endl;
    cout << "Exact CL = 0" << endl;

    cout << "Validation file: "
         << filename << endl;

    cout << "===================================="
         << endl;
}


// ------------------------------------------------------------
// Output Parameters
// ------------------------------------------------------------

vector<double> r;
vector<double> theta;


vector<double> PGauss;
vector<double> PCubic;
vector<double> PWedn;
vector<double> dPGaussX;
vector<double> dPGaussY;
vector<double> dPCubicX;
vector<double> dPCubicY;
vector<double> dPWednX;
vector<double> dPWednY;
vector<double> drhodtGauss;
vector<double> drhodtCubic;
vector<double> drhodtWedn;
vector<double> hp;


vector<double> x;
vector<double> y;

vector<double> u;
vector<double> v;

vector<double> xnew;
vector<double> ynew;

vector<double> unew;
vector<double> vnew;

vector<double> rho;
vector<double> rhonew;

vector<double> drhodtexact;
vector<double> pressureexact;

vector<double> pressure;
vector<double> pressurenew;

vector<double> drhodt;

vector<double> dudt;
vector<double> dvdt;

vector<double> xghost;
vector<double> yghost;
vector<double> rhoghost;
vector<double> drhoghostX;
vector<double> drhoghostY;

vector<double> xhalf;
vector<double> yhalf;
vector<double> uhalf;
vector<double> vhalf;
vector<double> rhohalf;

vector<double> drhodthalf;
vector<double> pressurehalf;
vector<double> dudthalf;
vector<double> dvdthalf;

vector<double> mass;
vector<double> Vi;

vector<double> vorticity;













// -----------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------
// Main Program
// -----------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------

int main()
{
    double qAlow = 0.0;
    double qAhigh = 2.0;
    double maxiter = 1000;
    double tolerance = 1.0e-6;
    double qA = 1.0;
    
    for (int iter = 0; iter < maxiter; ++iter)
    {
        double qAmid = 0.5*(qAlow+qAhigh);
        double error = radialSum(qAmid, Nr) - target;

        if (abs(error) < tolerance)
        {
            qA = qAmid;
            break;
        }
        if (radialSum(qAmid,Nr) < target)
        {
            qAlow = qAmid;
        }
        else
        {
            qAhigh = qAmid;
        }

        qA = qAmid;
        
    }

    double Apolar = log(qA);
    double Bpolar = dr0/(qA-1);

    cout << "qA = " << qA << endl;
    cout << "Apolar = " << Apolar << endl;
    cout << "Bpolar = " << Bpolar << endl;
    cout << "================================" << endl;
    cout << "" << endl;




    
    int NTheta =  round(2*PI*r0/dr0);
    double dTheta0 = 2*PI/NTheta;


    

    // ------------------------------------------------------------
    // Generate folder
    // ------------------------------------------------------------
    auto start = std::chrono::high_resolution_clock::now();

    cout<<"enter the kernel to be used (gaussian, cubic, wendland): ";
        cin>>kernel;

    string foldername = 
                "EWCSPHCylinder_dr0_" + to_string(dr0) + "_Nr_" + to_string(Nr) + "_Ntheta_" + to_string(NTheta) + "_Rr0_" + to_string(R/r0); 
            //system(("mkdir -p " + foldername).c_str());
            filesystem::create_directories(foldername);

        // Remove previous CSV output from this exact simulation folder.
        for (const auto& entry :
            std::filesystem::directory_iterator(foldername))
            {
                if (entry.is_regular_file() &&
                    entry.path().extension() == ".csv")
                {
                    std::filesystem::remove(entry.path());
                }
            }

        cout << "Previous CSV files deleted from: "
            << foldername << '\n';




    // ------------------------------------------------------------
    // Generate particles
    // ------------------------------------------------------------
    // ------------------------------------------------------------
    // Boundary particles
    // ------------------------------------------------------------

    int i = 0;
    double rp = r0;
    for (int j = 0; j < NTheta; ++j)
    {
        double thetap = j*dTheta0;
        double xp = rp*cos(thetap);
        double yp = rp*sin(thetap);

        r.push_back(rp);
        theta.push_back(thetap);

        x.push_back(xp);
        y.push_back(yp);


        cout << "ID: " << i
            << "  Type: Boundary"
            << "  r: " << rp
            << "  theta: " << thetap
            << "  x: " << xp
            << "  y: " << yp
            << endl;

        i++;

    }

    int Nboundary = i;

    // ------------------------------------------------------------
    // Fluid particles
    // ------------------------------------------------------------
    
    
    for (int ir = 1; ir <= NrTotal; ++ir)
    {
        double rp = Bpolar*(pow(qA,ir)-1)+r0;

        for (int j = 0; j < NTheta; ++j)
        {
            double thetap = j*dTheta0;
            double xp = rp*cos(thetap);
            double yp = rp*sin(thetap);

            r.push_back(rp);
            theta.push_back(thetap);

            x.push_back(xp);
            y.push_back(yp);

            

            cout << "ID: " << i
                << "  Type: "
                << ((ir <= Nr) ? "Fluid" : "Buffer")
                << "  r: " << rp
                << "  theta: " << thetap
                << "  x: " << xp
                << "  y: " << yp
                << endl;

        i++;

        }

    }

    int Nparticles = i;
    int NfluidEnd = (Nr + 1)*NTheta;
    int Nfluid = NfluidEnd - Nboundary;
    int NbufferParticles = Nparticles - NfluidEnd;

    cout << endl;
    cout << "Ntheta     = " << NTheta << endl;
    cout << "Nboundary  = " << Nboundary << endl;
    cout << "Nfluid     = " << Nfluid << endl;
    cout << "Nbufferparticles  = " << NbufferParticles << endl;
    cout << "Nparticles = " << Nparticles << endl;


    vector<double> hlist = {2};
    int Nh = hlist.size();

    for (int m = 0; m < Nh; m++)
    {
        
        
        double coefh = hlist[m];

        double L2PGauss = 0.0;
        double L2PCubic = 0.0;
        double L2PWedn = 0.0;

        double L2drhodtGauss = 0.0;
        double L2drhodtCubic = 0.0;
        double L2drhodtWedn = 0.0;

        double L2normPGauss = 0.0;
        double L2normPCubic = 0.0;
        double L2normPWedn = 0.0;

        double L2normdrhodtGauss = 0.0;
        double L2normdrhodtCubic = 0.0;
        double L2normdrhodtWedn = 0.0;

        


        u.resize(Nparticles);
        v.resize(Nparticles);
        rho.resize(Nparticles);

        PGauss.resize(Nparticles);
        PCubic.resize(Nparticles);
        PWedn.resize(Nparticles);

        dPGaussX.resize(Nparticles);
        dPGaussY.resize(Nparticles);

        dPCubicX.resize(Nparticles);
        dPCubicY.resize(Nparticles);

        dPWednX.resize(Nparticles);
        dPWednY.resize(Nparticles);

        drhodtGauss.resize(Nparticles);
        drhodtCubic.resize(Nparticles);
        drhodtWedn.resize(Nparticles);

        mass.resize(Nparticles);
        Vi.resize(Nparticles);

        hp.resize(Nparticles);

        pressure.resize(Nparticles);

        drhodt.resize(Nparticles);
        dudt.resize(Nparticles);
        dvdt.resize(Nparticles);

        drhodtexact.resize(Nparticles);
        pressureexact.resize(Nparticles);

        xhalf.resize(Nparticles);
        yhalf.resize(Nparticles);
        rhohalf.resize(Nparticles);
        uhalf.resize(Nparticles);
        vhalf.resize(Nparticles);
        drhodthalf.resize(Nparticles);
        pressurehalf.resize(Nparticles);
        dudthalf.resize(Nparticles);
        dvdthalf.resize(Nparticles);

        xnew.resize(Nparticles);
        ynew.resize(Nparticles);

        unew.resize(Nparticles);
        vnew.resize(Nparticles);

        rhonew.resize(Nparticles);

        xghost.resize(Nboundary);
        yghost.resize(Nboundary);

        rhoghost.resize(Nboundary);
        drhoghostX.resize(Nboundary);
        drhoghostY.resize(Nboundary);

        vorticity.resize(Nparticles);



// --------------------------------------------------
// Intial conditions
// 
// --------------------------------------------------
        //Slip Initial Conditions
        for (int i = 0; i < Nparticles; ++i)
        {
            double ri = r[i];
            double ct = cos(theta[i]);
            double st = sin(theta[i]);
            double a2_over_r2 = (r0*r0)/(ri*ri);

            double ur = inletVelocity*(1.0 - a2_over_r2)*ct;
            double ut = -inletVelocity*(1.0 + a2_over_r2)*st;

            u[i] = ur*ct - ut*st;
            v[i] = ur*st + ut*ct;

            double pInitial = 0.5*rho0*
                (inletVelocity*inletVelocity - u[i]*u[i] - v[i]*v[i]);

            rho[i] = rho0*pow(1.0 + pInitial/B, 1.0/gammaEOS);

            drhodtexact[i] = -rho[i]*velcoefX;

            pressureexact[i] = pInitial;
        }

        //No-slip Initial Conditions
        /*for (int i = 0; i < Nparticles; i++)
        {
            rho[i] = rho0;

            if (i < Nboundary)
            {
                // Cylinder
                u[i] = 0.0;
                v[i] = 0.0;
            }
            else
            {
                // Physical fluid + buffer
                u[i] = inletVelocity;
                v[i] = 0.0;
            }

            drhodtexact[i] =
                -rho[i]*velcoefX;

            pressureexact[i] =
                rho0*g*(waterheight-y[i]);
        }*/



        for (int i = 0; i < Nboundary; ++i)
        {
            // Corresponding particle on first fluid ring
            int j = i + NTheta;

            // Ghost point mirrored across the wall interface
            double dg = r[j] - r[i];

            double nx = cos(theta[i]);
            double ny = sin(theta[i]);

            xghost[i] = x[i] + dg*nx;
            yghost[i] = y[i] + dg*ny;
        }

        for (int i = 0; i < Nboundary; ++i)
        {
            applyCylinderFreeSlip(i, NTheta, theta, u, v);
        }
            



// ----------------------------------------------
// Calculate Vi, mass and h for ALL particles
// ----------------------------------------------

        for (int i = 0; i < Nparticles; i++)
        {
            int ir = i/NTheta;

            double dri;
            double rin, rout;

            if (ir == 0)
            {
                dri =
                    0.5*(r[i + NTheta] - r[i]);

                rin = r0;
                rout =
                    0.5*(r[i] + r[i + NTheta]);
            }
            else if (ir == NrTotal)
            {
                dri =
                    0.5*(r[i] - r[i - NTheta]);

                rin =
                    0.5*(r[i] + r[i - NTheta]);

                rout =
                    r[i] + 0.5*(r[i] - r[i - NTheta]);
            }
            else
            {
                dri =
                    0.5*(r[i + NTheta]
                    - r[i - NTheta]);

                rin =
                    0.5*(r[i] + r[i - NTheta]);

                rout =
                    0.5*(r[i] + r[i + NTheta]);
            }

            Vi[i] =
                0.5*(rout*rout-rin*rin)*dTheta0;

            mass[i] =
                rho[i]*Vi[i];

            double drthetai =
                r[i]*dTheta0;

            double dpi =
                sqrt(dri*drthetai);

            hp[i] =
                coefh*dpi;
        }


// ----------------------------------------------
// kernel consistency check
// ----------------------------------------------

        for (int i = Nboundary; i < NfluidEnd; i++)
        {
            PGauss[i] = 0.0;
            PCubic[i] = 0.0;
            PWedn[i] = 0.0;

            dPGaussX[i] = 0.0;
            dPGaussY[i] = 0.0; 

            dPCubicX[i] = 0.0;
            dPCubicY[i] = 0.0;

            dPWednX[i] = 0.0;
            dPWednY[i] = 0.0;

            drhodtGauss[i] = 0.0;
            drhodtCubic[i] = 0.0;
            drhodtWedn[i] = 0.0;
          

            for (int j = 0 ; j < Nparticles; j++)

            {

                double Vj = Vi[j];
                double mj = mass[j];

                double dx = x[i]-x[j];
                double dy = y[i]-y[j]; 

                double r2 = dx*dx+dy*dy;

                double s = sqrt(r2);
                double q = s/hp[i];              

                double du = u[i]-u[j];
                double dv = v[i]-v[j];

                double dirX = 0.0;
                double dirY = 0.0;


                if (r2>0.0)
                {
                    dirX = dx/s;
                    dirY = dy/s;
                }

                KernelResult result;

                result = gaussian(q, hp[i], dirX, dirY);
                PGauss[i] += result.Weight*Vj;
                dPGaussX[i] += result.dWeightX*Vj;
                dPGaussY[i] += result.dWeightY*Vj;
                drhodtGauss[i] += mj*(du*result.dWeightX+dv*result.dWeightY);
                
                result = cubicSpline(q, hp[i], dirX, dirY);
                PCubic[i] += result.Weight*Vj;
                dPCubicX[i] += result.dWeightX*Vj;
                dPCubicY[i] += result.dWeightY*Vj;
                drhodtCubic[i] += mj*(du*result.dWeightX+dv*result.dWeightY);

                result = Wendland(q, hp[i], dirX, dirY);
                PWedn[i] += result.Weight*Vj;
                dPWednX[i] += result.dWeightX*Vj;
                dPWednY[i] += result.dWeightY*Vj;
                drhodtWedn[i] += mj*(du*result.dWeightX+dv*result.dWeightY);

        
            }

            L2PGauss += (PGauss[i]-1.0)*(PGauss[i]-1.0);
            L2PCubic += (PCubic[i]-1.0)*(PCubic[i]-1.0);
            L2PWedn += (PWedn[i]-1.0)*(PWedn[i]-1.0); 
            L2drhodtGauss += (drhodtGauss[i]-drhodtexact[i])*(drhodtGauss[i]-drhodtexact[i]);
            L2drhodtCubic += (drhodtCubic[i]-drhodtexact[i])*(drhodtCubic[i]-drhodtexact[i]);
            L2drhodtWedn += (drhodtWedn[i]-drhodtexact[i])*(drhodtWedn[i]-drhodtexact[i]);
            
        }
        L2normPGauss = sqrt(L2PGauss/(Nfluid));
        L2normPCubic = sqrt(L2PCubic/(Nfluid));
        L2normPWedn = sqrt(L2PWedn/(Nfluid));

        L2normdrhodtGauss = sqrt(L2drhodtGauss/(Nfluid));
        L2normdrhodtCubic = sqrt(L2drhodtCubic/(Nfluid));
        L2normdrhodtWedn = sqrt(L2drhodtWedn/(Nfluid)); 
       
        

        

        // ------------------------------------------------------------
        // Save particle positions
        // ------------------------------------------------------------

        string filename =
            foldername  +  "/particles" + "_hdp_"+ to_string(coefh)+".csv";

        ofstream file(filename);


        // ------------------------------------------------------------
        // Simulation / particle information
        // ------------------------------------------------------------

        file << "dr0" << "," << dr0 << endl;
        file << "r0" << "," << r0 << endl;
        file << "R" << "," << R << endl;
        file << "R/r0" << "," << R/r0 << endl;
        file << "Nr" << "," << Nr << endl;
        file << "Ntheta" << "," << NTheta << endl;
        file << "Nboundary" << "," << Nboundary << endl;
        file << "Nfluid" << "," << Nfluid << endl;
        file << "Nparticles" << "," << Nparticles << endl;
        file << "Apolar" << "," << Apolar << endl;
        file << "Bpolar" << "," << Bpolar << endl;

        file << endl;

        file << "L2normPGauss" << "," << "L2normPCubic" << "," << "L2normPWedn" << "," << "L2normdrhodtGauss" << "," << "L2normdrhodtCubic" << "," << "L2normdrhodtWedn" << endl;
        file << L2normPGauss << "," << L2normPCubic << "," << L2normPWedn << "," << L2normdrhodtGauss << "," << L2normdrhodtCubic << "," << L2normdrhodtWedn << endl;
        file << endl;

        


        // ------------------------------------------------------------
        // Particle data
        // ------------------------------------------------------------

        file << "ID,r,theta,x,y,Type,,PGauss,PCubic,PWed,,dPGaussX,dPGaussY,,dPCubicX,dPCubicY,,dPWednX,dPWednY,,drhodtGauss,drhodtCubic,drhodtWedn" << endl;

        for (int i = 0; i < Nparticles; ++i)
        {
            string Type;

            if (i < Nboundary)
            {
                Type = "boundary";
            }
            else if (i < NfluidEnd)
            {
                Type = "fluid";
            }
            else
            {
                Type = "buffer";
            }

            file << i << ","
                << r[i] << ","
                << theta[i] << ","
                << x[i] << ","
                << y[i] << ","
                << Type;
            if (i >= Nboundary && i < NfluidEnd)
            {
                file << ",,"
                    << PGauss[i] << ","
                    << PCubic[i] << ","
                    << PWedn[i] << ",,"
                    << dPGaussX[i] << ","
                    << dPGaussY[i] << ",,"
                    << dPCubicX[i] << ","
                    << dPCubicY[i] << ",,"
                    << dPWednX[i] << ","
                    << dPWednY[i] << ",,"
                    << drhodtGauss[i] << ","
                    << drhodtCubic[i] << ","
                    << drhodtWedn[i];
            }
            else
            {
                file << ",,,,,,,,,,,,,,,,";
            }

            file << endl;
        }
        file.close();

        cout << endl;
        cout << "Particle file saved to: "
            << filename
            << endl;

        validatePotentialFlowForce(
            theta,
            Nboundary,
            r0,
            dTheta0,
            rho0,
            inletVelocity,
            Dcylinder,
            foldername);

// -----------------------------------------------------------------------------------------------------------------------------
// start time loop
// -----------------------------------------------------------------------------------------------------------------------------

        double t = 0.0;
        int n = 0;
        double nextOutputTime = 0.0;

        ofstream forceFile(
            foldername
            + "/CylinderForces_hdp_"
            + to_string(coefh)
            + ".csv");

        forceFile
            << "t,"
            << "FxPressure,FyPressure,"
            << "FxViscous,FyViscous,"
            << "Fx,Fy,"
            << "CD,CL"
            << endl;

        while (t < Totaltime)
        {


// -----------------------------------------------------------------------------------------------------------------------------
// mDBC Interpolation at n
// -----------------------------------------------------------------------------------------------------------------------------

            for (int i = 0; i < Nboundary; ++i)
            {
                double rhoBoundary;

                interpolateMDBC(
                    i,
                    Nboundary,
                    NfluidEnd,
                    NTheta,

                    x,
                    y,
                    rho,

                    xghost,
                    yghost,

                    Vi,
                    hp,

                    kernel,

                    rhoghost[i],
                    drhoghostX[i],
                    drhoghostY[i],
                    rhoBoundary);

                rho[i] = rhoBoundary;

                applyCylinderFreeSlip(i, NTheta, theta, u, v);

                // Stationary no-slip cylinder
                //u[i] = 0.0;
                //v[i] = 0.0;
            }

            applyRiemannBufferBoundary(
                u,
                v,
                rho,
                theta,
                NTheta,
                Nr,
                NrTotal);


// -----------------------------------------------------------------------------------------------------------------------------
// Compute Pressure and Force at n
// -----------------------------------------------------------------------------------------------------------------------------        
            
            for (int i = 0; i < Nparticles; i++)
            {
                pressure[i] = B*(pow(rho[i]/rho0,gammaEOS)-1.0);                         
            }


            double FxPressure = 0.0;
            double FyPressure = 0.0;

            double FxViscous = 0.0;
            double FyViscous = 0.0;

            // Dynamic viscosity
            double mu = rho0*viscosity;

            for (int i = 0; i < Nboundary; ++i)
            {
                double nx = cos(theta[i]);
                double ny = sin(theta[i]);

                double tx = -ny;
                double ty =  nx;

                // Arc length represented by boundary particle
                double ds = r0*dTheta0;

                // -----------------------------------------
                // Pressure force
                // -----------------------------------------

                FxPressure +=
                    -pressure[i]*nx*ds;

                FyPressure +=
                    -pressure[i]*ny*ds;


                // -----------------------------------------
                // Approximate viscous wall shear
                // using first physical-fluid ring
                // -----------------------------------------

                int j = i + NTheta;

                double dr =
                    r[j] - r[i];

                double utFluid =
                    u[j]*tx
                    + v[j]*ty;

                // Stationary cylinder: ut_wall = 0
                /*double tauWall =
                    mu*utFluid/dr;

                FxViscous +=
                    tauWall*tx*ds;

                FyViscous +=
                    tauWall*ty*ds;*/
            }

            double Fx = FxPressure + FxViscous;
            double Fy = FyPressure + FyViscous;
            double forceReference = 0.5 * rho0 * inletVelocity * inletVelocity * Dcylinder;
            double CD = Fx/forceReference;
            double CL = Fy/forceReference;

            forceFile
                << t << ","
                << FxPressure << ","
                << FyPressure << ","
                << FxViscous << ","
                << FyViscous << ","
                << Fx << ","
                << Fy << ","
                << CD << ","
                << CL
                << endl;


// -----------------------------------------------------------------------------------------------------------------------------
// Compute Kinetic Energy
// -----------------------------------------------------------------------------------------------------------------------------
            double KE = 0.0;
            for (int i = Nboundary; i < NfluidEnd; i++)
            {
                KE += 0.5 * mass[i]* ((u[i]*u[i])+(v[i]*v[i]));
            }

// -----------------------------------------------------------------------------------------------------------------------------
// Compute Vorticity
// -----------------------------------------------------------------------------------------------------------------------------

            for (int i = Nboundary;
                i < NfluidEnd;
                ++i)
            {
                vorticity[i] = 0.0;

                for (int j = 0;
                    j < Nparticles;
                    ++j)
                {
                    double rx = x[i]-x[j];
                    double ry = y[i]-y[j];

                    double r2 =
                        rx*rx + ry*ry;

                    if (r2 < 1.0e-14)
                    {
                        continue;
                    }

                    if ((kernel == "cubic" ||
                        kernel == "wendland") &&
                        r2 > 4.0*hp[i]*hp[i])
                    {
                        continue;
                    }

                    double distance =
                        sqrt(r2);

                    double q =
                        distance/hp[i];

                    if (kernel == "gaussian" &&
                        q > 3.0)
                    {
                        continue;
                    }

                    double dirX =
                        rx/distance;

                    double dirY =
                        ry/distance;

                    KernelResult result;

                    if (kernel == "gaussian")
                        result =
                            gaussian(q,hp[i],dirX,dirY);

                    else if (kernel == "cubic")
                        result =
                            cubicSpline(q,hp[i],dirX,dirY);

                    else
                        result =
                            Wendland(q,hp[i],dirX,dirY);

                    double du =
                        u[j]-u[i];

                    double dv =
                        v[j]-v[i];

                    vorticity[i] +=
                        Vi[j]*
                        (
                            dv*result.dWeightX
                            -
                            du*result.dWeightY
                        );
                }
            }


// -----------------------------------------------------------------------------------------------------------------------------
// Fluid Continuity and Momentum at n
// -----------------------------------------------------------------------------------------------------------------------------
            //#pragma omp parallel for
            for (int i = Nboundary; i < NfluidEnd; i++)
            {
                //cummulative should be initialised with zero for each particle

                drhodt[i] = 0.0;
                dudt[i] = 0.0;
                dvdt[i] = 0.0;               

                for (int j = 0; j < Nparticles; j++)
                {
                    
                    accumulateInteraction(
                        i, j,
                        x, y, u, v, rho, pressure,
                        Vi[j], mass[j], hp[i], kernel,
                        drhodt[i], dudt[i], dvdt[i]);
                    
                }

                   
            }

            //Calculate Dt
            double dt = CFL*calculateDt(
                Nboundary,
                Nparticles,
                NfluidEnd,
                x,
                y,
                u,
                v,
                rho,
                pressure,
                dudt,
                dvdt,
                hp,
                c0,
                kernel);
            
                if (!isfinite(dt) || dt <= 0.0)
                {
                    cout << "Invalid timestep at t = " << t << endl;
                    return 1;
                }

                if (dt < dtMin)
                {
                    cout << "Timestep became too small"
                        << "  t = " << t
                        << "  dt = " << dt
                        << endl;

                    return 1;
                }


            dt = min(dt, dtMax);
            dt = min(dt, Totaltime - t);

// -----------------------------------------------------------------------------------------------------------------------------
// Compute error
// -----------------------------------------------------------------------------------------------------------------------------
            double L2drhodt = 0.0;
            double L2pressure = 0.0;           

            for (int i = Nboundary; i < NfluidEnd; i++)
            {
                //double errordrhodt = (drhodt[i]-drhodtexact[i])*(drhodt[i]-drhodtexact[i]);
                double errorpressure = (pressure[i]-pressureexact[i])*(pressure[i]-pressureexact[i]);
                
                //L2drhodt += errordrhodt;
                L2pressure += errorpressure;
            }
 
            double L2normdrhodt = sqrt(L2drhodt/(Nfluid));
            double L2normpressure = sqrt(L2pressure/Nfluid);


// -----------------------------------------------------------------------------------------------------------------------------
// Write output
// -----------------------------------------------------------------------------------------------------------------------------            
            if (t >= nextOutputTime - 1e-12)
            {
                cout << "t : " << t << endl;
                string filename = foldername + "/EWCSPH_hdp_" + to_string(coefh) + "_t_" + to_string(t) + ".csv";

                ofstream file(filename);
                file << "h/dp :" << "," << coefh << endl;
                file << "L2norm Pressure" << "," << "KE" << "," << "dr0" << ","  << "Nfluid" << endl;
                file << L2normpressure << "," << KE << "," << dr0 << "," << Nfluid <<endl;
                file << endl;
                file << "t" << "," << t << endl;
                file << "ID"<< ","<< "x"<< "," << "y" << ","<< ","<<"rho"<< "," << "drhodt" << "," << "pressure"<< "," << "," << "u" << "," << "v" << "," << "velocity" << "," << "," << "dudt" << "," << "dvdt" <<","<< "vorticity" <<"," << "type" << endl;
                for (int i = 0; i < Nparticles; i++)
                    {
                        if (i < Nboundary)
                        {
                            type = "boundary";
                        }
                        else if (i < NfluidEnd)
                        {
                            type = "fluid";
                        }
                        else
                        {
                            type = "buffer";
                        }
                        file << i << "," << x[i] << "," << y[i] << "," << "," << rho[i] << "," << drhodt[i] << "," << pressure[i] << "," << "," << u[i] << "," << v[i] << "," << sqrt(u[i]*u[i] + v[i]*v[i]) << "," << "," << dudt[i] << "," << dvdt[i] << ","<< ((i >= Nboundary && i < NfluidEnd)? vorticity[i]: 0.0) << "," << type << endl;
                    }
           
                file.close();
                
                nextOutputTime += 0.01;
            }

// -----------------------------------------------------------------------------------------------------------------------------
// Compute Time integration predictor for fluid to n+1/2
// -----------------------------------------------------------------------------------------------------------------------------

            for (int i = Nboundary;i < NfluidEnd;i++)
            {
                rhohalf[i] =
                    rho[i] + drhodt[i]*dt/2.0;

                uhalf[i] =
                    u[i] + dudt[i]*dt/2.0;

                vhalf[i] =
                    v[i] + dvdt[i]*dt/2.0;

                xhalf[i] = x[i];
                yhalf[i] = y[i];
            }

            for (int i = NfluidEnd;i < Nparticles;++i)
            {
                xhalf[i] = x[i];
                yhalf[i] = y[i];
            }

// -----------------------------------------------------------------------------------------------------------------------------
// Compute Time integration predictor for boundary to n+1/2
// -----------------------------------------------------------------------------------------------------------------------------

            for (int i = 0; i < Nboundary; i++)
            {

                xhalf[i] = x[i];
                yhalf[i] = y[i];


                double rhoBoundary;

                interpolateMDBC(
                    i,
                    Nboundary,
                    NfluidEnd,
                    NTheta,

                    xhalf,
                    yhalf,
                    rhohalf,

                    xghost,
                    yghost,

                    Vi,
                    hp,

                    kernel,

                    rhoghost[i],
                    drhoghostX[i],
                    drhoghostY[i],
                    rhoBoundary);

                rhohalf[i] = rhoBoundary;

                applyCylinderFreeSlip(i, NTheta, theta, uhalf, vhalf);

                // No-slip Cylinder
                //uhalf[i] = 0.0;
                //vhalf[i] = 0.0;
            }

            applyRiemannBufferBoundary(
                uhalf,
                vhalf,
                rhohalf,
                theta,
                NTheta,
                Nr,
                NrTotal);


// -----------------------------------------------------------------------------------------------------------------------------
// Compute Pressure at n+1/2
// -----------------------------------------------------------------------------------------------------------------------------        
            
            for (int i = 0; i < Nparticles; i++)
            {
                pressurehalf[i] = B*(pow(rhohalf[i]/rho0,gammaEOS)-1.0);                         
            }

// -----------------------------------------------------------------------------------------------------------------------------
// Fluid Continuity and Momentum at n+1/2
// -----------------------------------------------------------------------------------------------------------------------------
            //#pragma omp parallel for
            for (int i = Nboundary; i < NfluidEnd; i++)
            {
                //cummulative should be initialised with zero for each particle

                drhodthalf[i] = 0.0;
                dudthalf[i] = 0.0;
                dvdthalf[i] = 0.0;               

                for (int j = 0; j < Nparticles; j++)
                {
                    
                    accumulateInteraction(
                        i, j,
                        xhalf, yhalf,
                        uhalf, vhalf,
                        rhohalf, pressurehalf,
                        Vi[j], mass[j], hp[i], kernel,
                        drhodthalf[i], dudthalf[i], dvdthalf[i]);
                    
                }
            }


// -----------------------------------------------------------------------------------------------------------------------------
// Compute Time integration corrector for fluid to n+1/2
// -----------------------------------------------------------------------------------------------------------------------------

            for (int i = Nboundary; i < NfluidEnd; i++)
            {
                rhohalf[i] = rho[i]+drhodthalf[i]*dt/2;
                uhalf[i] = u[i]+dudthalf[i]*dt/2;
                vhalf[i] = v[i]+(dvdthalf[i])*dt/2;
                xhalf[i] = x[i];
                yhalf[i] = y[i];
            }

            for (int i = NfluidEnd;i < Nparticles;++i)
            {
                xhalf[i] = x[i];
                yhalf[i] = y[i];
            }



// -----------------------------------------------------------------------------------------------------------------------------
// Compute Time integration corrector for boundary to n+1/2
// -----------------------------------------------------------------------------------------------------------------------------

            for (int i = 0; i < Nboundary; ++i)
            {
                xhalf[i] = x[i];
                yhalf[i] = y[i];

                double rhoBoundary;

                interpolateMDBC(
                    i,
                    Nboundary,
                    NfluidEnd,
                    NTheta,

                    xhalf,
                    yhalf,
                    rhohalf,

                    xghost,
                    yghost,

                    Vi,
                    hp,

                    kernel,

                    rhoghost[i],
                    drhoghostX[i],
                    drhoghostY[i],
                    rhoBoundary);

                rhohalf[i] = rhoBoundary;

                applyCylinderFreeSlip(i, NTheta, theta, uhalf, vhalf);

                // No-slip Cylinder
                //uhalf[i] = 0.0;
                //vhalf[i] = 0.0;
            }

            applyRiemannBufferBoundary(
                uhalf,
                vhalf,
                rhohalf,
                theta,
                NTheta,
                Nr,
                NrTotal);


// -----------------------------------------------------------------------------------------------------------------------------
// Compute Time integration to n+1
// -----------------------------------------------------------------------------------------------------------------------------



// For Fluid          
            for (int i = Nboundary; i < NfluidEnd; i++)
            {
                rhonew[i] = 2*rhohalf[i]-rho[i];
                unew[i] = 2*uhalf[i]-u[i];
                vnew[i] = 2*vhalf[i]-v[i];
                xnew[i] = xhalf[i];
                ynew[i] = yhalf[i];
            }

            for (int i = NfluidEnd;i < Nparticles; ++i)
            {
                xnew[i] = x[i];
                ynew[i] = y[i];
            }


// For Boundary
            for (int i = 0; i < Nboundary; i++)
            {
                
                // No-slip Cylinder
                //unew[i] = 0.0;
                //vnew[i] = 0.0;

                applyCylinderFreeSlip(i, NTheta, theta, unew, vnew);
                xnew[i] = xhalf[i];
                ynew[i] = yhalf[i];

                double rhoBoundary;

                interpolateMDBC(
                    i,
                    Nboundary,
                    NfluidEnd,
                    NTheta,

                    xnew,
                    ynew,
                    rhonew,

                    xghost,
                    yghost,

                    Vi,
                    hp,

                    kernel,

                    rhoghost[i],
                    drhoghostX[i],
                    drhoghostY[i],
                    rhoBoundary);

                rhonew[i] = rhoBoundary;


            }

            applyRiemannBufferBoundary(
                unew,
                vnew,
                rhonew,
                theta,
                NTheta,
                Nr,
                NrTotal);




           
// -----------------------------------------------------------------------------------------------------------------------------
// next loop for all particles
// -----------------------------------------------------------------------------------------------------------------------------
            
            for (int i = 0; i < Nparticles; i++)
            {
                x[i] = xnew[i];
                y[i] = ynew[i];
                u[i] = unew[i];
                v[i] = vnew[i];
                rho[i] = rhonew[i];

                
                
            }                
                   
 
        t += dt;
        n++;
        
        }
        
        forceFile.close();


    }
    

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    cout << endl;
    cout << "==============================" << endl;
    cout << "Total simulation time = "
        << elapsed.count()
        << " seconds" << endl;
    cout << "==============================" << endl;

}
