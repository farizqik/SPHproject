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
const double Totaltime = 20.0;

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

const double inletVelocity = 0.5;  // m/s

const double c0 = 10.0*sqrt(g*(waterheight));

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

const int Nr = 31;

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
    int Nparticles,
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

        for (int i = Nboundary; i < Nparticles; ++i)
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
    
    
    for (int ir = 1; ir <= Nr; ++ir)
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
                << "  Type: Fluid"
                << "  r: " << rp
                << "  theta: " << thetap
                << "  x: " << xp
                << "  y: " << yp
                << endl;

        i++;

        }

    }

    int Nparticles = i;
    int Nfluid = Nparticles-Nboundary;

    cout << endl;
    cout << "Ntheta     = " << NTheta << endl;
    cout << "Nboundary  = " << Nboundary << endl;
    cout << "Nfluid     = " << Nfluid << endl;
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



// --------------------------------------------------
// Intial conditions
// 
// --------------------------------------------------
        for (int i = 0; i < Nparticles; i++)
        {
            u[i] = 0.0;
            v[i] = 0.0;
            rho[i] = rho0;
            if (i / NTheta == Nr && x[i] < 0) 
            {
                u[i] = inletVelocity;
                v[i] = 0.0;
            }
            drhodtexact[i] = -rho[i]*velcoefX;
            pressureexact[i] = rho0*g*(waterheight-y[i]);
        }

// ----------------------------------------------
// Calculate Vi, mass and h for ALL particles
//
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
                rout = 0.5*(r[i] + r[i + NTheta]);
            }
            else if (ir == Nr)
            {
                dri =
                    0.5*(r[i] - r[i - NTheta]);

                rin = 0.5*(r[i] + r[i - NTheta]);
                rout = R;
            }
            else
            {
                dri =
                    0.5*(r[i + NTheta]
                    - r[i - NTheta]);

                rin = 0.5*(r[i] + r[i - NTheta]);
                rout = 0.5*(r[i] + r[i + NTheta]);
            }


            Vi[i] = 0.5*(rout*rout - rin*rin)*dTheta0;

            mass[i] = rho[i]*Vi[i];


            double drthetai = r[i]*dTheta0;

            double dpi = sqrt(dri*drthetai);

            hp[i] = coefh*dpi;

        }


// ----------------------------------------------
// kernel consistency check
// ----------------------------------------------

        for (int i = 0; i < Nparticles; i++)
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
        L2normPGauss = sqrt(L2PGauss/(Nparticles));
        L2normPCubic = sqrt(L2PCubic/(Nparticles));
        L2normPWedn = sqrt(L2PWedn/(Nparticles));

        L2normdrhodtGauss = sqrt(L2drhodtGauss/(Nparticles));
        L2normdrhodtCubic = sqrt(L2drhodtCubic/(Nparticles));
        L2normdrhodtWedn = sqrt(L2drhodtWedn/(Nparticles)); 
       
        

        

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
            else
            {
                Type = "fluid";
            }

            file << i << ","
                << r[i] << ","
                << theta[i] << ","
                << x[i] << ","
                << y[i] << ","
                << Type << ",,"
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
                << drhodtWedn[i] << ","
                << endl;
        }

        file.close();

        cout << endl;
        cout << "Particle file saved to: "
            << filename
            << endl;

        // -----------------------------------------------------------------------------------------------------------------------------
        // start time loop
        // -----------------------------------------------------------------------------------------------------------------------------

        double t = 0.0;
        int n = 0;
        double nextOutputTime = 0.0;

        while (t < Totaltime)
        {

// -----------------------------------------------------------------------------------------------------------------------------
// Compute Pressure at n
// -----------------------------------------------------------------------------------------------------------------------------        
            
            for (int i = 0; i < Nparticles; i++)
            {
                pressure[i] = B*(pow(rho[i]/rho0,gammaEOS)-1.0);                         
            }
// -----------------------------------------------------------------------------------------------------------------------------
// Compute Kinetic Energy
// -----------------------------------------------------------------------------------------------------------------------------
            double KE = 0.0;
            for (int i = Nboundary; i < Nparticles; i++)
            {
                KE += 0.5 * mass[i]* ((u[i]*u[i])+(v[i]*v[i]));
            }
// -----------------------------------------------------------------------------------------------------------------------------
// Fluid Continuity and Momentum at n
// -----------------------------------------------------------------------------------------------------------------------------
            //#pragma omp parallel for
            for (int i = 0; i < Nparticles; i++)
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

            for (int i = Nboundary; i < Nparticles; i++)
            {
                double errordrhodt = (drhodt[i]-drhodtexact[i])*(drhodt[i]-drhodtexact[i]);
                double errorpressure = (pressure[i]-pressureexact[i])*(pressure[i]-pressureexact[i]);
                
                L2drhodt += errordrhodt;
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
                file << "L2norm Pressure" << "," << "KE" << "," << "dr0" << ","  << "Nparticles" << endl;
                file << L2normpressure << "," << KE << "," << dr0 << "," << Nfluid <<endl;
                file << endl;
                file << "t" << "," << t << endl;
                file << "ID"<< ","<< "x"<< "," << "y" << ","<< ","<<"rho"<< "," << "drhodt" << "," << "pressure"<< "," << "," << "u" << "," << "v" << "," << "," << "dudt" << "," << "dvdt" <<"," << "type" << endl;
                for (int i = 0; i < Nparticles; i++)
                    {
                        if (i < Nboundary)
                        {
                            type = "boundary";
                        }
                        else
                        {
                            type = "fluid";
                        }
                        file << i << "," << x[i] << "," << y[i] << "," << "," << rho[i] << "," << drhodt[i] << "," << pressure[i] << "," << "," << u[i] << "," << v[i] << "," << "," << dudt[i] << "," << dvdt[i] << "," << type << endl;
                    }
           
                file.close();
                
                nextOutputTime += 0.01;
            }

// -----------------------------------------------------------------------------------------------------------------------------
// Compute Time integration predictor for fluid to n+1/2
// -----------------------------------------------------------------------------------------------------------------------------

            for (int i = Nboundary; i < Nparticles; i++)
            {
                rhohalf[i] = rho[i]+drhodt[i]*dt/2;
                uhalf[i] = u[i]+dudt[i]*dt/2;
                vhalf[i] = v[i]+(dvdt[i])*dt/2;
                if (i / NTheta == Nr && x[i] < 0)
                {
                    uhalf[i] = u[i];
                    vhalf[i] = v[i];
                }
                xhalf[i] = x[i];
                yhalf[i] = y[i];
            }

// -----------------------------------------------------------------------------------------------------------------------------
// Compute Time integration predictor for boundary to n+1/2
// -----------------------------------------------------------------------------------------------------------------------------

            for (int i = 0; i < Nboundary; i++)
            {
                rhohalf[i] = rho[i]+drhodt[i]*dt/2;
                uhalf[i] = 0.0;
                vhalf[i] = 0.0;
                xhalf[i] = x[i];
                yhalf[i] = y[i];
            }


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
            for (int i = 0; i < Nparticles; i++)
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

            for (int i = Nboundary; i < Nparticles; i++)
            {
                rhohalf[i] = rho[i]+drhodthalf[i]*dt/2;
                uhalf[i] = u[i]+dudthalf[i]*dt/2;
                vhalf[i] = v[i]+(dvdthalf[i])*dt/2;
                if (i / NTheta == Nr && x[i] < 0)
                {
                    uhalf[i] = u[i];
                    vhalf[i] = v[i];
                }
                xhalf[i] = x[i];
                yhalf[i] = y[i];
            }

// -----------------------------------------------------------------------------------------------------------------------------
// Compute Time integration corrector for boundary to n+1/2
// -----------------------------------------------------------------------------------------------------------------------------

            for (int i = 0; i < Nboundary; i++)
            {
                rhohalf[i] = rho[i]+drhodthalf[i]*dt/2;
                uhalf[i] = 0.0;
                vhalf[i] = 0.0;
                xhalf[i] = x[i];
                yhalf[i] = y[i];
            }


// -----------------------------------------------------------------------------------------------------------------------------
// Compute Time integration to n+1
// -----------------------------------------------------------------------------------------------------------------------------

// For Boundary
            for (int i = 0; i < Nboundary; i++)
            {
                rhonew[i] = 2*rhohalf[i]-rho[i];
                unew[i] = 0.0;
                vnew[i] = 0.0;
                xnew[i] = xhalf[i];
                ynew[i] = yhalf[i];
            }

// For Fluid          
            for (int i = Nboundary; i < Nparticles; i++)
            {
                rhonew[i] = 2*rhohalf[i]-rho[i];
                unew[i] = 2*uhalf[i]-u[i];
                vnew[i] = 2*vhalf[i]-v[i];
                xnew[i] = xhalf[i];
                ynew[i] = yhalf[i];
            }


           
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
