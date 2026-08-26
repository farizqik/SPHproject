#include<iostream>
#include<cmath>
#include <string>
#include <iomanip>
#include <fstream>
#include <chrono>
#include <filesystem>
#include <sstream>
#include <vector>


using namespace std;

string type; 

// --------------------------------------------
// Timestep
// --------------------------------------------

const double dt = 0.0001; 
const double Totaltime= 20.0;
const int Nt = Totaltime / dt;


// ------------------------------------------------------------
// Geometry
// ------------------------------------------------------------

const double tanklength = 25.0;
const double tankheight = 10.0;

const double waterlength = 15.0;
const double freeboard = 2.0;
const double waterheight = tankheight-freeboard;

const double dp = 0.5;

const double boundthick = dp * 3;

// ------------------------------------------------------------
// Particle properties
// ------------------------------------------------------------

// parameters         

vector<double> x;
vector<double> y;

vector<double> xnew;
vector<double> ynew;

vector<double> u;
vector<double> v;

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


const double PI = 3.14159265358979323846;
const double g = 9.81;
const double rho0 = 1000.0;
const int boundpart = 0; 
const double c0 = 10.0*sqrt(g*(waterheight));
const double mass = rho0*dp*dp;

const double gamma = 7.0;
const double B = c0*c0*rho0/gamma;

const double alphaAV = 0.01;


//set deltadifussion = 0.0 if it's not considered
const double deltadifussion = 0.1;





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
KernelResult gaussian(double q, double h, double dirx, double diry)
{
    KernelResult result;
    double alpha = 1.0 / (PI * h * h);

    result.Weight = alpha*exp(-q*q);
    result.dWeightX = -2.0*alpha*q*exp(-q*q)/h * dirx;
    result.dWeightY = -2.0*alpha*q*exp(-q*q)/h * diry;

    return result;
}


// -----------------------------------------------------------------------------------------------------------------------------
// Cubic kernel
// -----------------------------------------------------------------------------------------------------------------------------
KernelResult cubicSpline(double q, double h, double dirx, double diry)
{
    double alpha = 10.0 / (7.0 * PI * h * h);
    KernelResult result;
    if (q<=1.0)
    {
        result.Weight = alpha*(1.0-1.5*pow(q,2)+0.75*pow(q,3));
        result.dWeightX = alpha*(-3.0*q+2.25*pow(q,2))/h * dirx;
        result.dWeightY = alpha*(-3.0*q+2.25*pow(q,2))/h * diry;
    }
    else if (q<=2.0)
    {
        result.Weight = alpha*0.25*pow(2.0-q,3);
        result.dWeightX = -0.75*alpha*pow(2.0-q,2)/h * dirx;
        result.dWeightY = -0.75*alpha*pow(2.0-q,2)/h * diry;
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
KernelResult Wendland(double q, double h, double dirx, double diry)
{
    double alpha = 7.0 / (4.0 * PI * h * h);
    KernelResult result;
    if (0.0<=q && q<=2.0)
    {
        result.Weight = alpha*pow(1.0-0.5*q,4)*(2.0*q+1.0);
        result.dWeightX = alpha*(4*pow(1.0-0.5*q,3)*(-0.5)*(2.0*q+1.0)+pow(1.0-0.5*q,4)*2.0)/h * dirx;
        result.dWeightY = alpha*(4*pow(1.0-0.5*q,3)*(-0.5)*(2.0*q+1.0)+pow(1.0-0.5*q,4)*2.0)/h * diry;
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
// Function to solve linear equations A*x = b
// --------------------------------------------------
bool solveLinear(double A[][3], double b[], double X[], int n)
{
    const double tolerance = 1e-12;
    // Forward elimination
    for (int k = 0; k < n - 1; k++)
    {
        if (abs(A[k][k]) < tolerance)
        {
            return false;
        }
        for (int i = k + 1; i < n; i++)
        {
            
            double factor = A[i][k] / A[k][k];

            for (int j = k; j < n; j++)
            {
                A[i][j] = A[i][j] - factor * A[k][j];
            }

            b[i] = b[i] - factor * b[k];
        }
    }

    // Back substitution
    for (int i = n - 1; i >= 0; i--)
    {
        if (abs(A[i][i]) < tolerance)
        {
            return false;
        }

        double sum = b[i];

        for (int j = i + 1; j < n; j++)
        {
            sum = sum - A[i][j] * X[j];
        }

        X[i] = sum / A[i][i];

        if (!isfinite(X[i]))
        {
            return false;
        }
    }
    return true;
}



// --------------------------------------------------
// Robust 3x3 linear solver with partial pivoting.
// Returns false when A is singular / nearly singular.
// --------------------------------------------------
bool solveLinear2(double A[][3], double b[], double X[], int n)
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




// -----------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------
// Main Program
// -----------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------

int main ()
{
    // --------------------------------------------
    // Generate particles
    // --------------------------------------------

    int i = 0;

    //Bottom Boundary Particles

    for (double yp = -0.5*dp; yp > -boundthick; yp -= dp)
    {
        for (double xp = -boundthick+0.5*dp; xp < tanklength + boundthick + 0.5*dp; xp += dp)
        {

            cout << "ID: " << i
             << "  Type: Boundary"
             << "  x: " << xp
             << "  y: " << yp
             << endl;

            i++;

        }
    }

    // Left boundary
    for (double xp = -0.5*dp; xp > -boundthick; xp -= dp)
    {
        for (double yp = 0.5*dp; yp < tankheight; yp += dp)
        {

            cout << "ID: " << i
             << "  Type: Boundary"
             << "  x: " << xp
             << "  y: " << yp
             << endl;

            i++;
        }
    }

    // Right boundary
    for (double xp = tanklength + 0.5*dp;
         xp < tanklength + boundthick;
         xp += dp)
    {
        for (double yp = 0.5*dp; yp < tankheight; yp += dp)
        {

            cout << "ID: " << i
             << "  Type: Boundary"
             << "  x: " << xp
             << "  y: " << yp
             << endl;


            i++;
            
        }
    }

    int Nboundary = i;

    for (double xp = 0.5*dp; xp < waterlength; xp += dp)
    {
        for (double yp = 0.5*dp; yp < waterheight; yp += dp)
        {

            cout << "ID: " << i
             << "  Type: Fluid"
             << "  x: " << xp
             << "  y: " << yp
             << endl;

            i++;

        }
    }

    // Total number of particles
    int Nparticles = i;

    int Nfluid = Nparticles - Nboundary;

    cout << endl;
    cout << "Boundary particles = " << Nboundary << endl;
    cout << "Fluid particles    = " << Nfluid << endl;
    cout << "Total particles    = " << Nparticles << endl;
    



    string kernel; 
    cout<<"enter the kernel to be used (gaussian, cubic, wendland): ";
    cin>>kernel;

    auto start = std::chrono::high_resolution_clock::now();


    const int Nh = 1;
    double hlist[Nh] = {2*dp};

    for (int m = 0; m < Nh; m++)
    {
        double h = hlist[m];
        cout << "h/dp : " << h/dp << endl;
        string foldername = 
            "Stillwater_dp_" + to_string(dp) + "_h_" + to_string(h) + "_Nparticles_" + to_string(Nparticles) + "_" + kernel; 
        //system(("mkdir -p " + foldername).c_str());
        std::filesystem::create_directories(foldername);

        x.resize(Nparticles);
        y.resize(Nparticles);

        xnew.resize(Nparticles);
        ynew.resize(Nparticles);

        u.resize(Nparticles);
        v.resize(Nparticles);

        unew.resize(Nparticles);
        vnew.resize(Nparticles);

        rho.resize(Nparticles);
        rhonew.resize(Nparticles);

        pressure.resize(Nparticles);
        pressurenew.resize(Nparticles);
        pressureexact.resize(Nparticles);

        drhodt.resize(Nparticles);
        drhodtexact.resize(Nparticles);

        dudt.resize(Nparticles);
        dvdt.resize(Nparticles);

        xghost.resize(Nboundary);
        yghost.resize(Nboundary);
        rhoghost.resize(Nboundary);
        drhoghostX.resize(Nboundary);
        drhoghostY.resize(Nboundary);

        xhalf.resize(Nparticles);
        yhalf.resize(Nparticles);
        rhohalf.resize(Nparticles);
        uhalf.resize(Nparticles);
        vhalf.resize(Nparticles);
        drhodthalf.resize(Nparticles);
        pressurehalf.resize(Nparticles);
        dudthalf.resize(Nparticles);
        dvdthalf.resize(Nparticles);

      

// -----------------------------------------------------------------------------------------------------------------------------
// Initial condition
// -----------------------------------------------------------------------------------------------------------------------------

    int i = 0;

    //Bottom Boundary Particles

    for (double yp = -0.5*dp; yp > -boundthick; yp -= dp)
    {
        for (double xp = -boundthick+0.5*dp; xp < tanklength + boundthick + 0.5*dp; xp += dp)
        {

            x[i] = xp;
            y[i] = yp;

            u[i] = 0.0;
            v[i] = 0.0;

            if (xp < 0)
            {
                xghost[i] = -x[i];
                yghost[i] = -y[i];
            }

            else if (xp > tanklength)
            {
                xghost[i] = 2*tanklength-x[i];
                yghost[i] = -y[i];
            }

            else
            {
                xghost[i] = x[i];
                yghost[i] = -y[i];
            }
            
            rho[i] = rho0;

            i++;

        }
    }

    // Left boundary
    for (double xp = -0.5*dp; xp > -boundthick; xp -= dp)
    {
        for (double yp = 0.5*dp; yp < tankheight; yp += dp)
        {

            x[i] = xp;
            y[i] = yp;

            xghost[i] = -x[i];
            yghost[i] = y[i];

            u[i] = 0.0;
            v[i] = 0.0;

            rho[i] = rho0;
            
            i++;
        }
    }

    // Right boundary
    for (double xp = tanklength + 0.5*dp;
         xp < tanklength + boundthick;
         xp += dp)
    {
        for (double yp = 0.5*dp; yp < tankheight; yp += dp)
        {

            x[i] = xp;
            y[i] = yp;

            xghost[i] = 2*tanklength-x[i];
            yghost[i] = y[i];

            u[i] = 0.0;
            v[i] = 0.0;

            rho[i] = rho0;
            
            i++;
            
        }
    }


    for (double xp = 0.5*dp; xp < waterlength; xp += dp)
    {
        for (double yp = 0.5*dp; yp < waterheight; yp += dp)
        {

            x[i] = xp;
            y[i] = yp;

            u[i] = 0.0;
            v[i] = 0.0;
            
            //rho[i] = rho0;

            // Exact hydrostatic initial density for the Tait EOS used here.
            // It satisfies dp/dy = -rho*g in the continuum.
            double depth = max(0.0, waterheight - y[i]);
            rho[i] = rho0 * pow(1.0 + (gamma - 1.0)*g*depth/(c0*c0), 1.0/(gamma - 1.0));
            
            i++;

        }
    }

    for (int i = 0; i < Nparticles; i++)
    {
        drhodtexact[i] = -rho[i]*(0.0);
        pressureexact[i] = rho0*g*(waterheight-y[i]);                                
    }
     




        
// -----------------------------------------------------------------------------------------------------------------------------
// start time loop
// -----------------------------------------------------------------------------------------------------------------------------

        for (int n=0;n<Nt;n++)
        {

            double t = n*dt;

            
// -----------------------------------------------------------------------------------------------------------------------------
// mDBC Interpolation at n
// -----------------------------------------------------------------------------------------------------------------------------

            for (int i = 0; i < Nboundary; i++)
            {
                //cummulative should be initialised with zero for each particle

                drhodt[i] = 0.0;
                dudt[i] = 0.0;
                dvdt[i] = 0.0;

                const int n = 3;

                double A[n][n]={0};
                double b[n]={0};
                double X[n];

                int Nneighbor = 0;

                // Shepard-filter fallback used when the mDBC correction matrix
                // does not have sufficient / reliable fluid support.
                double shepardNumerator = 0.0;
                double shepardDenominator = 0.0;

                                

                for (int j = Nboundary; j < Nparticles; j++)
                {

                    double rx = xghost[i]-x[j];
                    double ry = yghost[i]-y[j];

                    double r = sqrt(rx*rx+ry*ry);
                    double q = r/h;

                    double dirx = 0.0;
                    double diry = 0.0;


                    if (r>0.0)
                    {
                        dirx = rx/r;
                        diry = ry/r;
                    }



                    
                    KernelResult result;

                    if (kernel == "gaussian")
                    {
                        result = gaussian(q, h, dirx, diry);
                    }
                    else if (kernel == "cubic")
                    {
                        result = cubicSpline(q, h, dirx, diry);
                    }
                    else if (kernel == "wendland")
                    {
                        result = Wendland(q, h, dirx, diry);
                    }
                    else
                    {
                        cout << "Invalid kernel choice. Please choose 'gaussian', 'cubic', or 'wendland'." << endl;
                        return 1; // Exit the program with an error code
                    }

                    if (abs(result.Weight) < 1e-14)
                    {
                        continue;
                    }

                    Nneighbor++;

                                      
                    
                    double Vj = mass/rho[j];

                    shepardNumerator += rho[j] * result.Weight * Vj;
                    shepardDenominator += result.Weight * Vj;


                    double dA[n][n] = 
                    {
                        {result.Weight*Vj, result.Weight*Vj*(-rx), result.Weight*Vj*(-ry)},
                        {result.dWeightX*Vj, result.dWeightX*Vj*(-rx), result.dWeightX*Vj*(-ry)},
                        {result.dWeightY*Vj, result.dWeightY*Vj*(-rx), result.dWeightY*Vj*(-ry)}
                    };
                    
                    double db[n] = 
                    {
                        result.Weight*mass,
                        result.dWeightX*mass,
                        result.dWeightY*mass
                    };

                    

                    for (int row = 0; row < 3; row++)
                    {
                        for (int col = 0; col < 3; col++)
                        {
                            A[row][col]+=dA[row][col];
                            
                        }
                        b[row]+=db[row];
                    }


                }

                // The mDBC formulation becomes unreliable with very low ghost support.
                // Use at least four fluid neighbours and reject singular/nearly-singular solves.
                const int minMdbcNeighbors = 4;
                double drytolerance = 1e-12;
                double supporttolerance = 0.4;


                bool hasfluid = shepardDenominator > drytolerance;
                bool goodsupport = shepardDenominator > supporttolerance;
                
                bool solved = false;


                if (goodsupport && Nneighbor >= minMdbcNeighbors)
                {
                    solved = solveLinear2(A, b, X, n);
                }

                if (!hasfluid)
                {
                    rhoghost[i] = rho0;
                    drhoghostX[i] = 0.0;
                    drhoghostY[i] = 0.0;
                }

                else if (solved)
                {
                    rhoghost[i] = X[0];
                    drhoghostX[i] = X[1];
                    drhoghostY[i] = X[2];
                }
                else
                {                    
                    rhoghost[i] = shepardNumerator / shepardDenominator;
                    drhoghostX[i] = 0.0;
                    drhoghostY[i] = 0.0;
                }


            rho[i] = rhoghost[i]+drhoghostX[i]*(x[i]-xghost[i])+drhoghostY[i]*(y[i]-yghost[i]);    
            }
            

// -----------------------------------------------------------------------------------------------------------------------------
// Compute Pressure at n
// -----------------------------------------------------------------------------------------------------------------------------        
            
            for (int i = 0; i < Nparticles; i++)
            {
                pressure[i] = B*(pow(rho[i]/rho0,gamma)-1.0);                            
            }

// -----------------------------------------------------------------------------------------------------------------------------
// Compute Kinetic Energy
// -----------------------------------------------------------------------------------------------------------------------------
            double KE = 0.0;
            for (int i = 0; i < Nparticles; i++)
            {
                KE += 0.5 * mass* ((u[i]*u[i])+(v[i]*v[i]));
            }


// -----------------------------------------------------------------------------------------------------------------------------
// Fluid Continuity and Momentum at n
// -----------------------------------------------------------------------------------------------------------------------------

            for (int i = Nboundary; i < Nparticles; i++)
            {
                //cummulative should be initialised with zero for each particle

                drhodt[i] = 0.0;
                dudt[i] = 0.0;
                dvdt[i] = 0.0;
                

                for (int j = 0; j < Nparticles; j++)
                {

                    double rx = x[i]-x[j];
                    double ry = y[i]-y[j];

                    double r2 = rx*rx+ry*ry;

                    if (r2 < 1e-14)
                    {
                        continue;
                    }
                    
                    double r = sqrt(r2);
                    double q = r/h;

                    double dirx = rx/r;
                    double diry = ry/r;

                    double du = u[i]-u[j];
                    double dv = v[i]-v[j];

                    double vijrij = du*rx+dv*ry;

                    double Piij = 0.0;

                    

                    //artificial viscosity
                    if (vijrij<0.0)
                    {
                        double muij = h*vijrij/(r*r+0.01*h*h);
                        double rhoij = 0.5*(rho[i]+rho[j]);
                        Piij = -alphaAV*c0*muij/rhoij;
                    }


                    
                    KernelResult result;

                    if (kernel == "gaussian")
                    {
                        result = gaussian(q, h, dirx, diry);
                    }
                    else if (kernel == "cubic")
                    {
                        result = cubicSpline(q, h, dirx, diry);
                    }
                    else if (kernel == "wendland")
                    {
                        result = Wendland(q, h, dirx, diry);
                    }
                    else
                    {
                        cout << "Invalid kernel choice. Please choose 'gaussian', 'cubic', or 'wendland'." << endl;
                        return 1; // Exit the program with an error code
                    }

                    double difussionX = 2*(rho[i]-rho[j])*(rx/r2);
                    double difussionY = 2*(rho[i]-rho[j])*(ry/r2);


                    drhodt[i] += ((du*result.dWeightX+dv*result.dWeightY)*mass)+(deltadifussion*h*c0*(mass/rho[j])*(difussionX*result.dWeightX+difussionY*result.dWeightY));
                    
                    
                    //dudt[i] += -mass*((pressure[i]/(rho[i]*rho[i]))+(pressure[j]/(rho[j]*rho[j]))+Piij)*result.dWeightX;
                    //dvdt[i] += -mass*((pressure[i]/(rho[i]*rho[i]))+(pressure[j]/(rho[j]*rho[j]))+Piij)*result.dWeightY;
                    
                    dudt[i] += -mass*(((pressure[i]+pressure[j])/(rho[i]*rho[j]))+Piij)*result.dWeightX;
                    dvdt[i] += -mass*(((pressure[i]+pressure[j])/(rho[i]*rho[j]))+Piij)*result.dWeightY;
                }

            }

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
            if (n % int(0.01/dt) == 0)
            {
                cout << "t : " << t << endl;
                string filename = foldername + "/Stillwater_hdp_" + to_string(h/dp) + "_t_" + to_string(t) + ".csv";

                ofstream file(filename);
                file << "h/dp :" << "," << h/dp << endl;
                file << "L2norm Pressure" << "," << "KE" << "," << "dp" << ","  << "Nparticles" << endl;
                file << L2normpressure << "," << KE << "," << dp << "," << Nfluid <<endl;
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
            }

// -----------------------------------------------------------------------------------------------------------------------------
// Compute Time integration for fluid to n+1/2
// -----------------------------------------------------------------------------------------------------------------------------

            for (int i = Nboundary; i < Nparticles; i++)
            {
                rhohalf[i] = rho[i]+drhodt[i]*dt/2;
                uhalf[i] = u[i]+dudt[i]*dt/2;
                vhalf[i] = v[i]+(dvdt[i]-g)*dt/2;
                xhalf[i] = x[i]+u[i]*dt/2;
                yhalf[i] = y[i]+v[i]*dt/2;
            }

// -----------------------------------------------------------------------------------------------------------------------------
// Compute Time integration for boundary to n+1/2
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
// mDBC Interpolation to n+1/2
// -----------------------------------------------------------------------------------------------------------------------------

            for (int i = 0; i < Nboundary; i++)
            {
                //cummulative should be initialised with zero for each particle

                drhodthalf[i] = 0.0;
                dudthalf[i] = 0.0;
                dvdthalf[i] = 0.0;

                const int n = 3;

                double A[n][n]={0};
                double b[n]={0};
                double X[n];

                int Nneighbor = 0;

                // Shepard-filter fallback used when the mDBC correction matrix
                // does not have sufficient / reliable fluid support.
                double shepardNumerator = 0.0;
                double shepardDenominator = 0.0;

                                

                for (int j = Nboundary; j < Nparticles; j++)
                {

                    double rx = xghost[i]-xhalf[j];
                    double ry = yghost[i]-yhalf[j];

                    double r = sqrt(rx*rx+ry*ry);
                    double q = r/h;

                    double dirx = 0.0;
                    double diry = 0.0;


                    if (r>0.0)
                    {
                        dirx = rx/r;
                        diry = ry/r;
                    }



                    
                    KernelResult result;

                    if (kernel == "gaussian")
                    {
                        result = gaussian(q, h, dirx, diry);
                    }
                    else if (kernel == "cubic")
                    {
                        result = cubicSpline(q, h, dirx, diry);
                    }
                    else if (kernel == "wendland")
                    {
                        result = Wendland(q, h, dirx, diry);
                    }
                    else
                    {
                        cout << "Invalid kernel choice. Please choose 'gaussian', 'cubic', or 'wendland'." << endl;
                        return 1; // Exit the program with an error code
                    }

                    if (abs(result.Weight) < 1e-14)
                    {
                        continue;
                    }

                    Nneighbor++;

                                      
                    
                    double Vj = mass/rhohalf[j];

                    shepardNumerator += rhohalf[j] * result.Weight * Vj;
                    shepardDenominator += result.Weight * Vj;


                    double dA[n][n] = 
                    {
                        {result.Weight*Vj, result.Weight*Vj*(-rx), result.Weight*Vj*(-ry)},
                        {result.dWeightX*Vj, result.dWeightX*Vj*(-rx), result.dWeightX*Vj*(-ry)},
                        {result.dWeightY*Vj, result.dWeightY*Vj*(-rx), result.dWeightY*Vj*(-ry)}
                    };
                    
                    double db[n] = 
                    {
                        result.Weight*mass,
                        result.dWeightX*mass,
                        result.dWeightY*mass
                    };

                    

                    for (int row = 0; row < 3; row++)
                    {
                        for (int col = 0; col < 3; col++)
                        {
                            A[row][col]+=dA[row][col];
                            
                        }
                        b[row]+=db[row];
                    }


                }

                // The mDBC formulation becomes unreliable with very low ghost support.
                // Use at least four fluid neighbours and reject singular/nearly-singular solves.
                const int minMdbcNeighbors = 4;
                double drytolerance = 1e-12;
                double supporttolerance = 0.4;


                bool hasfluid = shepardDenominator > drytolerance;
                bool goodsupport = shepardDenominator > supporttolerance;
                
                bool solved = false;


                if (goodsupport && Nneighbor >= minMdbcNeighbors)
                {
                    solved = solveLinear2(A, b, X, n);
                }

                if (!hasfluid)
                {
                    rhoghost[i] = rho0;
                    drhoghostX[i] = 0.0;
                    drhoghostY[i] = 0.0;
                }

                else if (solved)
                {
                    rhoghost[i] = X[0];
                    drhoghostX[i] = X[1];
                    drhoghostY[i] = X[2];
                }
                else
                {                    
                    rhoghost[i] = shepardNumerator / shepardDenominator;
                    drhoghostX[i] = 0.0;
                    drhoghostY[i] = 0.0;
                }


            rhohalf[i] = rhoghost[i]+drhoghostX[i]*(x[i]-xghost[i])+drhoghostY[i]*(yhalf[i]-yghost[i]);    
            }




// -----------------------------------------------------------------------------------------------------------------------------
// Compute Pressure at n+1/2
// -----------------------------------------------------------------------------------------------------------------------------        
            
            for (int i = 0; i < Nparticles; i++)
            {
                pressurehalf[i] = B*(pow(rhohalf[i]/rho0,gamma)-1.0);                            
            }

// -----------------------------------------------------------------------------------------------------------------------------
// Fluid Continuity and Momentum at n+1/2
// -----------------------------------------------------------------------------------------------------------------------------

            for (int i = Nboundary; i < Nparticles; i++)
            {
                //cummulative should be initialised with zero for each particle

                drhodthalf[i] = 0.0;
                dudthalf[i] = 0.0;
                dvdthalf[i] = 0.0;
                

                for (int j = 0; j < Nparticles; j++)
                {

                    double rx = xhalf[i]-xhalf[j];
                    double ry = yhalf[i]-yhalf[j];

                    double r2 = rx*rx+ry*ry;

                    if (r2 < 1e-14)
                    {
                        continue;
                    }
                    
                    double r = sqrt(r2);
                    double q = r/h;

                    double dirx = rx/r;
                    double diry = ry/r;

                    double du = uhalf[i]-uhalf[j];
                    double dv = vhalf[i]-vhalf[j];

                    double vijrij = du*rx+dv*ry;

                    double Piij = 0.0;

                    

                    //artificial viscosity
                    if (vijrij<0.0)
                    {
                        double muij = h*vijrij/(r*r+0.01*h*h);
                        double rhoij = 0.5*(rhohalf[i]+rhohalf[j]);
                        Piij = -alphaAV*c0*muij/rhoij;
                    }


                    
                    KernelResult result;

                    if (kernel == "gaussian")
                    {
                        result = gaussian(q, h, dirx, diry);
                    }
                    else if (kernel == "cubic")
                    {
                        result = cubicSpline(q, h, dirx, diry);
                    }
                    else if (kernel == "wendland")
                    {
                        result = Wendland(q, h, dirx, diry);
                    }
                    else
                    {
                        cout << "Invalid kernel choice. Please choose 'gaussian', 'cubic', or 'wendland'." << endl;
                        return 1; // Exit the program with an error code
                    }

                    double difussionX = 2*(rhohalf[i]-rhohalf[j])*(rx/r2);
                    double difussionY = 2*(rhohalf[i]-rhohalf[j])*(ry/r2);


                    drhodthalf[i] += ((du*result.dWeightX+dv*result.dWeightY)*mass)+(deltadifussion*h*c0*(mass/rhohalf[j])*(difussionX*result.dWeightX+difussionY*result.dWeightY));
                    
                    
                    //dudt[i] += -mass*((pressure[i]/(rho[i]*rho[i]))+(pressure[j]/(rho[j]*rho[j]))+Piij)*result.dWeightX;
                    //dvdt[i] += -mass*((pressure[i]/(rho[i]*rho[i]))+(pressure[j]/(rho[j]*rho[j]))+Piij)*result.dWeightY;
                    
                    dudthalf[i] += -mass*(((pressurehalf[i]+pressurehalf[j])/(rhohalf[i]*rhohalf[j]))+Piij)*result.dWeightX;
                    dvdthalf[i] += -mass*(((pressurehalf[i]+pressurehalf[j])/(rhohalf[i]*rhohalf[j]))+Piij)*result.dWeightY;
                }


            }

// -----------------------------------------------------------------------------------------------------------------------------
// Compute Time integration to n+1
// -----------------------------------------------------------------------------------------------------------------------------

// For Boundary
            for (int i = 0; i < Nboundary; i++)
            {
                rhonew[i] = rho[i]+drhodthalf[i]*dt;
                unew[i] = 0.0;
                vnew[i] = 0.0;
                xnew[i] = xhalf[i];
                ynew[i] = yhalf[i];
            }

// For Fluid          
            for (int i = Nboundary; i < Nparticles; i++)
            {
                rhonew[i] = rho[i]+drhodthalf[i]*dt;
                unew[i] = u[i]+dudthalf[i]*dt;
                vnew[i] = v[i]+(dvdthalf[i]-g)*dt;
                xnew[i] = x[i]+uhalf[i]*dt;
                ynew[i] = y[i]+vhalf[i]*dt;
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