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

string Type; 

// --------------------------------------------
// Timestep
// --------------------------------------------

const double dt = 0.001; 
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

vector<double> r;
vector<double> z;

vector<double> rnew;
vector<double> znew;

vector<double> u_r;
vector<double> u_z;

vector<double> u_rnew;
vector<double> u_znew;

vector<double> rho;
vector<double> rhonew;

vector<double> drhodtexact;
vector<double> pressureexact;

vector<double> mass;

vector<double> pressure;
vector<double> pressurenew;

vector<double> drhodt;

vector<double> du_rdt;
vector<double> du_zdt;

vector<double> rghost;
vector<double> zghost;
vector<double> rhoghost;
vector<double> drhoghostR;
vector<double> drhoghostZ;

vector<double> rhalf;
vector<double> zhalf;
vector<double> u_rhalf;
vector<double> u_zhalf;
vector<double> rhohalf;

vector<double> drhodthalf;
vector<double> pressurehalf;
vector<double> du_rdthalf;
vector<double> du_zdthalf;


const double PI = 3.14159265358979323846;
const double g = 9.81;
const double rho0 = 1000.0;
const int boundpart = 0; 
const double c0 = 10.0*sqrt(g*(waterheight));


const double gammaEOS = 7.0;
const double B = c0*c0*rho0/gammaEOS;

const double alphaAV = 0.01;


//set deltadifussion = 0.0 if it's not considered
const double deltadifussion = 0.1;





// -----------------------------------------------------------------------------------------------------------------------------
// Kernel functions
// -----------------------------------------------------------------------------------------------------------------------------

struct KernelResult {
    double Weight;
    double dWeightR;
    double dWeightZ;
};


// -----------------------------------------------------------------------------------------------------------------------------
// Gaussian kernel
// -----------------------------------------------------------------------------------------------------------------------------
KernelResult gaussian(double q, double h, double dirR, double dirZ)
{
    KernelResult result;
    double alpha = 1.0 / (PI * h * h);

    result.Weight = alpha*exp(-q*q);
    result.dWeightR = -2.0*alpha*q*exp(-q*q)/h * dirR;
    result.dWeightZ = -2.0*alpha*q*exp(-q*q)/h * dirZ;

    return result;
}


// -----------------------------------------------------------------------------------------------------------------------------
// Cubic kernel
// -----------------------------------------------------------------------------------------------------------------------------
KernelResult cubicSpline(double q, double h, double dirR, double dirZ)
{
    double alpha = 10.0 / (7.0 * PI * h * h);
    KernelResult result;
    if (q<=1.0)
    {
        result.Weight = alpha*(1.0-1.5*pow(q,2)+0.75*pow(q,3));
        result.dWeightR = alpha*(-3.0*q+2.25*pow(q,2))/h * dirR;
        result.dWeightZ = alpha*(-3.0*q+2.25*pow(q,2))/h * dirZ;
    }
    else if (q<=2.0)
    {
        result.Weight = alpha*0.25*pow(2.0-q,3);
        result.dWeightR = -0.75*alpha*pow(2.0-q,2)/h * dirR;
        result.dWeightZ = -0.75*alpha*pow(2.0-q,2)/h * dirZ;
    }
    else
    {
        result.Weight = 0.0;
        result.dWeightR = 0.0;
        result.dWeightZ = 0.0;
    }

    return result;
}

// -----------------------------------------------------------------------------------------------------------------------------
// Wendland kernel
// -----------------------------------------------------------------------------------------------------------------------------
KernelResult Wendland(double q, double h, double dirR, double dirZ)
{
    double alpha = 7.0 / (4.0 * PI * h * h);
    KernelResult result;
    if (0.0<=q && q<=2.0)
    {
        result.Weight = alpha*pow(1.0-0.5*q,4)*(2.0*q+1.0);
        result.dWeightR = alpha*(4*pow(1.0-0.5*q,3)*(-0.5)*(2.0*q+1.0)+pow(1.0-0.5*q,4)*2.0)/h * dirR;
        result.dWeightZ = alpha*(4*pow(1.0-0.5*q,3)*(-0.5)*(2.0*q+1.0)+pow(1.0-0.5*q,4)*2.0)/h * dirZ;
    }
    else
    {
        result.Weight = 0.0;
        result.dWeightR = 0.0;
        result.dWeightZ = 0.0;
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

    for (double zp = -0.5*dp; zp > -boundthick; zp -= dp)
    {
        for (double rp = -boundthick+0.5*dp; rp < tanklength + boundthick + 0.5*dp; rp += dp)
        {

            cout << "ID: " << i
             << "  Type: Boundary"
             << "  x: " << rp
             << "  y: " << zp
             << endl;

            i++;

        }
    }

    // Left boundary
    for (double rp = -0.5*dp; rp > -boundthick; rp -= dp)
    {
        for (double zp = 0.5*dp; zp < tankheight; zp += dp)
        {

            cout << "ID: " << i
             << "  Type: Boundary"
             << "  r: " << rp
             << "  z: " << zp
             << endl;

            i++;
        }
    }

    // Right boundary
    for (double rp = tanklength + 0.5*dp;
         rp < tanklength + boundthick;
         rp += dp)
    {
        for (double zp = 0.5*dp; zp < tankheight; zp += dp)
        {

            cout << "ID: " << i
             << "  Type: Boundary"
             << "  x: " << rp
             << "  y: " << zp
             << endl;


            i++;
            
        }
    }

    int Nboundary = i;

    for (double rp = 0.5*dp; rp < waterlength; rp += dp)
    {
        for (double zp = 0.5*dp; zp < waterheight; zp += dp)
        {

            cout << "ID: " << i
             << "  Type: Fluid"
             << "  x: " << rp
             << "  y: " << zp
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
            "WCSPHpolar_dp_" + to_string(dp) + "_h_" + to_string(h) + "_Nparticles_" + to_string(Nparticles) + "_" + kernel; 
        system(("mkdir -p " + foldername).c_str());
        //std::filesystem::create_directories(foldername);

        r.resize(Nparticles);
        z.resize(Nparticles);

        rnew.resize(Nparticles);
        znew.resize(Nparticles);

        u_r.resize(Nparticles);
        u_z.resize(Nparticles);

        u_rnew.resize(Nparticles);
        u_znew.resize(Nparticles);

        rho.resize(Nparticles);
        rhonew.resize(Nparticles);

        pressure.resize(Nparticles);
        pressurenew.resize(Nparticles);
        pressureexact.resize(Nparticles);

        drhodt.resize(Nparticles);
        drhodtexact.resize(Nparticles);

        du_rdt.resize(Nparticles);
        du_zdt.resize(Nparticles);

        rghost.resize(Nboundary);
        zghost.resize(Nboundary);
        rhoghost.resize(Nboundary);
        drhoghostR.resize(Nboundary);
        drhoghostZ.resize(Nboundary);

        rhalf.resize(Nparticles);
        zhalf.resize(Nparticles);
        rhohalf.resize(Nparticles);
        u_rhalf.resize(Nparticles);
        u_zhalf.resize(Nparticles);
        drhodthalf.resize(Nparticles);
        pressurehalf.resize(Nparticles);
        du_rdthalf.resize(Nparticles);
        du_zdthalf.resize(Nparticles);

        mass.resize(Nparticles);

      

// -----------------------------------------------------------------------------------------------------------------------------
// Initial condition
// -----------------------------------------------------------------------------------------------------------------------------
    
    int i = 0;

    //Bottom Boundary Particles

    for (double zp = -0.5*dp; zp > -boundthick; zp -= dp)
    {
        for (double rp = -boundthick+0.5*dp; rp < tanklength + boundthick + 0.5*dp; rp += dp)
        {

            r[i] = rp;
            z[i] = zp;

            u_r[i] = 0.0;
            u_z[i] = 0.0;

            if (rp < 0)
            {
                rghost[i] = -r[i];
                zghost[i] = -z[i];
            }

            else if (rp > tanklength)
            {
                rghost[i] = 2*tanklength-r[i];
                zghost[i] = -z[i];
            }

            else
            {
                rghost[i] = r[i];
                zghost[i] = -z[i];
            }
            
            rho[i] = rho0;

            i++;

        }
    }

    // Left boundary
    for (double rp = -0.5*dp; rp > -boundthick; rp -= dp)
    {
        for (double zp = 0.5*dp; zp < tankheight; zp += dp)
        {

            r[i] = rp;
            z[i] = zp;

            rghost[i] = -r[i];
            zghost[i] = z[i];

            u_r[i] = 0.0;
            u_z[i] = 0.0;

            rho[i] = rho0;
            
            i++;
        }
    }

    // Right boundary
    for (double rp = tanklength + 0.5*dp;
         rp < tanklength + boundthick;
         rp += dp)
    {
        for (double zp = 0.5*dp; zp < tankheight; zp += dp)
        {

            r[i] = rp;
            z[i] = zp;

            rghost[i] = 2*tanklength-r[i];
            zghost[i] = z[i];

            u_r[i] = 0.0;
            u_z[i] = 0.0;

            rho[i] = rho0;
            
            i++;
            
        }
    }


    for (double rp = 0.5*dp; rp < waterlength; rp += dp)
    {
        for (double zp = 0.5*dp; zp < waterheight; zp += dp)
        {

            r[i] = rp;
            z[i] = zp;

            u_r[i] = 0.0;
            u_z[i] = 0.0;
            
            //rho[i] = rho0;

            // Exact hydrostatic initial density for the Tait EOS used here.
            // It satisfies dp/dy = -rho*g in the continuum.
            double depth = max(0.0, waterheight - z[i]);
            rho[i] = rho0 * pow(1.0 + (gammaEOS - 1.0)*g*depth/(c0*c0), 1.0/(gammaEOS - 1.0));
            
            i++;

        }
    }

    for (int i = 0; i < Nparticles; i++)
    {
        mass[i] = 2.0*PI*r[i]*rho[i]*dp*dp;
        drhodtexact[i] = -rho[i]*(0.0);
        pressureexact[i] = rho0*g*(waterheight-z[i]);                                
    }

     




        
// -----------------------------------------------------------------------------------------------------------------------------
// start time loop
// -----------------------------------------------------------------------------------------------------------------------------

        for (int n=0;n<Nt+1;n++)
        {

            double t = n*dt;

            
// -----------------------------------------------------------------------------------------------------------------------------
// mDBC Interpolation at n
// -----------------------------------------------------------------------------------------------------------------------------
            //#pragma omp parallel for
            for (int i = 0; i < Nboundary; i++)
            {
                //cummulative should be initialised with zero for each particle

                drhodt[i] = 0.0;
                du_rdt[i] = 0.0;
                du_zdt[i] = 0.0;

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

                    double sr = rghost[i]-r[j];
                    double sz = zghost[i]-z[j];

                    double s2 = sr*sr+sz*sz;

                    if ((kernel == "cubic" || kernel == "wendland") && s2 > 4*h*h)
                    {
                        continue;
                    }

                    if (s2 < 1e-14)
                    {
                        continue;
                    }
                    
                    double ds = sqrt(s2);
                    double q = ds/h;

                    double dirR = sr/ds;
                    double dirZ = sz/ds;

                    
                    KernelResult result;

                    if (kernel == "gaussian")
                    {
                        result = gaussian(q, h, dirR, dirZ);
                    }
                    else if (kernel == "cubic")
                    {
                        result = cubicSpline(q, h, dirR, dirZ);
                    }
                    else if (kernel == "wendland")
                    {
                        result = Wendland(q, h, dirR, dirZ);
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

                                      
                    
                    double Aj = mass[j]/(2*PI*r[j]*rho[j]);

                    shepardNumerator += rho[j] * result.Weight * Aj;
                    shepardDenominator += result.Weight * Aj;


                    double dA[n][n] = 
                    {
                        {result.Weight*Aj, result.Weight*Aj*(-sr), result.Weight*Aj*(-sz)},
                        {result.dWeightR*Aj, result.dWeightR*Aj*(-sr), result.dWeightR*Aj*(-sz)},
                        {result.dWeightZ*Aj, result.dWeightZ*Aj*(-sr), result.dWeightZ*Aj*(-sz)}
                    };
                    
                    double db[n] = 
                    {
                        result.Weight*rho[j]*Aj,
                        result.dWeightR*rho[j]*Aj,
                        result.dWeightZ*rho[j]*Aj
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
                    solved = solveLinear(A, b, X, n);
                }

                if (!hasfluid)
                {
                    rhoghost[i] = rho0;
                    drhoghostR[i] = 0.0;
                    drhoghostZ[i] = 0.0;
                }

                else if (solved)
                {
                    rhoghost[i] = X[0];
                    drhoghostR[i] = X[1];
                    drhoghostZ[i] = X[2];
                }
                else
                {                    
                    rhoghost[i] = shepardNumerator / shepardDenominator;
                    drhoghostR[i] = 0.0;
                    drhoghostZ[i] = 0.0;
                }


            rho[i] = rhoghost[i]+drhoghostR[i]*(r[i]-rghost[i])+drhoghostZ[i]*(z[i]-zghost[i]);    
            }
            

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
                if (z[i]<-boundthick)
                {
                    continue;
                }
                KE += 0.5 * mass[i]* ((u_r[i]*u_r[i])+(u_z[i]*u_z[i]));
            }


// -----------------------------------------------------------------------------------------------------------------------------
// Fluid Continuity and Momentum at n
// -----------------------------------------------------------------------------------------------------------------------------
            //#pragma omp parallel for
            for (int i = Nboundary; i < Nparticles; i++)
            {
                //cummulative should be initialised with zero for each particle

                drhodt[i] = 0.0;
                du_rdt[i] = 0.0;
                du_zdt[i] = 0.0;
                

                for (int j = 0; j < Nparticles; j++)
                {
                    
                    double sr = r[i]-r[j];
                    double sz = z[i]-z[j];

                    double s2 = sr*sr+sz*sz;

                    if ((kernel == "cubic" || kernel == "wendland") && s2 > 4*h*h)
                    {
                        continue;
                    }

                    if (s2 < 1e-14)
                    {
                        continue;
                    }
                    
                    double ds = sqrt(s2);
                    double q = ds/h;

                    double dirR = sr/ds;
                    double dirZ = sz/ds;

                    double du_r = u_r[i]-u_r[j];
                    double du_z = u_z[i]-u_z[j];

                    double vijrij = du_r* sr+du_z* sz;

                    double Piij = 0.0;

                    

                    //artificial viscosity
                    if (vijrij<0.0)
                    {
                        double muij = h*vijrij/(ds*ds+0.01*h*h);
                        double rhoij = 0.5*(rho[i]+rho[j]);
                        Piij = -alphaAV*c0*muij/rhoij;
                    }


                    
                    KernelResult result;

                    if (kernel == "gaussian")
                    {
                        result = gaussian(q, h, dirR, dirZ);
                    }
                    else if (kernel == "cubic")
                    {
                        result = cubicSpline(q, h, dirR, dirZ);
                    }
                    else if (kernel == "wendland")
                    {
                        result = Wendland(q, h, dirR, dirZ);
                    }
                    else
                    {
                        cout << "Invalid kernel choice. Please choose 'gaussian', 'cubic', or 'wendland'." << endl;
                        return 1; // Exit the program with an error code
                    }

                    double difussionR = 2*(rho[i]-rho[j])*(sr/s2);
                    double difussionZ = 2*(rho[i]-rho[j])*(sz/s2);


                    drhodt[i] += (1.0/(2*PI))*((mass[j]/r[j])*(du_r*result.dWeightR+du_z*result.dWeightZ));                 
                    du_rdt[i] -= 2*PI*mass[j]*(((pressure[i]*r[i]+pressure[j]*r[j])/(2*PI*r[i]*rho[i]*2*PI*r[j]*rho[j]))+(Piij/(2*PI)))*result.dWeightR;
                    du_zdt[i] -= 2*PI*mass[j]*(((pressure[i]*r[i]+pressure[j]*r[j])/(2*PI*r[i]*rho[i]*2*PI*r[j]*rho[j]))+(Piij/(2*PI)))*result.dWeightZ;
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
                string filename = foldername + "/WCSPHpolar_hdp_" + to_string(h/dp) + "_t_" + to_string(t) + ".csv";

                ofstream file(filename);
                file << "h/dp :" << "," << h/dp << endl;
                file << "L2norm Pressure" << "," << "KE" << "," << "dp" << ","  << "Nparticles" << endl;
                file << L2normpressure << "," << KE << "," << dp << "," << Nfluid <<endl;
                file << endl;
                file << "t" << "," << t << endl;
                file << "ID"<< ","<< "x"<< "," << "y" << ","<< ","<<"rho"<< "," << "drhodt" << "," << "pressure"<< "," << "," << "u" << "," << "v" << "," << "," << "du_rdt" << "," << "du_zdt" <<"," << "Type" << endl;
                for (int i = 0; i < Nparticles; i++)
                    {
                        if (i < Nboundary)
                        {
                            Type = "boundary";
                        }
                        else
                        {
                            Type = "fluid";
                        }
                        file << i << "," << r[i] << "," << z[i] << "," << "," << rho[i] << "," << drhodt[i] << "," << pressure[i] << "," << "," << u_r[i] << "," << u_z[i] << "," << "," << du_rdt[i] << "," << du_zdt[i] << "," << Type << endl;
                    }
           
                file.close();
            }

// -----------------------------------------------------------------------------------------------------------------------------
// Compute Time integration predictor for fluid to n+1/2
// -----------------------------------------------------------------------------------------------------------------------------

            for (int i = Nboundary; i < Nparticles; i++)
            {
                rhohalf[i] = rho[i]+(drhodt[i]-(rho[i]*u_r[i]/r[i]))*dt/2;
                u_rhalf[i] = u_r[i]+(du_rdt[i]+(pressure[i]/(rho[i]*r[i])))*dt/2;
                u_zhalf[i] = u_z[i]+(du_zdt[i]-g)*dt/2;
                rhalf[i] = r[i]+u_r[i]*dt/2;
                zhalf[i] = z[i]+u_z[i]*dt/2;
            }

// -----------------------------------------------------------------------------------------------------------------------------
// Compute Time integration predictor for boundary to n+1/2
// -----------------------------------------------------------------------------------------------------------------------------

            for (int i = 0; i < Nboundary; i++)
            {

                u_rhalf[i] = 0.0;
                u_zhalf[i] = 0.0;
                rhalf[i] = r[i];
                zhalf[i] = z[i];
            }

            
// -----------------------------------------------------------------------------------------------------------------------------
// mDBC Interpolation to n+1/2
// -----------------------------------------------------------------------------------------------------------------------------
            //#pragma omp parallel for
            for (int i = 0; i < Nboundary; i++)
            {
                //cummulative should be initialised with zero for each particle

                drhodthalf[i] = 0.0;
                du_rdthalf[i] = 0.0;
                du_zdthalf[i] = 0.0;

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

                    double sr = rghost[i]-rhalf[j];
                    double sz = zghost[i]-zhalf[j];

                    double s2 = sr*sr+sz*sz;

                    if ((kernel == "cubic" || kernel == "wendland") && s2 > 4*h*h)
                    {
                        continue;
                    }

                    if (s2 < 1e-14)
                    {
                        continue;
                    }
                    
                    double ds = sqrt(s2);
                    double q = ds/h;

                    double dirR = sr/ds;
                    double dirZ = sz/ds;

                    
                    KernelResult result;

                    if (kernel == "gaussian")
                    {
                        result = gaussian(q, h, dirR, dirZ);
                    }
                    else if (kernel == "cubic")
                    {
                        result = cubicSpline(q, h, dirR, dirZ);
                    }
                    else if (kernel == "wendland")
                    {
                        result = Wendland(q, h, dirR, dirZ);
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

                                      
                    
                    double Aj = mass[j]/(2*PI*rhalf[j]*rhohalf[j]);

                    shepardNumerator += rhohalf[j] * result.Weight * Aj;
                    shepardDenominator += result.Weight * Aj;


                    double dA[n][n] = 
                    {
                        {result.Weight*Aj, result.Weight*Aj*(-sr), result.Weight*Aj*(-sz)},
                        {result.dWeightR*Aj, result.dWeightR*Aj*(-sr), result.dWeightR*Aj*(-sz)},
                        {result.dWeightZ*Aj, result.dWeightZ*Aj*(-sr), result.dWeightZ*Aj*(-sz)}
                    };
                    
                    double db[n] = 
                    {
                        result.Weight*rhohalf[j]*Aj,
                        result.dWeightR*rhohalf[j]*Aj,
                        result.dWeightZ*rhohalf[j]*Aj
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
                    solved = solveLinear(A, b, X, n);
                }

                if (!hasfluid)
                {
                    rhoghost[i] = rho0;
                    drhoghostR[i] = 0.0;
                    drhoghostZ[i] = 0.0;
                }

                else if (solved)
                {
                    rhoghost[i] = X[0];
                    drhoghostR[i] = X[1];
                    drhoghostZ[i] = X[2];
                }
                else
                {                    
                    rhoghost[i] = shepardNumerator / shepardDenominator;
                    drhoghostR[i] = 0.0;
                    drhoghostZ[i] = 0.0;
                }


            rhohalf[i] = rhoghost[i]+drhoghostR[i]*(rhalf[i]-rghost[i])+drhoghostZ[i]*(zhalf[i]-zghost[i]);    
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
            for (int i = Nboundary; i < Nparticles; i++)
            {
                //cummulative should be initialised with zero for each particle

                drhodthalf[i] = 0.0;
                du_rdthalf[i] = 0.0;
                du_zdthalf[i] = 0.0;
                

                for (int j = 0; j < Nparticles; j++)
                {

                    double sr = rhalf[i]-rhalf[j];
                    double sz = zhalf[i]-zhalf[j];

                    double s2 = sr*sr+sz*sz;

                    if ((kernel == "cubic" || kernel == "wendland") && s2 > 4*h*h)
                    {
                        continue;
                    }

                    if (s2 < 1e-14)
                    {
                        continue;
                    }
                    
                    double ds = sqrt(s2);
                    double q = ds/h;

                    double dirR = sr/ds;
                    double dirZ = sz/ds;

                    double du_r = u_rhalf[i]-u_rhalf[j];
                    double du_z = u_zhalf[i]-u_zhalf[j];

                    double vijrij = du_r* sr+du_z* sz;

                    double Piij = 0.0;

                    

                    //artificial viscosity
                    if (vijrij<0.0)
                    {
                        double muij = h*vijrij/(ds*ds+0.01*h*h);
                        double rhoij = 0.5*(rhohalf[i]+rhohalf[j]);
                        Piij = -alphaAV*c0*muij/rhoij;
                    }


                    
                    KernelResult result;

                    if (kernel == "gaussian")
                    {
                        result = gaussian(q, h, dirR, dirZ);
                    }
                    else if (kernel == "cubic")
                    {
                        result = cubicSpline(q, h, dirR, dirZ);
                    }
                    else if (kernel == "wendland")
                    {
                        result = Wendland(q, h, dirR, dirZ);
                    }
                    else
                    {
                        cout << "Invalid kernel choice. Please choose 'gaussian', 'cubic', or 'wendland'." << endl;
                        return 1; // Exit the program with an error code
                    }

                    double difussionR = 2*(rhohalf[i]-rhohalf[j])*(sr/s2);
                    double difussionZ = 2*(rhohalf[i]-rhohalf[j])*(sz/s2);


                    drhodthalf[i] += (1.0/(2*PI))*((mass[j]/rhalf[j])*(du_r*result.dWeightR+du_z*result.dWeightZ));                               
                    du_rdthalf[i] -= 2*PI*mass[j]*(((pressurehalf[i]*rhalf[i]+pressurehalf[j]*rhalf[j])/(2*PI*rhalf[i]*rhohalf[i]*2*PI*rhalf[j]*rhohalf[j]))+(Piij/(2*PI)))*result.dWeightR;
                    du_zdthalf[i] -= 2*PI*mass[j]*(((pressurehalf[i]*rhalf[i]+pressurehalf[j]*rhalf[j])/(2*PI*rhalf[i]*rhohalf[i]*2*PI*rhalf[j]*rhohalf[j]))+(Piij/(2*PI)))*result.dWeightZ;
                }



            }


// -----------------------------------------------------------------------------------------------------------------------------
// Compute Time integration corrector for fluid to n+1/2
// -----------------------------------------------------------------------------------------------------------------------------

            for (int i = Nboundary; i < Nparticles; i++)
            {
                rhohalf[i] = rho[i]+(drhodthalf[i]-(rhohalf[i]*u_rhalf[i]/rhalf[i]))*dt/2;
                u_rhalf[i] = u_r[i]+(du_rdthalf[i]+(pressurehalf[i]/(rhohalf[i]*rhalf[i])))*dt/2;
                u_zhalf[i] = u_z[i]+(du_zdthalf[i]-g)*dt/2;
                rhalf[i] = r[i]+u_rhalf[i]*dt/2;
                zhalf[i] = z[i]+u_zhalf[i]*dt/2;
            }

// -----------------------------------------------------------------------------------------------------------------------------
// Compute Time integration corrector for boundary to n+1/2
// -----------------------------------------------------------------------------------------------------------------------------

            for (int i = 0; i < Nboundary; i++)
            {

                u_rhalf[i] = 0.0;
                u_zhalf[i] = 0.0;
                rhalf[i] = r[i];
                zhalf[i] = z[i];
            }


// -----------------------------------------------------------------------------------------------------------------------------
// Compute Time integration to n+1
// -----------------------------------------------------------------------------------------------------------------------------

// For Boundary
            for (int i = 0; i < Nboundary; i++)
            {
                rhonew[i] = rhohalf[i];
                u_rnew[i] = 0.0;
                u_znew[i] = 0.0;
                rnew[i] = rhalf[i];
                znew[i] = zhalf[i];
            }

// For Fluid          
            for (int i = Nboundary; i < Nparticles; i++)
            {
                rhonew[i] = 2*rhohalf[i]-rho[i];
                u_rnew[i] = 2*u_rhalf[i]-u_r[i];
                u_znew[i] = 2*u_zhalf[i]-u_z[i];
                rnew[i] = 2*rhalf[i]-r[i];
                znew[i] = 2*zhalf[i]-z[i];
            }

           
// -----------------------------------------------------------------------------------------------------------------------------
// next loop for all particles
// -----------------------------------------------------------------------------------------------------------------------------
            
            for (int i = 0; i < Nparticles; i++)
            {
                r[i] = rnew[i];
                z[i] = znew[i];
                u_r[i] = u_rnew[i];
                u_z[i] = u_znew[i];
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