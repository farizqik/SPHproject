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

const double dt = 0.001; 
const double Totaltime= 1.0;
const int Nt = Totaltime / dt;


// ------------------------------------------------------------
// Geometry
// ------------------------------------------------------------

const double tanklength = 25.0;
const double tankheight = 10.0;

const double waterlength = 25.0;
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

vector<double> pressure;
vector<double> pressurenew;

vector<double> drhodt;

vector<double> dudt;
vector<double> dvdt;

vector<double> xghost;
vector<double> yghost;
vector<double> rhoghost;






const double PI = 3.14159265358979323846;
const double g = 9.81;
const double rho0 = 1000.0;
const int boundpart = 0; 
const double c0 = 10.0*sqrt(g*dp*(waterheight));
const double mass = rho0*dp*dp;

const double gamma = 7.0;
const double B = c0*c0*rho0/gamma;





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
void solveLinear(double A[][3], double b[], double X[], int n)
{
    // Forward elimination
    for (int k = 0; k < n - 1; k++)
    {
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
        double sum = b[i];

        for (int j = i + 1; j < n; j++)
        {
            sum = sum - A[i][j] * X[j];
        }

        X[i] = sum / A[i][i];
    }
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
        cout << "h : " << h << endl;
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

        drhodt.resize(Nparticles);
        drhodtexact.resize(Nparticles);

        dudt.resize(Nparticles);
        dvdt.resize(Nparticles);

        xghost.resize(Nboundary);
        yghost.resize(Nboundary);
        rhoghost.resize(Nboundary);




        double L2normdrhodt = 0.0;




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

            rho[i] = rho0;
            
            i++;

        }
    }

    
    /*for (int i = 0; i < Nboundary; i++)
    {
        cout << "ID: " << i
             << "  Type: Boundary"
             << "  x: " << x[i]
             << "  y: " << y[i]
             << "  xghost: " << xghost[i]
             << "  yghost: " << yghost[i]
             << endl;
    }
    */
     




        
// -----------------------------------------------------------------------------------------------------------------------------
// start time loop
// -----------------------------------------------------------------------------------------------------------------------------

        for (int n=0;n<Nt+1;n++)
        {
            double L2drhodt = 0.0;

            double t = n*dt;
            
// -----------------------------------------------------------------------------------------------------------------------------
// Compute Pressure
// -----------------------------------------------------------------------------------------------------------------------------        
            
            for (int i = 0; i < Nparticles; i++)
            {
                pressure[i] = B*(pow(rho[i]/rho0,gamma)-1.0);
                //drhodtexact [i] = -rho[i]*(VelcoefX+VelcoefY);                               
            }

// -----------------------------------------------------------------------------------------------------------------------------
// Compute Continuity and Momentum
// -----------------------------------------------------------------------------------------------------------------------------

// For Fluid
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

                    double r = sqrt(rx*rx+ry*ry);
                    double q = r/h;

                    double dirx = 0.0;
                    double diry = 0.0;

                    double du = u[i]-u[j];
                    double dv = v[i]-v[j];




                    if (abs(rx)>0.0)
                    {
                        dirx = rx/r;
                    }
                    
                    if (abs(ry)>0.0)
                    {
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

                    drhodt[i] += (du*result.dWeightX+dv*result.dWeightY)*mass;
                    dudt[i] += -mass*((pressure[i]/(rho[i]*rho[i]))+(pressure[j]/(rho[j]*rho[j])))*result.dWeightX;
                    dvdt[i] += -mass*((pressure[i]/(rho[i]*rho[i]))+(pressure[j]/(rho[j]*rho[j])))*result.dWeightY;                                             
                }

                                  
                

            }

// For Boundary

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

                                

                for (int j = 0; j < Nparticles; j++)
                {

                    double rx = xghost[i]-x[j];
                    double ry = yghost[i]-y[j];

                    double r = sqrt(rx*rx+ry*ry);
                    double q = r/h;

                    double dirx = 0.0;
                    double diry = 0.0;

                    double du = u[i]-u[j];
                    double dv = v[i]-v[j];




                    if (abs(rx)>0.0)
                    {
                        dirx = rx/r;
                    }
                    
                    if (abs(ry)>0.0)
                    {
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

                    
                    double dA[n][n] = 
                    {
                        {result.Weight*dp*dp, result.Weight*dp*dp*(-rx), result.Weight*dp*dp*(-ry)},
                        {result.dWeightX*dp*dp, result.dWeightX*dp*dp*(-rx), result.dWeightX*dp*dp*(-ry)},
                        {result.dWeightY*dp*dp, result.dWeightY*dp*dp*(-rx), result.dWeightY*dp*dp*(-ry)}
                    };
                    
                    double db[n] = 
                    {
                        {result.Weight*mass},
                        {result.dWeightX*mass},
                        {result.dWeightY*mass}
                    };

                    

                    for (int row = 0; row < 3; row++)
                    {
                        for (int col = 0; col < 3; col++)
                        {
                            A[row][col]+=dA[row][col];
                            
                        }
                        b[row]+=db[row];
                    }

                    drhodt[i] += (du*result.dWeightX+dv*result.dWeightY)*mass;                    
                    dudt[i] = 0;
                    dvdt[i] = 0;
                }

                
                /*cout << "\nBoundary particle i = " << i << endl;

                cout << "A =" << endl;

                for (int row = 0; row < 3; row++)
                {
                    for (int col = 0; col < 3; col++)
                    {
                        cout << A[row][col] << "\t";
                    }

                    cout << endl;
                }

                cout << "b = ";

                for (int row = 0; row < 3; row++)
                {
                    cout << b[row] << "\t";
                }

                cout << endl;*/


                solveLinear (A, b, X, n);

                for (int p =0; p < n; p++)
                {
                    cout << "X[" << p << "] = " << X[p] << endl;
                }
                    
                    


                      
                    
                


                 
            }
// -----------------------------------------------------------------------------------------------------------------------------
// Compute Time integration
// -----------------------------------------------------------------------------------------------------------------------------

// For Boundary
            for (int i = 0; i < Nboundary; i++)
            {
                rhonew[i] = rho[i]+drhodt[i]*dt;
                unew[i] = 0.0;
                vnew[i] = 0.0;
                xnew[i] = x[i];
                ynew[i] = y[i];
            }

// For Fluid          
            for (int i = Nboundary; i < Nparticles; i++)
            {
                rhonew[i] = rho[i]+drhodt[i]*dt;
                unew[i] = u[i]+dudt[i]*dt;
                vnew[i] = v[i]+(dvdt[i]-g)*dt;
                xnew[i] = x[i]+unew[i]*dt;
                ynew[i] = y[i]+vnew[i]*dt;
            }

// next loop for all particles 
            for (int i = 0; i < Nparticles; i++)
            {
                x[i] = xnew[i];
                y[i] = ynew[i];
                u[i] = unew[i];
                v[i] = vnew[i];
                rho[i] = rhonew [i];


                double errordrhodt = (pow(drhodt[i]-drhodtexact[i],2));
                L2drhodt += errordrhodt;
            }
               
                
            

            L2normdrhodt = sqrt(L2drhodt/(Nparticles));

// -----------------------------------------------------------------------------------------------------------------------------
// Write output
// -----------------------------------------------------------------------------------------------------------------------------
            if (n % int(0.01/dt) == 0)
            {
                cout << "t : " << t << endl;
                string filename = foldername + "/Stillwater_hdp_" + to_string(h/dp) + "_t_" + to_string(t) + ".csv";

                ofstream file(filename);
                file << "h/dp :" << "," << h/dp << endl;
                file << "L2normdrhodt" << "," << "dp" << ","  << "Nparticles" << endl;
                file << L2normdrhodt << "," << dp << "," << Nparticles <<endl;
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