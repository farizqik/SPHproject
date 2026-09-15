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

const double c0 = 10.0*sqrt(g*(waterheight));

const double gammaEOS = 7.0;
const double BEOS = c0*c0*rho0/gammaEOS;

const double viscosity = 1.0e-6;
const double tvis = 1.0;
const double Dcylinder = 0.1;
const double r0 = 0.5*Dcylinder;
const double R = 20*Dcylinder;
const double dr0 = sqrt(2.0*viscosity*tvis);

const int N = 20;




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

double target = (R-r0)/dr0;

double radialSum(double qA, int N)
{
    if (abs(qA - 1.0) < 1.0e-12)
    {
        return N;
    }

    return (pow(qA, N) - 1.0) / (qA - 1.0);
}

int main()
{
    double qAlow = 1.0;
    double qAhigh = 1.5;
    double maxiter = 1000;
    double tolerance = 1.0e-6;
    double qA = 1.0;
    
    for (int iter = 0; iter < 1000; ++iter)
    {
        double qAmid = 0.5*(qAlow+qAhigh);
        double error = radialSum(qAmid, N) - target;

        if (abs(error) < tolerance)
        {
            qA = qAmid;
            break;
        }
        if (radialSum(qAmid,N) < target)
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




    double dTheta = dr0/r0;
    double NTheta =  2*PI/dTheta;


    // ------------------------------------------------------------
    // Particle coordinates
    // ------------------------------------------------------------

    vector<double> r;
    vector<double> theta;

    vector<double> x;
    vector<double> y;

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
        double thetap = j*dTheta;
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


    for (int ir = 1; ir <= N; ++ir)
    {
        double rp = Bpolar*(pow(qA,ir)-1)+r0;

        for (int j = 0; j < NTheta; ++j)
        {
            double thetap = j*dTheta;
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

    string kernel; 
    cout<<"enter the kernel to be used (gaussian, cubic, wendland): ";
    cin>>kernel;

    string foldername = 
            "EWCSPHCylinder_dr0_" + to_string(dr0) + "_Nr_" + to_string(N) + "_Ntheta_" + to_string(NTheta) + "_Rr0_" + to_string(R/r0)+ "_" + kernel; 
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
    // Save particle positions
    // ------------------------------------------------------------

    string filename =
        foldername + "/particles.csv";

    ofstream file(filename);


    // ------------------------------------------------------------
    // Simulation / particle information
    // ------------------------------------------------------------

    file << "dr0" << "," << dr0 << endl;
    file << "r0" << "," << r0 << endl;
    file << "R" << "," << R << endl;
    file << "R/r0" << "," << R/r0 << endl;
    file << "Nr" << "," << N << endl;
    file << "Ntheta" << "," << NTheta << endl;
    file << "Nboundary" << "," << Nboundary << endl;
    file << "Nfluid" << "," << Nfluid << endl;
    file << "Nparticles" << "," << Nparticles << endl;
    file << "Apolar" << "," << Apolar << endl;
    file << "Bpolar" << "," << Bpolar << endl;

    file << endl;


    // ------------------------------------------------------------
    // Particle data
    // ------------------------------------------------------------

    file << "ID,r,theta,x,y,Type" << endl;

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
            << Type
            << endl;
    }

    file.close();

    cout << endl;
    cout << "Particle file saved to: "
        << filename
        << endl;
    


    
}
