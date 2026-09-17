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
const double rho = 1000.0;

const double c0 = 10.0*sqrt(g*(waterheight));

const double gammaEOS = 7.0;
const double BEOS = c0*c0*rho/gammaEOS;

const double viscosity = 1.0e-6;
const double tvis = 1.0;
const double Dcylinder = 0.1;
const double r0 = 0.5*Dcylinder;
const double R = 20*Dcylinder;
const double dr0 = sqrt(2.0*viscosity*tvis);
double drhodtexact = 0;

const int Nr = 100;

string kernel; 




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


// ------------------------------------------------------------
// Particle coordinates
// ------------------------------------------------------------

vector<double> r;
vector<double> theta;

vector<double> x;
vector<double> y;

vector<double> u;
vector<double> v;
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
vector<double> hloc;




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
    double qAlow = 0.0;
    double qAhigh = 2.0;
    double maxiter = 1000;
    double tolerance = 1.0e-6;
    double qA = 1.0;
    
    for (int iter = 0; iter < 1000; ++iter)
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

    //cout<<"enter the kernel to be used (gaussian, cubic, wendland): ";
        //cin>>kernel;

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


    vector<double> hlist = {1.2, 1.5, 1.8, 2};
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


        for (int i = 0; i < Nparticles; i++)
        {
            u[i] = 1.0;
            v[i] = 0.0;
        }


        // kernel consistency check

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

            int ir = i/NTheta;

            double dri;
            if (ir == 0)
            {
                dri =
                    0.5*(r[i + NTheta] - r[i]);
            }
            else if (ir == Nr)
            {
                dri =
                    0.5*(r[i] - r[i - NTheta]);
            }
            else
            {
                dri =
                    0.5*(r[i + NTheta]
                    - r[i - NTheta]);
            }

            double drthetai = r[i]*dTheta0;

            double dpi = sqrt(dri*drthetai);

            double hp = coefh*dpi;
            

            for (int j = 0 ; j < Nparticles; j++)

            {
                int jr = j/NTheta;

                double drj;
                double rin, rout;

                if (jr == 0)
                {
                    rin = r0;
                    rout = 0.5*(r[j] + r[j + NTheta]);
                }
                else if (jr == Nr)
                {
                    rin = 0.5*(r[j] + r[j - NTheta]);
                    rout = R;
                }
                else
                {
                    rin = 0.5*(r[j] + r[j - NTheta]);
                    rout = 0.5*(r[j] + r[j + NTheta]);
                }

                double Vj =
                    0.5*(rout*rout - rin*rin)*dTheta0;

                double mass =
                    rho*Vj;

                double dx = x[i]-x[j];
                double dy = y[i]-y[j];

                double dr = r[i]-r[j];
                double dtheta = theta[i]-theta[j]; 

                double r2 = dx*dx+dy*dy;

                double s = sqrt(r2);
                double q = s/hp;              

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

                result = gaussian(q, hp, dirX, dirY);
                PGauss[i] += result.Weight*Vj;
                dPGaussX[i] += result.dWeightX*Vj;
                dPGaussY[i] += result.dWeightY*Vj;
                drhodtGauss[i] += mass*(du*result.dWeightX+dv*result.dWeightY);
                
                result = cubicSpline(q, hp, dirX, dirY);
                PCubic[i] += result.Weight*Vj;
                dPCubicX[i] += result.dWeightX*Vj;
                dPCubicY[i] += result.dWeightY*Vj;
                drhodtCubic[i] += mass*(du*result.dWeightX+dv*result.dWeightY);

                result = Wendland(q, hp, dirX, dirY);
                PWedn[i] += result.Weight*Vj;
                dPWednX[i] += result.dWeightX*Vj;
                dPWednY[i] += result.dWeightY*Vj;
                drhodtWedn[i] += mass*(du*result.dWeightX+dv*result.dWeightY);

        
            }

            L2PGauss += (PGauss[i]-1.0)*(PGauss[i]-1.0);
            L2PCubic += (PCubic[i]-1.0)*(PCubic[i]-1.0);
            L2PWedn += (PWedn[i]-1.0)*(PWedn[i]-1.0); 
            L2drhodtGauss += (drhodtGauss[i]-drhodtexact)*(drhodtGauss[i]-drhodtexact);
            L2drhodtCubic += (drhodtCubic[i]-drhodtexact)*(drhodtCubic[i]-drhodtexact);
            L2drhodtWedn += (drhodtWedn[i]-drhodtexact)*(drhodtWedn[i]-drhodtexact);
            
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
