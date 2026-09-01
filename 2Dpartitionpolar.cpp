#include<iostream>
#include<cmath>
#include <string>
#include <iomanip>
#include <fstream>
#include <chrono>
#include <vector>


using namespace std;


const double VelcoefR = 0.1;
const double VelcoefZ = 0.0;

const double PI = 3.14159265358979323846;

const double rho = 1000.0;

const double dp = 0.1;

const double width = 5.0;
const double height = 5.0;


// ------------------------------------------------------------
// Particle properties
// ------------------------------------------------------------
vector<double> r;
vector<double> z;
vector<double> u_r;
vector<double> u_z;
vector<double> PGauss;
vector<double> PCubic;
vector<double> PWedn;
vector<double> dPGaussR;
vector<double> dPGaussZ;
vector<double> dPCubicR;
vector<double> dPCubicZ;
vector<double> dPWednR;
vector<double> dPWednZ;
vector<double> drhodtGauss;
vector<double> drhodtCubic;
vector<double> drhodtWedn;


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


int main ()
{
    int i = 0;
    for (double zp = 0.5*dp; zp < height; zp+=dp)
    {
        for (double rp = 0.5*dp; rp < width; rp+=dp)
        {
            cout << "ID: " << i            
             << "  r: " << rp
             << "  z: " << zp
             << "  Type: Fluid"
             << endl;

            i++;   
        }
    }
    
    
    // Total number of particles
    int Nparticles = i;
    cout << "Total particles    = " << Nparticles << endl;


    


    auto start = std::chrono::high_resolution_clock::now();
    const int Nh = 4;
    double hlist[Nh] = {0.12, 0.15, 0.18, 0.2};   

    for (int m = 0; m < Nh; m++)
    {
        double h = hlist[m];
        cout << "h : " << h << endl;

        // parameters

        r.resize(Nparticles);
        z.resize(Nparticles);

        u_r.resize(Nparticles);
        u_z.resize(Nparticles);

        PGauss.resize(Nparticles);
        PCubic.resize(Nparticles);
        PWedn.resize(Nparticles);

        dPGaussR.resize(Nparticles);
        dPGaussZ.resize(Nparticles);

        dPCubicR.resize(Nparticles);
        dPCubicZ.resize(Nparticles);

        dPWednR.resize(Nparticles);
        dPWednZ.resize(Nparticles);

        drhodtGauss.resize(Nparticles);
        drhodtCubic.resize(Nparticles);
        drhodtWedn.resize(Nparticles);

        

        double drhodtexact = -rho*(2*VelcoefR+VelcoefZ);


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
        
        string kernel;

        int i = 0;   
        for (double zp = 0.5*dp; zp < height; zp+=dp)
        {
            for (double rp = 0.5*dp; rp < width; rp+=dp)
            {
                r[i] =  rp;
                z[i] =  zp;
                u_r[i] =  VelcoefR*r[i];
                u_z[i] =  VelcoefZ*z[i];

                i++;   
            }
        }
        
        for (int i = 0; i < Nparticles; i++)

        {
            //cummulative should be initialised with zero for each particle
            PGauss[i] = 0.0;
            PCubic[i] = 0.0;
            PWedn[i] = 0.0;
            dPGaussR[i] = 0.0;
            dPGaussZ[i] = 0.0; 
            dPCubicR[i] = 0.0;
            dPCubicZ[i] = 0.0;
            dPWednR[i] = 0.0;
            dPWednZ[i] = 0.0;
            drhodtGauss[i] = 0.0;
            drhodtCubic[i] = 0.0;
            drhodtWedn[i] = 0.0;


            for (int j = 0 ; j < Nparticles; j++)

            {

                double sr = r[i]-r[j];
                double sz = z[i]-z[j];

                double r2 = sr*sr+sz*sz;

                double s = sqrt(r2);
                double q = s/h;              

                double du_r = u_r[i]-u_r[j];
                double du_z = u_z[i]-u_z[j];

                double dirR = 0.0;
                double dirZ = 0.0;
                double mass = 2*PI*r[j]*dp*dp*rho;

                if (r2>0.0)
                {
                    dirR = sr/s;
                    dirZ = sz/s;
                }

                KernelResult result;

                result = gaussian(q, h, dirR, dirZ);
                PGauss[i] += result.Weight*dp*dp;
                dPGaussR[i] += result.dWeightR*dp*dp;
                dPGaussZ[i] += result.dWeightZ*dp*dp;
                drhodtGauss[i] += (1.0/(2*PI))*((mass/r[j])*(du_r*result.dWeightR+du_z*result.dWeightZ));
                
                result = cubicSpline(q, h, dirR, dirZ);
                PCubic[i] += result.Weight*dp*dp;
                dPCubicR[i] += result.dWeightR*dp*dp;
                dPCubicZ[i] += result.dWeightZ*dp*dp;
                drhodtCubic[i] += (1.0/(2*PI))*((mass/r[j])*(du_r*result.dWeightR+du_z*result.dWeightZ));

                result = Wendland(q, h, dirR, dirZ);
                PWedn[i] += result.Weight*dp*dp;
                dPWednR[i] += result.dWeightR*dp*dp;
                dPWednZ[i] += result.dWeightZ*dp*dp;
                drhodtWedn[i] += (1.0/(2*PI))*((mass/r[j])*(du_r*result.dWeightR+du_z*result.dWeightZ));

        
            }
            drhodtGauss[i] = -(rho*u_r[i]/r[i])+drhodtGauss[i];
            drhodtCubic[i] = -(rho*u_r[i]/r[i])+drhodtCubic[i];
            drhodtWedn[i] = -(rho*u_r[i]/r[i])+drhodtWedn[i];

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


        string filename = "2DPartitionpolar_hdp_" + to_string(h/dp) + ".csv";

        ofstream file(filename);
        file << "h/dp :" << "," << h/dp << endl;
        file << "L2normPGauss" << "," << "L2normPCubic" << "," << "L2normPWedn" << "," << "L2normdrhodtGauss" << "," << "L2normdrhodtCubic" << "," << "L2normdrhodtWedn" << endl;
        file << L2normPGauss << "," << L2normPCubic << "," << L2normPWedn << "," << L2normdrhodtGauss << "," << L2normdrhodtCubic << "," << L2normdrhodtWedn << endl;
        file << endl;
        file << "id"<< "," <<"r"<< "," << "z" << ","<< "," << "PGauss" << "," << "PCubic" << "," << "PWedn" << "," << ","<< "dPGaussR" << "," << "dPGaussZ" << "," << ","<< "dPCubicR" << "," << "dPCubicZ" << ","<< "," << "dPWednR" << "," << "dPWednZ" << "," << "," << "drhodtGauss" << "," << "drhodtCubic" << "," << "drhodtWedn" << endl;
        for (int i = 0; i < Nparticles; i++)
            {
                file << i << "," << r[i] << "," << z[i] << "," << "," << PGauss[i] << "," << PCubic[i] << "," << PWedn[i] << "," << "," << dPGaussR[i] << "," << dPGaussZ[i] << "," << "," << dPCubicR[i] << "," << dPCubicZ[i] << "," <<"," << dPWednR[i] << "," << dPWednZ[i]<< "," << "," << drhodtGauss[i] << "," << drhodtCubic[i] << "," << drhodtWedn[i] << endl;
            }
        file << endl;
        file << endl;
        file << "#######################################################################################################################################################################################################################################################################################################################################################################################################################################################################################" << endl;
        file << endl;
        file << endl;
        file.close();

    

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
    
        
