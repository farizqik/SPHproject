#include<iostream>
#include<cmath>
#include <string>
#include <iomanip>
#include <fstream>
#include <chrono>


using namespace std;


const double PI = 3.14159265358979323846;

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


int main ()
{
    auto start = std::chrono::high_resolution_clock::now();
    const int Nh = 4;
    double hlist[Nh] = {0.12, 0.15, 0.18, 0.2};   

    for (int m = 0; m < Nh; m++)
    {
        double h = hlist[m];
        cout << "h : " << h << endl;

        // parameters
        
        const double dx = 0.1;
        const double dy = 0.1;
        const int Ncol = 50;
        const int Nrow = 50;
        const double VelcoefX = 0.1;
        const double VelcoefY = 0.0;

        const double PI = 3.14159265358979323846;

        const double rho = 1000.0;

        double drhodtexact = -rho*(VelcoefX+VelcoefY);

        double x[Ncol][Nrow];
        double y[Ncol][Nrow];

        double u[Ncol][Nrow];
        double v[Ncol][Nrow];

        double PGauss[Ncol][Nrow];
        double PCubic[Ncol][Nrow];
        double PWedn[Ncol][Nrow];

        double dPGaussX[Ncol][Nrow];
        double dPGaussY[Ncol][Nrow];

        double dPCubicX[Ncol][Nrow];
        double dPCubicY[Ncol][Nrow];

        double dPWednX[Ncol][Nrow];
        double dPWednY[Ncol][Nrow];

        double drhodtGauss[Ncol][Nrow];
        double drhodtCubic[Ncol][Nrow];
        double drhodtWedn[Ncol][Nrow];

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
        
        for (int i = 0; i < Ncol; i++)
        for (int j = 0; j < Nrow; j++)
        {
            x[i][j] =  i * dx;
            y[i][j] =  j * dy;
            u[i][j] =  VelcoefX*x[i][j];
            v[i][j] =  VelcoefY*y[i][j];
        }

        
        for (int i = 0; i < Ncol; i++)
        for (int j = 0; j < Nrow; j++)
        {
            //cummulative should be initialised with zero for each particle
            PGauss[i][j] = 0.0;
            PCubic[i][j] = 0.0;
            PWedn[i][j] = 0.0;
            dPGaussX[i][j] = 0.0;
            dPGaussY[i][j] = 0.0; 
            dPCubicX[i][j] = 0.0;
            dPCubicY[i][j] = 0.0;
            dPWednX[i][j] = 0.0;
            dPWednY[i][j] = 0.0;
            drhodtGauss[i][j] = 0.0;
            drhodtCubic[i][j] = 0.0;
            drhodtWedn[i][j] = 0.0;


            for (int k = 0 ; k < Ncol; k++)
            for (int l = 0 ; l < Nrow; l++)
            {

                double rx = x[i][j]-x[k][l];
                double ry = y[i][j]-y[k][l];

                double r2 = rx*rx+ry*ry;

                double r = sqrt(r2);
                double q = r/h;              

                double du = u[i][j]-u[k][l];
                double dv = v[i][j]-v[k][l];

                double dirx = 0.0;
                double diry = 0.0;
                if (r2>0.0)
                {
                    dirx = rx/r;
                    diry = ry/r;
                }

                KernelResult result;

                result = gaussian(q, h, dirx, diry);
                PGauss[i][j] += result.Weight*dx*dy;
                dPGaussX[i][j] += result.dWeightX*dx*dy;
                dPGaussY[i][j] += result.dWeightY*dx*dy;
                drhodtGauss[i][j] += (du*result.dWeightX+dv*result.dWeightY)*dx*dy*rho;
                
                result = cubicSpline(q, h, dirx, diry);
                PCubic[i][j] += result.Weight*dx*dy;
                dPCubicX[i][j] += result.dWeightX*dx*dy;
                dPCubicY[i][j] += result.dWeightY*dx*dy;
                drhodtCubic[i][j] += (du*result.dWeightX+dv*result.dWeightY)*dx*dy*rho;

                result = Wendland(q, h, dirx, diry);
                PWedn[i][j] += result.Weight*dx*dy;
                dPWednX[i][j] += result.dWeightX*dx*dy;
                dPWednY[i][j] += result.dWeightY*dx*dy;
                drhodtWedn[i][j] += (du*result.dWeightX+dv*result.dWeightY)*dx*dy*rho;

        
            }
            L2PGauss += (PGauss[i][j]-1.0)*(PGauss[i][j]-1.0);
            L2PCubic += (PCubic[i][j]-1.0)*(PCubic[i][j]-1.0);
            L2PWedn += (PWedn[i][j]-1.0)*(PWedn[i][j]-1.0); 
            L2drhodtGauss += (drhodtGauss[i][j]-drhodtexact)*(drhodtGauss[i][j]-drhodtexact);
            L2drhodtCubic += (drhodtCubic[i][j]-drhodtexact)*(drhodtCubic[i][j]-drhodtexact);
            L2drhodtWedn += (drhodtWedn[i][j]-drhodtexact)*(drhodtWedn[i][j]-drhodtexact);
            
        }
        L2normPGauss = sqrt(L2PGauss/(Ncol*Nrow));
        L2normPCubic = sqrt(L2PCubic/(Ncol*Nrow));
        L2normPWedn = sqrt(L2PWedn/(Ncol*Nrow));

        L2normdrhodtGauss = sqrt(L2drhodtGauss/(Ncol*Nrow));
        L2normdrhodtCubic = sqrt(L2drhodtCubic/(Ncol*Nrow));
        L2normdrhodtWedn = sqrt(L2drhodtWedn/(Ncol*Nrow));


        string filename = "2DPartition_hdp_" + to_string(h/dx) + ".csv";

        ofstream file(filename);
        file << "h/dx :" << "," << h/dx << endl;
        file << "h/dy :" << "," << h/dy << endl;
        file << "L2normPGauss" << "," << "L2normPCubic" << "," << "L2normPWedn" << "," << "L2normdrhodtGauss" << "," << "L2normdrhodtCubic" << "," << "L2normdrhodtWedn" << endl;
        file << L2normPGauss << "," << L2normPCubic << "," << L2normPWedn << "," << L2normdrhodtGauss << "," << L2normdrhodtCubic << "," << L2normdrhodtWedn << endl;
        file << endl;
        file << "x"<< "," << "y" << ","<< "," << "PGauss" << "," << "PCubic" << "," << "PWedn" << "," << ","<< "dPGaussX" << "," << "dPGaussY" << "," << ","<< "dPCubicX" << "," << "dPCubicY" << ","<< "," << "dPWednX" << "," << "dPWednY" << "," << "," << "drhodtGauss" << "," << "drhodtCubic" << "," << "drhodtWedn" << endl;
        for (int i = 0; i < Ncol; i++)
        for (int j = 0; j < Nrow; j++)
            {
                file << x[i][j] << "," << y[i][j] << "," << "," << PGauss[i][j] << "," << PCubic[i][j] << "," << PWedn[i][j] << "," << "," << dPGaussX[i][j] << "," << dPGaussY[i][j] << "," << "," << dPCubicX[i][j] << "," << dPCubicY[i][j] << "," <<"," << dPWednX[i][j] << "," << dPWednY[i][j]<< "," << "," << drhodtGauss[i][j] << "," << drhodtCubic[i][j] << "," << drhodtWedn[i][j] << endl;
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
    
        
