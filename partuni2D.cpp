#include<iostream>
#include<cmath>
#include <string>
#include <iomanip>
#include <fstream>


using namespace std;


int main ()
{
    
    // parameters
    
    const double h = 0.2;
    const double dx = 0.1;
    const double dy = 0.1;
    /*const int N = (4*h/dx)+1;*/
    const int Ncol = 101;
    const int Nrow = 101;


    const double PI = 3.14159265358979323846;
    double alphagauss = 1.0 / (PI * h * h);
    double alphacubic = 10.0 / (7.0 * PI * h * h);
    double alphawedn = 7.0 / (4.0 * PI * h * h);
    double x[Ncol][Nrow];
    double y[Ncol][Nrow];
    double PGauss[Ncol][Nrow];
    double PCubic[Ncol][Nrow];
    double PWedn[Ncol][Nrow];
    double dPGauss[Ncol][Nrow];
    double dPCubic[Ncol][Nrow];
    double dPWedn[Ncol][Nrow];
    string kernel;

    for (int i = 0; i < Ncol; i++)
    for (int j = 0; j < Nrow; j++)
    {
        x[i][j] =  i * dx;
        y[i][j] =  j * dy;
    }

    for (int i = 0; i < Ncol; i++)
    for (int j = 0; j < Nrow; j++)
    {
        PGauss[i][j] = 0.0;
        PCubic[i][j] = 0.0;
        PWedn[i][j] = 0.0;
        dPGauss[i][j] = 0.0;
        dPCubic[i][j] = 0.0;
        dPWedn[i][j] = 0.0;

        for (int k = 0 ; k < Ncol; k++)
        for (int l = 0 ; l < Nrow; l++)
        {
            double q = sqrt(pow(x[i][j]-x[k][l],2)+pow(y[i][j]-y[k][l],2))/h;
            double dirx=(x[i][j]-x[k][l])/(abs(x[i][j]-x[k][l])+1e-10);
            double diry=(y[i][j]-y[k][l])/(abs(y[i][j]-y[k][l])+1e-10);

            double WGauss = 0.0;
            double WCubic = 0.0;
            double WWedn = 0.0;
            double dWGauss = 0.0;
            double dWCubic = 0.0;
            double dWWedn = 0.0;


            WGauss = alphagauss*exp(-q*q);
            dWGauss = -2.0*alphagauss*q*exp(-q*q)/h;

            
            if (q<=1.0)
            {
                WCubic = alphacubic*(1.0-1.5*pow(q,2)+0.75*pow(q,3));
                dWCubic = alphacubic*(-3.0*q+2.25*pow(q,2))/h;           
            }
            else if (q<=2.0)
            {
                WCubic = alphacubic*0.25*pow(2.0-q,3);
                dWCubic = -0.75*alphacubic*pow(2.0-q,2)/h;
            }
            else
            {
                WCubic = 0.0;
                dWCubic = 0.0;
            }


            if (0.0<=q && q<=2.0)
            {
                WWedn = alphawedn*pow(1.0-0.5*q,4)*(2.0*q+1.0);
                dWWedn = alphawedn*(4*pow(1.0-0.5*q,3)*(-0.5)*(2.0*q+1.0)+pow(1.0-0.5*q,4)*2.0)/h;
            }
            else
            {
                WWedn = 0.0;
                dWWedn = 0.0;
            }

            PGauss[i][j] += WGauss*dx*dy;
            PCubic[i][j] += WCubic*dx*dy;
            PWedn[i][j] += WWedn*dx*dy;
            dPGauss[i][j] += dWGauss*dx*dy*dirx*diry;
            dPCubic[i][j] += dWCubic*dx*dy*dirx*diry;
            dPWedn[i][j] += dWWedn*dx*dy*dirx*diry;
        }
    }

    for (int i = 0; i < Ncol; i++)
    for (int j = 0; j < Nrow; j++)
        {
            cout << "x: " << x[i][j] << " y: " << y[i][j] << " PGauss: " << PGauss[i][j] << " PCubic: " << PCubic[i][j] << " PWedn: " << PWedn[i][j] << " dPGauss: " << dPGauss[i][j] << " dPCubic: " << dPCubic[i][j] << " dPWedn: " << dPWedn[i][j] <<endl;
        }

    ofstream file("partuni2D.xls");

    file << "x"<< "\t" << "y" << "\t" << "PGauss" << "\t" << "PCubic" << "\t" << "PWedn" << "\t" << "dPGauss" << "\t" << "dPCubic" << "\t" << "dPWedn" <<  endl;
    for (int i = 0; i < Ncol; i++)
    for (int j = 0; j < Nrow; j++)
        {
            file << x[i][j] << "\t" << y[i][j] << "\t" << PGauss[i][j] << "\t" << PCubic[i][j] << "\t" << PWedn[i][j] <<"\t" << dPGauss[i][j] << "\t" << dPCubic[i][j] << "\t" << dPWedn[i][j] << endl;
        }
    file.close();
}
    
        
