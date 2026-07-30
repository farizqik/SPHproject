#include<iostream>
#include<cmath>
#include <string>
#include <iomanip>
#include <fstream>


using namespace std;


int main ()
{
    
    // parameters
    
    const double h = 1.0;
    const double dx = 0.05;
    const double dy = 0.05;
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

        for (int k = 0 ; k < Ncol; k++)
        for (int l = 0 ; l < Nrow; l++)
        {
            double q = sqrt(pow(x[i][j]-x[k][l],2)+pow(y[i][j]-y[k][l],2))/h;
            double WGauss = 0.0;
            double WCubic = 0.0;
            double WWedn = 0.0;
            WGauss = alphagauss*exp(-q*q);
            
            if (q<=1.0)
            {
                WCubic = alphacubic*(1.0-1.5*pow(q,2)+0.75*pow(q,3));           
            }
            else if (q<=2.0)
            {
                WCubic = alphacubic*0.25*pow(2.0-q,3);
            }
            else
            {
                WCubic = 0.0;
            }


            if (0.0<=q && q<=2.0)
            {
                WWedn = alphawedn*pow(1.0-0.5*q,4)*(2.0*q+1.0);
            }
            else
            {
                WWedn = 0.0;
            }

            PGauss[i][j] += WGauss*dx*dy;
            PCubic[i][j] += WCubic*dx*dy;
            PWedn[i][j] += WWedn*dx*dy;
        }
    }

    for (int i = 0; i < Ncol; i++)
    for (int j = 0; j < Nrow; j++)
        {
            cout << "x: " << x[i][j] << " y: " << y[i][j] << " PGauss: " << PGauss[i][j] << " PCubic: " << PCubic[i][j] << " PWedn: " << PWedn[i][j] << endl;
        }

    ofstream file("partuni2D.xls");

    file << "x"<< "\t" << "y" << "\t" << "PGauss" << "\t" << "PCubic" << "\t" << "PWedn" << endl;
    for (int i = 0; i < Ncol; i++)
    for (int j = 0; j < Nrow; j++)
        {
            file << x[i][j] << "\t" << y[i][j] << "\t" << PGauss[i][j] << "\t" << PCubic[i][j] << "\t" << PWedn[i][j] << endl;
        }
    file.close();
}
    
        
