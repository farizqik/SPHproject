#include<iostream>
#include<cmath>
#include <string>
#include <iomanip>
#include <fstream>


using namespace std;


int main ()
{
    
    // parameters
    const int N = 401;
    double h = 1.0;
    double dx = 4*h/(N-1);
    const double PI = 3.14159265358979323846;
    double alphagauss = 1.0 / (sqrt(PI) * h);
    double alphacubic = 2.0 / (3.0 * h);
    double alphawedn = 3.0 / (4.0*h);;
    double x[N];
    double PGauss[N];
    double PCubic[N];
    double PWedn[N];
    string kernel;


    for (int i = 0; i < N; i++)
    {
        x[i] =  i * dx;
    }

    for (int i = 0; i < N; i++)
    {
        PGauss[i] = 0.0;
        PCubic[i] = 0.0;
        PWedn[i] = 0.0;
    
        for (int j = 0; j < N; j++)
        {
            double q = abs(x[i]-x[j])/h;
            double WGauss = 0.0;
            double WCubic = 0.0;
            double WWedn = 0.0;
            
            WGauss = alphagauss*exp(-q*q);

            if (q<=1.0)
            {
                WCubic = alphacubic*(1.0-1.5*q*q+0.75*q*q*q);
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
            
                
            PGauss[i]+=WGauss*dx;
            PCubic[i]+=WCubic*dx;
            PWedn[i]+=WWedn*dx;
        }
    }
                   

    for (int i = 0; i < N; i++)
    {
        cout<<"x = "<<x[i]<<" ,PGauss_i = "<<PGauss[i]<<" ,PCubic_i = "<<PCubic[i]<<" ,PWedn_i = "<<PWedn[i]<< endl;
    }

    ofstream file("partuni1D.xls");

    file<<"x"<<"\t"<<"PGauss"<<"\t"<<"PCubic"<<"\t"<<"PWedn"<<endl;
    for (int i = 0; i<N; i++)
    {
        file<<x[i]<<"\t"<<PGauss[i]<<"\t"<<PCubic[i]<<"\t"<<PWedn[i]<<endl;
    }
    file.close();

}

