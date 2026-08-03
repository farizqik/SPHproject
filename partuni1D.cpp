#include<iostream>

#include <fstream>


using namespace std;


int main ()
{
    
    // parameters
    
    
    double h = 0.2;
    const double dx = 0.05;
    /*const int N = (4*h/dx)+1;*/
    const int N = 401;
    
    const double PI = 3.14159265358979323846;
    double alphagauss = 1.0 / (sqrt(PI) * h);
    double alphacubic = 2.0 / (3.0 * h);
    double alphawedn = 3.0 / (4.0*h);;
    double x[N];
    double PGauss[N];
    double PCubic[N];
    double PWedn[N];
    double dPGauss[N];
    double dPCubic[N]; 
    double dPWedn[N];


    for (int i = 0; i < N; i++)
    {
        x[i] =  i * dx;
    }

    for (int i = 0; i < N; i++)
    {
        PGauss[i] = 0.0;
        PCubic[i] = 0.0;
        PWedn[i] = 0.0;
        dPGauss[i] = 0.0;
        dPCubic[i] = 0.0;
        dPWedn[i] = 0.0;

        for (int j = 0; j < N; j++)
        {
            double q = abs(x[i]-x[j])/h;
            double dir=(x[i]-x[j])/(abs(x[i]-x[j])+1e-10);

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
                WCubic = alphacubic*(1.0-1.5*q*q+0.75*q*q*q);
                dWCubic = alphacubic*(-3.0*q+2.25*pow(q,2))/h;
            }
            else if (q<=2.0)
            {
                WCubic = alphacubic*0.25*pow(2.0-q,3);
                dWCubic = -alphacubic*0.75*pow(2.0-q,2)/h;
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
            
                
            PGauss[i]+=WGauss*dx;
            PCubic[i]+=WCubic*dx;
            PWedn[i]+=WWedn*dx;
            dPGauss[i]+=dWGauss*dx*dir;
            dPCubic[i]+=dWCubic*dx*dir;
            dPWedn[i]+=dWWedn*dx*dir;
        }
    }
                   

    for (int i = 0; i < N; i++)
    {
        cout<<"x = "<<x[i]<<" ,PGauss_i = "<<PGauss[i]<<" ,PCubic_i = "<<PCubic[i]<<" ,PWedn_i = "<<PWedn[i]<<" ,dPGauss_i = "<<dPGauss[i]<<" ,dPCubic_i = "<<dPCubic[i]<<" ,dPWedn_i = "<<dPWedn[i]<< endl;
    }

    ofstream file("partuni1D.xls");

    file<<"x"<<"\t"<<"PGauss"<<"\t"<<"PCubic"<<"\t"<<"PWedn"<<"\t"<<"dPGauss"<<"\t"<<"dPCubic"<<"\t"<<"dPWedn"<<endl;
    for (int i = 0; i<N; i++)
    {
        file<<x[i]<<"\t"<<PGauss[i]<<"\t"<<PCubic[i]<<"\t"<<PWedn[i]<<"\t"<<dPGauss[i]<<"\t"<<dPCubic[i]<<"\t"<<dPWedn[i]<<endl;
    }
    file.close();

}

