#include<iostream>
#include<cmath>
#include <string>
#include <iomanip>
#include <fstream>
#include <chrono>


using namespace std;


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
        
        /*const double h = 0.12;*/
        const double dx = 0.1;
        const double dy = 0.1;
        /*const int N = (4*h/dx)+1;*/
        const int Ncol = 60;
        const int Nrow = 60;


        const double PI = 3.14159265358979323846;
        double alphagauss = 1.0 / (PI * h * h);
        double alphacubic = 10.0 / (7.0 * PI * h * h);
        double alphawedn = 7.0 / (4.0 * PI * h * h);
        double x[Ncol][Nrow];
        double y[Ncol][Nrow];
        double PGauss[Ncol][Nrow];
        double PCubic[Ncol][Nrow];
        double PWedn[Ncol][Nrow];
        double dPGaussX[Ncol][Nrow];
        double dPGaussY[Ncol][Nrow];
        double dPCubicX[Ncol][Nrow];
        double dPCubicY[Ncol][Nrow];
        double dPWednX[Ncol][Nrow];
        double dPWednY[Ncol][Nrow];

        double L2PGauss = 0.0;
        double L2PCubic = 0.0;
        double L2PWedn = 0.0;

        double L2normPGauss = 0.0;
        double L2normPCubic = 0.0;
        double L2normPWedn = 0.0;
        
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
            dPGaussX[i][j] = 0.0;
            dPGaussY[i][j] = 0.0;
            dPCubicX[i][j] = 0.0;
            dPCubicY[i][j] = 0.0;
            dPWednX[i][j] = 0.0;
            dPWednY[i][j] = 0.0;

            for (int k = 0 ; k < Ncol; k++)
            for (int l = 0 ; l < Nrow; l++)
            {
                double q = sqrt(pow(x[i][j]-x[k][l],2)+pow(y[i][j]-y[k][l],2))/h;
                double r = sqrt(pow(x[i][j]-x[k][l],2)+pow(y[i][j]-y[k][l],2));
                double rx = x[i][j]-x[k][l];
                double ry = y[i][j]-y[k][l];
                double dirx = 0.0;
                double diry = 0.0;
                if (abs(rx)>0.0)
                {
                    dirx = rx/r;
                }
                
                if (abs(ry)>0.0)
                {
                    diry = ry/r;
                }


                double WGauss = 0.0;
                double WCubic = 0.0;
                double WWedn = 0.0;
                double dWGaussX = 0.0;
                double dWGaussY = 0.0;
                double dWCubicX = 0.0;
                double dWCubicY = 0.0;
                double dWWednX = 0.0;
                double dWWednY = 0.0;


                WGauss = alphagauss*exp(-q*q);
                dWGaussX = -2.0*alphagauss*q*exp(-q*q)/h * dirx;
                dWGaussY = -2.0*alphagauss*q*exp(-q*q)/h * diry;

                
                if (q<=1.0)
                {
                    WCubic = alphacubic*(1.0-1.5*pow(q,2)+0.75*pow(q,3));
                    dWCubicX = alphacubic*(-3.0*q+2.25*pow(q,2))/h * dirx;
                    dWCubicY = alphacubic*(-3.0*q+2.25*pow(q,2))/h * diry;
                }
                else if (q<=2.0)
                {
                    WCubic = alphacubic*0.25*pow(2.0-q,3);
                    dWCubicX = -0.75*alphacubic*pow(2.0-q,2)/h * dirx;
                    dWCubicY = -0.75*alphacubic*pow(2.0-q,2)/h * diry;
                }
                else
                {
                    WCubic = 0.0;
                    dWCubicX = 0.0;
                    dWCubicY = 0.0;
                }


                if (0.0<=q && q<=2.0)
                {
                    WWedn = alphawedn*pow(1.0-0.5*q,4)*(2.0*q+1.0);
                    dWWednX = alphawedn*(4*pow(1.0-0.5*q,3)*(-0.5)*(2.0*q+1.0)+pow(1.0-0.5*q,4)*2.0)/h * dirx;
                    dWWednY = alphawedn*(4*pow(1.0-0.5*q,3)*(-0.5)*(2.0*q+1.0)+pow(1.0-0.5*q,4)*2.0)/h * diry;
                }
                else
                {
                    WWedn = 0.0;
                    dWWednX = 0.0;
                    dWWednY = 0.0;
                }

                PGauss[i][j] += WGauss*dx*dy;
                PCubic[i][j] += WCubic*dx*dy;
                PWedn[i][j] += WWedn*dx*dy;
                dPGaussX[i][j] += dWGaussX*dx*dy;
                dPGaussY[i][j] += dWGaussY*dx*dy;
                dPCubicX[i][j] += dWCubicX*dx*dy;
                dPCubicY[i][j] += dWCubicY*dx*dy;
                dPWednX[i][j] += dWWednX*dx*dy;
                dPWednY[i][j] += dWWednY*dx*dy;

                
        
            }
            double errorPGauss = (pow(PGauss[i][j]-1.0,2));
            double errorPCubic = (pow(PCubic[i][j]-1.0,2));
            double errorPWedn = (pow(PWedn[i][j]-1.0,2));
            L2PGauss += errorPGauss;
            L2PCubic += errorPCubic;
            L2PWedn += errorPWedn;    
                
            
        }
        L2normPGauss = sqrt(L2PGauss/(Ncol*Nrow));
        L2normPCubic = sqrt(L2PCubic/(Ncol*Nrow));
        L2normPWedn = sqrt(L2PWedn/(Ncol*Nrow));



        cout << "x"<< "\t" << "y" << "\t" << "PGauss" << "\t" << "PCubic" << "\t" << "PWedn" << "\t" << "dPGaussX" << "\t" << "dPGaussY" << "\t" << "dPCubicX" << "\t" << "dPCubicY" << "\t" << "dPWednX" << "\t" << "dPWednY" << "\t" << "L2normPGauss" << "\t" << "L2normPCubic" << "\t" << "L2normPWedn" << endl;
        for (int i = 0; i < Ncol; i++)
        for (int j = 0; j < Nrow; j++)
            {
                cout <<  x[i][j] << "\t" << y[i][j] << "\t" << PGauss[i][j] << "\t" << PCubic[i][j] << "\t" << PWedn[i][j] << "\t" << dPGaussX[i][j] << "\t" << dPGaussY[i][j] << "\t" << dPCubicX[i][j] << "\t" << dPCubicY[i][j] << "\t" << dPWednX[i][j] << "\t" << dPWednY[i][j] << "\t" << L2normPGauss << "\t" << L2normPCubic << "\t" << L2normPWedn << endl;
            }


        string filename = "partuni2D_dxh_" + to_string(dx/h) + ".csv";

        ofstream file(filename);
        file << "dx/h :" << "," << dx/h << endl;
        file << "dy/h :" << "," << dy/h << endl;
        file << "L2normPGauss" << "," << "L2normPCubic" << "," << "L2normPWedn" << endl;
        file << L2normPGauss << "," << L2normPCubic << "," << L2normPWedn << endl;
        file << endl;
        file << "x"<< "," << "y" << ","<< "," << "PGauss" << "," << "PCubic" << "," << "PWedn" << "," << ","<< "dPGaussX" << "," << "dPGaussY" << "," << ","<< "dPCubicX" << "," << "dPCubicY" << ","<< "," << "dPWednX" << "," << "dPWednY" <<  endl;
        for (int i = 0; i < Ncol; i++)
        for (int j = 0; j < Nrow; j++)
            {
                file << x[i][j] << "," << y[i][j] << "," << "," << PGauss[i][j] << "," << PCubic[i][j] << "," << PWedn[i][j] << "," << "," << dPGaussX[i][j] << "," << dPGaussY[i][j] << "," << "," << dPCubicX[i][j] << "," << dPCubicY[i][j] << "," <<"," << dPWednX[i][j] << "," << dPWednY[i][j] << endl;
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
    
        
