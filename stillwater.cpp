#include<iostream>
#include<cmath>
#include <string>
#include <iomanip>
#include <fstream>
#include <chrono>
#include <filesystem>
#include <sstream>


using namespace std;

const double dt = 0.0001; 
const double Totaltime= 1.0;
const int Nt = Totaltime / dt;


const double dx = 0.5;
const double dy = 0.5;
const int Ncol = 50;
const int Nrow = 50;
const double VelcoefX = 0.0;
const double VelcoefY = 0.0;

const double PI = 3.14159265358979323846;
const double g = 9.81;
const double rho0 = 1000.0;
const int boundpart = 3; 
const double c0 = 10.0*sqrt(g*dy*(Nrow-1-boundpart));
const double mass = rho0*dx*dy;

const double gamma = 7.0;
const double B = c0*c0*rho0/gamma;





struct KernelResult {
    double Weight;
    double dWeightX;
    double dWeightY;
};


// --------------------------------------------
// Gaussian kernel
// --------------------------------------------
KernelResult gaussian(double q, double h, double dirx, double diry)
{
    KernelResult result;
    double alpha = 1.0 / (PI * h * h);

    result.Weight = alpha*exp(-q*q);
    result.dWeightX = -2.0*alpha*q*exp(-q*q)/h * dirx;
    result.dWeightY = -2.0*alpha*q*exp(-q*q)/h * diry;

    return result;
}


// --------------------------------------------
// Cubic kernel
// --------------------------------------------
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



// --------------------------------------------
// Wendland kernel
// --------------------------------------------
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
        




// --------------------------------------------
// Main Program
// --------------------------------------------

int main ()
{
    string kernel; 
    cout<<"enter the kernel to be used (gaussian, cubic, wendland): ";
    cin>>kernel;

    auto start = std::chrono::high_resolution_clock::now();


    const int Nh = 1;
    double hlist[Nh] = {2*dx}; 

    for (int m = 0; m < Nh; m++)
    {
        double h = hlist[m];
        cout << "h : " << h << endl;
        string foldername = 
            "Stillwater_dx_" + to_string(dx) + "_h_" + to_string(h) + "_Ncol_" + to_string(Ncol) + "_Nrow_" + to_string(Nrow)+ "_" + kernel; 
        //system(("mkdir -p " + foldername).c_str());
        std::filesystem::create_directories(foldername);

        // parameters         

        double x[Ncol][Nrow];
        double y[Ncol][Nrow];

        double xnew[Ncol][Nrow];
        double ynew[Ncol][Nrow];

        double u[Ncol][Nrow];
        double v[Ncol][Nrow];

        double unew[Ncol][Nrow];
        double vnew[Ncol][Nrow];

        double rho[Ncol][Nrow];
        double rhonew[Ncol][Nrow];

        double drhodtexact[Ncol][Nrow];

        double pressure[Ncol][Nrow];
        double pressurenew[Ncol][Nrow];

        double drhodt[Ncol][Nrow];

        double dudt[Ncol][Nrow];
        double dvdt[Ncol][Nrow];

        

        double L2normdrhodt = 0.0;


        for (int n=0;n<Nt+1;n++)
        {
            double L2drhodt = 0.0;

            double t = n*dt;
            
            
            for (int i = 0; i < Ncol; i++)
            for (int j = 0; j < Nrow; j++)
            {
                if (t==0.0)
                {
                    x[i][j] = (i-boundpart)*dx;
                    y[i][j] = (j-boundpart)*dy;
                    u[i][j] = 0.0;
                    v[i][j] = 0.0;
                    //rho [i][j] = rho0*pow(1.0+((rho0*g*(dy*(Nrow-1-boundpart)-y[i][j]))/B), 1.0/gamma);
                    rho [i][j]=rho0;
                                     
                }
                else
                {
                    x[i][j] = xnew[i][j];
                    y[i][j] = ynew[i][j];
                    u[i][j] = unew[i][j];
                    v[i][j] = vnew[i][j];
                    rho [i][j] = rhonew[i][j];
                    
                                       
                }
                pressure[i][j] = B*(pow(rho[i][j]/rho0,gamma)-1.0);
                drhodtexact [i][j] = -rho[i][j]*(VelcoefX+VelcoefY);                
                
            }

            
            for (int i = 0; i < Ncol; i++)
            for (int j = 0; j < Nrow; j++)
            {
                //cummulative should be initialised with zero for each particle

                drhodt[i][j] = 0.0;
                dudt[i][j] = 0.0;
                dvdt[i][j] = 0.0;

                bool boundary = ( i<=boundpart-1 || i >= Ncol -boundpart || j<=boundpart-1 );
                

                for (int k = 0 ; k < Ncol; k++)
                for (int l = 0 ; l < Nrow; l++)
                {

                    double rx = x[i][j]-x[k][l];
                    double ry = y[i][j]-y[k][l];

                    double r = sqrt(rx*rx+ry*ry);
                    double q = r/h;

                    double dirx = 0.0;
                    double diry = 0.0;

                    double du = u[i][j]-u[k][l];
                    double dv = v[i][j]-v[k][l];

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

                    

                    
                    drhodt[i][j] += (du*result.dWeightX+dv*result.dWeightY)*mass;

                    if (!boundary)

                    {
                        dudt[i][j] += -mass*((pressure[i][j]/(rho[i][j]*rho[i][j]))+(pressure[k][l]/(rho[k][l]*rho[k][l])))*result.dWeightX;
                        dvdt[i][j] += -mass*((pressure[i][j]/(rho[i][j]*rho[i][j]))+(pressure[k][l]/(rho[k][l]*rho[k][l])))*result.dWeightY;
                        
                        
                    }
                }

                if (boundary)
                {
                    rhonew[i][j] = rho[i][j]+drhodt[i][j]*dt;
                    unew[i][j] = 0.0;
                    vnew[i][j] = 0.0;
                    xnew[i][j] = x[i][j];
                    ynew[i][j] = y[i][j];
                }
                else
                {

                rhonew[i][j] = rho[i][j]+drhodt[i][j]*dt;
                unew[i][j] = u[i][j]+dudt[i][j]*dt;
                vnew[i][j] = v[i][j]+(dvdt[i][j]-g)*dt;
                xnew[i][j] = x[i][j]+unew[i][j]*dt;
                ynew[i][j] = y[i][j]+vnew[i][j]*dt;

                }

                double errordrhodt = (pow(drhodt[i][j]-drhodtexact[i][j],2));

                L2drhodt += errordrhodt;
                
            }

            L2normdrhodt = sqrt(L2drhodt/(Ncol*Nrow));


            if (n % int(0.01/dt) == 0)
            {
                cout << "t : " << t << endl;
                string filename = foldername + "/Stillwater_dxh_" + to_string(dx/h) + "_t_" + to_string(t) + ".csv";

                ofstream file(filename);
                file << "dx/h :" << "," << dx/h << endl;
                file << "dy/h :" << "," << dy/h << endl;
                file << "L2normdrhodt" << "," << "dx" << "," << "dy" << "," << "Ncol" << "," << "Nrow" << endl;
                file << L2normdrhodt << "," << dx << "," << dy << "," << Ncol << "," << Nrow << endl;
                file << endl;
                file << "t" << "," << t << endl;
                file << "x"<< "," << "y" << ","<< ","<<"rho"<< "," << "drhodt" << "," << "pressure"<< "," << "," << "u" << "," << "v" << "," << "," << "dudt" << "," << "dvdt" << endl;
                for (int i = 0; i < Ncol; i++)
                for (int j = 0; j < Nrow; j++)
                    {
                        file << x[i][j] << "," << y[i][j] << "," << "," << rho[i][j] << "," << drhodt[i][j] << "," << pressure[i][j] << "," << "," << u[i][j] << "," << v[i][j] << "," << "," << dudt[i][j] << "," << dvdt[i][j] << endl;
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