#include<iostream>
#include<cmath>

using namespace std;


int main ()
{
    
    // parameters
    const int N = 101;
    double h = 0.5;
    double dx = 0.01;
    const double PI = 3.14159265358979323846;
    double alphagauss = 1.0 / (sqrt(PI) * h);
    double alphacubic = 2.0 / (3.0 * h);
    double x[N];
    double q [N][N];
    string kernel;
    double  W[N][N];
    double P[N];



    cout<<"Kernel = "<<endl;
    cin>>kernel;

    for (int i = 0; i < N; i++)
    {
        x[i] = i*dx;
        for (int j = 0; j < N; j++)
        {
            x[j] = j*dx;
            q[i][j] = abs(x[i]-x[j])/h;
        }
    }

    if (kernel == "Gaussian")
    {
        for (int i = 0; i < N; i++)
        {
            P[i]=0.0;
            for (int j = 0; j < N; j++)
            {
                W[i][j] = alphagauss*exp(-q[i][j]*q[i][j]);
                
                P[i]+=W[i][j]*dx;
            }
        }
    }
    if (kernel == "Cubic")
    {
        for (int i = 0; i < N; i++)
        {
            P[i]=0.0;
            for (int j = 0; j < N; j++)
            {
                if (q[i][j]<=1.0)
                {
                    W[i][j] = alphacubic*(1.0-1.5*pow(q[i][j], 2)+0.75*pow(q[i][j], 3));
                }
                else if (q[i][j]>1.0 && q[i][j]<=2.0)
                {
                    W[i][j] = alphacubic*0.25*pow((2-q[i][j]), 3);
                }
                else
                {
                    W[i][j] = 0.0;
                }
                
                P[i]+=W[i][j]*dx;
            }
        }
    }



    /*for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    {
        cout<<"distances = "<< abs(x[i]-x[j])<<",     q_ij = "<<q[i][j]<<",     x_i = "<<x[i]<<",     x_j = "<<x[j]<<",    W_ij = " << W[i][j]
            <<endl;

    }*/

    for (int i = 0; i < N; i++)
    {
        cout<<"P_i = "<<P[i]<<endl;
    }

}



