#include<iostream>
#include<cmath>

using namespace std;


int main ()
{
    
    // parameters
    const int N = 10;
    double h = 1.0;
    double dx = 1.0;
    double alpha = 2.0/(3.0*h);
    double x[N];
    double q [N][N];
    string kernel;
    double  W[N][N];



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
            for (int j = 0; j < N; j++)
            {
                W[i][j] = alpha*exp(-q[i][j]*q[i][j]);
            }
        }
    }
    if (kernel == "Cubic")
    {
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                if (q[i][j]<=1.0)
                {
                    W[i][j] = alpha*(1.0-1.5*pow(q[i][j], 2)+0.75*pow(q[i][j], 3));
                }
                if (q[i][j]>1.0 && q[i][j]<=2.0)
                {
                    W[i][j] = alpha*0.25*pow((2-q[i][j]), 3);
                }
                else
                {
                    W[i][j] = 0.0;
                }
            }
        }
    }



    for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    {
        cout<<"distances = "<< abs(x[i]-x[j])<<",     q = "<<q[i][j]<<",     x_i = "<<x[i]<<",     x_j = "<<x[j]<<",    W_ij = " << W[i][j]
            <<endl;

    }
}


