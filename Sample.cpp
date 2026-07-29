#include <iostream>
#include <cmath>

using namespace std;

double kernel(double h, double distance)
{
    double q = distance / h;
    double alpha = 2.0 / (3.0 * h);

    if (q <= 1.0)
    {
        return alpha *
               (1.0 - 1.5 * q * q
               + 0.75 * q * q * q);
    }
    else if (q <= 2.0)
    {
        return alpha *
               0.25 * pow(2.0 - q, 3);
    }
    else
    {
        return 0.0;
    }
}

int main()
{
    const int N = 21;

    double x[N];
    double P[N];

    double dx = 0.1;
    double h = 0.2;

    // Create particle positions
    for (int i = 0; i < N; i++)
    {
        x[i] = i * dx;
        P[i] = 0.0;
    }

    // Partition-of-unity summation
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            double distance = abs(x[i] - x[j]);

            double Wij = kernel(h, distance);

            P[i] = P[i] + Wij * dx;
        }
    }

    // Display results
    for (int i = 0; i < N; i++)
    {
        cout << "Particle " << i
             << ", x = " << x[i]
             << ", partition sum = " << P[i]
             << endl;
    }

    return 0;
}