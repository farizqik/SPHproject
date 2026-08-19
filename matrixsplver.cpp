#include <iostream>

using namespace std;

int main()
{
    const int n = 3;

    double A[n][n] =
    {
        { 2,  1, -1},
        {-3, -1,  2},
        {-2,  1,  2}
    };

    double b[n] = {8, -11, -3};
    double x[n];

    // Forward elimination
    for (int k = 0; k < n - 1; k++)
    {
        for (int i = k + 1; i < n; i++)
        {
            double factor = A[i][k] / A[k][k];

            for (int j = k; j < n; j++)
            {
                A[i][j] = A[i][j] - factor * A[k][j];
            }

            b[i] = b[i] - factor * b[k];
        }
    }

    // Back substitution
    for (int i = n - 1; i >= 0; i--)
    {
        double sum = b[i];

        for (int j = i + 1; j < n; j++)
        {
            sum = sum - A[i][j] * x[j];
        }

        x[i] = sum / A[i][i];
    }

    // Print result
    for (int i = 0; i < n; i++)
    {
        cout << "x[" << i << "] = " << x[i] << endl;
    }

    return 0;
}