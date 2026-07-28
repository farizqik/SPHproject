#include <iostream>
#include <iomanip>

using namespace std;

double p = 0.50;
double q = p*4;

// Function to be integrated
double partuni(double x)
{
    return x * x ;
}

int main()
{
    // Lower integration limit
    double a = 0.0;

    // Upper integration limit
    double b = 1.0;

    // Number of small intervals
    int n = 1000;

    // Width of each interval
    double h = (b - a) / n;

    // Variable for storing the total area
    double integral = 0.0;

    // Numerical integration loop
    for (int i = 0; i < n; i++)
    {
        // Left point of each interval
        double x = a + i * h;

        // Area of one rectangle
        double rectangleArea = partuni(x) * h * q;

        // Add rectangle area to the total
        integral = integral + rectangleArea;
    }

    // Display more decimal places
    cout << fixed << setprecision(6);

    cout << "Lower limit, a = " << a << endl;
    cout << "Upper limit, b = " << b << endl;
    cout << "Number of intervals, n = " << n << endl;
    cout << "Interval width, h = " << h << endl;
    cout << "Approximate integral = " << integral << endl;

    return 0;
}