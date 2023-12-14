#include<iostream>
#include<cmath>


using namespace std;
double func_var8(double x) 
{
    return pow(log(x)*log(x) + 1, -1);
}


void trapezoid(double a, double b, int n) 
{
    double h = (b - a) / n;
    double integral =(func_var8(a) + func_var8(b))/2;

    for (int i = 1; i < n; i++) 
    {
        integral += func_var8(a + i * h);
    }

    integral *= h;
    cout << "Trapezoidal rule: " << integral << endl;
}

void quad_simpson(double a, double b, int n) {
    double h = (b - a) / n;
    double integral = func_var8(a) + func_var8(b);

    for (int i = 1; i < n; i++) {
        if (i % 2 != 0)
            integral += 4 * func_var8(a + i * h);
        else
            integral += 2 * func_var8(a + i * h);
    }

    integral *= h / 3;
    cout << "Quadrature Simpson's rule: " << integral << endl;
}

double cube_simpson(double x, double y)
{
    return pow(x * x / (1 + y * y), 2);
}
void cube_simpson(double a, double b, double c, double d, int m, int n)
{
    double hx = (b - a) / (2 * n);
    double hy = (d - c) / m;
    double integral = 0;

    for (int i = 0; i < 2 * n; i += 2)
    {
        for (int j = 0; j < m; j++)
        {
            integral += cube_simpson(a + i * hx, c + j * hy);
            integral += 4 * cube_simpson(a + (i + 1) * hx, c + j * hy);
            integral += cube_simpson(a + (i + 2) * hx, c + j * hy);
            integral += 4 * cube_simpson(a + i * hx, c + (j + 1) * hy);
            integral += 16 * cube_simpson(a + (i + 1) * hx, c + (j + 1) * hy);
            integral += 4 * cube_simpson(a + (i + 2) * hx, c + (j + 1) * hy);
            integral += cube_simpson(a + i * hx, c + (j + 2) * hy);
            integral += 4 * cube_simpson(a + (i + 1) * hx, c + (j + 2) * hy);
            integral += cube_simpson(a + (i + 2) * hx, c + (j + 2) * hy);
        }
    }

    integral *= (hx * hy / 9);
    cout << "Cube Simpson's rule: " << integral << endl;
}

int main() 
{
    const double accuracy = 0.00001;
    double a = 1.0;
    double b = 2.835;
    double c = 1.0;
    double d = 2.0;

    int trapezoid_n = 100; 
    int quad_simpson_n = 100;
    int cube_simpson_m = 100; 
    int cube_simpson_n = 100; 

    trapezoid(a, b, trapezoid_n);
    quad_simpson(a, b, quad_simpson_n);
    cube_simpson(a, b, c, d, cube_simpson_m, cube_simpson_n);

    return 0;
}
