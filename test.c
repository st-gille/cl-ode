#include <math.h>
#include <stdio.h>

#include "newton.h"

void f(const double *x, double *res)
{
    *res = 1 - *x * exp(*x);
    res[1] = 2 - x[1] * x[1];
}

void df(const double *x, double **res)
{
    res[0][0] = - *x * exp(*x) - exp(*x);
    res[1][1] = - 2 * x[1];
    res[0][1] = res[1][0] = 0.0;
}
int impl()

{
    size_t dim = 2;
    double x[] = {0.5, 1.0};
    double r[] = {1.5, 0.5};
    if(newtons_method(dim, x, f, df) == FAILURE)
        printf("failed!\n");
    else
    {
        f(x, r);
        printf("f(x) = ");
        print_vector(dim, r);
    }
}

void f2(const double *y, const double *x, double *res)
{
    *res = 1 - *x * exp(*y);
}

void df2(const double *y, const double *x, double **res)
{
    **res = - *x * exp(*y);
}
int impl2()
{
    double y[] = {0.6};
    double x[] = {1.5};
    double r[] = {1.5};
    print_vector(1, y);
    if(newtons_method_implicit(1, y, x, f2, df2) == FAILURE)
        printf("failed!\n");
    else
    {

        printf("y = ");
        print_vector(1, y);
        printf("x = ");
        print_vector(1, x);
        f2(y, x, r);
        printf("f(x, y) = ");
        print_vector(1, r);
    }
}


int main()
{
    impl();
    //impl2();
    return 1;
}
