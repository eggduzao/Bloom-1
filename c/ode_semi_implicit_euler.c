#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define order 2

void problem(const double *x, double *y, double *dy)
{
    const double omega = 1.F;
    dy[0] = y[1];
    dy[1] = -omega * omega * y[0];
}

void exact_solution(const double *x, double *y)
{
    y[0] = cos(x[0]);
    y[1] = -sin(x[0]);
}

void semi_implicit_euler_step(double dx, double *x, double *y, double *dy)
{
    int o;

    problem(x, y, dy);
    y[0] += dx * dy[0];

    problem(x, y, dy);

    for (o = 1; o < order; o++)
        y[o] += dx * dy[o];
    *x += dx;
}

double semi_implicit_euler(double dx, double x0, double x_max, double *y,
                           char save_to_file)
{
    double dy[order];

    FILE *fp = NULL;
    if (save_to_file)
    {
        fp = fopen("semi_implicit_euler.csv", "w+");
        if (fp == NULL)
        {
            perror("Error! ");
            return -1;
        }
    }

    clock_t t1 = clock();
    double x = x0;
    do
    {
        if (save_to_file && fp)
            fprintf(fp, "%.4g,%.4g,%.4g\n", x, y[0], y[1]);
        semi_implicit_euler_step(dx, &x, y, dy);
        x += dx;
    } while (x <= x_max);

    clock_t t2 = clock();

    if (save_to_file && fp)
        fclose(fp);

    return (double)(t2 - t1) / CLOCKS_PER_SEC;
}

