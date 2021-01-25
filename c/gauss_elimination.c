#include <math.h>
#include <stdio.h>

#define ARRAY_SIZE 20

void display(float a[ARRAY_SIZE][ARRAY_SIZE], int n)
{
    int i, j;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j <= n; j++)
        {
            printf("%.2f \t", a[i][j]);
        }
        printf("\n");
    }
}

float interchange(float m[ARRAY_SIZE][ARRAY_SIZE], int i, int n)
{
    float tmp[ARRAY_SIZE][ARRAY_SIZE];
    float max = fabs(m[i][i]);
    int j, k = i;

    for (j = i; j < n; j++)
    {
        if (max < fabs(m[j][i]))
        {
            max = fabs(m[j][i]);
            k = j;
        }
    }
    for (j = 0; j <= n; j++)
    {
        tmp[i][j] = m[i][j];
        m[i][j] = m[k][j];
        m[k][j] = tmp[i][j];
    }
    return m[ARRAY_SIZE - 1][ARRAY_SIZE - 1];
}
float eliminate(float m[ARRAY_SIZE][ARRAY_SIZE], int i, int n)
{
    float tmp;
    int k = 1, l, j;
    for (j = i; j < n - 1; j++)
    {
        tmp = -((m[i + k][i]) / (m[i][i]));
        for (l = 0; l <= n; l++)
        {
            m[i + k][l] = (m[i + k][l]) + (m[i][l] * tmp);
        }
        k++;
    }
    return m[ARRAY_SIZE - 1][ARRAY_SIZE - 1];
}

