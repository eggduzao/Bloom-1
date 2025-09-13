
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

void euclidean_distance_matrix(int N) {
    double** coords = malloc(N * sizeof(double*));
    double** dist_matrix = malloc(N * sizeof(double*));
    for (int i = 0; i < N; i++) {
        coords[i] = malloc(3 * sizeof(double));
        dist_matrix[i] = malloc(N * sizeof(double));
        for (int j = 0; j < 3; j++) {
            coords[i][j] = (double)rand() / RAND_MAX;
        }
    }

    clock_t start = clock();

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            double dx = coords[i][0] - coords[j][0];
            double dy = coords[i][1] - coords[j][1];
            double dz = coords[i][2] - coords[j][2];
            dist_matrix[i][j] = sqrt(dx * dx + dy * dy + dz * dz);
        }
    }

    clock_t end = clock();
    double time_spent = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Time: %.2f seconds\n", time_spent);

    for (int i = 0; i < N; i++) {
        free(coords[i]);
        free(dist_matrix[i]);
    }
    free(coords);
    free(dist_matrix);
}

int main() {
    euclidean_distance_matrix(100000);
    return 0;
}

// Compilation:
// gcc -O3 pure_c_algorithm.c -o dist_c -lm
// ./dist_c

