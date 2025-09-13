
import time
import numpy as np

def euclidean_distance_matrix(N=1000):
    coords = np.random.rand(N, 3)
    dist_matrix = np.zeros((N, N))

    start = time.perf_counter()

    for i in range(N):
        for j in range(N):
            dx = coords[i][0] - coords[j][0]
            dy = coords[i][1] - coords[j][1]
            dz = coords[i][2] - coords[j][2]
            dist_matrix[i][j] = (dx * dx + dy * dy + dz * dz) ** 0.5

    end = time.perf_counter()
    print(f"Time: {end - start:.2f} seconds")

    return dist_matrix

if __name__ == "__main__":
    euclidean_distance_matrix(50000)

