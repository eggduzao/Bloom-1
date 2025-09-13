
import numpy as np
cimport numpy as np
cimport cython
from libc.time cimport clock, CLOCKS_PER_SEC, clock_t

@cython.boundscheck(False)
@cython.wraparound(False)
def euclidean_distance_matrix(int N=1000):
    cdef int i, j
    cdef double dx, dy, dz

    cdef double[:, ::1] coords = np.random.rand(N, 3)
    cdef double[:, ::1] dist_matrix = np.zeros((N, N), dtype=np.float64)

    cdef clock_t start = clock()

    for i in range(N):
        for j in range(N):
            dx = coords[i, 0] - coords[j, 0]
            dy = coords[i, 1] - coords[j, 1]
            dz = coords[i, 2] - coords[j, 2]
            dist_matrix[i, j] = (dx * dx + dy * dy + dz * dz) ** 0.5

    cdef clock_t end = clock()
    print(f"Time: {(end - start) / <double>CLOCKS_PER_SEC:.2f} seconds")

    return np.asarray(dist_matrix)

