import numpy as np
from time import time

n = 300
x = np.ones((n, n, n),  dtype="float64")


def row(mat):
    out = 0
    for i in range(n):
        sub = mat[i, :, :]
        for i in range(n):
            out += sub[i, :]
        for i in range(n):
            out += sub[:, i]
    return out


def col(mat):
    out = 0
    for i in range(n):
        sub = mat[:, i, :]
        for i in range(n):
            out += sub[i, :]
        for i in range(n):
            out += sub[:, i]

    return out


def col2(mat):
    out = 0
    for i in range(n):
        sub = mat[:, :, i]
        for i in range(n):
            out += sub[i, :]
        for i in range(n):
            out += sub[:, i]

    return out


p = 10
t = time()
for i in range(p):
    s = row(x)
print(time()-t)


t = time()
for i in range(p):
    s = col(x)
print(time()-t)


t = time()
for i in range(p):
    s = col2(x)
print(time()-t)
