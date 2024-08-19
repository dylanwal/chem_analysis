
import numpy as np

data = np.loadtxt(r"C:\Users\nicep\Desktop\dynamic_poly\data\DW2-15\DW2-15-SEC-RI.csv", delimiter=',', skiprows=2)

x = data[:, 0]
y = data[:, 1::2]
np.savez(r"C:\Users\nicep\Desktop\dynamic_poly\data\DW2-15\DW2-15-SEC-RI.npz", x=x, y=y)
