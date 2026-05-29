import numpy as np
import matplotlib.pyplot as plt
import h5py

ib_start, ib_end = 0, 110
f = h5py.File("WFN.h5", 'r')
Eig= f['mf_header/kpoints/el'][:]
print(Eig.shape)

eigenvalues = Eig[0]
eigenvalues = eigenvalues*13.605698 # convert from Hartree to eV
eigenvalues = eigenvalues.transpose()
print(eigenvalues.shape)

print(eigenvalues[:, 0])

efermi = np.max(eigenvalues[46-1, :])
print(efermi)