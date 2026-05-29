import numpy as np
import sisl
from scipy.linalg import eig, eigh


def dotprod(coeff1, coeff2):
    if len(coeff1) == 1:
        dp = np.dot(np.conj(coeff1[0,:]), coeff2[0,:])
    else:
        dp_up = np.dot(np.conj(coeff1[0,:]), coeff2[0,:])
        dp_dn = np.dot(np.conj(coeff1[1,:]), coeff2[1,:])
        dp = dp_up + dp_dn
    return dp



hsx = sisl.io.siesta.hsxSileSiesta("WSe2.HSX")

H = hsx.read_hamiltonian()
# print(hsx.read_geometry())


k = [0.66666667, 0.33333333, 0.0]
Hk = np.asarray(H.Hk(k=k, format="array"))
Sk = np.asarray(H.Sk(k=k, format="array"))

kx = [0.66667667, 0.33333833, 0.0]
Hkx = np.asarray(H.Hk(k=kx, format="array"))
Skx = np.asarray(H.Sk(k=kx, format="array"))

ky = [0.66666667, 0.33334199, 0.0]
Hky = np.asarray(H.Hk(k=ky, format="array"))
Sky = np.asarray(H.Sk(k=ky, format="array"))

Ek, ck = eigh(Hk, Sk)
"""print([f"{(x-3.621):.3f}" for x in Ek])
#eigenvalues might have an offset equal to the fermi level in the .dat file"""
Ekx, ckx = eigh(Hkx, Skx)
Eky, cky = eigh(Hky, Sky)


for i in range(44, 48):
    print(Ek[i])

# checking generalized inner product
"""for n in range(48):
    print(np.conj(ck[:, n]).T @ Sk @ ck[:, n+1])"""


for n in range (44, 48):
    phase_x = np.conj(ck[:, n]).T @ (Sk) @ ckx[:, n]
    phase_y = np.conj(ck[:, n]).T @ (Sk) @ cky[:, n]
    phase_x = phase_x/np.abs(phase_x)
    phase_y = phase_y/np.abs(phase_y)

    #print(phase_x); print(phase_y)
    #old comment: got reasonably accurate phases (up to the -7 decimal digit)

    blat = 1.00755 #2pi/alat in XV file
    deltakx = 1e-5*blat
    deltaky = 1e-5*blat

    du_dkx = (np.conj(phase_x)*ckx[:, n] - ck[:, n])/deltakx
    du_dky = (np.conj(phase_y)*cky[:, n] - ck[:, n])/deltaky

    #berry_1 = np.dot(np.conj(du_dkx), du_dky)
    #berry_2 = np.dot(np.conj(du_dky), du_dkx)
    berry_1 = np.conj(du_dkx).T @ (Sk) @ du_dky
    berry_2 = np.conj(du_dky).T @ (Sk) @ du_dkx  
    print("Berry curvature of state:", n, -1j*(berry_1 - berry_2)*0.52918**2)

    rhs_x = Hk @ du_dkx
    rhs_x2 = Ek[n] * ((Sk) @ du_dkx) # probably need to include this in the generalized case?
    rhs_x = rhs_x - rhs_x2 #(H-E)*du_dkx
    rhs_y = Hk @ du_dky
    rhs_y2 = Ek[n] * ((Sk) @ du_dky)
    rhs_y = rhs_y - rhs_y2 #(H-E)*du_dky


    fin1 = du_dkx[np.newaxis, :]
    fin2 = du_dky[np.newaxis, :]
    rhs1 = rhs_x[:, np.newaxis]
    rhs2 = rhs_y[:, np.newaxis]

    swich_xy = fin1.conj() @ rhs2
    swich_yx = fin2.conj() @ rhs1

    swich = swich_xy - swich_yx
    print("m_z:", 1j*swich[0][0]/(2*13.605698))


