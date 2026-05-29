import numpy as np
import sisl
from scipy.linalg import eig, eigh


hsx = sisl.io.siesta.hsxSileSiesta("WSe2.HSX")
H = hsx.read_hamiltonian()

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
Ekx, ckx = eigh(Hkx, Skx)
Eky, cky = eigh(Hky, Sky)


Hdk = np.asarray(H.dHk(k = k, format = "array"))
Sdk = np.asarray(H.dSk(k = k, format = "array"))
dHdk_x = Hdk[0]
dSdk_x = Sdk[0]
dHdk_y = Hdk[1]
dSdk_y = Sdk[1]


nbands = len(Ek)
def align_phase(u_ref, u, S):
    phase = np.angle(np.conj(u_ref) @ (S @ u))
    return u * np.exp(-1j * phase)

for n in range(nbands):
    ckx[:, n] = align_phase(ck[:, n], ckx[:, n], Sk)
    cky[:, n] = align_phase(ck[:, n], cky[:, n], Sk)

du_dkx = np.zeros_like(ck, dtype=complex)
du_dky = np.zeros_like(ck, dtype=complex)

for n in range(nbands):
    En = Ek[n]
    un = ck[:, n]

    du_n = np.zeros_like(un, dtype=complex)
    du_n2 = np.zeros_like(un, dtype=complex)

    for m in range(nbands):
        if m == n:
            continue

        Em = Ek[m]
        um = ck[:, m]

        # matrix element
        Mmn = np.conj(um) @ ((dHdk_x - En * dSdk_x) @ un)
        Mmn2 = np.conj(um) @ ((dHdk_y - En * dSdk_y) @ un)

        if abs(En - Em) > 1e-5:
            du_n += um * (Mmn / (En - Em))
            du_n2 += um * (Mmn2 / (En - Em))
        else:
            du_n = 0
            du_n2 = 0

    du_dkx[:, n] = du_n
    du_dky[:, n] = du_n2


print(du_dkx[:, 47])
print(du_dky[:, 47])


"""for n in range(nbands):
    val = np.conj(ck[:, n]) @ (Sk @ du_dkx[:, n])
    print(n, val)  # should be ~0 (gauge condition)"""

for n in range(44, 48):
    En = Ek[n]
    # operator (H - E_n S)
    H_eff = Hk - En * Sk

    dukx = du_dkx[:, n]
    duky = du_dky[:, n]

    # matrix element: <∂kx u | (H - E S) | ∂ky u>
    val = np.conj(dukx) @ (H_eff @ duky)

    m_z = -np.imag(val)
    print("m_z: " + str(m_z))