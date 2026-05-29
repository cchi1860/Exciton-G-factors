import numpy as np
import h5py
import math

def dotprod(coeff1, coeff2):
    # Check spinor
    if len(coeff1) == 1:
        dp = np.dot(np.conj(coeff1[0,:]), coeff2[0,:])
    else:
        dp_up = np.dot(np.conj(coeff1[0,:]), coeff2[0,:])
        dp_dn = np.dot(np.conj(coeff1[1,:]), coeff2[1,:])
        dp = dp_up + dp_dn
    return dp

gvecs = np.load('Gvecs_k1.npy')
gvecs2 = np.load('Gvecs_k2.npy')
gvecs3 = np.load('Gvecs_k3.npy')
gvecs4 = np.load('Gvecs_k4.npy')
gvecs5 = np.load('Gvecs_k5.npy')
unk = np.load('Unkg_k1.npy')
unk2 = np.load('Unkg_k2.npy')
unk3 = np.load('Unkg_k3.npy')
unk4 = np.load('Unkg_k4.npy')
unk5 = np.load('Unkg_k5.npy')

sum = np.sum(np.dot(unk[0, :], np.conj(unk[0, :])))
print(sum)

unk = unk.reshape(110, 2, unk.shape[1]//2)
unk2 = unk2.reshape(110, 2, unk2.shape[1]//2)
unk3 = unk3.reshape(110, 2, unk3.shape[1]//2)
unk4 = unk4.reshape(110, 2, unk4.shape[1]//2)
unk5 = unk5.reshape(110, 2, unk5.shape[1]//2)
c_kb = np.array([unk, unk2, unk5])


print("Coefficients shape: " + str(c_kb.shape))
print("G-vectors shape: " + str(gvecs.shape))
# firat half is all of the same spin

with open("enk.txt", "r") as f:
    el = [line.strip().split() for line in f]

blat = 1.00755 #2pi/alat in XV file
deltaky = 1e-5*blat
deltakx = 1e-5*blat

print(dotprod(c_kb[0,0], c_kb[1,0]), np.sum(np.abs(c_kb[0,0])**2))
print(dotprod(c_kb[0,44], c_kb[1,44]))

phase_y = dotprod(c_kb[0,0], c_kb[1,0])
phase_y = phase_y/np.abs(phase_y)
print("phase: " + str(phase_y))

for i in range(0, 50):
	print(abs(dotprod(c_kb[0,i], c_kb[1,i])))

"""print("Eigs for k-point")
for i in range(80):
     number = float(el[0][i])
     print(f"{number:.3f}" + " for point " + str(i+1))"""

"""
for i in range(44, 48):
	phase_y = dotprod(c_kb[0,i], c_kb[2,i])
	phase_y = phase_y/np.abs(phase_y)
	phase_x = dotprod(c_kb[0,i], c_kb[1,i])
	phase_x = phase_x/np.abs(phase_x)
	du_dkx = (np.conj(phase_x)*c_kb[1,i] - c_kb[0,i])/deltakx
	du_dky = (np.conj(phase_y)*c_kb[2,i] - c_kb[0,i])/deltaky
	berry_1 = dotprod(du_dkx, du_dky)
	berry_2 = dotprod(du_dky, du_dkx)

	print("Berry curvature of state:", i, -1j*(berry_1 - berry_2)*.52918**2)
	#print(el[0,3,i+24]*13.605698, el[0,1,25]*13.605698)

	swich_xy = 0j
	swich_yx = 0j
	# Check for the completeness of basis
	# by setting H-E = 1
	swich_xy_berry = 0j
	swich_yx_berry = 0j
	for n in range(105): #can change this value
		dp1 = dotprod(du_dkx, c_kb[0,n])
		dp2 = dotprod(c_kb[0,n], du_dky)
		ediff = (float(el[0][n]) - float(el[0][i]))/13.605698
		swich_xy += ediff*dp1*dp2
		swich_xy_berry += 1*dp1*dp2
    
		dp1 = dotprod(du_dky, c_kb[0,n])
		dp2 = dotprod(c_kb[0,n], du_dkx)
		swich_yx += ediff*dp1*dp2
		swich_yx_berry += 1*dp1*dp2
		
	swich = swich_yx - swich_xy
	swich_berry = swich_yx_berry - swich_xy_berry
	print("m_z:", -1j*swich/2) #, "berry:", -1j*swich_berry*(bohr2A**2))"""