import numpy as np
import h5py

def dunk_dkDotdunpk(unk,unkpdk, unpk, dk):
	phase = dotprod(unk,unkpdk)
	duv_dkx = (phase*unkpdk - unk)/dk
	
def check_G_mismatch(gv1, gv2):
	print("ng:", len(gv1))
	for ig in range(len(gv1)):
		if np.linalg.norm( gv1[ig] - gv2[ig]) > 1e-6:
			print("Mismatch:", gv1[ig], gv2[ig])

def dotprod(coeff1, coeff2):
    add = 0.0
    # Check spinor
    if len(coeff1) == 1:
        dp = np.dot(np.conj(coeff1[0,:]), coeff2[0,:])
    else:
        dp_up = np.dot(np.conj(coeff1[0,:]), coeff2[0,:])
        dp_dn = np.dot(np.conj(coeff1[1,:]), coeff2[1,:])
        dp = dp_up + dp_dn
    #for ig in range(len(coeff1)):
    #    add = add + np.conj(coeff1[ig]) * coeff2[ig]
    #print("Up:",dp_up, "Dn:", dp_dn)
    return dp

def get_coeffs(ik, ib, gvecs, ngk, coeffs):
    # coeffs_x = []
    # gvecs_x = []
    kstart = np.sum(ngk[:ik])
    kend = kstart + ngk[ik]
    coeffs_x = coeffs[ib, :, kstart:kend, 0] + 1j*coeffs[ib, :, kstart:kend, 1]
    gvecs_x = gvecs[kstart:kend]
    # for ig in range(kstart, kend):
    #    coeffs_x.append(np.complex(coeffs[ib,0,ig,0], coeffs[ib,0,ig,1]))
    #    gvecs_x.append(gvecs[ig])
    return coeffs_x, gvecs_x

def get_coeffs_kbarr(ik_arr, ib_arr, gvecs, ngk, coeffs):
    coeffs_kbarr = []
    gvecs_kbarr = []
    for ik in ik_arr:
        tmp_c = []
        tmp_g = []
        for ib in ib_arr:
            coeffs_ikib, gvecs_ikib = get_coeffs(ik, ib, gvecs, ngk, coeffs)
            tmp_c.append(coeffs_ikib)
            tmp_g.append(gvecs_ikib)
        coeffs_kbarr.append(tmp_c)
        gvecs_kbarr.append(tmp_g)
    return np.array(coeffs_kbarr), np.array(gvecs_kbarr)

f = h5py.File("WFN.h5", 'r')

coeffs = f['wfns/coeffs'][:]
ngk = f['mf_header/kpoints/ngk'][:]
rk= f['mf_header/kpoints/rk'][:]
el= f['mf_header/kpoints/el'][:]
nrk= f['mf_header/kpoints/nrk'][()]
alat = f['mf_header/crystal/alat'][()]
blat = f['mf_header/crystal/blat'][()]
bvec= f['mf_header/crystal/bvec'][:]
gvecs = f['wfns/gvecs'][:]

# working on the K point
"""print("K-point:", rk[0])
print("K-point:", rk[1])"""
iks = [0, 1, 3]
# ik = 0 is the K point
# ik = 1 is the K + dkx
# ik = 3 is the K + dky

num_bands = 350
ibs = range(num_bands)
c_kb, gv = get_coeffs_kbarr(iks, ibs, gvecs, ngk, coeffs)
print(c_kb.shape)
print(gv.shape)

# c_kb and gv are now the main arrays
#----------------------------------------

# Check dotprod
print(dotprod(c_kb[0,0], c_kb[1,0]))
print(dotprod(c_kb[0,44], c_kb[1,44]))

phase_y = dotprod(c_kb[0,0], c_kb[1,0])
phase_y = phase_y/np.abs(phase_y)
print("phase: " + str(phase_y))
print("du with phase: " + str(np.conj(phase_y)*c_kb[1,0] - c_kb[0,0]))

for i in range(0, 50):
    print(abs(dotprod(c_kb[0,i], c_kb[1,i])))

"""print(el.shape)
print("Eigs", el[0,0,0:80]*13.605698)"""


"""
# Find du_val/dky:
# e^-itheta u_n(k+dky) - u_n k
deltaky = 1e-5*blat#*2*np.pi
deltakx = 1e-5*blat#*2*np.pi
bohr2A = 0.52918
Ryd2ev = 13.605698
#HARTREE UNITS!
e = 1
bohr2A = 0.52918
mu_B = 1/2.
hcross = 1

for i in range(44, 48):
	phase_y = dotprod(c_kb[0,i], c_kb[2,i])
	phase_y = phase_y/np.abs(phase_y)
	phase_x = dotprod(c_kb[0,i], c_kb[1,i])
	phase_x = phase_x/np.abs(phase_x)
	du_dkx = (np.conj(phase_x)*c_kb[1,i] - c_kb[0,i])/deltakx
	du_dky = (np.conj(phase_y)*c_kb[2,i] - c_kb[0,i])/deltaky
	berry_1 = dotprod(du_dkx, du_dky)
	berry_2 = dotprod(du_dky, du_dkx)

	print(phase_x); print(phase_y)
	print("Berry curvature of state:", i, -1j*(berry_1 - berry_2)*bohr2A**2)
	#print(el[0,3,i+24]*13.605698, el[0,1,25]*13.605698)

	swich_xy = 0j
	swich_yx = 0j
	# Check for the completeness of basis
	# by setting H-E = 1
	swich_xy_berry = 0j
	swich_yx_berry = 0j
	for n in range(num_bands):
		dp1 = dotprod(du_dkx, c_kb[0,n])
		dp2 = dotprod(c_kb[0,n], du_dky)
		ediff = (el[0,0,n] - el[0,0,i])
		swich_xy += ediff*dp1*dp2
		swich_xy_berry += 1*dp1*dp2
    
		dp1 = dotprod(du_dky, c_kb[0,n])
		dp2 = dotprod(c_kb[0,n], du_dkx)
		swich_yx += ediff*dp1*dp2
		swich_yx_berry += 1*dp1*dp2
		
	swich = swich_yx - swich_xy
	swich_berry = swich_yx_berry - swich_xy_berry
	pref = -1j*e*(bohr2A**2)/(2*hcross*mu_B)
	#print("Eig of K in Ryd:", el[0, 0, i])
	print("m_z:", -1j*swich/2) #, "berry:", -1j*swich_berry*(bohr2A**2))"""