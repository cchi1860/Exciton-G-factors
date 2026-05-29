In this folder, we use Siesta instead of Quantum Espresso to calculate the wavefunction and eigenvalues in the atomic orbital basis.

In Siesta2Bgw.save, we convert siesta to BGW and obtain u_nk and g_vector arrays, in the basis of plane waves instead of atomic orbitals. This is done using siesta2bgw.inp.
In the main folder, we stay in the atomic orbital basis and use siesta_ham.py to extract the Hamiltonian and overlap matrix of the necessary k-points. Then, we calculate the g-factors.
