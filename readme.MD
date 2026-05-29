Repository for calculating exciton g-factors, uploaded May 29, 2026. See folders for extra readme's, and the g-factor results PDF.
The three methods are: QE calculation -> WFN.h5 -> g-factors; Siesta calculation -> siesta2BGW -> g-factors; Siesta calculation -> extract H -> g-factors.
WFN.h5 files were not uploaded for sake of size.

Current issue:
The third method may not be as accurate as possible, as we want to calculate the correct overlap integral <u_nk | u_nk+dk> between k-points, yet the vectors are under different S_k, S_k+dk overlap matrices. We want to calculate this as done on page 22 of https://wiki.physics.udel.edu/wiki_qttg/images/7/79/SIESTA%3Dmethod_for_ab_initio_order_N_simulations.pdf. This may require modifying Siesta source code or modifying Siesta's wannier90 implementation, which computes the overlap integrals.


Future goal: solving g-factors for 4x4 or other supercells, using Siesta and extracting the Hamiltonian directly.
