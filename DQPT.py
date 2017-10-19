#Jonathan Curtis
#This code will compute the matrix exponential of the driven Jaynes-Cummings Hamiltonian
#In particular, we are looking for the < vac| exp(-isH)|vac> elements where |vac> is the undriven vacuum and H has a quenched drive 
#We will exploit the fast numerical calculation to compute this in the complex plane and then study the behavior 

import numpy as np
import qutip as qt

def H(cutoff, drive, coupling):
	"""
	This will be the Hamiltonian 
		H = drive*(a+a^\dagger) + coupling*( a*\sigma^\dagger + a^\dagger*\sigma)
	with photon Hilbert space truncated to cutoff photons.
	We take the parameters appearing to be real numbers with units of frequency 
	"""

	a = qt.tensor(	qt.eye(2), 	qt.destroy(cutoff))
	s = qt.tensor(	qt.destroy(2),	qt.eye(cutoff))

	h = drive*( a + a.dag() ) + coupling*(a*s.dag() + s*a.dag() )

	return h 

def propagator(t, H):
	"""Computes the propagator exp(-itH) extended to complex time """

	return (-np. 
