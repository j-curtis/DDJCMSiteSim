#Jonathan Curtis
#This program will use qutip to compute the propagator for the driven Jaynes-Cummings model 
#We will assume the atom, cavity, and drive are all resonant with each other and we will go to the rotating frame (check this to see if it is ok)

import numpy as np
import qutip as qt

def H(N,x,y):
	"""
	The Hamiltonian in the rotating frame
	We use unitless couplings of 
		x = gt/hbar
		y = Et/hbar
	where t is the time (taken real) hbar is Planks constant
	g is the atomic-field coupling and E is the classical drive on the cavity
	This yields
		H = x(a\sigma^\dagger + a^\dagger \sigma) + y(a+a^\dagger)
	N is the cutoff on the cavity Hilbert space
	"""

	a = qt.tensor( qt.qeye(2), qt.destroy(N) )
	s = qt.tensor( qt.destroy(2), qt.qeye(N) )

	h = x*(a*s.dag() + s*a.dag() ) + y*(a + a.dag() )
	
	return h, a, s

def U(N,x,y):
	"""
	Matrix exponential of the Hamiltonian
	
	U = exp(-i*H(N,x,y) )
	"""
	h,a,s = H(N,x,y)

	u = ( (-1.j) *h ).expm() 

	return u

def main():
	
	N = 100
	vac = qt.tensor(qt.basis(0), qt.basis(N) )

	x = 1.0
	y = .25 

	h,a,s = H(N,x,y)

	u = U(N,x,y)

	print 

	

if __name__ == "__main__":
	main()

	 
