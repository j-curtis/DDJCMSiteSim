#Jonathan Curtis
#This program will use qutip to compute the propagator for the driven Jaynes-Cummings model 
#We will assume the atom, cavity, and drive are all resonant with each other and we will go to the rotating frame (check this to see if it is ok)

import numpy as np
import qutip as qt

def Propagator(N,x,y):
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
	It also returns the matrix exponential of e^(-i H(x,y) )
	"""

	a = qt.tensor( qt.qeye(2), qt.destroy(N) )
	s = qt.tensor( qt.destroy(2), qt.qeye(N) )

	h = x*(a*s.dag() + s*a.dag() ) + y*(a + a.dag() )
	u = ( -1.j * h).expm()	

	return h, a, s, u


def main():
	
	N = 200
	vac = qt.tensor(qt.basis(2), qt.basis(N) )

	ratio = .25	#ratio of y/x = E/g = .5 control-parameter
			#critical point occurs at ratio = .5
	
	numSweep = 50
	OTOC = np.zeros(shape = numSweep)

	scale = np.linspace( 0.0, 10.0, numSweep)	#x = scale sets the overall time scale
						#this is the parameter we sweep to look for DQPT

	index = 0
	for x in scale:
		y = ratio*x
		h,a,s,u = Propagator(N,x,y)

		#Compute OTOC
		#OTOC of interest is <vac | ( [ a(t) , a(0)^\dagger ] )^2 |vac> 
		at = u.dag() * a * u 
		commutator = ( at*a.dag() - (a.dag() )*at)
		commSquare = commutator**2

		OTOC[index] += qt.expect( commSquare, vac)
	
		print( OTOC[index] )
		index+= 1
	

if __name__ == "__main__":
	main()

	 
