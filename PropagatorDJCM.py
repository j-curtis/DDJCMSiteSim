#Jonathan Curtis
#This program will use qutip to compute the propagator for the driven Jaynes-Cummings model 
#We will assume the atom, cavity, and drive are all resonant with each other and we will go to the rotating frame (check this to see if it is ok)

import numpy as np
import qutip as qt
from matplotlib import pyplot as plt

def propagator(size_p,drive_p,time_p):
	"""
	The Hamiltonian in the rotating frame
	We use unitless drive parameter of 
		d = E/g
	g is the atomic-field coupling which we set to 1.0 
	E is the classical drive on the cavity
	t is the unitless time 
	
	This yields
		H = a\sigma^\dagger + a^\dagger \sigma + w(a+a^\dagger)
	
	N is the cutoff on the cavity Hilbert space
	"""

	a = qt.tensor( qt.qeye(2), qt.destroy(size_p) )
	s = qt.tensor( qt.destroy(2), qt.qeye(size_p) )

	h = a*s.dag() + s*a.dag() + drive_p*(a + a.dag() )
	u = ( -1.j *time_p* h).expm()	

	return u

def ppCommutator(size_p,drive_p,time_p):
    	"""
    	Computes the time-dependent p-p commutator operator as 
    	-i*[p(t),p(0)] 
    	with p = i(a-a^+)/root2
    	"""
    	U = propagator(size_p,drive_p,time_p)
    	a = qt.tensor(qt.qeye(2),qt.destroy(size_p))
    
    	return 0.5j*( U.dag()*( a-a.dag() )*U*( a-a.dag() ) - ( a-a.dag() )*U.dag()*( a-a.dag() )*U )


def ppCEigenval(size_p,drive_p,time_p):
	"""
	Returns an array of the eigenvalues of the pp commutator
	"""
	return ppCommutator(size_p,drive_p,time_p).eigenenergies()

def main():
	
	cutoffList = [100]		#different Hilbert space cutoffs we try, to extract some finite-size scaling
	numCutoff = len(cutoffList)		#number in the Hilbert space cutoff list

	dList = [0.0,.25,.5,.75]		#different drive strengths 
						#d = 0 is the free model (obviously integrable)
						#0<d<.5 is discrete spectrum (possibly integrable via hidden symmetry, exactly solvable but matrix elements extend throughout hilbert space)
						#d = .5 is critical point
						#d > .5 is continous spectrum, should be "isomorphic" to inverted oscillator, dynamical symemtry group solution still conceivable?
	numd = len(dList)			#number of d parameters we try

	numTimes = 4				#number of time-slices we calculate for 
	maxTime = 50				#largest time we compute out to

	times = np.linspace(0.0,maxTime,numTimes)	#time array


	"""Let us now compute the time dependent commutator and use this to determine the eigenvalue spectrum"""

	eigenvalues = np.array( [[[ ppCEigenval(N,d,t) for t in times ] for d in dList] for N in cutoffList] )

	###Save data to numpy file 
	np.save("cutoffList.npy", cutoffList)
	np.save("dList.npy", dList)
	np.save("times.npy", times)
	np.save("eigenvalues.npy", eigenvalues)
	

if __name__ == "__main__":
	main()

	 
