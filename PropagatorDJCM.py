#Jonathan Curtis
#This program will use qutip to compute the propagator for the driven Jaynes-Cummings model 
#We will assume the atom, cavity, and drive are all resonant with each other and we will go to the rotating frame (check this to see if it is ok)

import numpy as np
import qutip as qt
from matplotlib import pyplot as plt

def Propagator(N,d,t):
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

	a = qt.tensor( qt.qeye(2), qt.destroy(N) )
	s = qt.tensor( qt.destroy(2), qt.qeye(N) )

	h = a*s.dag() + s*a.dag() + d*(a + a.dag() )
	u = ( -1.j *t* h).expm()	

	return u


def main():
	
	cutoffList = [100,200,300]		#different Hilbert space cutoffs we try, to extract some finite-size scaling
	numCutoff = len(cutoffList)		#number in the Hilbert space cutoff list

	dList = [0.0,.25,.5,.75]		#different drive strengths 
						#d = 0 is the free model (obviously integrable)
						#0<d<.5 is discrete spectrum (possibly integrable via hidden symmetry, exactly solvable but matrix elements extend throughout hilbert space)
						#d = .5 is critical point
						#d > .5 is continous spectrum, should be "isomorphic" to inverted oscillator, dynamical symemtry group solution still conceivable?
	numd = len(dList)			#number of d parameters we try

	numTimes = 200				#number of time-slices we calculate for 
	maxTime = 100				#largest time we compute out to

	times = np.linspace(0.0,maxTime,numTimes)	#time array


	"""
	We are interested in the behavior of the OTOC for the various parameter values given above 
	We will store this in a large array with dimension 3 and lengths 
	(numCutoff x numd x numTimes)
	"""
	
	OTOC = np.zeros( shape = (numCutoff,numd,numTimes) )

	"""
	We will compute the OTOC for the state |vac> defined as 
	|0>_cav \otimes |down>_spin
	It is the vacuum of the d = 0 model

	The operators we commute are 
	a(t) and a^\dagger(0) 

	and the OTOC is defined as -<vac| [a(t), a^+(0) ]^2|vac >
	"""
	
	for iCutoff in np.arange(numCutoff):
		N = cutoffList[iCutoff]
		vac = qt.tensor( qt.basis(2,0), qt.basis(N,0) )
	
		a = qt.tensor( qt.qeye(2) , qt.destroy(N) )

		for idrive in np.arange(numd):
			d = xList[idrive]
		
			for iTime in np.arange(numTimes):
				t = times[iTime]
		
				U = Propagator(N,d,t)

				OTOOperator = -1.0 *( (U.dag() *a*U *a.dag() - a.dag() *U.dag() * a *U) **2 )		

				OTOC[iCutoff,idrive,iTime] += np.real( qt.expect( OTOOperator, vac) )


	###Save data to numpy file 
	np.save("cutoffList.npy", cutoffList)
	np.save("dList.npy", dList)
	np.save("times.npy", times)
	np.save("OTOC.npy", OTOC)
	

if __name__ == "__main__":
	main()

	 
