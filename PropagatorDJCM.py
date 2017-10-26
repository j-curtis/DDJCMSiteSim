#Jonathan Curtis
#This program will use qutip to compute the propagator for the driven Jaynes-Cummings model 
#We will assume the atom, cavity, and drive are all resonant with each other and we will go to the rotating frame (check this to see if it is ok)

import numpy as np
import qutip as qt
from matplotlib import pyplot as plt

def Propagator(N,x,y,t):
	"""
	The Hamiltonian in the rotating frame
	We use unitless couplings of 
		x = g*(time-scale)/hbar
		y = E*(time-scale)/hbar
	time-scale is some characteristic time scale (so that physical time is unitless) 
	hbar is Planks constant
	g is the atomic-field coupling 
	E is the classical drive on the cavity
	t is the unitless time 
	
	This yields
		H = x(a\sigma^\dagger + a^\dagger \sigma) + y(a+a^\dagger)
	
	N is the cutoff on the cavity Hilbert space
	"""

	a = qt.tensor( qt.qeye(2), qt.destroy(N) )
	s = qt.tensor( qt.destroy(2), qt.qeye(N) )

	h = x*(a*s.dag() + s*a.dag() ) + y*(a + a.dag() )
	u = ( -1.j *t* h).expm()	

	return u


def main():
	
	cutoffList = [25,50]		#different Hilbert space cutoffs we try, to extract some finite-size scaling
	numCutoff = len(cutoffList)		#number in the Hilbert space cutoff list

	xList = [1.0]				#different atom-cavity couplings we try
	numX = len(xList)			#number of different atom-cavity couplings

	yList = [0.0,.25,.5,.75]	#different drive strengths 
						#y/x = 0 is the free model (obviously integrable)
						#0<y/x<.5 is discrete spectrum (possibly integrable via hidden symmetry, exactly solvable but matrix elements extend throughout hilbert space)
						#y/x = .5 is critical point
						#y > .5 is continous spectrum, should be "isomorphic" to inverted oscillator, dynamical symemtry group solution still conceivable?
	numY = len(yList)			#number of y parameters we try

	numTimes = 25				#number of time-slices we calculate for 
	maxTime = 100				#largest time we compute out to

	times = np.linspace(0.0,maxTime,numTimes)	#time array


	"""
	We are interested in the behavior of the OTOC for the various parameter values given above 
	We will store this in a large array with dimension 4 and lengths 
	(numCutoff x numX x numY x numTimes)
	"""
	
	OTOC = np.zeros( shape = (numCutoff,numX,numY,numTimes) )

	"""
	We will compute the OTOC for the state |vac> defined as 
	|0>_cav \otimes |down>_spin
	It is the vacuum of the y= 0 model

	The operators we commute are 
	a(t) and a^\dagger(0) 

	and the OTOC is defined as -<vac| [a(t), a^+(0) ]^2|vac >
	"""
	
	for iCutoff in np.arange(numCutoff):
		N = cutoffList[iCutoff]
		vac = qt.tensor( qt.basis(2,0), qt.basis(N,0) )
	
		a = qt.tensor( qt.qeye(0) , qt.destroy(N) )

		for iX in np.arange(numX):
			x = xList[iX]
			
			for iY in np.arange(numY):
				y = yList[iY]
		
				for iTime in np.arange(numTimes):
					t = times[iTime]
			
					U = Propagator(N,x,y,t)
	
					OTOOperator = -1.0 *( (U.dag() *a*U *a.dag() - a.dag() *U.dag() * a *U) **2 )		

					OTOC[iCutoff,iX,iY,iTime] += np.real( qt.expect( OTOOperator, vac) )
	

	###Save data to file 
	np.savetxt("cutoffList.txt", cutoffList)
	np.savetxt("xList.txt", xList)
	np.savetxt("yList.txt", yList)
	np.savetxt("times.txt", times)
	np.savetxt("OTOC.txt", OTOC)
	

if __name__ == "__main__":
	main()

	 
