#Jonathan Curtis
#This code will compute the steady state density matrix for the driven dissipative Jaynes Cummings model 
#It will then be used to study the photon number scaling. Hopefully we can learn something that perhaps Mike's numerics were not uncovering

import numpy as np
import qutip as qt

class rhoSS:
	"""Class for computing steady state density matrix
	Accepts as parameters: 
		system size (cutoff)
		decay rate (cavity photon loss rate)
		drive strength (cavity coherent pump field)
	
	Takes zero detuning and sets cavity/atom coupling to one
	In other words,

	H = a^+ sig + sig^+ a + drive(a+a^+) 

	Loss is 
	k(a rho a^+ -.5 {a^+a ,rho)
	"""
	
	def __init__(self,cutoff_p,drive_p,decay_p):
		"""Computes the steady state"""

		self.cutoff = cutoff_p
		self.drive = drive_p
		self.decay = decay_p
	
		self.a = qt.tensor(qt.destroy(self.cutoff), qt.qeye(2) )
		self.sig = qt.tensor(qt.qeye(self.cutoff), qt.destroy(2) )

		self.rho = qt.steadystate( self.a*self.sig.dag() + self.a.dag()*self.sig + self.drive*(self.a + self.a.dag() ) , [np.sqrt(self.decay)*self.a] )
		
		self.cavNum = qt.expect( self.a.dag()*self.a, self.rho)

def main():
	N = 100
	drive = .2 
	decay = .01

	rho = rhoSS(N,drive,decay)

	print( rho.cavNum)

if __name__ == "__main__":
	main()	
