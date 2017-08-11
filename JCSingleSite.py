#Jonathan Curtis
#This code will run numerical solutions of the driven-dissipative Jaynes-Cummings master equation 

import numpy as np
import qutip as qt

class DDJCMSim:
	"""Class for running single site of the Jaynes-Cummings model"""
	
	def __init__(self,cutoff_param):
		"""Initializes an instance of the simulation and sets parameters to defaults"""
		self.cutoff = cutoff_param 	#Max photon count in Hilbert space 		

		self.detuning = 0	#We assume cavity and atom are resonant with each other. This is cav_freq - drive_freq
		self.drive = 0		#External drive strength. We drive the field X quadrature

		self.coupling = 1	#Jaynes-Cummings coupling between atom and cavity in RWA
	
		self.decays = [0.5,0.0]	#The decay rates are [0] = cavity SE, [1] = atom SE
		
		self.maxTime = 40	#Max time integrated to
		self.numTimeSteps = 400	#Number of time steps integrated for 
		
		self.a = qt.tensor( qt.qeye(2),qt.destroy(self.cutoff))	#Lowering operator for the cavity in the Hilbert space H = Hspin x Hcav
		self.sigma = [ qt.tensor( qt.sigmax(), qt.qeye(self.cutoff) ), qt.tensor( qt.sigmay(), qt.qeye(self.cutoff) ), qt.tensor(qt.sigmaz(),qt.qeye(self.cutoff))]	#Pauli algebra for the spin
		self.sigmin = 0.5*(self.sigma[0] -0.j*self.sigma[1] )	#atomic lowering operator |0><1|
		
		self.collapse = [np.sqrt( self.decays[0])*self.a, np.sqrt(self.decays[1])*self.sigmin ]	#collapse operators for the master equation
		self.measure = [self.a, self.a.dag()*self.a ]	#operators we want to measure the expected values for 
		self.numMeasure = len(self.measure)	#The number of operators we are measuring

		self.initKet = qt.tensor( qt.basis(2,0),qt.basis(self.cutoff,0) )	#initializes the state as a pure state in the vacuum 

	def genHamiltonian(self):
		"""Generates the Hamiltonian from the current parameters"""
		self.H = self.detuning*(self.a.dag()*self.a + 0.5*self.sigma[2] )+self.drive*(self.a + self.a.dag() ) + self.coupling*(self.a*self.sigmin.dag() + self.a.dag()*self.sigmin )	#Driven Jaynes-Cummings model in the rotating frame 
		
	def genTimes(self):
		"""Generates an array of times given the current parameters"""
		self.times = np.linspace(0,self.maxTime,self.numTimeSteps) 

	def holsteinPrimakoff(self,hp_cutoff):
		"""Re-runs the simulation but this time uses a Holstein-Primakoff Hamiltoninan"""
		self.cutoffHP = hp_cutoff
		self.aHP = qt.tensor(qt.eye(self.cutoffHP),qt.destroy(self.cutoff))	#Cavity boson 
		self.bHP = qt.tensor(qt.destroy(self.cutoffHP),qt.qeye(self.cutoff))	#Holstein Primakoff boson

		self.collapseHP = [np.sqrt(self.decays[0])*self.a]	#A the moment we do not allow for spontaneous emission in Holstein Primakoff
		self.measure = [self.aHP,self.aHP.dag()*self.aHP,self.bHP,self.bHP.dag()*self.bHP ]	#Measure population and coherence for each

		self.initKetHP = qt.tensor(qt.bsis(self.cutoffHP,0),qt.basis(self,cutoff,0))	#We start at the critical point values of the unentangled fields
		self.HHP = self.coupling/2.0 *( self.aHP.dag() + self.aHP)*(self.bHP.dag()*self.bHP - 1) + drive*(self.aHP.dag()+self.aHP) + coupling/2.0*(1.j)*(self.aHP.dag() - self.aHP)*(self.bHP.dag()+self.bHP) 	#This is the form of the HP Hamiltonian to quadratic order with no detuning

		self.rhoHP = qt.mesolve(self.HHP,self.initKetHP,self.times,self.collapseHP,self.measureHP)

	def runSim(self):
		"""Runs the QuTip mesolve method and computes the density matrix"""
		self.rho = qt.mesolve(self.H,self.initKet,self.times,self.collapse,self.measure)
	

def main():
	
	sim1 = DDJCMSim(300)	

	sim1.coupling = 0

	sim1.genHamiltonian()
	sim1.genTimes()

	sim1.runSim()


if __name__=="__main__":
	main()



