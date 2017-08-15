#Jonathan Curtis
#This program will study the photon blockade breakdown phase transition numerically 

#First we will fix the system size (Hilbert space cutoff) and sweep the drive through the critical point
#We will always take the atomic spontaneous emission to be zero
#We will take all frequencies to be resonant with the drive 
#We will also take the dissipation to be sufficiently small so that the dissipation doesn't play a large role 

#We will then increase the system size and repeat  

import JCSingleSite as JC
import numpy as np


#We will set the coupling strength to unity (which sets the overall energy scale)
#We then will take the decay rate to be very small (roughly 1/hilbert space size)
#We will then sweep the drive strength across the critical point at .5 

def main():
	
	#DRIVE PARAMETERS
	min_drive = 0.0	#Sweep drive stength from this
	max_drive = 1.0	#To this
	num_drive = 8	#With this many steps

	drives = np.linspace(min_drive,max_drive,num_drive)	#array of the drives to use

	#TIME PARAMETERS
	max_time = 100	#Run until this time
	num_time = 400	#With this many steps

	#SYSTEM SIZE PARAMETERS
	min_size = 100	#Minimum Hilbert space size
	max_size = 400	#Maximum Hilber space size
	num_size = 4	#Number of system size points

	sizes = np.linspace(min_size,max_size,num_size,dtype=np.int16)		#array of system sizes to use

	#DECAY RATE
	decay = 1.0/max_size 	#Decay rate of the cavity (we hold fixed at the smallest value it will need to be)

	print("MAXTIME=",max_time)
	print("NUMTIME=",num_time)
	print("DECAY=",decay)
	print("\n")
	
	for s in sizes:
		for d in drives:
		
			sim = JC.DDJCMSim(s)
			sim.drive = d
			
			sim.maxtime = max_time
			sim.numtimesteps = num_time

			sim.genHamiltonian()
			sim.genTimes()

			sim.runSim()

			print("DRIVE=",d)
			print("SIZE=",s)
			print("COHERENCE=",sim.rho.expect[0][-1])
			print("PHOTON NUMBER=",sim.rho.expect[1][-1])
			print("\n")
			
	

main()




