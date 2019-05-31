"""
init.py

Startup script to run model.

Model was developed using NetPyNE, a python package to facilitate the development, 
parallel simulation and analysis of biological neuronal networks using the NEURON simulator.

www.neurosimlab.org/netpyne

Usage:
    python init.py # Run simulation, optionally plot a raster

MPI usage:
    mpiexec -n 4 nrniv -python -mpi init.py

Contributors: salvadordura@gmail.com
"""

from netpyne import sim
import HHNet, IzhiNet, HybridNet


###############################################################################
# Sequence of commands to run full model
###############################################################################
def runModel():
    sim.initialize(                
        simConfig = HHNet.simConfig, 
        netParams = HHNet.netParams)  
    sim.net.createPops()                  # instantiate network populations
    sim.net.createCells()                 # instantiate network cells based on defined populations
    sim.net.connectCells()                # create connections between cells based on params
    sim.setupRecording()                  # setup variables to record for each cell (spikes, V traces, etc)
    sim.runSim()                          # run parallel Neuron simulation  
    sim.gatherData()                      # gather spiking data and cell info from each node
    sim.saveData()                        # save params, cell info and sim output to file (pickle,mat,txt,etc)
    sim.analysis.plotData()               # plot spike raster

runModel()                              # execute sequence of commands to run full model
