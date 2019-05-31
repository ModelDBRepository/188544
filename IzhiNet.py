
"""
params.py 

netParams is an object containing a set of network parameters using a standardized structure

simConfig is an object containing a set of simulation configurations using a standardized structure

"""

from netpyne import specs

netParams = specs.NetParams()  # object to store sets of network parameters
simConfig = specs.SimConfig()  # object to store sets of simulation configurations


###############################################################################
#
# MPI HH TUTORIAL PARAMS
#
###############################################################################

###############################################################################
# NETWORK PARAMETERS
###############################################################################

# Population parameters
netParams.popParams['PYR'] = {'cellModel': 'Izhi2007b', 'cellType': 'PYR', 'numCells': 500} # add dict with params for this pop 
netParams.popParams['background'] = {'cellModel': 'NetStim', 'rate': 10, 'noise': 0.5, 'source': 'random'}  # background inputs


# PYR cell properties
netParams.importCellParams(label='PYR_rule', conds= {'cellType': 'PYR', 'cellModel': 'Izhi2007b'},
    fileName='izhi2007Wrapper.py', cellName='IzhiCell',  cellArgs={'type':'RS'})  


# Synaptic mechanism parameters
netParams.synMechParams['exc'] = {'mod': 'ExpSyn', 'tau': 0.1, 'e': 0}


# Connectivity parameters
netParams.connParams['PYR->PYR'] = {
    'preConds': {'popLabel': 'PYR'}, 'postConds': {'popLabel': 'PYR'},
    'weight': 0.0005,                    # weight of each connection 
    'delay': '0.2+gauss(13.0, 1.4)',    # delay function (min=0.2, mean=13, var=1.4)
    'threshold': 10,                    # threshold
    'convergence': 'uniform(0, 10)'}   # convergence function (num of presyn conns per postsyn)      

netParams.connParams['bg->PYR'] = {
    'preConds': {'popLabel': 'background'}, 'postConds': {'cellType': 'PYR'}, # background -> PYR
    'weight': 0.5,  
    'synMech': 'exc',
    'delay': 5}


###############################################################################
# SIMULATION PARAMETERS
###############################################################################

# Simulation parameters
simConfig.duration = 5*1e3 # Duration of the simulation, in ms
simConfig.dt = 0.025 # Internal integration timestep to use
simConfig.verbose = False  # show detailed messages 
simConfig.timing = True  # record timing
simConfig.cache_efficient = True  # use CVode cache_efficient option to optimize load when running on many cores

# Recording 
simConfig.recordStim = False  # record spikes of cell stims
simConfig.recordStep = 0.1 # Step size in ms to save data (eg. V traces, LFP, etc)

# Saving
simConfig.filename = 'IzhiNet'  # Set file output name
simConfig.saveFileStep = 1000 # step size in ms to save data to disk
simConfig.savePickle = False # Whether or not to write spikes etc. to a .mat file
simConfig.saveJson = False # Whether or not to write spikes etc. to a .mat file
simConfig.saveMat = False # Whether or not to write spikes etc. to a .mat file


# Analysis and plotting 
simConfig.analysis['plotRaster'] = True # Whether or not to plot a raster