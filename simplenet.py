from matplotlib import pyplot
import random
from datetime import datetime
import pickle
from neuron import h, gui


class Cell(object):
    def __init__(self):
        self.synlist = []
        self.createSections()
        self.buildTopology()
        self.defineGeometry()
        self.defineBiophysics()
        self.createSynapses()
        self.nclist = []

    def createSections(self):
        pass

    def buildTopology(self):
        pass

    def defineGeometry(self):
        """Set the 3D geometry of the cell."""
        self.soma.L = 18.8
        self.soma.diam = 18.8
        self.soma.Ra = 123.0

    def defineBiophysics(self):
        pass

    def createSynapses(self):
        """Add an exponentially decaying synapse """
        syn = h.ExpSyn(self.soma(0.5))
        syn.tau = 2
        syn.e = 0
        self.synlist.append(syn) # synlist is defined in Cell

    def associateGid (self):
        pc.set_gid2node(self.gid, idhost)
        nc = h.NetCon(self.soma(0.5)._ref_v, None, sec=self.soma)
        nc.threshold = 10
        pc.cell(self.gid, nc)
        del nc # discard netcon


    def createNetcon(self, thresh=10):
        """ created netcon to record spikes """
        nc = h.NetCon(self.soma(0.5)._ref_v, None, sec = self.soma)
        nc.threshold = thresh
        return nc

    def createStim(self, number=1e9, start=1, noise=0.5, rate=50, weight=1, delay=5):
        self.stim = h.NetStim()
        self.stim.number = number
        self.stim.start = start
        self.stim.noise = noise
        self.stim.interval = 1000.0/rate
        self.ncstim = h.NetCon(self.stim, self.synlist[0])
        self.ncstim.delay = delay
        self.ncstim.weight[0] = noise # NetCon weight is a vector.

    def connect2Source(self, sourceCell, thresh=10):
        """Make a new NetCon with the source cell's membrane
        potential at the soma as the source (i.e. the spike detector)
        onto the target (i.e. a synapse on this cell)."""
        nc = h.NetCon(sourceCell.soma(1)._ref_v, self.synlist[0], sec = sourceCell.soma)
        nc.threshold = thresh
        return nc

    def setRecording(self):
        """Set soma, dendrite, and time recording vectors on the cell. """
        self.soma_v_vec = h.Vector()   # Membrane potential vector at soma
        self.tVec = h.Vector()        # Time stamp vector
        self.soma_v_vec.record(self.soma(0.5)._ref_v)
        self.tVec.record(h._ref_t)

    def plotTraces(self):
        """Plot the recorded traces"""
        pyplot.figure(figsize=(8,4)) # Default figsize is (8,6)
        somaPlot = pyplot.plot(self.tVec, self.soma_v_vec, color='black')
        pyplot.legend(somaPlot, ['soma'])
        pyplot.xlabel('time (ms)')
        pyplot.ylabel('mV')
        pyplot.title('Cell %d voltage trace'%(self.gid))
        pyplot.show()
        #pyplot.savefig('traces')


class HHCell(Cell): 
    """HH cell: A soma with active channels""" 
    def createSections(self):
        """Create the sections of the cell."""
        self.soma = h.Section(name='soma', cell=self)
    
    def defineBiophysics(self):
        """Assign the membrane properties across the cell."""
        # Insert active Hodgkin-Huxley current in the soma
        self.soma.insert('hh')
        self.soma.gnabar_hh = 0.12  # Sodium conductance in S/cm2
        self.soma.gkbar_hh = 0.036  # Potassium conductance in S/cm2
        self.soma.gl_hh = 0.003    # Leak conductance in S/cm2
        self.soma.el_hh = -70       # Reversal potential in mV


class Net:
    """Creates Network of N neurons (using parallelContext)
    Connectivity and stimulation params provided as arguments
    Also ncludes methods to gather and plot spikes
    """
    def __init__(self, N=5, cellType=HHCell, connParams={}, stimParams={}):

        """
        N: Number of cells.
        cellType: class of cell type
        connParams: dict of connectivity params
        stimParams: dict of stimulation params

        """
        self.cellType = cellType
        self.N = N                      # number of cells
        self.connParams = connParams    # connectivity params
        self.stimParams = stimParams    # backgroudn stim params
        self.cells = []                 # List of Cell objects in the net
        self.nclist = []                # List of NetCon in the net
        self.tVec = h.Vector()         # spike time of all cells
        self.idVec = h.Vector()        # cell ids of spike times
        self.createNet()  # Actually build the net
    

    def createNet(self):
        """Create, layout, and connect N cells."""
        self.setGids() #### set global ids (gids), used to connect cells
        self.createCells()
        self.connectCells() 
        self.createStims()


    def setGids(self):
        self.gidList = []
        #### Round-robin counting. Each host as an id from 0 to nhost - 1.
        for i in range(idhost, self.N, nhost):
            self.gidList.append(i)

    
    def createCells(self):
        """Create and layout cells (in this host) in the network."""
        self.cells = []

        for i in self.gidList: #### Loop over cells in this node/host
            cell = self.cellType() # dynamically create cell object 
            self.cells.append(cell)  # add cell object to net cell list
            cell.gid = i # assign gid (can be any unique integer)
            cell.associateGid() # associated gid to each cell
            pc.spike_record(cell.gid, self.tVec, self.idVec) # Record spikes of this cell
            
            print 'Created cell %d on host %d out of %d'%(i, idhost, nhost) 

    def connectCells(self):
        """Connect cells"""
        connType = self.connParams['type']
        if connType == 'rand':
            weight = self.connParams['weight']
            delayMean = self.connParams['delayMean']
            delayVar = self.connParams['delayVar']
            delayMin = self.connParams['delayMin']
            maxConnsPerCell = self.connParams['maxConnsPerCell']
            self.nclist = []

            ## create random delays
            random.seed(randSeed)  # Reset random number generator  
            randDelays = [max(delayMin, random.gauss(delayMean, delayVar)) for pre in range(maxConnsPerCell*self.N)] # select random delays based on mean and var params    

            ## loop over postsyn gids in this host
            for postCell in self.cells:  
                preGids = [gid for gid in self.gidList if gid != postCell.gid] # get list of presyn cell gids (omit post to prevent self connection)
                randPreGids = random.sample(preGids, random.randint(0, min(maxConnsPerCell, len(preGids)))) 
                for preGid in randPreGids: # for each presyn cell
                    nc = pc.gid_connect(preGid, postCell.synlist[0]) # create NetCon by associating pre gid to post synapse
                    nc.weight[0] = weight
                    nc.delay = randDelays.pop()
                    nc.threshold = 10
                    self.nclist.append((preGid,postCell.gid,nc))
                    postCell.nclist.append((preGid,nc))
                    
                    print 'Created conn between pregid %d and postgid %d on host %d'%(preGid, postCell.gid, idhost) 


    def createStims(self):
        """Connect a spiking generator to the first cell to get
        the network going."""
        #### Only continue if the first cell is not on this host
        for cell in self.cells:
            cell.createStim(
                noise=self.stimParams['noise'], 
                rate=self.stimParams['rate'], 
                weight=self.stimParams['weight'],
                delay=self.stimParams['delay'])
            print 'Created stim on cell %d on host %d'%(cell.gid, idhost) 

    def gatherSpikes(self):
        """Gather spikes from all nodes/hosts"""
        if idhost==0: print 'Gathering spikes ...'
        data = [None]*nhost
        data[0] = {'tVec': self.tVec, 'idVec': self.idVec}
        pc.barrier()
        gather=pc.py_alltoall(data)
        pc.barrier()
        self.tVecAll = [] 
        self.idVecAll = [] 
        if idhost==0:
            for d in gather:
                self.tVecAll.extend(list(d['tVec']))
                self.idVecAll.extend(list(d['idVec']))

    def plotRaster(self):

        print 'Plotting raster ...'
        pyplot.figure()
        pyplot.scatter(self.tVecAll,self.idVecAll,marker=".",s=1,color='blue')
        pyplot.xlabel('Time (ms)')
        pyplot.ylabel('Cell ID')
        pyplot.title('Raster Plot of Network with 500 HH Cells')
        pyplot.xlim(0,max(self.tVecAll))
        pyplot.ylim(0,self.N)
        pyplot.show()
        #pyplot.savefig('raster')

    def saveData(self):
        print 'Savind data ...'
        dataSave = {'N': self.N, 'connParams': self.connParams, 'stimParams': self.stimParams, 'tVec': self.tVecAll, 'idVec': self.idVec}
        with open('output.pkl', 'wb') as f:
            pickle.dump(dataSave, f)

#### New ParallelContext object 
pc = h.ParallelContext()
pc.set_maxstep(10)
idhost = int(pc.id())
nhost = int(pc.nhost())

# set randomizer seed
randSeed = 1

# create network
net = Net(N=500, cellType=HHCell,
    connParams={'type': 'rand', 'weight': 0.004, 'delayMean': 13.0, 'delayVar': 1.4, 'delayMin': 0.2, 'maxConnsPerCell': 20}, 
    stimParams={'rate': 50, 'noise': 0.5, 'weight': 50, 'delay':5}) 

# set voltage recording for cell 0 in net
net.cells[0].setRecording()

# run sim and gather spikes
h.stdinit()
h.dt = 0.025
duration = 1000.0
if idhost==0: 
    print 'Running sim...'
    startTime = datetime.now() # store sim start time
pc.psolve(duration)  # actually run the sim in all nodes
if idhost==0: 
    runTime = (datetime.now() - startTime).total_seconds()  # calculate run time
    print "Run time for %d sec sim = %.2f sec"%(int(duration/1000.0), runTime)
net.gatherSpikes()  # gather spikes from all nodes onto master node

# plot net raster, save net data and plot cell 0 traces 
if idhost==0:
    net.plotRaster()
    net.saveData()
    net.cells[0].plotTraces()

pc.done()
