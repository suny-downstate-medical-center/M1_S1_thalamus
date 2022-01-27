"""
cfg.py 

Simulation configuration for S1 model (using NetPyNE)
This file has sim configs as well as specification for parameterized values in netParams.py 

Contributors: salvadordura@gmail.com, fernandodasilvaborges@gmail.com, joaovvitor@gmail.com
"""

from netpyne import specs
import pickle

cfg = specs.SimConfig()  

#------------------------------------------------------------------------------
#
# SIMULATION CONFIGURATION
#
#------------------------------------------------------------------------------

cfg.simType='M1THS1'
cfg.coreneuron = False

#------------------------------------------------------------------------------
# Run parameters
#------------------------------------------------------------------------------
cfg.duration = 1.0*1e2
cfg.dt = 0.025
cfg.seeds = {'conn': 4321, 'stim': 1234, 'loc': 4321} 
cfg.hParams = {'celsius': 34, 'v_init': -80}  
cfg.verbose = 0
cfg.createNEURONObj = True
cfg.createPyStruct = True
cfg.connRandomSecFromList = False  # set to false for reproducibility 
cfg.cvode_active = False
cfg.cvode_atol = 1e-6
cfg.cache_efficient = True
cfg.printRunTime = 0.1
cfg.pt3dRelativeToCellLocation = True
cfg.oneSynPerNetcon = False  # only affects conns not in subconnParams; produces identical results

cfg.includeParamsLabel = False
cfg.printPopAvgRates = [0, cfg.duration]

cfg.checkErrors = False
cfg.checkErrorsVerbose = False

cfg.rand123GlobalIndex = None

#------------------------------------------------------------------------------
#
# Cells
#
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------  
# S1 Cells
#------------------------------------------------------------------------------  
# Load 55 Morphological Names and Cell pop numbers -> L1:6 L23:10 L4:12 L5:13 L6:14
# Load 207 Morpho-electrical Names used to import the cells from 'cell_data/' -> L1:14 L23:43 L4:46 L5:52 L6:52

cfg.poptypeNumberS1 = 55 # max 55
cfg.celltypeNumberS1 = 207 # max 207

# TO DEBUG - import and simulate only the Cell soma (to study only the Net)
cfg.reducedtestS1 = False    

with open('../info/anatomy/S1-cells-distributions-Mouse.txt') as mtype_file:
    mtype_content = mtype_file.read()       

cfg.popNumberS1 = {}
cfg.cellNumberS1 = {} 
cfg.popLabelS1 = {} 
popParamS1 = []
cellParamS1 = []
cfg.meParamLabelsS1 = {} 
cfg.popLabelElS1 = {} 
cfg.cellLabelS1 = {}

for line in mtype_content.split('\n')[:-1]:
    cellname, mtype, etype, n, m = line.split()
    metype = mtype + '_' + etype[0:3]
    cfg.cellNumberS1[metype] = int(n)
    cfg.popLabelS1[metype] = mtype
    cfg.popNumberS1[mtype] = int(m)
    cfg.cellLabelS1[metype] = cellname

    if mtype not in popParamS1:
        popParamS1.append(mtype)
        cfg.popLabelElS1[mtype] = [] 
               
    cfg.popLabelElS1[mtype].append(metype)
    
    cellParamS1.append(mtype + '_' + etype[0:3])
    
cfg.S1pops = popParamS1[0:55]
cfg.S1cells = cellParamS1[0:207] # pop used in conns 

#used to convengence rules S1-Th
cfg.popNumberS1['ss_RTN_o'] = 382
cfg.popNumberS1['ss_RTN_m'] = 382
cfg.popNumberS1['ss_RTN_i'] = 765
cfg.popNumberS1['VPL_sTC'] = 656
cfg.popNumberS1['VPM_sTC'] = 839
cfg.popNumberS1['POm_sTC_s1'] = 685

#------------------------------------------------------------------------------
# Recording 
#------------------------------------------------------------------------------
allpops = ['NGF1', 	'IT2', 	'PV2', 	 'SOM2',  'VIP2', 	'NGF2',
           'IT4', 	'PV4', 	'SOM4',  'VIP4',  'NGF4',
           'IT5A', 	'PV5A', 'SOM5A', 'VIP5A', 'NGF5A',
           'IT5B', 	'PT5B', 'PV5B',  'SOM5B', 'VIP5B',	'NGF5B',
           'IT6',	'CT6',	'PV6',	 'SOM6',  'VIP6',	'NGF6',
			'VPL_sTC', 	'VPM_sTC', 		'POm_sTC_s1', 
			'VL_sTC', 	'VM_sTC_m1', 	'POm_sTC_m1',
			'mt_RTN',	'ss_RTN_o',		'ss_RTN_m',		'ss_RTN_i'
	   ]

for metype in cfg.S1cells:
	allpops.append(metype)

cfg.cellsrec = 1
if cfg.cellsrec == 0:  cfg.recordCells = ['all'] # record all cells
elif cfg.cellsrec == 1: cfg.recordCells = [(pop,0) for pop in allpops] # record one cell of each pop
elif cfg.cellsrec == 2: cfg.recordCells = [('IT2',10), ('IT5A',10), ('PT5B',10), ('PV5B',10), ('SOM5B',10)] # record selected cells
elif cfg.cellsrec == 3: cfg.recordCells = [(pop,50) for pop in ['IT5A', 'PT5B']]+[('PT5B',x) for x in [393,579,19,104]] #,214,1138,799]] # record selected cells # record selected cells
elif cfg.cellsrec == 4: cfg.recordCells = [(pop,50) for pop in ['IT2', 'IT4', 'IT5A', 'PT5B']] \
										+ [('IT5A',x) for x in [393,447,579,19,104]] \
										+ [('PT5B',x) for x in [393,447,579,19,104,214,1138,979,799]] # record selected cells
elif cfg.cellsrec == 5: 
	cfg.recordCells = [] # record 5 cells of each pop
	for pop in allpops:
		for x in range(5):
			cfg.recordCells.append((pop,x))

cfg.recordTraces = {'V_soma': {'sec':'soma', 'loc':0.5, 'var':'v'}}  ## Dict with traces to record
cfg.recordStim = False			
cfg.recordTime = False  		
cfg.recordStep = 0.1    

#------------------------------------------------------------------------------
# Saving
#------------------------------------------------------------------------------
cfg.simLabel = 'v0_batch0'
cfg.saveFolder = '../data/'+cfg.simLabel
# cfg.filename =                	## Set file output name
cfg.savePickle = True	        	## Save pkl file
cfg.saveJson = False           	## Save json file
cfg.saveDataInclude = ['simData', 'simConfig', 'netParams', 'net'] ## , 'simConfig', 'netParams'
cfg.backupCfgFile = None 		##  
cfg.gatherOnlySimData = False	##  
cfg.saveCellSecs = False			
cfg.saveCellConns = True	

#------------------------------------------------------------------------------
# Analysis and plotting 
#------------------------------------------------------------------------------
with open('../cells/popColors.pkl', 'rb') as fileObj: popColors = pickle.load(fileObj)['popColors']

cfg.analysis['plotRaster'] = {	'include': allpops, 			'orderBy': ['pop', 'y'], 
								'timeRange': [0,cfg.duration], 	'saveFig': True, 			'showFig': False, 
								'labels': 'overlay', 			'popRates': True, 			'orderInverse': True, 
								'popColors': popColors, 		'figSize': (24,20), 		'lw': 0.3, 
								'markerSize':2, 				'marker': '.', 				'dpi': 300} 

cfg.analysis['plotTraces'] = {'include': cfg.recordCells, 'timeRange': [0,cfg.duration], 'ylim': [-100,55],'overlay': False, 'oneFigPer': 'cell', 'figSize': (10,4), 'saveFig': True, 'showFig': False} 

cfg.analysis['plot2Dnet']   = {'include': allpops, 'saveFig': True, 'showConns': False, 'figSize': (24,24), 'fontSize':16}   # Plot 2D cells xy


allpopsplotConn = [	'NGF1', 'IT2', 'SOM2', 'PV2', 'VIP2', 'NGF2', 
				'IT4', 'SOM4', 'PV4', 'VIP4', 'NGF4', 'IT5A', 
				'SOM5A', 'PV5A', 'VIP5A', 'NGF5A', 'IT5B', 
				'PT5B', 'SOM5B', 'PV5B', 'VIP5B', 'NGF5B', 
				'IT6', 'CT6', 'SOM6', 'PV6', 'VIP6', 'NGF6',	
				'VPL_sTC', 	'VPM_sTC', 		'POm_sTC_s1', 
				'VL_sTC', 	'VM_sTC_m1', 	'POm_sTC_m1',
				'mt_RTN', 	'ss_RTN_o', 	'ss_RTN_m', 	'ss_RTN_i',
				'L23_PC_cAD', 'L4_PC_cAD', 'L4_SS_cAD', 'L4_SP_cAD', 
             	'L5_TTPC1_cAD', 'L5_TTPC2_cAD', 'L5_STPC_cAD', 'L5_UTPC_cAD',
             	'L6_TPC_L1_cAD', 'L6_TPC_L4_cAD', 'L6_BPC_cAD', 'L6_IPC_cAD', 'L6_UTPC_cAD']

cfg.analysis['plotConn'] = {'includePre': allpopsplotConn, 
							'includePost': allpopsplotConn, 
							'feature': 'probability', 
							'figSize': (20,20), 
							'groupBy': 'pop', \
 							# 'graphType': 'bar', 
							'synOrConn': 'conn', 
							'synMech': None, 
							'saveData': True, 
							'saveFig': True, 
							'showFig': False
							}
#------------------------------------------------------------------------------
# Synapses
#------------------------------------------------------------------------------

cfg.intraThalamicGain = 1.0
cfg.addTopographicalConn = 1

# Used for NetStim and GroupNetStim implementation 
cfg.synWeightFractionEE = [0.5, 0.5] # E->E AMPA to NMDA ratio
cfg.synWeightFractionEI = [0.5, 0.5] # E->I AMPA to NMDA ratio
cfg.synWeightFractionIE = [0.9, 0.1]  # SOM -> E GABAASlow to GABAB ratio (update this) 
cfg.synWeightFractionII = [0.9, 0.1]  # SOM -> E GABAASlow to GABAB ratio (update this)
cfg.synWeightFractionSOME = [0.9, 0.1] # SOM -> E GABAASlow to GABAB ratio
cfg.synWeightFractionNGF = [0.5, 0.5] # NGF GABAA to GABAB ratio

cfg.synsperconn = {'HH_full': 5, 'HH_reduced': 1, 'HH_simple': 1}

cfg.addSynMechs = True
cfg.distributeSynsUniformly = True



#------------------------------------------------------------------------------
# Connectivity
#------------------------------------------------------------------------------

# S1 Spontaneous synapses + background - data from Rat
cfg.addStimSynS1 = True
cfg.rateStimES1= 9.0
cfg.rateStimIS1 = 9.0
#---------------------
## S1->S1
cfg.connect_S1_S1 = True
#------------------------------------------------------------------------------
## Th->S1
cfg.connect_Th_S1 = True
cfg.TC_S1 = {}
cfg.TC_S1['VPL_sTC'] = True
cfg.TC_S1['VPM_sTC'] = True
cfg.TC_S1['POm_sTC_s1'] = True
#------------------------------------------------------------------------------
## S1->Th 
cfg.connect_S1_Th = True

cfg.connect_S1_RTN = True
cfg.convergence_S1_RTN         = 30.0  # dist_2D<R
cfg.connWeight_S1_RTN       = 0.500

cfg.connect_S1_TC = True
cfg.convergence_S1_TC         = 30.0  # dist_2D<R
cfg.connWeight_S1_TC       = 0.250

#------------------------------------------------------------------------------
## Cortico-cortical connections - 2022_01_21
cfg.connect_M1_S1 = 1
cfg.connect_S1_M1 = 1

#------------------------------------------------------------------------------
# Network 
#------------------------------------------------------------------------------
cfg.singleCellPops = 0  # Create pops with 1 single cell (to debug)
cfg.weightNorm = 1  # use weight normalization
cfg.weightNormThreshold = 4.0  # weight normalization factor threshold


# ===== Actively changing parameters ===== #

cfg.M1_pops=[	'NGF1', 'IT2', 'SOM2', 'PV2', 'VIP2', 'NGF2', 
				'IT4', 'SOM4', 'PV4', 'VIP4', 'NGF4', 'IT5A', 
				'SOM5A', 'PV5A', 'VIP5A', 'NGF5A', 'IT5B', 
				'PT5B', 'SOM5B', 'PV5B', 'VIP5B', 'NGF5B', 
				'IT6', 'CT6', 'SOM6', 'PV6', 'VIP6', 'NGF6',]

cfg.S1_pops= cfg.S1cells
				
cfg.Th_pops=[	
				'VPL_sTC', 	'VPM_sTC', 		'POm_sTC_s1', 
				'VL_sTC', 	'VM_sTC_m1', 	'POm_sTC_m1',
				'mt_RTN', 	'ss_RTN_o', 	'ss_RTN_m', 	'ss_RTN_i']

# ----- Network Parameters ------ #
cfg.removeM1 = False # removes M1 pops
cfg.removeS1 = False # removes M1 pops
cfg.removeTh = False # removes Th pops
cfg.scaleDensity = 0.1 # 1.0

cfg.addThalSs=1
cfg.addThalMt=1

# ----- Network Connections ----- #
cfg.addConn 				= True
cfg.addSubConn 				= True
cfg.addLongConn 			= True
cfg.connectThalamusNetwork 	= True

# Connections under cfg.connectThalamusNetwork
cfg.connect_RTN_RTN     = 1
cfg.connect_TC_RTN      = 1
cfg.connect_RTN_TC      = 1
cfg.connect_TC_CTX      = 1
cfg.connect_CTX_TC      = 1
# ------------------------------- #

cfg.allowConnsWithWeight0 = True
cfg.allowSelfConns = False
cfg.scale = 1.0
cfg.sizeY = 1350.0
cfg.sizeX = 300.0
cfg.sizeZ = 300.0
cfg.correctBorderThreshold = 150.0 # M1
cfg.sizeYS1 = 1378.8 # resized to 1350.0=M1

cfg.L5BrecurrentFactor = 1.0
cfg.ITinterFactor = 1.0
cfg.strengthFactor = 1.0

cfg.EEGain = 0.530860873959182
cfg.EIGain = 1.0
cfg.IEGain = 1.0
cfg.IIGain = 1.0

## E->I by target cell type
cfg.EICellTypeGain= {'PV': 2.588295268601415, 'SOM': 0.6568380849927258, 'VIP': 1.4582025338644486, 'NGF': 3.355557614291127}

# t_allpops joao parameters
cfg.yConnFactor             = 10
cfg.connProb_RTN_RTN        = 0.5
cfg.connProb_TC_RTN         = 1.0
cfg.connProb_RTN_TC         = 1.0
cfg.connWeight_RTN_RTN      = 1.0
cfg.connWeight_TC_RTN       = 1.5
cfg.connWeight_RTN_TC       = 0.5

cfg.divergenceHO = 10

cfg.connLenghtConst = 200

cfg.sRE_model = '../cells/sRE_jv_00.json'
cfg.sTC_model = '../cells/sTC_jv_00.json'

cfg.addBicuculline = False

#------------------------------------------------------------------------------
## I->E/I layer weights (L2/3+4, L5, L6)
cfg.IEweights = [0.5175411466399648, 0.7434834613857577, 1.0101817500320014]
cfg.IIweights = [1.449601171855032, 0.7831317900654744, 1.141724408254077]
#------------------------------------------------------------------------------
# Subcellular distribution
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Long range inputs
#------------------------------------------------------------------------------
cfg.numCellsLong = int(1000 * cfg.scaleDensity) # num of cells per population
cfg.noiseLong = 0.0 #1.0  # firing rate random noise
cfg.delayLong = 5.0  # (ms)
factor = 1
cfg.weightLong = {	
					'POm_sTC_m1': 0.0*factor, 
					'VL_sTC': 0.0*factor, 
					'S1': 0.5*factor, 
					'S2': 0.5*factor, 
					'cM1': 0.5*factor, 
					'M2': 0.5*factor, 
					'OC': 0.5*factor}  # corresponds to unitary connection somatic EPSP (mV)
cfg.startLong = 0  # start at 0 ms
cfg.ratesLong = {	'POm_sTC_m1': [0,0.005], 
					'VL_sTC': [0,0.005], 
					'S1': [0,5], 
					'S2': [0,5], 
					'cM1': [0,5], 
					'M2': [0,5], 
					'OC': [0,5]}

cfg.weightLong_thalM1=0.5
cfg.weightLong_M1thal=0.5

cfg.ratesLong_thalM1 = {	
							'POm_sTC_m1': [0,0.005], 
							'VL_sTC': [0,0.005], 
							}

## input pulses
cfg.addPulses = 0
cfg.pulse = {'pop': 'None', 'start': 1000, 'end': 1100, 'rate': 20, 'noise': 0.8}
cfg.pulse2 = {'pop': 'None', 'start': 1000, 'end': 1200, 'rate': 20, 'noise': 0.5, 'duration': None}


#------------------------------------------------------------------------------
# Current inputs 
#------------------------------------------------------------------------------
cfg.addIClamp = 0

cfg.IClamp1 = {'pop': 'POm_sTC_m1', 'sec': 'soma', 'loc': 0.5, 'start': 0, 'dur': 1000, 'amp': 0.50}
# cfg.IClamp1 = {'pop': 'IT5B', 'sec': 'soma', 'loc': 0.5, 'start': 0, 'dur': 1000, 'amp': 0.50}


#------------------------------------------------------------------------------
# NetStim inputs 
#------------------------------------------------------------------------------
cfg.addNetStim = 0

 			   ## pop, sec, loc, synMech, start, interval, noise, number, weight, delay 
# cfg.NetStim1 = {'pop': 'IT2', 'sec': 'soma', 'loc': 0.5, 'synMech': ['AMPA','NMDA'], 'synMechWeightFactor': cfg.synWeightFractionEE,
# 				'start': 500, 'interval': 50.0, 'noise': 0.2, 'number': 1000.0/50.0, 'weight': 10.0, 'delay': 1}
cfg.NetStim1 = {'pop': 'IT2', 'ynorm':[0,1], 'sec': 'soma', 'loc': 0.5, 'synMech': ['AMPA'], 'synMechWeightFactor': [1.0],
				'start': 500, 'interval': 1000.0/60.0, 'noise': 0.0, 'number': 60.0, 'weight': 30.0, 'delay': 0}
