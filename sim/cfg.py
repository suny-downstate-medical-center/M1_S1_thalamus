"""
cfg.py 

Simulation configuration for M1 model (using NetPyNE)

Contributors: salvadordura@gmail.com
"""

from netpyne import specs
import pickle

cfg = specs.SimConfig()  

#------------------------------------------------------------------------------
#
# SIMULATION CONFIGURATION
#
#------------------------------------------------------------------------------

cfg.simType='network_setup'

#------------------------------------------------------------------------------
# Run parameters
#------------------------------------------------------------------------------
cfg.duration = 0.1 
# cfg.duration = 0.4*1e3 
# cfg.duration = 2.0*1e3 
# cfg.duration = 5.0*1e3 
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
cfg.printRunTime = 0.01
# cfg.printRunTime = 0.1
cfg.printSynsAfterRule = False
cfg.pt3dRelativeToCellLocation = True
cfg.oneSynPerNetcon = False  # only affects conns not in subconnParams; produces identical results

cfg.includeParamsLabel = False
cfg.printPopAvgRates = [0, cfg.duration]

cfg.checkErrors = False
cfg.checkErrorsVerbose = False

cfg.rand123GlobalIndex = None

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

cfg.cellsrec = 1
if cfg.cellsrec == 0:  cfg.recordCells = ['all'] # record all cells
elif cfg.cellsrec == 1: cfg.recordCells = [(pop,0) for pop in allpops] # record one cell of each pop
elif cfg.cellsrec == 2: cfg.recordCells = [('IT2',10), ('IT5A',10), ('PT5B',10), ('PV5B',10), ('SOM5B',10)] # record selected cells
elif cfg.cellsrec == 3: cfg.recordCells = [(pop,50) for pop in ['IT5A', 'PT5B']]+[('PT5B',x) for x in [393,579,19,104]] #,214,1138,799]] # record selected cells # record selected cells
elif cfg.cellsrec == 4: cfg.recordCells = [(pop,50) for pop in ['IT2', 'IT4', 'IT5A', 'PT5B']] \
										+ [('IT5A',x) for x in [393,447,579,19,104]] \
										+ [('PT5B',x) for x in [393,447,579,19,104,214,1138,979,799]] # record selected cells
cfg.recordTraces = {'V_soma': {'sec':'soma', 'loc':0.5, 'var':'v'}}#,
					# 'V_apic_23': {'sec':'apic_23', 'loc':0.5, 'var':'v', 'conds':{'pop': 'PT5B'}},
					# 'V_apic_26': {'sec':'apic_26', 'loc':0.5, 'var':'v', 'conds':{'pop': 'PT5B'}},
					# 'V_dend_5': {'sec':'dend_5', 'loc':0.5, 'var':'v', 'conds':{'pop': 'PT5B'}}}
					#'I_AMPA_Adend2': {'sec':'Adend2', 'loc':0.5, 'synMech': 'AMPA', 'var': 'i'}}

#cfg.recordLFP = [[150, y, 150] for y in range(200,1300,100)]

cfg.recordStim = False
cfg.recordTime = False  
cfg.recordStep = cfg.dt


#------------------------------------------------------------------------------
# Saving
#------------------------------------------------------------------------------
simNum='sim_09'
simFlag='fullModel_density_100_conn_00'
simLabel='tm1_joao_12_29_2021_'+simNum
cfg.simLabel = 'init_sims_joao/'+simLabel+'/'+simNum+simFlag
cfg.saveFolder = '../data'
cfg.savePickle = True
cfg.saveJson = False
cfg.saveDataInclude = ['simData', 'simConfig', 'netParams', 'net']
cfg.backupCfgFile = None #['cfg.py', 'backupcfg/'] 
cfg.gatherOnlySimData = False
# cfg.saveCellSecs = 1
# cfg.saveCellConns = 1
cfg.saveCellSecs = 0
cfg.saveCellConns = 1
cfg.compactConnFormat = 0

#------------------------------------------------------------------------------
# Analysis and plotting 
#------------------------------------------------------------------------------
with open('../cells/popColors.pkl', 'rb') as fileObj: popColors = pickle.load(fileObj)['popColors']
# cfg.analysis['plotRaster'] = {'include': allpops, 'orderBy': ['pop', 'y'], 'timeRange': [0,cfg.duration], 'saveFig': True, 'showFig': False, 'labels': 'overlay', 'popRates': True, 'orderInverse': True, 'popColors': popColors, 'figSize': (12,10), 'lw': 0.3, 'markerSize':3, 'marker': '.', 'dpi': 300} 
cfg.analysis['plotRaster'] = {	'include': allpops, 			'orderBy': ['pop', 'y'], 
								'timeRange': [0,cfg.duration], 	'saveFig': True, 			'showFig': False, 
								'labels': 'overlay', 			'popRates': True, 			'orderInverse': True, 
								'popColors': popColors, 		'figSize': (12,10), 		'lw': 0.3, 
								'markerSize':2, 				'marker': '.', 				'dpi': 300} 


# cfg.analysis['plotSpikeHist'] = {'include': ['IT2','IT4','IT5A','IT5B','PT5B','IT6','CT6'], 'timeRange': [1000,6000], 'yaxis':'rate', 'binSize':5, 'graphType':'bar',
#  								'saveFig': True, 'showFig': False, 'popColors': popColors, 'figSize': (10,4), 'dpi': 300} 

# cfg.analysis['plotLFP'] = {'plots': ['spectrogram'], 'figSize': (6,10), 'timeRange': [1000,6000], 'NFFT': 256*20, 'noverlap': 128*20, 'nperseg': 132*20, 
# 							'saveFig': True, 'showFig':False} 


# cfg.analysis['plotTraces'] = {'include': cfg.recordCells, 'timeRange': [0,cfg.duration], 'ylim': [-100,55],'overlay': False, 'oneFigPer': 'cell', 'figSize': (10,4), 'saveFig': True, 'showFig': False} 

#cfg.analysis['plotShape'] = {'includePre': ['all'], 'includePost': [('PT5B',100)], 'cvar':'numSyns','saveFig': True, 'showFig': False, 'includeAxon': False}
#cfg.analysis['plotConn'] = {'include': ['allCells']}
# cfg.analysis['calculateDisynaptic'] = True

cfg.analysis['plotConn'] = {'includePre': allpops, 
							'includePost': allpops, 
							'feature': 'probability', 
							'figSize': (20,20), 
							'groupBy': 'pop', \
 							# 'graphType': 'bar', 
							'synOrConn': 'conn', 
							'synMech': None, 
							'saveData': None, 
							'saveFig': 1, 
							'showFig': 0
							}

# cfg.analysis['plotConn'] = {'includePre': allpops, 
# 							'includePost': allpops, 
# 							'feature': 'probability', 
# 							'figSize': (20,20), 
# 							'groupBy': 'pop', \
#  							# 'graphType': 'bar', 
# 							'synOrConn': 'conn', 
# 							'synMech': None, 
# 							'saveData': None, 
# 							'saveFig': 1, 
# 							'showFig': 0
# 							}

# cfg.analysis['plotConn'] = {'includePre': allpops, 
# 							'includePost': allpops, 
# 							'feature': 'weight', 
# 							'figSize': (20,20), 
# 							'groupBy': 'pop', \
#  							# 'graphType': 'bar', 
# 							'synOrConn': 'conn', 
# 							'synMech': None, 
# 							'saveData': None, 
# 							'saveFig': 1, 
# 							'showFig': 0
# 							}

# cfg.analysis['plot2Dnet']   = { 'saveFig': True, 'showConns':True, 'figSize': (20,15)}   # Plot 2D net cells and connections


#------------------------------------------------------------------------------
# Cells
#------------------------------------------------------------------------------
cfg.cellmod =  {'IT2': 'HH_reduced',
				'IT4': 'HH_reduced',
				'IT5A': 'HH_full',
				'IT5B': 'HH_reduced',
				'PT5B': 'HH_full',
				'IT6': 'HH_reduced',
				'CT6': 'HH_reduced'}

cfg.ihModel = 'migliore'  # ih model
cfg.ihGbar = 1.0  # multiplicative factor for ih gbar in PT cells - Joao 2021-12-14 - Restored the value from the M1 master branch
# cfg.ihGbar = 0.75  # multiplicative factor for ih gbar in PT cells //// REVERTED CHANGE - Github commit comment: actionpotential  actionpotential, 9 months ago   (March 11th, 2021 12:17pm) - Inserted parameters from best trial (Trial_188)
cfg.ihGbarZD = None # multiplicative factor for ih gbar in PT cells
cfg.ihGbarBasal = 1.0 # 0.1 # multiplicative factor for ih gbar in PT cells
cfg.ihlkc = 0.2 # ih leak param (used in Migliore)
cfg.ihlkcBasal = 1.0
cfg.ihlkcBelowSoma = 0.01
cfg.ihlke = -86  # ih leak param (used in Migliore)
cfg.ihSlope = 14*2

cfg.removeNa = False  # simulate TTX; set gnabar=0s
cfg.somaNa = 5
cfg.dendNa = 0.3
cfg.axonNa = 7
cfg.axonRa = 0.005

cfg.gpas = 0.5  # multiplicative factor for pas g in PT cells
cfg.epas = 0.9  # multiplicative factor for pas e in PT cells

#------------------------------------------------------------------------------
# Synapses
#------------------------------------------------------------------------------

cfg.addIntraThalamicConn = 0  
cfg.intraThalamicGain = 1.0
cfg.addTopographicalConn = 1
cfg.removeWeightNorm = 0 

# Used for NetStim and GroupNetStim implementation 
cfg.synWeightFractionEE = [0.5, 0.5] # E->E AMPA to NMDA ratio
cfg.synWeightFractionEI = [0.5, 0.5] # E->I AMPA to NMDA ratio
cfg.synWeightFractionIE = [0.9, 0.1]  # SOM -> E GABAASlow to GABAB ratio (update this) 
cfg.synWeightFractionII = [0.9, 0.1]  # SOM -> E GABAASlow to GABAB ratio (update this)
cfg.synWeightFractionSOME = [0.9, 0.1] # SOM -> E GABAASlow to GABAB ratio
cfg.synWeightFractionNGF = [0.5, 0.5] # NGF GABAA to GABAB ratio

cfg.synsperconn = {'HH_full': 5, 'HH_reduced': 1, 'HH_simple': 1}
cfg.AMPATau2Factor = 1.0

cfg.addSynMechs = True
cfg.distributeSynsUniformly = True

#------------------------------------------------------------------------------
# Network 
#------------------------------------------------------------------------------
cfg.singleCellPops = 0  # Create pops with 1 single cell (to debug)
cfg.singlePop = 0 
cfg.weightNorm = 1  # use weight normalization
cfg.weightNormThreshold = 4.0  # weight normalization factor threshold


# ===== Actively changing parameters ===== #

cfg.M1_pops=[	'NGF1', 'IT2', 'SOM2', 'PV2', 'VIP2', 'NGF2', 
				'IT4', 'SOM4', 'PV4', 'VIP4', 'NGF4', 'IT5A', 
				'SOM5A', 'PV5A', 'VIP5A', 'NGF5A', 'IT5B', 
				'PT5B', 'SOM5B', 'PV5B', 'VIP5B', 'NGF5B', 
				'IT6', 'CT6', 'SOM6', 'PV6', 'VIP6', 'NGF6',]

cfg.S1_pops=[	
				# ADD S1 POPS HERE
			]
				
cfg.Th_pops=[	
				'VPL_sTC', 	'VPM_sTC', 		'POm_sTC_s1', 
				'VL_sTC', 	'VM_sTC_m1', 	'POm_sTC_m1',
				'mt_RTN', 	'ss_RTN_o', 	'ss_RTN_m', 	'ss_RTN_i']

# ----- Network Parameters ------ #
cfg.removeM1=0 # removes M1 pops
cfg.removeS1=0 # removes M1 pops
cfg.removeTh=0 # removes Th pops
cfg.scaleDensity = 0.02 # 1.0
# cfg.scaleDensity = 0.5 # 1.0
# cfg.scaleDensity = 1.0 # 1.0

cfg.addThalSs=1
cfg.addThalMt=1

# ----- Network Connections ----- #
cfg.addConn 				= 1
cfg.addSubConn 				= 1
cfg.addLongConn 			= 1
cfg.connectThalamusNetwork 	= 1

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
cfg.correctBorderThreshold = 150.0

cfg.L5BrecurrentFactor = 1.0
cfg.ITinterFactor = 1.0
cfg.strengthFactor = 1.0

cfg.EEGain = 0.530860873959182
cfg.EIGain = 1.0
cfg.IEGain = 1.0
cfg.IIGain = 1.0

## E->I by target cell type
cfg.EICellTypeGain= {'PV': 2.588295268601415, 'SOM': 0.6568380849927258, 'VIP': 1.4582025338644486, 'NGF': 3.355557614291127}

cfg.IEdisynapticBias = None  # increase prob of I->Ey conns if Ex->I and Ex->Ey exist 

# # Don paramters
# cfg.yConnFactor             = 10
# cfg.connProb_RTN_RTN        = 1.0
# cfg.connWeight_RTN_RTN      = 2.0
# cfg.connProb_TC_RTN         = 0.75
# cfg.connWeight_TC_RTN       = 1.5
# cfg.connProb_RTN_TC         = 0.75
# cfg.connWeight_RTN_TC       = 0.25

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

cfg.nothing=False
cfg.addBicuculline = False

#------------------------------------------------------------------------------
## (deprecated) E->I gains 
cfg.EPVGain 	= 1.0
cfg.ESOMGain 	= 1.0

#------------------------------------------------------------------------------
## (deprecated) I->E gains
cfg.PVEGain 	= 1.0
cfg.SOMEGain 	= 1.0

#------------------------------------------------------------------------------
## (deprecated) I->I gains
cfg.PVSOMGain 	= None #0.25
cfg.SOMPVGain 	= None #0.25
cfg.PVPVGain 	= None # 0.75
cfg.SOMSOMGain 	= None #0.75

#------------------------------------------------------------------------------
## I->E/I layer weights (L2/3+4, L5, L6)
cfg.IEweights = [0.5175411466399648, 0.7434834613857577, 1.0101817500320014]
cfg.IIweights = [1.449601171855032, 0.7831317900654744, 1.141724408254077]

cfg.IPTGain = 1.0
cfg.IFullGain = 1.0  # deprecated

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
					# 'POm_sTC_m1': 0.5*factor, 
					# 'VL_sTC': 0.5*factor, 
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
# cfg.weightLong = {'POm_sTC_m1': 0.5*factor, 'VL_sTC': 0.5*factor, 'S1': 0.5*factor, 'S2': 0.5*factor, 'cM1': 0.5*factor, 'M2': 0.5*factor, 'OC': 0.5*factor}  # corresponds to unitary connection somatic EPSP (mV)
# cfg.startLong = 0  # start at 0 ms
# cfg.ratesLong = {'POm_sTC_m1': [0,0.005], 'VL_sTC': [0,0.005], 'S1': [0,5], 'S2': [0,5], 'cM1': [0,5], 'M2': [0,5], 'OC': [0,5]}

# cfg.weightLong_thalM1=0.01
cfg.weightLong_thalM1=0.5

cfg.weightLong_M1thal=0.5

# cfg.weightLong_thalM1 = {	
# 							# 'POm_sTC_m1': 0.0*factor, 
# 							# 'VL_sTC': 0.0*factor, 
# 							'POm_sTC_m1': 0.5*factor, 
# 							'VL_sTC': 0.5*factor, 
							# }

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
