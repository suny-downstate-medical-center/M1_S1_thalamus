"""
batch.py 

Batch simulation for M1 model using NetPyNE

Contributors: salvadordura@gmail.com
"""
from netpyne.batch import Batch
from netpyne import specs
import numpy as np

# ----------------------------------------------------------------------------------------------
# Joao - M1+thalamus batch sims
# ----------------------------------------------------------------------------------------------
def projectionWeights():
    
    # ----- Batch Variables ------ #
    params = specs.ODict()
    
    params['weightLong_thalM1']=[0.1,0.25,0.5,0.75,1.0,1.5,1.75]
    params['weightLong_M1thal']=[0.1,0.25,0.5,0.75,1.0,1.5,1.75]
    
    # ----- Simulation Parameters ------ #
    initCfg = {}
    initCfg['duration']         = 0.1
    initCfg['savePickle']       = True
    initCfg['saveJson']         = False
    initCfg['saveCellSecs']     = False
    initCfg['saveCellConns']    = True

    # ----- Population Parameters ------ #
    initCfg['M1_pops']  =[	'NGF1', 'IT2', 'SOM2', 'PV2', 'VIP2', 'NGF2', 'IT4', 'SOM4', 'PV4', 'VIP4', 'NGF4', 'IT5A', 
                            'SOM5A', 'PV5A', 'VIP5A', 'NGF5A', 'IT5B', 'PT5B', 'SOM5B', 'PV5B', 'VIP5B', 'NGF5B', 
                            'IT6', 'CT6', 'SOM6', 'PV6', 'VIP6', 'NGF6',]
    initCfg['S1_pops']  =[	
                            # ADD S1 POPS HERE
                            ]
    initCfg['Th_pops']  =[	
                            'VPL_sTC', 	'VPM_sTC', 		'POm_sTC_s1', 
                            'VL_sTC', 	'VM_sTC_m1', 	'POm_sTC_m1',
                            'mt_RTN', 	'ss_RTN_o', 	'ss_RTN_m', 	'ss_RTN_i']

    # ----- Network Parameters ------ #
    initCfg['removeM1']         = 0 # removes M1 pops
    initCfg['removeS1']         = 0 # removes M1 pops
    initCfg['removeTh']         = 0 # removes Th pops
    
    initCfg['scaleDensity']     = 0.02 # 1.0
    # initCfg['scaleDensity']     = 1.0 # 1.0

    initCfg['addThalSs']        =1
    initCfg['addThalMt']        =1

    # ----- Network Connections ----- #
    initCfg['addConn']   		        = 1
    initCfg['addSubConn']    	        = 1
    initCfg['addLongConn']   	        = 1
    initCfg['connectThalamusNetwork']   = 1

    # Connections under cfg.connectThalamusNetwork
    initCfg['connect_RTN_RTN']  = 1
    initCfg['connect_TC_RTN']   = 1
    initCfg['connect_RTN_TC']   = 1
    initCfg['connect_TC_CTX']   = 1
    initCfg['connect_CTX_TC']   = 1
    # ------------------------------- #

    groupedParams = []

    b = Batch(params=params, netParamsFile='netParams.py', cfgFile='cfg.py', initCfg=initCfg, groupedParams=groupedParams)

    return b

# ----------------------------------------------------------------------------------------------
# Weight Normalization Exc
# ----------------------------------------------------------------------------------------------
def weightNorm(pops=['IT2', 'IT4', 'IT5A', 'IT5B', 'PT5B', 'IT6', 'CT6', 'PV2', 'SOM2'], 
    segs = None, allSegs = True, rule = 'IT2_reduced', weights=list(np.arange(0.01, 0.2, 0.01)/100.0)):

    # Add params
    from cfg_cell import cfg
    from netParams_cell import netParams

    excludeSegs = ['axon']
    if not segs:
        secs = []
        locs = []
        for secName,sec in netParams.cellParams[rule]['secs'].items():
            if secName not in excludeSegs:
                if allSegs:
                    nseg = sec['geom']['nseg']
                    for iseg in range(nseg):
                        secs.append(secName) 
                        locs.append((iseg+1)*(1.0/(nseg+1)))
                else:
                    secs.append(secName) 
                    locs.append(0.5)

    params = specs.ODict()
    params[('NetStim1', 'pop')] = pops
    params[('NetStim1', 'sec')] = secs
    params[('NetStim1', 'loc')] = locs
    params[('NetStim1', 'weight')] = weights

    groupedParams = [('NetStim1', 'sec'), ('NetStim1', 'loc')] 

    initCfg = {}
    initCfg['duration'] = 1.0*1e3
    initCfg[('analysis','plotTraces','timeRange')] = [0, 1000]
    initCfg['weightNorm'] = False
    initCfg['stimSubConn'] = False
    initCfg['addNetStim'] = True
    initCfg[('NetStim1', 'synMech')] = ['AMPA','NMDA']
    initCfg[('NetStim1','synMechWeightFactor')] = [0.5,0.5]
    initCfg[('NetStim1', 'start')] = 700
    initCfg[('NetStim1', 'interval')] = 1000
    initCfg[('NetStim1','ynorm')] = [0.0, 1.0]

    initCfg[('NetStim1', 'noise')] = 0
    initCfg[('NetStim1', 'number')] = 1
    initCfg[('NetStim1', 'delay')] = 1
    #initCfg[('GroupNetStimW1', 'pop')] = 'None'
    initCfg[('NetStim1', 'delay')] = 1
    initCfg['addIClamp'] = 0
    
    b = Batch(params=params, netParamsFile='netParams_cell.py', cfgFile='cfg_cell.py', initCfg=initCfg, groupedParams=groupedParams)

    return b

# ----------------------------------------------------------------------------------------------
# Exc-Inh balance
# ----------------------------------------------------------------------------------------------
def EIbalance():
    params = specs.ODict()

    params[('ratesLong', 'TPO', 1)] = [3, 4] #[2, 4]#[2,4,2,2,4,2,4,4]
    params[('ratesLong', 'TVL', 1)] = [1, 1] #[1, 2]#[2,4,2,2,4,2,4,4]
    params[('ratesLong', 'S1', 1)] =  [3, 4] #[2, 4] #[2,2,4,2,4,4,2,4]
    params[('ratesLong', 'S2', 1)] =  [3, 4] #[2, 4] #[2,2,4,2,4,4,2,4]
    params[('ratesLong', 'cM1', 1)] = [1, 1] #[1, 2]#[2,2,2,4,2,4,4,4]
    params[('ratesLong', 'M2', 1)] =  [1, 1] #[1, 2] #[2,2,2,4,2,4,4,4]
        
    params['EEGain'] = [0.5,0.6] #[0.6, 0.8] #[0.8, 1.0]

    # # L2/3+4
    params[('IEweights',0)] =  [0.8, 1.0, 1.2]
    params[('IIweights',0)] =  [0.8, 1.0, 1.2]   
    # L5
    params[('IEweights',1)] = [0.8, 1.0, 1.2]   
    params[('IIweights',1)] =  [0.8, 1.0, 1.2]
    # L6
    # params[('IEweights',2)] =  [0.8, 1.0, 1.2]  
    # params[('IIweights',2)] =  [0.8, 1.0, 1.2]

    params['ihGbar'] = [0.0, 0.25, 1.0]

    groupedParams = [('ratesLong', 'TPO', 1), ('ratesLong', 'TVL', 1),
                    ('ratesLong', 'S1', 1), ('ratesLong', 'S2', 1),
                    ('ratesLong', 'cM1', 1), ('ratesLong', 'M2', 1)]# ['IEGain','IIGain'] #'EEGain', 'EPVGain', 'ESOMGain', 'PVEGain', 'SOMEGain', 'PVIGain', 'SOMIGain']

    # initial config
    initCfg = {}
    initCfg['duration'] = 2.0*1e3
    initCfg['ihModel'] = 'migliore'  # ih model

    initCfg['ihGbarBasal'] = 1.0 # multiplicative factor for ih gbar in PT cells
    initCfg['ihlkc'] = 0.2 # ih leak param (used in Migliore)
    initCfg['ihLkcBasal'] = 1.0 # multiplicative factor for ih lk in PT cells
    initCfg['ihLkcBelowSoma'] = 0.01 # multiplicative factor for ih lk in PT cells
    initCfg['ihlke'] = -86  # ih leak param (used in Migliore)
    initCfg['ihSlope'] = 28  # ih leak param (used in Migliore)

    initCfg['somaNa'] = 5.0  # somatic Na conduct
    initCfg['dendNa'] = 0.3  # dendritic Na conduct (reduced to avoid dend spikes) 
    initCfg['axonNa'] = 7   # axon Na conduct (increased to compensate) 
    initCfg['axonRa'] = 0.005
    initCfg['gpas'] = 0.5
    initCfg['epas'] = 0.9

    initCfg[('pulse', 'pop')] = 'S2'
    initCfg[('pulse', 'rate')] = 10.0
    initCfg[('pulse', 'start')] = 1000.0
    initCfg[('pulse', 'end')] = 1100.0
    initCfg[('pulse', 'noise')] = 0.8

    initCfg['IEdisynapticBias'] = None

    initCfg['weightNormThreshold'] = 4.0
    initCfg['IEGain'] = 1.0
    initCfg['IIGain'] = 1.0
    initCfg['IPTGain'] = 1.0
    initCfg['IIweights'] =  [1.0, 1.0, 1.0]

    initCfg['saveCellSecs'] = False
    initCfg['saveCellConns'] = False
    
    b = Batch(params=params, groupedParams=groupedParams, initCfg=initCfg)

    return b


# ----------------------------------------------------------------------------------------------
# Exc-Inh balance
# ----------------------------------------------------------------------------------------------
def longBalance():
    params = specs.ODict()

    params[('ratesLong', 'TPO', 1)] = [2,4]
    params[('ratesLong', 'TVL', 1)] = [2,4]
    params[('ratesLong', 'S1', 1)] = [2,4]
    params[('ratesLong', 'S2', 1)] = [2,4]
    params[('ratesLong', 'cM1', 1)] = [2,4]
    params[('ratesLong', 'M2', 1)] = [2,4]
    params[('ratesLong', 'OC', 1)] = [2,4]

    # 
    params['IEweights'] = [[0.8,0.8,0.8], [1.0,1.0,1.0], [1.2,1.2,1.2]]
    params['IIweights'] =  [[0.8,0.8,0.80], [1.0, 1.0, 1.0], [1.2,1.2,1.2]]

    params['ihGbar'] = [0.25, 1.0]

    groupedParams = []

    # initial config
    initCfg = {}
    initCfg['duration'] = 2.0*1e3
    initCfg['ihModel'] = 'migliore'  # ih model

    initCfg['ihGbarBasal'] = 1.0 # multiplicative factor for ih gbar in PT cells
    initCfg['ihlkc'] = 0.2 # ih leak param (used in Migliore)
    initCfg['ihLkcBasal'] = 1.0 # multiplicative factor for ih lk in PT cells
    initCfg['ihLkcBelowSoma'] = 0.01 # multiplicative factor for ih lk in PT cells
    initCfg['ihlke'] = -86  # ih leak param (used in Migliore)
    initCfg['ihSlope'] = 28  # ih leak param (used in Migliore)

    initCfg['somaNa'] = 5.0  # somatic Na conduct
    initCfg['dendNa'] = 0.3  # dendritic Na conduct (reduced to avoid dend spikes) 
    initCfg['axonNa'] = 7   # axon Na conduct (increased to compensate) 
    initCfg['axonRa'] = 0.005
    initCfg['gpas'] = 0.5
    initCfg['epas'] = 0.9

    initCfg[('pulse', 'pop')] = 'S2'
    initCfg[('pulse', 'rate')] = 10.0
    initCfg[('pulse', 'start')] = 1000.0
    initCfg[('pulse', 'end')] = 1100.0
    initCfg[('pulse', 'noise')] = 0.8

    initCfg['IEdisynapticBias'] = None

    initCfg['weightNormThreshold'] = 4.0
    initCfg['IEGain'] = 1.0
    initCfg['IIGain'] = 1.0
    initCfg['IPTGain'] = 1.0

    initCfg['saveCellSecs'] = False
    initCfg['saveCellConns'] = False
    
    b = Batch(params=params, groupedParams=groupedParams, initCfg=initCfg)

    return b

# ----------------------------------------------------------------------------------------------
# Long-range pop stimulation
# ----------------------------------------------------------------------------------------------
def longPopStims():
    params = specs.ODict()
    
    params['ihGbar'] = [0.25, 1.0] # [0.2, 0.25, 0.3, 1.0]
    params[('seeds', 'conn')] = [4321] #+(17*i) for i in range(5)]
    params[('seeds', 'stim')] = [1234] #+(17*i) for i in range(5)]

    params[('pulse', 'pop')] = ['None'] #, 'TPO', 'TVL', 'S2', 'M2'] #, 'OC'] # 'S1','cM1',
    #params[('pulse', 'end')] = [1100, 1500]

    groupedParams = []

    # initial config
    initCfg = {}
    initCfg['duration'] = 2.*1e3 #2.5*1e3
    initCfg['ihModel'] = 'migliore'  # ih model

    initCfg['ihGbarBasal'] = 1.0 # multiplicative factor for ih gbar in PT cells
    initCfg['ihlkc'] = 0.2 # ih leak param (used in Migliore)
    initCfg['ihLkcBasal'] = 1.0 # multiplicative factor for ih lk in PT cells
    initCfg['ihLkcBelowSoma'] = 0.01 # multiplicative factor for ih lk in PT cells
    initCfg['ihlke'] = -86  # ih leak param (used in Migliore)
    initCfg['ihSlope'] = 28  # ih leak param (used in Migliore)

    initCfg['somaNa'] = 5.0  # somatic Na conduct
    initCfg['dendNa'] = 0.3  # dendritic Na conduct (reduced to avoid dend spikes) 
    initCfg['axonNa'] = 7   # axon Na conduct (increased to compensate) 
    initCfg['axonRa'] = 0.005
    initCfg['gpas'] = 0.5
    initCfg['epas'] = 0.9

    #initCfg[('pulse', 'pop')] = 'None'
    initCfg[('pulse', 'rate')] = 10.0
    initCfg[('pulse', 'start')] = 1000.0
    initCfg[('pulse', 'end')] = 1100.0
    initCfg[('pulse', 'noise')] = 0.8

    initCfg['IEdisynapticBias'] = None

    initCfg['weightNormThreshold'] = 4.0
    initCfg['EEGain'] = 0.5 
    initCfg['IEGain'] = 1.0
    initCfg['IIGain'] = 1.0
    initCfg['IPTGain'] = 1.0

    initCfg[('ratesLong', 'TPO', 1)] = 5 	
    initCfg[('ratesLong', 'TVL', 1)] = 2.5
    initCfg[('ratesLong', 'S1', 1)] = 5
    initCfg[('ratesLong', 'S2', 1)] = 5 
    initCfg[('ratesLong', 'cM1', 1)] = 2.5
    initCfg[('ratesLong', 'M2', 1)] = 2.5
    initCfg[('ratesLong', 'OC', 1)] = 5	

    # # L2/3+4
    initCfg[('IEweights',0)] =  0.8
    initCfg[('IIweights',0)] =  1.2 
    # L5
    initCfg[('IEweights',1)] = 0.8   
    initCfg[('IIweights',1)] = 1.0
    # L6
    initCfg[('IEweights',2)] =  1.0  
    initCfg[('IIweights',2)] =  1.0

    initCfg['saveCellSecs'] = False
    initCfg['saveCellConns'] = False

    groupedParams = [] #('IEweights',0), ('IIweights',0), ('IEweights',1), ('IIweights',1), ('IEweights',2), ('IIweights',2)]

    b = Batch(params=params, initCfg=initCfg, groupedParams=groupedParams)

    return b

# ----------------------------------------------------------------------------------------------
# Simultaenous long-range pop stimulations
# ----------------------------------------------------------------------------------------------
def simultLongPopStims():
    params = specs.ODict()
    
    params[('pulse', 'pop')] = ['TPO', 'M2', 'TVL', 'S2', 'S2', 'M2', 'TVL', 'TPO']
    params[('pulse2', 'pop')] = ['M2', 'TPO', 'S2', 'TVL', 'M2', 'S2', 'TPO', 'TVL']
    params[('pulse2', 'start')] = list(np.arange(1500, 2020, 20))
    params['ihGbar'] = [0.25, 1.0]


    # initial config
    initCfg = {}
    initCfg['duration'] = 3.0*1e3
    initCfg['ihModel'] = 'migliore'  # ih model

    initCfg['ihGbarBasal'] = 1.0 # multiplicative factor for ih gbar in PT cells
    initCfg['ihlkc'] = 0.2 # ih leak param (used in Migliore)
    initCfg['ihLkcBasal'] = 1.0 # multiplicative factor for ih lk in PT cells
    initCfg['ihLkcBelowSoma'] = 0.01 # multiplicative factor for ih lk in PT cells
    initCfg['ihlke'] = -86  # ih leak param (used in Migliore)
    initCfg['ihSlope'] = 28  # ih leak param (used in Migliore)

    initCfg['somaNa'] = 5.0  # somatic Na conduct
    initCfg['dendNa'] = 0.3  # dendritic Na conduct (reduced to avoid dend spikes) 
    initCfg['axonNa'] = 7   # axon Na conduct (increased to compensate) 
    initCfg['axonRa'] = 0.005
    initCfg['gpas'] = 0.5
    initCfg['epas'] = 0.9

    #initCfg[('pulse', 'pop')] = 'None'
    initCfg[('pulse', 'rate')] = 10.0
    initCfg[('pulse', 'start')] = 1500.0
    initCfg[('pulse', 'end')] = 1700.0
    initCfg[('pulse', 'noise')] = 0.8

    #initCfg[('pulse2', 'start')] = 1500.0
    initCfg[('pulse2', 'rate')] = 10.0
    initCfg[('pulse2', 'duration')] = 200.0
    initCfg[('pulse2', 'noise')] = 0.8


    initCfg['IEdisynapticBias'] = None

    initCfg['weightNormThreshold'] = 4.0
    initCfg['EEGain'] = 0.5 
    initCfg['IEGain'] = 1.0
    initCfg['IIGain'] = 1.0
    initCfg['IPTGain'] = 1.0

    initCfg[('ratesLong', 'TPO', 1)] = 5 	
    initCfg[('ratesLong', 'TVL', 1)] = 2.5
    initCfg[('ratesLong', 'S1', 1)] = 5
    initCfg[('ratesLong', 'S2', 1)] = 5 
    initCfg[('ratesLong', 'cM1', 1)] = 2.5
    initCfg[('ratesLong', 'M2', 1)] = 2.5
    initCfg[('ratesLong', 'OC', 1)] = 5	

    # # L2/3+4
    initCfg[('IEweights',0)] =  0.8
    initCfg[('IIweights',0)] =  1.2 
    # L5
    initCfg[('IEweights',1)] = 0.8   
    initCfg[('IIweights',1)] = 1.0
    # L6
    initCfg[('IEweights',2)] =  1.0  
    initCfg[('IIweights',2)] =  1.0

    initCfg['saveCellSecs'] = False
    initCfg['saveCellConns'] = False

    groupedParams = [('pulse', 'pop'),('pulse2', 'pop')] 
    b = Batch(params=params, initCfg=initCfg, groupedParams=groupedParams)

    return b



# ----------------------------------------------------------------------------------------------
# Recorded stimulation
# ----------------------------------------------------------------------------------------------
def recordedLongPopStims():
    params = specs.ODict()
    
    high = 'cells/ssc-3_spikes.json'
    low  = 'cells/ssc-3_lowrate_spikes.json'
    low2 = 'cells/ssc-3_lowrate2_spikes.json'


    # 1) normal, 2) S2high+lowbkg, 3) S2low+bkg0.1, 4) S2low2+bkg0.1, 5) S2low2+M2low+bkg0.1, 6) S2low, 
    # 7) S2high, 8) S1high, 9) S1low, 10) M2low, 11) M2high
    params[('ratesLong','S2')] =  [[0,2]]#, high,    low,     low2,	 low2,		high, 	low,   [0,2], [0,2], [0,2], [0,2]]
    params[('ratesLong','S1')] =  [[0,2]]#, [0,0.1], [0,0.1], [0,0.1], [0,0.1],	[0,2], 	[0,2], high,  low,   [0,2], [0,2]]
    params[('ratesLong','M2')] =  [[0,2]]#, [0,0.1], [0,0.1], [0,0.1], low, 		[0,2], 	[0,2], [0,2], [0,2], high, 	low]
    params[('ratesLong','TPO')] = [[0,4]]#, [0,0.1], [0,0.1], [0,0.1], [0,0.1],	[0,4],	[0,4], [0,4], [0,4], [0,4], [0,4]]
    params[('ratesLong','TVL')] = [[0,4]]#, [0,0.1], [0,0.1], [0,0.1], [0,0.1],	[0,4],	[0,4], [0,4], [0,4], [0,4], [0,4]]
    params[('ratesLong','cM1')] = [[0,4]]#, [0,0.1], [0,0.1], [0,0.1], [0,0.1],	[0,4],	[0,4], [0,4], [0,4], [0,4], [0,4]]
    params[('ratesLong','OC')] =  [[0,2]]#, [0,0.1], [0,0.1], [0,0.1], [0,0.1],	[0,2], 	[0,2], [0,2], [0,2], [0,2], [0,2]]
    #params['ihGbar'] = [0.3, 0.4, 0.5, 1.0]
    params['ihGbar'] = [0.3] #, 1.0]
    
    # initial config
    initCfg = {}

    initCfg['duration'] = 6.0*1e3
    initCfg['ihModel'] = 'migliore'  # ih model

    initCfg['ihGbarBasal'] = 1.0 # multiplicative factor for ih gbar in PT cells
    initCfg['ihlkc'] = 0.2 # ih leak param (used in Migliore)
    initCfg['ihLkcBasal'] = 1.0 # multiplicative factor for ih lk in PT cells
    initCfg['ihLkcBelowSoma'] = 0.01 # multiplicative factor for ih lk in PT cells
    initCfg['ihlke'] = -86  # ih leak param (used in Migliore)
    initCfg['ihSlope'] = 28  # ih leak param (used in Migliore)

    initCfg['somaNa'] = 5.0
    initCfg['dendNa'] = 0.3  # dendritic Na conduct (reduced to avoid dend spikes) 
    initCfg['axonNa'] = 7   # axon Na conduct (increased to compensate) 
    initCfg['axonRa'] = 0.005
    initCfg['gpas'] = 0.5
    initCfg['epas'] = 0.9

    initCfg[('analysis','plotRaster','timeRange')] = [500, 5500]

    initCfg['weightNormThreshold'] = 4.0

    initCfg['saveCellSecs'] = False
    initCfg['saveCellConns'] = False

    initCfg['IEGain'] = 1.0
    initCfg['IIGain'] = 1.0
    initCfg['IEdisynapticBias'] = None

    # 1101 222

    # # L2/3+4
    initCfg[('IEweights',0)] = 1.2
    initCfg[('IIweights',0)] =  1.0  
    # L5
    initCfg[('IEweights',1)] = 1.2
    initCfg[('IIweights',1)] = 1.0
    # L6
    initCfg[('IEweights',2)] =  1.2  
    initCfg[('IIweights',2)] =  1.0

    # groupedParams = [('ratesLong','S2'), ('ratesLong','S1'), ('ratesLong','M2'), 
    # 				('ratesLong','TPO'), ('ratesLong','TVL'), ('ratesLong','cM1'), ('ratesLong','OC')]
    groupedParams = []

    b = Batch(params=params, initCfg=initCfg, groupedParams=groupedParams)

    return b



# ----------------------------------------------------------------------------------------------
# Frequency stimulation
# ----------------------------------------------------------------------------------------------
def freqStims():
    params = specs.ODict()

    params[('NetStim1', 'interval')] = [1000.0/f for f in [4,8,12,16,20,24,28,32,36,40]]
    params[('NetStim1', 'number')] = [f for f in [4,8,12,16,20,24,28,32,36,40]]	
    params[('NetStim1', 'start')] = [500, 550]
    params['ihGbar'] = [0.5, 1.0]
    params[('NetStim1', 'ynorm', 1)] = [0.15+x*(0.31-0.12) for x in [0.1, 0.2, 0.3]]  # 10, 20, 30% of cells; L23 NCD = 0.12 - 0.31


    # initial config
    initCfg = {}
    initCfg['addNetStim'] = True
    initCfg[('NetStim1', 'pop')] = 'IT2'
    initCfg[('NetStim1', 'ynorm', 0)] = 0.15
    initCfg[('NetStim1', 'weight')] = 30.0	
    initCfg[('NetStim1', 'noise')] = 0.01	

    initCfg['duration'] = 2.0*1e3
    initCfg['ihModel'] = 'migliore'  # ih model

    initCfg['ihGbarBasal'] = 1.0 # multiplicative factor for ih gbar in PT cells
    initCfg['ihlkc'] = 0.2 # ih leak param (used in Migliore)
    initCfg['ihLkcBasal'] = 1.0 # multiplicative factor for ih lk in PT cells
    initCfg['ihLkcBelowSoma'] = 0.01 # multiplicative factor for ih lk in PT cells
    initCfg['ihlke'] = -86  # ih leak param (used in Migliore)
    initCfg['ihSlope'] = 28  # ih leak param (used in Migliore)

    initCfg['somaNa'] = 5.0
    initCfg['dendNa'] = 0.3  # dendritic Na conduct (reduced to avoid dend spikes) 
    initCfg['axonNa'] = 7   # axon Na conduct (increased to compensate) 
    initCfg['axonRa'] = 0.005
    initCfg['gpas'] = 0.5
    initCfg['epas'] = 0.9

    initCfg['weightNormThreshold'] = 4.0

    initCfg['saveCellSecs'] = False
    initCfg['saveCellConns'] = False

    initCfg['IEGain'] = 1.0
    initCfg['IIGain'] = 1.0
    initCfg['IEdisynapticBias'] = None


    # 1101 222
    initCfg[('ratesLong', 'TPO', 1)] = 4
    initCfg[('ratesLong', 'TVL', 1)] = 4
    initCfg[('ratesLong', 'S1', 1)] = 2
    initCfg[('ratesLong', 'cM1', 1)] = 4

    # # L2/3+4
    initCfg[('IEweights',0)] = 1.2
    initCfg[('IIweights',0)] =  1.0  
    # L5
    initCfg[('IEweights',1)] = 1.2
    initCfg[('IIweights',1)] = 1.0
    # L6
    initCfg[('IEweights',2)] =  1.2  
    initCfg[('IIweights',2)] =  1.0
    initCfg[('IIweights',2)] =  0.8

    groupedParams = [('NetStim1', 'interval'), ('NetStim1', 'number')] 

    b = Batch(params=params, initCfg=initCfg, groupedParams=groupedParams)

    return b

# ----------------------------------------------------------------------------------------------
# Local pop stimulation
# ----------------------------------------------------------------------------------------------
def localPopStims():
    params = specs.ODict()

    params['ihGbar'] = [0.0, 1.0, 2.0]
    params[('NetStim1', 'pop')] = ['IT2','IT4','IT5A','IT5B','PT5B','IT6','CT6']
    params[('NetStim1', 'interval')] = [1000.0/20.0, 1000.0/30.0]

    b = Batch(params=params)

    grouped = []

    for p in b.params:
        if p['label'] in grouped: 
            p['group'] = True

    return b


# ----------------------------------------------------------------------------------------------
# EPSPs via NetStim
# ----------------------------------------------------------------------------------------------
def EPSPs():
    params = specs.ODict()

    params['groupWeight'] = [x*0.05 for x in np.arange(1, 8, 1)]
    params['ihGbar'] = [0.0, 1.0]
 
    
    # initial config
    initCfg = {}
    initCfg['duration'] = 0.5*1e3
    initCfg['addIClamp'] = False
    initCfg['addNetStim'] = True
    initCfg[('GroupNetStimW1', 'pop')] = 'PT5B'
    initCfg[('analysis','plotTraces','timeRange')] = [0, 500]
    initCfg['excTau2Factor'] = 2.0
    initCfg['weightNorm'] = True
    initCfg['stimSubConn'] = False
    initCfg['ihGbarZD'] = None

    groupedParams = [] 

    b = Batch(params=params, netParamsFile='netParams_cell.py', cfgFile='cfg_cell.py', initCfg=initCfg, groupedParams=groupedParams)

    return b


# ----------------------------------------------------------------------------------------------
# f-I curve
# ----------------------------------------------------------------------------------------------
def fIcurve():
    params = specs.ODict()

    params[('IClamp1', 'pop')] = ['IT2', 'IT4', 'IT5A', 'IT5B', 'PT5B', 'IT6', 'CT6', 'PV2', 'SOM2']
    params[('IClamp1', 'amp')] = list(np.arange(0.0, 6.5, 0.5)/10.0) 
    #params['ihGbar'] = [0.0, 1.0, 2.0]
    # params['axonNa'] = [5, 6, 7, 8] 
    # params['gpas'] = [0.6, 0.65, 0.70, 0.75] 
    # params['epas'] = [1.0, 1.05] 
    # params['ihLkcBasal'] = [0.0, 0.01, 0.1, 0.5, 1.0] 

    # initial config
    initCfg = {}
    initCfg['duration'] = 1.5*1e3
    initCfg['addIClamp'] = True
    initCfg['addNetStim'] = False
    initCfg['weightNorm'] = True
    initCfg[('IClamp1','sec')] = 'soma'
    initCfg[('IClamp1','loc')] = 0.5
    initCfg[('IClamp1','start')] = 500
    initCfg[('IClamp1','dur')] = 1000
    initCfg[('analysis','plotTraces','timeRange')] = [0, 1500]

    groupedParams = [] 

    b = Batch(params=params, netParamsFile='netParams_cell.py', cfgFile='cfg_cell.py', initCfg=initCfg, groupedParams=groupedParams)

    return b


# ----------------------------------------------------------------------------------------------
# Custom
# ----------------------------------------------------------------------------------------------
def custom():
    params = specs.ODict()
    # long-range inputs
    params[('weightLong', 'TPO')] = [0.5896891377382606] 
    params[('weightLong', 'TVL')] = [0.42705597272266954] 
    params[('weightLong', 'S1')] =  [0.5360314931416853] 
    params[('weightLong', 'S2')] =  [0.4390764895101832] 
    params[('weightLong', 'cM1')] = [0.36758768725488966] 
    params[('weightLong', 'M2')] =  [0.2516229090853658] 
    params[('weightLong', 'OC')] =  [0.6813698683166248]	

    # EEgain
    params['EEGain'] = [0.5157775394221152] 

    # IEgain
    ## L2/3+4
    params[('IEweights',0)] =  [0.7041397442431576]
    ## L5
    params[('IEweights',1)] = [1.1737318529671519] #[0.8, 1.0]   
    ## L6
    params[('IEweights',2)] =  [1.4290869898071543] # [0.8, 1.0]  

    # IIGain
    params['IIGain'] = [1.0]


    groupedParams = [('weightLong', 'TPO'), 
                    ('weightLong', 'TVL'), 
                    ('weightLong', 'S1'), 
                    ('weightLong', 'S2'), 
                    ('weightLong', 'cM1'), 
                    ('weightLong', 'M2'), 
                    ('weightLong', 'OC')] 
    
    # --------------------------------------------------------
    # initial config
    initCfg = {}
    initCfg['duration'] = 1500
    initCfg['printPopAvgRates'] = [500, 1500] 
    initCfg['dt'] = 0.025

    initCfg['scaleDensity'] = 1.0

    # cell params
    initCfg['ihGbar'] = 0.75  # ih (for quiet/sponti condition)
    initCfg['ihModel'] = 'migliore'  # ih model
    initCfg['ihGbarBasal'] = 1.0 # multiplicative factor for ih gbar in PT cells
    initCfg['ihlkc'] = 0.2 # ih leak param (used in Migliore)
    initCfg['ihLkcBasal'] = 1.0 # multiplicative factor for ih lk in PT cells
    initCfg['ihLkcBelowSoma'] = 0.01 # multiplicative factor for ih lk in PT cells
    initCfg['ihlke'] = -86  # ih leak param (used in Migliore)
    initCfg['ihSlope'] = 28  # ih leak param (used in Migliore)

    initCfg['somaNa'] = 5.0  # somatic Na conduct
    initCfg['dendNa'] = 0.3  # dendritic Na conduct (reduced to avoid dend spikes) 
    initCfg['axonNa'] = 7   # axon Na conduct (increased to compensate) 
    initCfg['axonRa'] = 0.005
    initCfg['gpas'] = 0.5
    initCfg['epas'] = 0.9
    
    # long-range input params
    initCfg['numCellsLong'] = 1000
    initCfg[('pulse', 'pop')] = 'None'
    initCfg[('pulse', 'start')] = 1000.0
    initCfg[('pulse', 'end')] = 1100.0
    initCfg[('pulse', 'noise')] = 0.8

    # conn params
    initCfg['IEdisynapticBias'] = None

    initCfg['weightNormThreshold'] = 4.0
    initCfg['IEGain'] = 1.0
    initCfg['IPTGain'] = 1.0
    initCfg['IIweights'] = [1.4915732417052807, 0.8355817935676599, 1.4072445490624554]

    # plotting and saving params
    initCfg[('analysis','plotRaster','timeRange')] = initCfg['printPopAvgRates']
    initCfg[('analysis', 'plotTraces', 'timeRange')] = initCfg['printPopAvgRates']
    
    initCfg[('analysis', 'plotTraces', 'oneFigPer')] = 'trace'

    initCfg['saveCellSecs'] = False
    initCfg['saveCellConns'] = False


    b = Batch(params=params, netParamsFile='netParams.py', cfgFile='cfg.py', initCfg=initCfg, groupedParams=groupedParams)

    return b


# ----------------------------------------------------------------------------------------------
# Evol
# ----------------------------------------------------------------------------------------------
def evolRates():
    # --------------------------------------------------------
    # parameters
    params = specs.ODict()

    # long-range inputs
    params[('weightLong', 'TPO')] = [0.25, 0.75] 
    params[('weightLong', 'TVL')] = [0.25, 0.75] 
    params[('weightLong', 'S1')] =  [0.25, 0.75] 
    params[('weightLong', 'S2')] =  [0.25, 0.75] 
    params[('weightLong', 'cM1')] = [0.25, 0.75] 
    params[('weightLong', 'M2')] =  [0.25, 0.75] 
    params[('weightLong', 'OC')] =  [0.25, 0.75]	

    # EEgain
    params['EEGain'] = [0.5, 1.5] 

    # IEgain
    ## L2/3+4
    params[('IEweights',0)] =  [0.5, 1.5]
    ## L5
    params[('IEweights',1)] = [0.5, 1.5] #[0.8, 1.0]   
    ## L6
    params[('IEweights',2)] =  [0.5, 1.5] # [0.8, 1.0]  

    # IIGain
    ## L2/3+4
    params[('IIweights',0)] =  [0.5, 1.5]
    ## L5
    params[('IIweights',1)] = [0.5, 1.5] #[0.8, 1.0]   
    ## L6
    params[('IIweights',2)] =  [0.5, 1.5] # [0.8, 1.0]  

    groupedParams = []

    # --------------------------------------------------------
    # initial config
    initCfg = {}
    initCfg['duration'] = 1500
    initCfg['printPopAvgRates'] = [500, 1500] 
    initCfg['dt'] = 0.025

    initCfg['scaleDensity'] = 1.0

    # cell params
    initCfg['ihGbar'] = 0.75  # ih (for quiet/sponti condition)
    initCfg['ihModel'] = 'migliore'  # ih model
    initCfg['ihGbarBasal'] = 1.0 # multiplicative factor for ih gbar in PT cells
    initCfg['ihlkc'] = 0.2 # ih leak param (used in Migliore)
    initCfg['ihLkcBasal'] = 1.0 # multiplicative factor for ih lk in PT cells
    initCfg['ihLkcBelowSoma'] = 0.01 # multiplicative factor for ih lk in PT cells
    initCfg['ihlke'] = -86  # ih leak param (used in Migliore)
    initCfg['ihSlope'] = 28  # ih leak param (used in Migliore)

    initCfg['somaNa'] = 5.0  # somatic Na conduct
    initCfg['dendNa'] = 0.3  # dendritic Na conduct (reduced to avoid dend spikes) 
    initCfg['axonNa'] = 7   # axon Na conduct (increased to compensate) 
    initCfg['axonRa'] = 0.005
    initCfg['gpas'] = 0.5
    initCfg['epas'] = 0.9
    
    # long-range input params
    initCfg['numCellsLong'] = 1000
    initCfg[('pulse', 'pop')] = 'None'
    initCfg[('pulse', 'start')] = 1000.0
    initCfg[('pulse', 'end')] = 1100.0
    initCfg[('pulse', 'noise')] = 0.8

    # conn params
    initCfg['IEdisynapticBias'] = None

    initCfg['weightNormThreshold'] = 4.0
    initCfg['IEGain'] = 1.0
    initCfg['IIGain'] = 1.0
    initCfg['IPTGain'] = 1.0

    # plotting and saving params
    initCfg[('analysis','plotRaster','timeRange')] = initCfg['printPopAvgRates']
    initCfg[('analysis','plotTraces','timeRange')] = initCfg['printPopAvgRates']

    initCfg['saveCellSecs'] = False
    initCfg['saveCellConns'] = False


    # --------------------------------------------------------
    # fitness function
    fitnessFuncArgs = {}
    pops = {}
    
    ## Exc pops
    Epops = ['IT2', 'IT4', 'IT5A', 'IT5B', 'PT5B', 'IT6', 'CT6']

    Etune = {'target': 5, 'width': 5, 'min': 0.5}
    for pop in Epops:
        pops[pop] = Etune
    
    ## Inh pops 
    Ipops = ['NGF1', 'PV2', 'SOM2', 'VIP2', 'NGF2',
            'PV4', 'SOM4', 'VIP4', 'NGF4',
            'PV5A', 'SOM5A','VIP5A','NGF5A',
            'PV5B', 'SOM5B','VIP5B','NGF5B',
            'PV6', 'SOM6', 'VIP6', 'NGF6']

    Itune = {'target': 10, 'width': 15, 'min': 0.25}
    for pop in Ipops:
        pops[pop] = Itune
    
    fitnessFuncArgs['pops'] = pops
    fitnessFuncArgs['maxFitness'] = 1000


    def fitnessFunc(simData, **kwargs):
        import numpy as np
        pops = kwargs['pops']
        maxFitness = kwargs['maxFitness']
        popFitness = [min(np.exp(abs(v['target'] - simData['popRates'][k])/v['width']), maxFitness) 
                if simData['popRates'][k] > v['min'] else maxFitness for k,v in pops.items()]
        fitness = np.mean(popFitness)

        popInfo = '; '.join(['%s rate=%.1f fit=%1.f'%(p, simData['popRates'][p], popFitness[i]) for i,p in enumerate(pops)])
        print('  '+popInfo)
        return fitness
    
    #from IPython import embed; embed()

    b = Batch(params=params, groupedParams=groupedParams, initCfg=initCfg)

    # Set evol alg configuration
    b.evolCfg = {
        'evolAlgorithm': 'custom',
        'fitnessFunc': fitnessFunc, # fitness expression (should read simData)
        'fitnessFuncArgs': fitnessFuncArgs,
        'pop_size': 100,
        'num_elites': 2,
        'mutation_rate': 0.5,
        'crossover': 0.5,
        'maximize': False, # maximize fitness function?
        'max_generations': 200,
        'time_sleep': 300, # 5min wait this time before checking again if sim is completed (for each generation)
        'maxiter_wait': 64, # (5h20) max number of times to check if sim is completed (for each generation)
        'defaultFitness': 1000, # set fitness value in case simulation time is over
        'scancelUser': 'ext_donald_doherty_actionpotenti'
    }


    return b

# ----------------------------------------------------------------------------------------------
# Adaptive Stochastic Descent (ASD)
# ----------------------------------------------------------------------------------------------
def optunaRates():

    # --------------------------------------------------------
    # parameters
    params = specs.ODict()

    # long-range inputs
    params[('weightLong', 'TPO')] = [0.25, 0.75]
    params[('weightLong', 'TVL')] = [0.25, 0.75]
    params[('weightLong', 'S1')] =  [0.25, 0.75]
    params[('weightLong', 'S2')] =  [0.25, 0.75]
    params[('weightLong', 'cM1')] = [0.25, 0.75]
    params[('weightLong', 'M2')] =  [0.25, 0.75]
    params[('weightLong', 'OC')] =  [0.25, 0.75]

    # EEgain
    params['EEGain'] = [0.5, 1.5]

    # IEgain
    ## L2/3+4
    params[('IEweights',0)] =  [0.5, 1.5]
    ## L5
    params[('IEweights',1)] = [0.5, 1.5] #[0.8, 1.0]   
    ## L6
    params[('IEweights',2)] =  [0.5, 1.5] # [0.8, 1.0]  

    # IIGain
    ## L2/3+4
    params[('IIweights',0)] =  [0.5, 1.5]
    ## L5
    params[('IIweights',1)] = [0.5, 1.5] #[0.8, 1.0]   
    ## L6
    params[('IIweights',2)] =  [0.5, 1.5] # [0.8, 1.0]  

    groupedParams = []

    # --------------------------------------------------------
    # initial config
    initCfg = {}
    initCfg['duration'] = 1500
    initCfg['printPopAvgRates'] = [[500, 750], [750, 1000], [1000, 1250], [1250, 1500]]
    initCfg['dt'] = 0.025

    initCfg['scaleDensity'] = 1.0

    # cell params
    initCfg['ihGbar'] = 0.75  # ih (for quiet/sponti condition)
    initCfg['ihModel'] = 'migliore'  # ih model
    initCfg['ihGbarBasal'] = 1.0 # multiplicative factor for ih gbar in PT cells
    initCfg['ihlkc'] = 0.2 # ih leak param (used in Migliore)
    initCfg['ihLkcBasal'] = 1.0 # multiplicative factor for ih lk in PT cells
    initCfg['ihLkcBelowSoma'] = 0.01 # multiplicative factor for ih lk in PT cells
    initCfg['ihlke'] = -86  # ih leak param (used in Migliore)
    initCfg['ihSlope'] = 28  # ih leak param (used in Migliore)

    initCfg['somaNa'] = 5.0  # somatic Na conduct
    initCfg['dendNa'] = 0.3  # dendritic Na conduct (reduced to avoid dend spikes) 
    initCfg['axonNa'] = 7   # axon Na conduct (increased to compensate) 
    initCfg['axonRa'] = 0.005
    initCfg['gpas'] = 0.5
    initCfg['epas'] = 0.9

    # long-range input params
    initCfg['numCellsLong'] = 1000
    initCfg[('pulse', 'pop')] = 'None'
    initCfg[('pulse', 'start')] = 1000.0
    initCfg[('pulse', 'end')] = 1100.0
    initCfg[('pulse', 'noise')] = 0.8

    # conn params
    initCfg['IEdisynapticBias'] = None

    initCfg['weightNormThreshold'] = 4.0
    initCfg['IEGain'] = 1.0
    initCfg['IIGain'] = 1.0
    initCfg['IPTGain'] = 1.0

    # plotting and saving params
    initCfg[('analysis','plotRaster','timeRange')] = [500,1500]
    initCfg[('analysis', 'plotTraces', 'timeRange')] = [500,1500]
    initCfg[('analysis', 'plotTraces', 'oneFigPer')] = 'trace'
    initCfg['recordLFP'] = None
    initCfg[('analysis', 'plotLFP')] = False

    initCfg['saveCellSecs'] = False
    initCfg['saveCellConns'] = False


    # --------------------------------------------------------
    # fitness function
    fitnessFuncArgs = {}
    pops = {}
    
    ## Exc pops
    Epops = ['IT2', 'IT4', 'IT5A', 'IT5B', 'PT5B', 'IT6', 'CT6']  # all layers

    Etune = {'target': 5, 'width': 20, 'min': 0.05}
    for pop in Epops:
        pops[pop] = Etune
    
    ## Inh pops 
    Ipops = ['NGF1',                            # L1
            'PV2', 'SOM2', 'VIP2', 'NGF2',      # L2/3
            'PV4', 'SOM4', 'VIP4', 'NGF4',      # L4
            'PV5A', 'SOM5A', 'VIP5A', 'NGF5A',  # L5A  
            'PV5B', 'SOM5B', 'VIP5B', 'NGF5B',  # L5B
            'PV6', 'SOM6', 'VIP6', 'NGF6']       # L6

    Itune = {'target': 10, 'width': 30, 'min': 0.05}
    for pop in Ipops:
        pops[pop] = Itune
    
    fitnessFuncArgs['pops'] = pops
    fitnessFuncArgs['maxFitness'] = 1000
    fitnessFuncArgs['tranges'] = initCfg['printPopAvgRates']


    def fitnessFunc(simData, **kwargs):
        import numpy as np
        pops = kwargs['pops']
        maxFitness = kwargs['maxFitness']
        tranges = kwargs['tranges']

        popFitnessAll = []

        for trange in tranges:
            popFitnessAll.append([min(np.exp(abs(v['target'] - simData['popRates'][k]['%d_%d'%(trange[0], trange[1])])/v['width']), maxFitness) 
                if simData['popRates'][k]['%d_%d'%(trange[0], trange[1])] > v['min'] else maxFitness for k, v in pops.items()])
        
        popFitness = np.mean(np.array(popFitnessAll), axis=0)
        
        fitness = np.mean(popFitness)

        popInfo = '; '.join(['%s rate=%.1f fit=%1.f' % (p, np.mean(list(simData['popRates'][p].values())), popFitness[i]) for i,p in enumerate(pops)])
        print('  ' + popInfo)

        return fitness
        
    # create Batch object with paramaters to modify, and specifying files to use
    b = Batch(params=params, netParamsFile='netParams.py', cfgFile='cfg.py', initCfg=initCfg, groupedParams=groupedParams)

    b.method = 'optuna'

    b.optimCfg = {
        'fitnessFunc': fitnessFunc, # fitness expression (should read simData)
        'fitnessFuncArgs': fitnessFuncArgs,
        'maxFitness': fitnessFuncArgs['maxFitness'],
        'maxiters':     1e6,    #    Maximum number of iterations (1 iteration = 1 function evaluation)
        'maxtime':      None,    #    Maximum time allowed, in seconds
        'maxiter_wait': 16,     # check every time step (time_sleep), 16 times
        'time_sleep': 900,     # 16 x 900 seconds = 4 hours
        'popsize': 1  # unused - run with mpi 
    }

    return b

# ----------------------------------------------------------------------------------------------
# Optuna optimization
# ----------------------------------------------------------------------------------------------
def optunaRatesCellTypes():
    # --------------------------------------------------------
    # parameters
    params = specs.ODict()

    # from prev
    import json
    with open('../data/v103_batch3/trial_718/trial_718_cfg.json', 'rb') as f:
        cfgLoad = json.load(f)['simConfig']

    params[('EICellTypeGain', 'PV')] = [0.1, 4.0]
    params[('EICellTypeGain', 'SOM')] = [0.1, 4.0]
    params[('EICellTypeGain', 'VIP')] = [0.1, 4.0]
    params[('EICellTypeGain', 'NGF')] = [0.1, 4.0]

    rangeG = 0.25
    minG = 0.5
    maxG = 1.5

    # EEgain
    params['EEGain'] = [max(cfgLoad['EEGain']-rangeG, minG), min(cfgLoad['EEGain']+rangeG, maxG)]

    # IEgain
    ## L2/3+4
    params[('IEweights',0)] =  [max(cfgLoad['IEweights'][0]-rangeG, minG), min(cfgLoad['IEweights'][0]+rangeG, maxG)]
    ## L5
    params[('IEweights',1)] = [max(cfgLoad['IEweights'][1]-rangeG, minG), min(cfgLoad['IEweights'][1]+rangeG, maxG)] #[0.8, 1.0]   
    ## L6
    params[('IEweights',2)] =  [max(cfgLoad['IEweights'][2]-rangeG, minG), min(cfgLoad['IEweights'][2]+rangeG, maxG)] # [0.8, 1.0]  

    # IIGain
    ## L2/3+4
    params[('IIweights',0)] =  [max(cfgLoad['IIweights'][0]-rangeG, minG), min(cfgLoad['IIweights'][0]+rangeG, maxG)]
    ## L5
    params[('IIweights',1)] = [max(cfgLoad['IIweights'][1]-rangeG, minG), min(cfgLoad['IIweights'][1]+rangeG, maxG)] #[0.8, 1.0]   
    ## L6
    params[('IIweights',2)] =  [max(cfgLoad['IIweights'][2]-rangeG, minG), min(cfgLoad['IIweights'][2]+rangeG, maxG)] # [0.8, 1.0]  

    groupedParams = []

    # --------------------------------------------------------
    # initial config
    initCfg = {}
    initCfg['duration'] = 1500
    initCfg['printPopAvgRates'] = [[500, 750], [750, 1000], [1000, 1250], [1250, 1500]]
    initCfg['dt'] = 0.025

    initCfg['scaleDensity'] = 1.0

    # cell params
    initCfg['ihGbar'] = 0.75  # ih (for quiet/sponti condition)
    initCfg['ihModel'] = 'migliore'  # ih model
    initCfg['ihGbarBasal'] = 1.0 # multiplicative factor for ih gbar in PT cells
    initCfg['ihlkc'] = 0.2 # ih leak param (used in Migliore)
    initCfg['ihLkcBasal'] = 1.0 # multiplicative factor for ih lk in PT cells
    initCfg['ihLkcBelowSoma'] = 0.01 # multiplicative factor for ih lk in PT cells
    initCfg['ihlke'] = -86  # ih leak param (used in Migliore)
    initCfg['ihSlope'] = 28  # ih leak param (used in Migliore)

    initCfg['somaNa'] = 5.0  # somatic Na conduct
    initCfg['dendNa'] = 0.3  # dendritic Na conduct (reduced to avoid dend spikes) 
    initCfg['axonNa'] = 7   # axon Na conduct (increased to compensate) 
    initCfg['axonRa'] = 0.005
    initCfg['gpas'] = 0.5
    initCfg['epas'] = 0.9

    # long-range input params
    initCfg['numCellsLong'] = 1000
    initCfg[('pulse', 'pop')] = 'None'
    initCfg[('pulse', 'start')] = 1000.0
    initCfg[('pulse', 'end')] = 1100.0
    initCfg[('pulse', 'noise')] = 0.8

    # long-range inputs
    params[('weightLong', 'TPO')] = [cfgLoad['weightLong']['TPO']]
    params[('weightLong', 'TVL')] = [cfgLoad['weightLong']['TVL']]
    params[('weightLong', 'S1')] =  [cfgLoad['weightLong']['S1']]
    params[('weightLong', 'S2')] =  [cfgLoad['weightLong']['S2']]
    params[('weightLong', 'cM1')] = [cfgLoad['weightLong']['cM1']]
    params[('weightLong', 'M2')] =  [cfgLoad['weightLong']['M2']]
    params[('weightLong', 'OC')] =  [cfgLoad['weightLong']['OC']]

    # conn params
    initCfg['IEdisynapticBias'] = None

    initCfg['weightNormThreshold'] = 4.0
    initCfg['IEGain'] = 1.0
    initCfg['IIGain'] = 1.0
    initCfg['IPTGain'] = 1.0

    # plotting and saving params
    initCfg[('analysis','plotRaster','timeRange')] = [500,1500]
    initCfg[('analysis', 'plotTraces', 'timeRange')] = [500,1500]
    initCfg[('analysis', 'plotTraces', 'oneFigPer')] = 'trace'
    initCfg['recordLFP'] = None
    initCfg[('analysis', 'plotLFP')] = False

    initCfg['saveCellSecs'] = False
    initCfg['saveCellConns'] = False

    # --------------------------------------------------------
    # fitness function
    fitnessFuncArgs = {}
    pops = {}

    ## Exc pops
    Epops = ['IT2', 'IT4', 'IT5A', 'IT5B', 'PT5B', 'IT6', 'CT6']  # all layers

    Etune = {'target': 5, 'width': 20, 'min': 0.05}
    for pop in Epops:
        pops[pop] = Etune

    ## Inh pops 
    Ipops = ['NGF1',                            # L1
            'PV2', 'SOM2', 'VIP2', 'NGF2',      # L2/3
            'PV4', 'SOM4', 'VIP4', 'NGF4',      # L4
            'PV5A', 'SOM5A', 'VIP5A', 'NGF5A',  # L5A  
            'PV5B', 'SOM5B', 'VIP5B', 'NGF5B',  # L5B
            'PV6', 'SOM6', 'VIP6', 'NGF6']       # L6

    Itune = {'target': 10, 'width': 30, 'min': 0.05}
    for pop in Ipops:
        pops[pop] = Itune

    fitnessFuncArgs['pops'] = pops
    fitnessFuncArgs['maxFitness'] = 1000
    fitnessFuncArgs['tranges'] = initCfg['printPopAvgRates']


    def fitnessFunc(simData, **kwargs):
        import numpy as np
        pops = kwargs['pops']
        maxFitness = kwargs['maxFitness']
        tranges = kwargs['tranges']

        popFitnessAll = []

        for trange in tranges:
            popFitnessAll.append([min(np.exp(abs(v['target'] - simData['popRates'][k]['%d_%d'%(trange[0], trange[1])])/v['width']), maxFitness)
                if simData['popRates'][k]['%d_%d'%(trange[0], trange[1])] > v['min'] else maxFitness for k, v in pops.items()])

        popFitness = np.mean(np.array(popFitnessAll), axis=0)

        fitness = np.mean(popFitness)

        popInfo = '; '.join(['%s rate=%.1f fit=%1.f' % (p, np.mean(list(simData['popRates'][p].values())), popFitness[i]) for i,p in enumerate(pops)])
        print('  ' + popInfo)

        return fitness

    # create Batch object with paramaters to modify, and specifying files to use
    b = Batch(params=params, netParamsFile='netParams.py', cfgFile='cfg.py', initCfg=initCfg, groupedParams=groupedParams)

    b.method = 'optuna'

    b.optimCfg = {
        'fitnessFunc': fitnessFunc, # fitness expression (should read simData)
        'fitnessFuncArgs': fitnessFuncArgs,
        'maxFitness': fitnessFuncArgs['maxFitness'],
        'maxiters':     1e6,    #    Maximum number of iterations (1 iteration = 1 function evaluation)
        'maxtime':      None,    #    Maximum time allowed, in seconds
        'maxiter_wait': 16,     # check every time step (time_sleep), 16 times
        'time_sleep': 900,     # 16 x 900 seconds = 4 hours
        'popsize': 1  # unused - run with mpi 
    }

    return b

# ----------------------------------------------------------------------------------------------
# Run configurations
# ----------------------------------------------------------------------------------------------
def setRunCfg(b, type='mpi_bulletin'):
    if type=='mpi_bulletin' or type=='mpi':
        b.runCfg = {'type': 'mpi_bulletin', 
            'script': 'init.py', 
            'skip': True}

    elif type=='mpi_direct':
        b.runCfg = {'type': 'mpi_direct',
            'cores': 4,
            'script': 'init_cell.py',
            'mpiCommand': 'mpirun',
            'skip': True}

    elif type=='hpc_torque':
        b.runCfg = {'type': 'hpc_torque',
             'script': 'init.py',
             'nodes': 3,
             'ppn': 8,
             'walltime': "12:00:00",
             'queueName': 'longerq',
             'sleepInterval': 5,
             'skip': True}

    elif type=='hpc_slurm_comet':
        b.runCfg = {'type': 'hpc_slurm', 
            'allocation': 'shs100', # bridges='ib4iflp', comet m1='shs100', comet nsg='csd403'
            #'reservation': 'salva1',
            'walltime': '6:00:00',
            'nodes': 4,
            'coresPerNode': 24,  # comet=24, bridges=28
            'email': 'donald.doherty@actionpotential.com',
            'folder': '/home/salvadord/m1/sim/',  # comet='/salvadord', bridges='/salvi82'
            'script': 'init.py', 
            'mpiCommand': 'ibrun', # comet='ibrun', bridges='mpirun'
            'skipCustom': '_raster.png'}

    elif type=='hpc_slurm_gcp':
        b.runCfg = {'type': 'hpc_slurm', 
            'allocation': 'default', # bridges='ib4iflp', comet m1='shs100', comet nsg='csd403', gcp='default'
            'walltime': '24:00:00', #'48:00:00',
            'nodes': 1,
            'coresPerNode': 80,  # comet=24, bridges=28, gcp=32
            'email': 'donald.doherty@actionpotential.com',
            'folder': '/home/ext_donald_doherty_actionpotenti/sim/',  # comet,gcp='/salvadord', bridges='/salvi82'
            'script': 'init.py', 
            'mpiCommand': 'mpirun', # comet='ibrun', bridges,gcp='mpirun' 
            'skipCustom': '_raster.png'}
            #'custom': '#SBATCH --exclude=compute[17-64000]'} # only use first 16 nodes (non-preemptible for long runs )
            # --nodelist=compute1


    elif type=='hpc_slurm_bridges':
        b.runCfg = {'type': 'hpc_slurm', 
            'allocation': 'ib4iflp', # bridges='ib4iflp', comet m1='shs100', comet nsg='csd403'
            'walltime': '06:00:00',
            'nodes': 2,
            'coresPerNode': 28,  # comet=24, bridges=28
            'email': 'donald.doherty@actionpotential.com',
            'folder': '/home/salvi82/m1/sim/',  # comet='/salvadord', bridges='/salvi82'
            'script': 'init.py', 
            'mpiCommand': 'mpirun', # comet='ibrun', bridges='mpirun'
            'skip': True}

# ----------------------------------------------------------------------------------------------
# Main code
# ----------------------------------------------------------------------------------------------

# if __name__ == '__main__': 
#     #b = custom() 
#     b = optunaRatesCellTypes()
#     b.batchLabel = 'v104_batch3'  
#     b.saveFolder = '../data/'+b.batchLabel
#     #b.method = 'grid'
#     setRunCfg(b, 'hpc_slurm_gcp')
#     b.run() # run batch

if __name__ == '__main__': 
    #b = custom() 
    b = projectionWeights()

    b.dataFolder    = '../data/batch_sims'
    simDate         = '2021_12_29'
    b.date          = simDate

    b.batchLabel    = simDate + '_fullModel_Th_M1_density002_batch00'
    b.saveFolder    = b.dataFolder + '/' +  b.batchLabel


    b.saveFolder = '../data/batch_sims_joao/'+b.batchLabel
    #b.method = 'grid'
    setRunCfg(b, 'mpi_bulletin')
    # setRunCfg(b, 'hpc_slurm_gcp')
    b.run() # run batch

"""  """

# ----------------------------------------------------------------------------------------------
# Git Commands
# ----------------------------------------------------------------------------------------------
'''

MOST USEFUL COMMANDS:

mpiexec -n 8 nrniv -python -mpi init.py 
mpiexec -n 80 nrniv -python -mpi init.py 

sudo git commit -a -m " "; sudo git push; sudo git pull


mpiexec -n 8 nrniv -python -mpi batch.py
mpiexec -n 80 nrniv -python -mpi batch.py


# ---       Git update                --- #

# check the status of the repository - new/modified/unstaged files
sudo git status

# LOCAL MACHINE
# update repository - running sim 
sudo git commit -a -m " "; sudo git push; sudo git pull

# update repository - custom modification
sudo git commit -a -m "insert modification description"; sudo git push; sudo git pull

# CLOUD MACHINE
sudo git pull

# FORCE PULL
git reset --hard HEAD
git pull


gsutil -m cp -r ~/Research/Models/NetPyNE/thalamus/data/new_sim gs://joao-bucket

gsutil -m cp -r /home/joao_moreira_downstate_edu/Research/Models/NetPyNE/thalamus/data/thalamus_allpops/2021_04_30/thalamus_allpops_thalamus_allpops_t_allpops_012 gs://joao-bucket/2021_04_30
gsutil -m cp -r /home/joao/Research/Models/NetPyNE/thalamus/data/thalamus_allpops/2021_04_30/thalamus_allpops_thalamus_allpops_t_allpops_012 gs://joao-bucket/2021_04_30

gsutil -m cp -r gs://joao-bucket/2021_04_30/thalamus_allpops_thalamus_allpops_t_allpops_012 /Users/joao/Research/Models/NetPyNE/thalamus/data/batch_sims/gcloud/2021_04_30




Permanently authenticating with Git repositories
    Run the following command to enable credential caching:

        $ git config credential.helper store
        $ git push https://github.com/owner/repo.git

        Username for 'https://github.com': <USERNAME>
        Password for 'https://USERNAME@github.com': <PASSWORD>
    
    You should also specify caching expire,

        git config --global credential.helper 'cache --timeout 7200'

    After enabling credential caching, it will be cached for 7200 seconds (2 hour).


'''
# ----------------------------------------------------------------------------------------------
# Bash Commands
# ----------------------------------------------------------------------------------------------
'''
# ---       Run Sim commands                --- #

# run single (init) simulation
mpiexec -n 8 nrniv -python -mpi init.py 
mpiexec -n 80 nrniv -python -mpi init.py 

screen -L mpiexec -n 8 nrniv -python -mpi init.py 
screen -L mpiexec -n 80 nrniv -python -mpi init.py 

# run batch simulation
mpiexec -n 6 nrniv -python -mpi batch.py

# ---       Run Sim commands (background)   --- #
# run batch simulation in the background
screen -L mpiexec -n 6 nrniv -python -mpi batch.py

#GCloud
screen -L mpiexec -n 6 nrniv -python -mpi batch.py 
# OBS: GCloud PC "joao-machine" has up to 80 cores

# detach
ctrl+A , D

# retrun to the process
screen -r

# list processes
screen -ls

# kill a process
screen -r <id of the process gotten from "screen -ls">

# Kill all screens
pkill screen

'''
# ----------------------------------------------------------------------------------------------
# Cluster Commands
# ----------------------------------------------------------------------------------------------
'''
# Queue alias
alias sq=‘squeue -u $USER; squeue -u $USER | wc -l’

# Check queue
sq

# Cancel Jobs
scancel -u $USER

# Inspect jobs
sdetail [job id]
'''