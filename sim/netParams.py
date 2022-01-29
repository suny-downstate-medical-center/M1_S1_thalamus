
"""
netParams.py 

High-level specifications for M1 + Th + S1 network model using NetPyNE

Contributors: salvadordura@gmail.com, fernandodasilvaborges@gmail.com, joaovvitor@gmail.com
"""

from netpyne import specs
import pickle, json
import os
import numpy as np
import pandas as pd


netParams = specs.NetParams()   # object of class NetParams to store the network parameters

netParams.version = 1

try:
    from __main__ import cfg  # import SimConfig object with params from parent module
except:
    from cfg import cfg

#------------------------------------------------------------------------------
#
# NETWORK PARAMETERS
#
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# General network parameters
#------------------------------------------------------------------------------
netParams.scale = cfg.scale # Scale factor for number of cells
netParams.sizeX = cfg.sizeX # x-dimension (horizontal length) size in um
netParams.sizeY = cfg.sizeY # y-dimension (vertical height or cortical depth) size in um
netParams.sizeZ = cfg.sizeZ # z-dimension (horizontal depth) size in um
netParams.shape = 'cylinder' # cylindrical (column-like) volume

#------------------------------------------------------------------------------
# General connectivity parameters
#------------------------------------------------------------------------------
netParams.scaleConnWeight = 1.0 # Connection weight scale factor (default if no model specified)
netParams.scaleConnWeightModels = {'HH_simple': 1.0, 'HH_reduced': 1.0, 'HH_full': 1.0} #scale conn weight factor for each cell model
netParams.scaleConnWeightNetStims = 1.0 # scale conn weight factor for NetStims
netParams.defaultThreshold = 0.0 # spike threshold, 10 mV is NetCon default, lower it for all cells
netParams.defaultDelay = 2.0 # default conn delay (ms)
netParams.propVelocity = 500.0 # propagation velocity (um/ms)
netParams.probLambda = 100.0  # length constant (lambda) for connection probability decay (um)
netParams.defineCellShapes = True  # convert stylized geoms to 3d points

netParams.defaultDelayS1 = 0.1 # default conn delay (ms)
netParams.propVelocityS1 = 300.0 #  300 μm/ms (Stuart et al., 1997)

#------------------------------------------------------------------------------
# Cell parameters
#------------------------------------------------------------------------------
cellModels  = ['HH_reduced', 'HH_full']
excTypes    = ['IT', 'CT', 'PT']
Etypes      = ['E_VL']
inhTypes    = ['PV', 'SOM', 'VIP', 'NGF']
Itypes      = ['I_VL','I_RTN']

ymin={ 
		'mt_RTN':       1688,
		'ss_RTN_o':     1688,
		'ss_RTN_m':     1766,
		'ss_RTN_i':     1844,
		'VPL_sTC':      2000,
		'VPM_sTC':      2156,
		'POm_sTC_m1':   2312,
		'POm_sTC_s1':   2312,
		'VL_sTC':       2000,
		'VM_sTC_m1':    2000
	}

ymax={ 
		'mt_RTN':       2000,
		'ss_RTN_o':     1766,
		'ss_RTN_m':     1844,
		'ss_RTN_i':     2000,
		'VPL_sTC':      2156,
		'VPM_sTC':      2312,
		'POm_sTC_m1':   2624,
		'POm_sTC_s1':   2624,
		'VL_sTC':       2468,
		'VM_sTC_m1':    2468
	}

t_layer = {
        'mt_RTN':       [ymin['mt_RTN'],        ymax['mt_RTN']],
		'ss_RTN_o':     [ymin['ss_RTN_o'],      ymax['ss_RTN_o']],
		'ss_RTN_m':     [ymin['ss_RTN_m'],      ymax['ss_RTN_m']],
		'ss_RTN_i':     [ymin['ss_RTN_i'],      ymax['ss_RTN_i']],
		'VPL_sTC':      [ymin['VPL_sTC'],       ymax['VPL_sTC']],
		'VPM_sTC':      [ymin['VPM_sTC'],       ymax['VPM_sTC']],
		'POm_sTC_m1':   [ymin['POm_sTC_m1'],    ymax['POm_sTC_m1']],
		'POm_sTC_s1':   [ymin['POm_sTC_s1'],    ymax['POm_sTC_s1']],
		'VL_sTC':       [ymin['VL_sTC'],        ymax['VL_sTC']],
		'VM_sTC_m1':    [ymin['VM_sTC_m1'],     ymax['VM_sTC_m1']],
	}  # normalized layer boundaries 

layer = {   '1':                [0.0,   0.1], 
            '2':                [0.1,   0.29], 
            '4':                [0.29,  0.37], 
            '5A':               [0.37,  0.47], 
            '24':               [0.1,   0.37], 
            '5B':               [0.47,  0.8], 
            '6':                [0.8,   1.0], 
            'longS1':           [2.2,   2.3], 
            'longS2':           [2.3,   2.4], 
            'longcM1':          [2.4,   2.5], 
            'longM2':           [2.5,   2.6], 
            'longOC':           [2.6,   2.7]
            }  # normalized layer boundaries


netParams.correctBorder = {'threshold': [cfg.correctBorderThreshold, cfg.correctBorderThreshold, cfg.correctBorderThreshold], 
                        'yborders': [layer['2'][0], layer['5A'][0], layer['6'][0], layer['6'][1]]}  # correct conn border effect

#------------------------------------------------------------------------------
# Thalamocortical Relay cell model 
netParams.loadCellParamsRule(label='sTC_cell', fileName=cfg.sTC_model)
netParams.cellParams['sTC_cell']['conds']={} 

# Reticular Nucleus cell model
netParams.loadCellParamsRule(label='sRE_cell', fileName=cfg.sRE_model)
netParams.cellParams['sRE_cell']['conds']={}

#------------------------------------------------------------------------------
## M1 rules 
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

cellParamLabels = ['IT2_reduced', 'IT4_reduced', 'IT5A_reduced', 'IT5B_reduced', 'PT5B_reduced',
    'IT6_reduced', 'CT6_reduced', 'SOM_reduced', 'IT5A_full']#  'PV_reduced', 'VIP_reduced', 'NGF_reduced','PT5B_full'] #  # list of cell rules to load from file
loadCellParams  = cellParamLabels
saveCellParams  = False #True

for ruleLabel in loadCellParams:
    netParams.loadCellParamsRule(label=ruleLabel, fileName='../cells/' + ruleLabel + '_cellParams.pkl')

#------------------------------------------------------------------------------
# Reduced cell model params (6-comp) 
reducedCells    = { # layer and cell type for reduced cell models
                    'IT2_reduced':  {'layer': '2',  'cname': 'CSTR6', 'carg': 'BS1578'}, 
                    'IT4_reduced':  {'layer': '4',  'cname': 'CSTR6', 'carg': 'BS1578'},
                    'IT5A_reduced': {'layer': '5A', 'cname': 'CSTR6', 'carg': 'BS1579'},
                    'IT5B_reduced': {'layer': '5B', 'cname': 'CSTR6', 'carg': 'BS1579'},
                    'PT5B_reduced': {'layer': '5B', 'cname': 'SPI6',  'carg':  None},
                    'IT6_reduced':  {'layer': '6',  'cname': 'CSTR6', 'carg': 'BS1579'},
                    'CT6_reduced':  {'layer': '6',  'cname': 'CSTR6', 'carg': 'BS1578'}}

reducedSecList  = { # section Lists for reduced cell model
                    'alldend':  ['Adend1', 'Adend2', 'Adend3', 'Bdend'],
                    'spiny':    ['Adend1', 'Adend2', 'Adend3', 'Bdend'],
                    'apicdend': ['Adend1', 'Adend2', 'Adend3'],
                    'perisom':  ['soma']}

for label, p in reducedCells.items():  # create cell rules that were not loaded 
    if label not in loadCellParams:
        cellRule = netParams.importCellParams(label=label, conds={'cellType': label[0:2], 'cellModel': 'HH_reduced', 'ynorm': layer[p['layer']]},
        fileName='../cells/'+p['cname']+'.py', cellName=p['cname'], cellArgs={'params': p['carg']} if p['carg'] else None)
        dendL = (layer[p['layer']][0]+(layer[p['layer']][1]-layer[p['layer']][0])/2.0) * cfg.sizeY  # adapt dend L based on layer
        for secName in ['Adend1', 'Adend2', 'Adend3', 'Bdend']: cellRule['secs'][secName]['geom']['L'] = dendL / 3.0  # update dend L
        for k,v in reducedSecList.items(): cellRule['secLists'][k] = v  # add secLists
        netParams.addCellParamsWeightNorm(label, '../conn/'+label+'_weightNorm.pkl', threshold=cfg.weightNormThreshold)  # add weightNorm

        # set 3d points
        offset, prevL = 0, 0
        somaL = netParams.cellParams[label]['secs']['soma']['geom']['L']
        for secName, sec in netParams.cellParams[label]['secs'].items():
            sec['geom']['pt3d'] = []
            if secName in ['soma', 'Adend1', 'Adend2', 'Adend3']:  # set 3d geom of soma and Adends
                sec['geom']['pt3d'].append([offset+0, prevL, 0, sec['geom']['diam']])
                prevL = float(prevL + sec['geom']['L'])
                sec['geom']['pt3d'].append([offset+0, prevL, 0, sec['geom']['diam']])
            if secName in ['Bdend']:  # set 3d geom of Bdend
                sec['geom']['pt3d'].append([offset+0, somaL, 0, sec['geom']['diam']])
                sec['geom']['pt3d'].append([offset+sec['geom']['L'], somaL, 0, sec['geom']['diam']])        
            if secName in ['axon']:  # set 3d geom of axon
                sec['geom']['pt3d'].append([offset+0, 0, 0, sec['geom']['diam']])
                sec['geom']['pt3d'].append([offset+0, -sec['geom']['L'], 0, sec['geom']['diam']])   

        if saveCellParams: netParams.saveCellParamsRule(label=label, fileName='../cells/'+label+'_cellParams.pkl')

#------------------------------------------------------------------------------
## PT5B full cell model params (700+ comps)
if 'PT5B_full' not in loadCellParams:
    ihMod2str = {'harnett': 1, 'kole': 2, 'migliore': 3}
    cellRule = netParams.importCellParams(label='PT5B_full', conds={'cellType': 'PT', 'cellModel': 'HH_full'},
      fileName='../cells/PTcell.hoc', cellName='PTcell', cellArgs=[ihMod2str[cfg.ihModel], cfg.ihSlope], somaAtOrigin=True)
    nonSpiny = ['apic_0', 'apic_1']
    netParams.addCellParamsSecList(label='PT5B_full', secListName='perisom', somaDist=[0, 50])  # sections within 50 um of soma
    netParams.addCellParamsSecList(label='PT5B_full', secListName='below_soma', somaDistY=[-600, 0])  # sections within 0-300 um of soma
    for sec in nonSpiny: cellRule['secLists']['perisom'].remove(sec)
    cellRule['secLists']['alldend'] = [sec for sec in cellRule.secs if ('dend' in sec or 'apic' in sec)] # basal+apical
    cellRule['secLists']['apicdend'] = [sec for sec in cellRule.secs if ('apic' in sec)] # apical
    cellRule['secLists']['spiny'] = [sec for sec in cellRule['secLists']['alldend'] if sec not in nonSpiny]
    # Adapt ih params based on cfg param
    for secName in cellRule['secs']:
        for mechName,mech in cellRule['secs'][secName]['mechs'].items():
            if mechName in ['ih','h','h15', 'hd']: 
                mech['gbar'] = [g*cfg.ihGbar for g in mech['gbar']] if isinstance(mech['gbar'],list) else mech['gbar']*cfg.ihGbar
                if cfg.ihModel == 'migliore':   
                    mech['clk'] = cfg.ihlkc  # migliore's shunt current factor
                    mech['elk'] = cfg.ihlke  # migliore's _S1shunt current reversal potential
                if secName.startswith('dend'): 
                    mech['gbar'] *= cfg.ihGbarBasal  # modify ih conductance in soma+basal dendrites
                    mech['clk'] *= cfg.ihlkcBasal  # modify ih conductance in soma+basal dendrites
                if secName in cellRule['secLists']['below_soma']: #secName.startswith('dend'): 
                    mech['clk'] *= cfg.ihlkcBelowSoma  # modify ih conductance in soma+basal dendrites
    # Reduce dend Na to avoid dend spikes (compensate properties by modifying axon params)
    for secName in cellRule['secLists']['alldend']:
        cellRule['secs'][secName]['mechs']['nax']['gbar'] = 0.0153130368342 * cfg.dendNa # 0.25 
    cellRule['secs']['soma']['mechs']['nax']['gbar'] = 0.0153130368342  * cfg.somaNa
    cellRule['secs']['axon']['mechs']['nax']['gbar'] = 0.0153130368342  * cfg.axonNa # 11  
    cellRule['secs']['axon']['geom']['Ra'] = 137.494564931 * cfg.axonRa # 0.005
    # Remove Na (TTX)
    if cfg.removeNa:
        for secName in cellRule['secs']: cellRule['secs'][secName]['mechs']['nax']['gbar'] = 0.0
    netParams.addCellParamsWeightNorm('PT5B_full', '../conn/PT5B_full_weightNorm.pkl', threshold=cfg.weightNormThreshold)  # load weight norm
    if saveCellParams: netParams.saveCellParamsRule(label='PT5B_full', fileName='../cells/PT5B_full_cellParams.pkl')

#------------------------------------------------------------------------------
## IT5A full cell model params (700+ comps)
if 'IT5A_full' not in loadCellParams:
    cellRule = netParams.importCellParams(label='IT5A_full', conds={'cellType': 'IT', 'cellModel': 'HH_full', 'ynorm': layer['5A']},
      fileName='./cells/ITcell.py', cellName='ITcell', cellArgs={'params': 'BS1579'}, somaAtOrigin=True)
    netParams.renameCellParamsSec(label='IT5A_full', oldSec='soma_0', newSec='soma')
    netParams.addCellParamsWeightNorm('IT5A_full', '../conn/IT_full_BS1579_weightNorm.pkl', threshold=cfg.weightNormThreshold) # add weightNorm before renaming soma_0
    netParams.addCellParamsSecList(label='IT5A_full', secListName='perisom', somaDist=[0, 50])  # sections within 50 um of soma
    cellRule['secLists']['alldend'] = [sec for sec in cellRule.secs if ('dend' in sec or 'apic' in sec)] # basal+apical
    cellRule['secLists']['apicdend'] = [sec for sec in cellRule.secs if ('apic' in sec)] # basal+apical
    cellRule['secLists']['spiny'] = [sec for sec in cellRule['secLists']['alldend'] if sec not in ['apic_0', 'apic_1']]
    if saveCellParams: netParams.saveCellParamsRule(label='IT5A_full', fileName='../cells/IT5A_full_cellParams.pkl')

#------------------------------------------------------------------------------
## PV cell params (3-comp)
if 'PV_reduced' not in loadCellParams:
    cellRule = netParams.importCellParams(label='PV_reduced', conds={'cellType':'PV', 'cellModel':'HH_reduced'}, 
      fileName='../cells/FS3.hoc', cellName='FScell1', cellInstance = True)
    cellRule['secLists']['spiny'] = ['soma', 'dend']
    netParams.addCellParamsWeightNorm('PV_reduced', '../conn/PV_reduced_weightNorm.pkl', threshold=cfg.weightNormThreshold)
    # cellRule['secs']['soma']['weightNorm'][0] *= 1.5
    if saveCellParams: netParams.saveCellParamsRule(label='PV_reduced', fileName='../cells/PV_reduced_cellParams.pkl')

#------------------------------------------------------------------------------
## SOM cell params (3-comp)
if 'SOM_reduced' not in loadCellParams:
    cellRule = netParams.importCellParams(label='SOM_reduced', conds={'cellType':'SOM', 'cellModel':'HH_reduced'}, 
      fileName='../cells/LTS3.hoc', cellName='LTScell1', cellInstance = True)
    cellRule['secLists']['spiny'] = ['soma', 'dend']
    netParams.addCellParamsWeightNorm('SOM_reduced', '../conn/SOM_reduced_weightNorm.pkl', threshold=cfg.weightNormThreshold)
    if saveCellParams: netParams.saveCellParamsRule(label='SOM_reduced', fileName='../cells/SOM_reduced_cellParams.pkl')
    
#------------------------------------------------------------------------------
## VIP cell params (5-comp)
if 'VIP_reduced' not in loadCellParams:
    cellRule = netParams.importCellParams(label='VIP_reduced', conds={'cellType': 'VIP', 'cellModel': 'HH_reduced'}, fileName='../cells/vipcr_cell.hoc',         cellName='VIPCRCell_EDITED', importSynMechs = True)
    cellRule['secLists']['spiny'] = ['soma', 'rad1', 'rad2', 'ori1', 'ori2']
    netParams.addCellParamsWeightNorm('VIP_reduced', '../conn/VIP_reduced_weightNorm.pkl', threshold=cfg.weightNormThreshold)
    if saveCellParams: netParams.saveCellParamsRule(label='VIP_reduced', fileName='../cells/VIP_reduced_cellParams.pkl')

#------------------------------------------------------------------------------
## NGF cell params (1-comp)
if 'NGF_reduced' not in loadCellParams:
    cellRule = netParams.importCellParams(label='NGF_reduced', conds={'cellType': 'NGF', 'cellModel': 'HH_reduced'}, fileName='../cells/ngf_cell.hoc', cellName='ngfcell', importSynMechs = True)
    cellRule['secLists']['spiny'] = ['soma', 'dend']
    netParams.addCellParamsWeightNorm('NGF_reduced', '../conn/NGF_reduced_weightNorm.pkl', threshold=cfg.weightNormThreshold)
    cellRule['secs']['soma']['weightNorm'][0] *= 1.5
    cellRule['secs']['soma']['weightNorm'][0] *= 1.5
    if saveCellParams: netParams.saveCellParamsRule(label='NGF_reduced', fileName='../cells/NGF_reduced_cellParams.pkl')

#------------------------------------------------------------------------------
# Population parameters
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
## load densities
with open('../cells/cellDensity.pkl', 'rb') as fileObj:         density = pickle.load(fileObj)['density']
density         = {k: [x * cfg.scaleDensity for x in v] for k,v in density.items()} # Scale densities 

with open('../cells/thal_cellDensity.pkl', 'rb') as fileObj:    thal_density = pickle.load(fileObj)['density']
thal_density    = {k: [x * cfg.scaleDensity for x in v] for k,v in thal_density.items()} # Scale densities 

# # BLUE BRAIN DATA - cell density - BBP(https://bbp.epfl.ch/nexus/cell-atlas/)
# VL  = 67711.5
# VM  = 68816
# VPL = 60916.1
# VPM = 76151.1
# POm = 62144.8
# RTN = 69417.7

zmin=0
zmax=300

xmin={
    'VPL_sTC':      300,
	'VPM_sTC':      300,
	'POm_sTC_m1':   300,
	'POm_sTC_s1':   300,
	'VL_sTC':       150,
	'VM_sTC_m1':    0,  
	'mt_RTN':       0, 
	'ss_RTN_o':     300,
	'ss_RTN_m':     300,
	'ss_RTN_i':     300,
	}

xmax={
	'VPL_sTC':      600,
	'VPM_sTC':      600,
	'POm_sTC_m1':   600,
	'POm_sTC_s1':   600,
	'VL_sTC':       300,
	'VM_sTC_m1':    150,
	'mt_RTN':       300,
	'ss_RTN_o':     600,
	'ss_RTN_m':     600,
	'ss_RTN_i':     600,
	}

## Local populations
### Layer 1:
netParams.popParams['NGF1']  =   {'cellModel': 'HH_reduced',       'cellType': 'NGF', 'ynormRange': layer['1'], 'density': density[('M1','nonVIP')][0]}

### Layer 2/3:
netParams.popParams['IT2']  =   {'cellModel': cfg.cellmod['IT2'],  'cellType': 'IT',  'ynormRange': layer['2'], 'density': density[('M1','E')][1]}
netParams.popParams['SOM2'] =   {'cellModel': 'HH_reduced',        'cellType': 'SOM', 'ynormRange': layer['2'], 'density': density[('M1','SOM')][1]}
netParams.popParams['PV2']  =   {'cellModel': 'HH_reduced',        'cellType': 'PV',  'ynormRange': layer['2'], 'density': density[('M1','PV')][1]}
netParams.popParams['VIP2']  =  {'cellModel': 'HH_reduced',        'cellType': 'VIP', 'ynormRange': layer['2'], 'density': density[('M1','VIP')][1]}
netParams.popParams['NGF2']  =  {'cellModel': 'HH_reduced',        'cellType': 'NGF', 'ynormRange': layer['2'], 'density': density[('M1','nonVIP')][1]}

### Layer 4:
netParams.popParams['IT4']  =   {'cellModel': cfg.cellmod['IT4'],  'cellType': 'IT',  'ynormRange': layer['4'], 'density': density[('M1','E')][2]}
netParams.popParams['SOM4'] =   {'cellModel': 'HH_reduced',        'cellType': 'SOM', 'ynormRange': layer['4'], 'density': density[('M1','SOM')][2]}
netParams.popParams['PV4']  =   {'cellModel': 'HH_reduced',        'cellType': 'PV',  'ynormRange': layer['4'], 'density': density[('M1','PV')][2]}
netParams.popParams['VIP4']  =  {'cellModel': 'HH_reduced',        'cellType': 'VIP', 'ynormRange': layer['4'], 'density': density[('M1','VIP')][2]}
netParams.popParams['NGF4']  =  {'cellModel': 'HH_reduced',        'cellType': 'NGF', 'ynormRange': layer['4'], 'density': density[('M1','nonVIP')][2]}

### Layer 5A:
netParams.popParams['IT5A'] =   {'cellModel': cfg.cellmod['IT5A'], 'cellType': 'IT',  'ynormRange': layer['5A'], 'density': density[('M1','E')][3]}
netParams.popParams['SOM5A'] =  {'cellModel': 'HH_reduced',        'cellType': 'SOM', 'ynormRange': layer['5A'], 'density': density[('M1','SOM')][3]}
netParams.popParams['PV5A']  =  {'cellModel': 'HH_reduced',        'cellType': 'PV',  'ynormRange': layer['5A'], 'density': density[('M1','PV')][3]}
netParams.popParams['VIP5A']  = {'cellModel': 'HH_reduced',        'cellType': 'VIP', 'ynormRange': layer['5A'], 'density': density[('M1','VIP')][3]}
netParams.popParams['NGF5A']  = {'cellModel': 'HH_reduced',        'cellType': 'NGF', 'ynormRange': layer['5A'], 'density': density[('M1','nonVIP')][3]}

### Layer 5B:
netParams.popParams['IT5B'] =   {'cellModel': cfg.cellmod['IT5B'], 'cellType': 'IT',  'ynormRange': layer['5B'], 'density': 0.5*density[('M1','E')][4]}
netParams.popParams['PT5B'] =   {'cellModel': cfg.cellmod['PT5B'], 'cellType': 'PT',  'ynormRange': layer['5B'], 'density': 0.5*density[('M1','E')][4]}
netParams.popParams['SOM5B'] =  {'cellModel': 'HH_reduced',        'cellType': 'SOM', 'ynormRange': layer['5B'], 'density': density[('M1','SOM')][4]}
netParams.popParams['PV5B']  =  {'cellModel': 'HH_reduced',        'cellType': 'PV',  'ynormRange': layer['5B'], 'density': density[('M1','PV')][4]}
netParams.popParams['VIP5B']  = {'cellModel': 'HH_reduced',        'cellType': 'VIP', 'ynormRange': layer['5B'], 'density': density[('M1','VIP')][4]}
netParams.popParams['NGF5B']  = {'cellModel': 'HH_reduced',        'cellType': 'NGF', 'ynormRange': layer['5B'], 'density': density[('M1','nonVIP')][4]}

### Layer 6:
netParams.popParams['IT6']  =   {'cellModel': cfg.cellmod['IT6'],  'cellType': 'IT',  'ynormRange': layer['6'], 'density': 0.5*density[('M1','E')][5]}
netParams.popParams['CT6']  =   {'cellModel': cfg.cellmod['CT6'],  'cellType': 'CT',  'ynormRange': layer['6'], 'density': 0.5*density[('M1','E')][5]}
netParams.popParams['SOM6'] =   {'cellModel': 'HH_reduced',        'cellType': 'SOM', 'ynormRange': layer['6'], 'density': density[('M1','SOM')][5]}
netParams.popParams['PV6']  =   {'cellModel': 'HH_reduced',        'cellType': 'PV',  'ynormRange': layer['6'], 'density': density[('M1','PV')][5]}
netParams.popParams['VIP6']  =  {'cellModel': 'HH_reduced',        'cellType': 'VIP', 'ynormRange': layer['6'], 'density': density[('M1','VIP')][1]}
netParams.popParams['NGF6']  =  {'cellModel': 'HH_reduced',        'cellType': 'NGF', 'ynormRange': layer['6'], 'density': density[('M1','nonVIP')][1]}

## THALAMIC POPULATIONS (from prev model) 

# excitatory 
if cfg.addThalSs:
    netParams.popParams['VPL_sTC']      =   {'cellModel': 'HH_full', 'cellType': 'sTC_cell',  'yRange': t_layer['VPL_sTC'],     'xRange':[xmin['VPL_sTC'],    xmax['VPL_sTC']],    'zRange':[zmin,zmax],  'density': thal_density[('thal','sTC')][2]}
    netParams.popParams['VPM_sTC']      =   {'cellModel': 'HH_full', 'cellType': 'sTC_cell',  'yRange': t_layer['VPM_sTC'],     'xRange':[xmin['VPM_sTC'],    xmax['VPM_sTC']],    'zRange':[zmin,zmax],  'density': thal_density[('thal','sTC')][3]}
    netParams.popParams['POm_sTC_s1']   =   {'cellModel': 'HH_full', 'cellType': 'sTC_cell',  'yRange': t_layer['POm_sTC_s1'],  'xRange':[xmin['POm_sTC_s1'], xmax['POm_sTC_s1']], 'zRange':[zmin,zmax],  'density': thal_density[('thal','sTC')][4]/2}
if cfg.addThalMt:
    netParams.popParams['VL_sTC']       =   {'cellModel': 'HH_full', 'cellType': 'sTC_cell',  'yRange': t_layer['VL_sTC'],      'xRange':[xmin['VL_sTC'],     xmax['VL_sTC']],     'zRange':[zmin,zmax],  'density': thal_density[('thal','sTC')][0]}
    netParams.popParams['VM_sTC_m1']    =   {'cellModel': 'HH_fu2ll', 'cellType': 'sTC_cell',  'yRange': t_layer['VM_sTC_m1'],   'xRange':[xmin['VM_sTC_m1'],  xmax['VM_sTC_m1']],  'zRange':[zmin,zmax],  'density': thal_density[('thal','sTC')][1]/2}
    netParams.popParams['POm_sTC_m1']   =   {'cellModel': 'HH_full', 'cellType': 'sTC_cell',  'yRange': t_layer['POm_sTC_m1'],  'xRange':[xmin['POm_sTC_m1'], xmax['POm_sTC_m1']], 'zRange':[zmin,zmax],  'density': thal_density[('thal','sTC')][4]/2}

# inhibitory - RTN
if cfg.addThalSs or cfg.addThalMt:
    netParams.popParams['mt_RTN']       =   {'cellModel': 'HH_full', 'cellType': 'sRE_cell',  'yRange': t_layer['mt_RTN'],       'xRange':[xmin['mt_RTN'],    xmax['mt_RTN']],     'zRange':[zmin,zmax],  'density': thal_density[('thal','sRE')][5]}
    netParams.popParams['ss_RTN_o']     =   {'cellModel': 'HH_full', 'cellType': 'sRE_cell',  'yRange': t_layer['ss_RTN_o'],     'xRange':[xmin['ss_RTN_o'],  xmax['ss_RTN_o']],   'zRange':[zmin,zmax],  'density': thal_density[('thal','sRE')][5]}   
    netParams.popParams['ss_RTN_m']     =   {'cellModel': 'HH_full', 'cellType': 'sRE_cell',  'yRange': t_layer['ss_RTN_m'],     'xRange':[xmin['ss_RTN_m'],  xmax['ss_RTN_m']],   'zRange':[zmin,zmax],  'density': thal_density[('thal','sRE')][5]}
    netParams.popParams['ss_RTN_i']     =   {'cellModel': 'HH_full', 'cellType': 'sRE_cell',  'yRange': t_layer['ss_RTN_i'],     'xRange':[xmin['ss_RTN_i'],  xmax['ss_RTN_i']],   'zRange':[zmin,zmax],  'density': thal_density[('thal','sRE')][5]}



#------------------------------------------------------------------------------
## S1 cells property rules
#------------------------------------------------------------------------------

cellModelsS1 = ['HH_full']
EpopsS1 = ['L23_PC', 'L4_PC', 'L4_SS', 'L4_SP', 
             'L5_TTPC1', 'L5_TTPC2', 'L5_STPC', 'L5_UTPC',
             'L6_TPC_L1', 'L6_TPC_L4', 'L6_BPC', 'L6_IPC', 'L6_UTPC']
IpopsS1 = []
for popName in cfg.S1pops:
    if popName not in EpopsS1:
        IpopsS1.append(popName)

layerS1 = {'1':[0.0, 0.089], '2': [0.089,0.159], '3': [0.159,0.286], '23': [0.089,0.286], '4':[0.286,0.421], 
'5': [0.421,0.684], '6': [0.684,1.0]}  # normalized layer boundaries

for cellName in cfg.S1cells:
	layernumber = cellName[1:2]
	if layernumber == '2':
		netParams.popParams[cellName] = {'cellType': cellName, 'cellModel': 'HH_full', 'ynormRange': layerS1['23'], 'xRange':[300.0,600.0], 
                                        'numCells': int(np.ceil(cfg.scaleDensity*cfg.cellNumberS1[cellName])), 'diversity': True}
	else:
		netParams.popParams[cellName] = {'cellType': cellName, 'cellModel': 'HH_full', 'ynormRange': layerS1[layernumber], 'xRange':[300.0,600.0], 
                                        'numCells': int(np.ceil(cfg.scaleDensity*cfg.cellNumberS1[cellName])), 'diversity': True}

for cellName in cfg.S1cells:
    
    if cfg.cellNumberS1[cellName] < 5:
        morphoNumbers = cfg.cellNumberS1[cellName]
    else:
        morphoNumbers = 5
    
    cellFraction = 1.0/morphoNumbers
    
    for morphoNumber in range(morphoNumbers):
        cellMe = cfg.cellLabelS1[cellName] + '_' + str(morphoNumber+1)
        
        netParams.loadCellParamsRule(label = cellMe, fileName = 'cellsS1/' + cellMe + '_cellParams.json')   
        cellRule = {'conds': {'cellType': cellName}, 'diversityFraction': cellFraction, 'secs': {}}  # cell rule dict
        cellRule['secs'] = netParams.cellParams[cellMe]['secs']     
        cellRule['conds'] = netParams.cellParams[cellMe]['conds']    
        cellRule['conds']['cellType'] = cellName
        cellRule['globals'] = netParams.cellParams[cellMe]['globals']       
        cellRule['secLists'] = netParams.cellParams[cellMe]['secLists']      
        cellRule['secLists']['spiny'] = {}
        cellRule['secLists']['spinyEE'] = {}
        nonSpiny = ['axon_0', 'axon_1']
        cellRule['secLists']['spiny'] = [sec for sec in cellRule['secLists']['all'] if sec not in nonSpiny]
        nonSpinyEE = ['axon_0', 'axon_1', 'soma']
        cellRule['secLists']['spinyEE'] = [sec for sec in cellRule['secLists']['all'] if sec not in nonSpinyEE]
        netParams.cellParams[cellMe] = cellRule   # add dict to list of cell params  

        #-----------------------------------------------------------------------------------#
        if cfg.reducedtestS1:
            cellRule = {'conds': {'cellType': cellName}, 'diversityFraction': cellFraction, 'secs': {}}  # cell rule dict
            cellRule['conds'] = netParams.cellParams[cellMe]['conds']    
            cellRule['secs'] = {}
            cellRule['secs']['soma'] = netParams.cellParams[cellMe]['secs']['soma']
            cellRule['secLists'] = {}
            cellRule['secLists']['spiny'] = ['soma']
            cellRule['secLists']['spinyEE'] = ['soma']
            cellRule['secLists']['all'] = ['soma']
            cellRule['secLists']['basal'] = ['soma']   
            cellRule['secLists']['apical'] = ['soma']    
            netParams.cellParams[cellMe] = cellRule   # add dict to list of cell params   
        #-----------------------------------------------------------------------------------#
       

#------------------------------------------------------------------------------
# load data from S1 conn pre-processing file 
#------------------------------------------------------------------------------
with open('connS1/conn.pkl', 'rb') as fileObj: connData = pickle.load(fileObj)

lmat = connData['lmat']
a0mat = connData['a0mat']
d0 = connData['d0']

a0e = connData['a0mat_exp']
l0e = connData['lmat_exp']
d0e = connData['d0_exp']

a0g = connData['a0mat_gauss']
x0g = connData['x0_gauss']
l0g = connData['lmat_gauss']
d0g = connData['d0_gauss']

dfinal = connData['dfinal']
pmat = {}
pmat[12.5] = connData['pmat12um']
pmat[25] = connData['pmat25um']
pmat[50] = connData['pmat50um']
pmat[75] = connData['pmat75um']
pmat[100] = connData['pmat100um']
pmat[125] = connData['pmat125um']
pmat[150] = connData['pmat150um']
pmat[175] = connData['pmat175um']
pmat[200] = connData['pmat200um'] #max value for d0=200

synperconnNumber = connData['synperconnNumber']
connNumber = connData['connNumber']
decay = connData['decay']
gsyn = connData['gsyn']
use = connData['use']

ConnTypesNumber = connData['ConnTypesNumber'] 
ConnTypes = connData['ConnTypes']   

connIEtype = connData['connIEtype']  
connEItype = connData['connEItype']
parameters_syn = connData['parameters_syn']

physColumnNames = []
syntypes = []
for name,syntype in parameters_syn:    
    if name not in physColumnNames:
        physColumnNames.append(name) 
    if syntype not in syntypes:
        syntypes.append(syntype)
        
dfS6 = pd.DataFrame(index=syntypes, columns=physColumnNames)
for syntype in syntypes:
    for name in physColumnNames:    
        dfS6.loc[syntype][name] = parameters_syn[name,syntype]

#------------------------------------------------------------------------------
# Synaptic mechanism parameters with STP
#------------------------------------------------------------------------------
#  mods from S1 BBP - deterministic version
for syntype in syntypes:
    if syntype > 50:  # Exc
        
        netParams.synMechParams['S1_EE_STP_Det_' + str(syntype)] = {'mod': 'DetAMPANMDA',
                                         'Use': dfS6['use'][syntype], # ± dfS6['useStd'][syntype]
                                         'Dep': dfS6['dep'][syntype], # ± dfS6['depStd'][syntype] 
                                         'Fac': dfS6['fac'][syntype], # ± dfS6['facStd'][syntype]
                                         'tau_d_AMPA': 1.74, # ± 0.18 ms
                                         'tau_r_AMPA': 0.2, 
                                         'tau_r_NMDA': 0.29,
                                         'tau_d_NMDA': 43,   
                                         'NMDA_ratio': 0.8, # ± 0.1 for EE -- experimentally measured for some path?
                                         'mg':1.0, #    0.5mM where exceptionally specified?                                                                
                                            }
        netParams.synMechParams['S1_EI_STP_Det_' + str(syntype)] = {'mod': 'DetAMPANMDA',
                                         'Use': dfS6['use'][syntype], # ± dfS6['useStd'][syntype]
                                         'Dep': dfS6['dep'][syntype], # ± dfS6['depStd'][syntype] 
                                         'Fac': dfS6['fac'][syntype], # ± dfS6['facStd'][syntype]
                                         'tau_d_AMPA': 1.74, # ± 0.18 ms
                                         'tau_r_AMPA': 0.2,
                                         'tau_r_NMDA': 0.29,
                                         'tau_d_NMDA': 43,   
                                         'NMDA_ratio': 0.4, # ± 0.1  for EI -- experimentally measured for some path?
                                         'mg':1.0, #    0.5mM where exceptionally specified?                                                                
                                            }
    else: # Inh
        
        netParams.synMechParams['S1_II_STP_Det_' + str(syntype)] = {'mod': 'DetGABAAB',
                                         'Use': dfS6['use'][syntype], # ± dfS6['useStd'][syntype]
                                         'Dep': dfS6['dep'][syntype], # ± dfS6['depStd'][syntype]  
                                         'Fac': dfS6['fac'][syntype], # ± dfS6['facStd'][syntype]
                                         'tau_d_GABAA': dfS6['decay'][syntype], # ± dfS6['decayStd'][syntype]
                                         'tau_r_GABAA': 0.2,   #rng.lognormal(0.2, 0.1) in synapses.hoc  
                                         'tau_d_GABAB': 260.9,
                                         'tau_r_GABAB': 3.5,
#                                          'GABAB_ratio': 1.0,  #=0(1):The ratio of GABAB to GABAA  ?          
                                            }
        
        netParams.synMechParams['S1_IE_STP_Det_' + str(syntype)] = {'mod': 'DetGABAAB',
                                         'Use': dfS6['use'][syntype], # ± dfS6['useStd'][syntype]
                                         'Dep': dfS6['dep'][syntype], # ± dfS6['depStd'][syntype]  
                                         'Fac': dfS6['fac'][syntype], # ± dfS6['facStd'][syntype]
                                         'tau_d_GABAA': dfS6['decay'][syntype], # ± dfS6['decayStd'][syntype]
                                         'tau_r_GABAA': 0.2,   #rng.lognormal(0.2, 0.1) in synapses.hoc  
                                         'tau_d_GABAB': 260.9,
                                         'tau_r_GABAB': 3.5,
#                                          'GABAB_ratio': 1.0,  #=0(1):The ratio of GABAB to GABAA   ?       
                                            }

# Th->S1 NEW
#E2 -> syn 134
netParams.synMechParams['TC:S1'] = {'mod': 'DetAMPANMDA',
                                          'Dep': 227.0,
                                          'Fac': 13.0,
                                          'Use': 0.72,
                                          'tau_r_AMPA': 0.2,
                                          'tau_d_AMPA': 1.74,
                                          'NMDA_ratio': 0.4,
                                          'tau_r_NMDA': 0.29,
                                          'tau_d_NMDA': 43.0}
# M1_>S1
#E2 -> syn 134
netParams.synMechParams['M1:S1'] = {'mod': 'DetAMPANMDA',
                                          'Dep': 227.0,
                                          'Fac': 13.0,
                                          'Use': 0.72,
                                          'tau_r_AMPA': 0.2,
                                          'tau_d_AMPA': 1.74,
                                          'NMDA_ratio': 0.4,
                                          'tau_r_NMDA': 0.29,
                                          'tau_d_NMDA': 43.0}


#------------------------------------------------------------------------------
# Synaptic mechanism parameters
#------------------------------------------------------------------------------
netParams.synMechParams['NMDA']             = {'mod': 'MyExp2SynNMDABB',    'tau1NMDA': 15,     'tau2NMDA': 150,                'e': 0}
netParams.synMechParams['AMPA']             = {'mod': 'MyExp2SynBB',        'tau1': 0.05,       'tau2': 5.3,                    'e': 0}
netParams.synMechParams['GABAB']            = {'mod': 'MyExp2SynBB',        'tau1': 3.5,        'tau2': 260.9,                  'e': -93} 
netParams.synMechParams['GABAA']            = {'mod': 'MyExp2SynBB',        'tau1': 0.07,       'tau2': 18.2,                   'e': -80}
netParams.synMechParams['GABAA_VIP']        = {'mod': 'MyExp2SynBB',        'tau1': 0.3,        'tau2': 6.4,                    'e': -80}  # Pi et al 2013
netParams.synMechParams['GABAASlow']        = {'mod': 'MyExp2SynBB',        'tau1': 2,          'tau2': 100,                    'e': -80}
netParams.synMechParams['GABAASlowSlow']    = {'mod': 'MyExp2SynBB',        'tau1': 200,        'tau2': 400,                    'e': -80}

ESynMech    = ['AMPA', 'NMDA']
SOMESynMech = ['GABAASlow','GABAB']
SOMISynMech = ['GABAASlow']
PVSynMech   = ['GABAA']
VIPSynMech  = ['GABAA_VIP']
NGFSynMech  = ['GABAA', 'GABAB']


# Spont and BG_S1
netParams.synMechParams['NMDA_S1'] = {'mod': 'MyExp2SynNMDABB', 'tau1NMDA': 0.29, 'tau2NMDA': 43, 'e': 0}
netParams.synMechParams['AMPA_S1'] = {'mod':'MyExp2SynBB', 'tau1': 0.2, 'tau2': 1.74, 'e': 0}
netParams.synMechParams['GABAA_S1'] = {'mod':'MyExp2SynBB', 'tau1': 0.2, 'tau2': 8.3, 'e': -80}
netParams.synMechParams['GABAB_S1'] = {'mod':'MyExp2SynBB', 'tau1': 3.5, 'tau2': 260.9, 'e': -93} 
ESynMech_S1 = ['AMPA_S1', 'NMDA_S1']
ISynMech_S1 = ['GABAA_S1', 'GABAB_S1']

# Th
netParams.synMechParams['NMDA_Th']             = {'mod': 'MyExp2SynNMDABB',    'tau1NMDA': 15, 'tau2NMDA': 150,                'e': 0}
netParams.synMechParams['AMPA_Th']             = {'mod': 'MyExp2SynBB',        'tau1': 0.05,   'tau2': 5.3, 'e': 0}
netParams.synMechParams['GABAB_Th']            = {'mod': 'MyExp2SynBB',        'tau1': 3.5,    'tau2': 260.9,                  'e': -93} 
netParams.synMechParams['GABAA_Th']            = {'mod': 'MyExp2SynBB',        'tau1': 0.07,   'tau2': 18.2,                   'e': -80}
ESynMech_Th    = ['AMPA_Th', 'NMDA_Th']
PVSynMech_Th   = ['GABAA_Th']
NGFSynMech_Th  = ['GABAA_Th', 'GABAB_Th']

#------------------------------------------------------------------------------
# S1 Local connectivity parameters 
#------------------------------------------------------------------------------
contA = 0

synWeightFactor = [1.0]

if cfg.connect_S1_S1:    
    for pre in IpopsS1+EpopsS1:
        for post in IpopsS1+EpopsS1:
            if float(connNumber[pre][post]) > 0:           
                # ------------------------------------------------------------------------------    
                #  2D distance prob rules
                # ------------------------------------------------------------------------------ 
                if int(float(d0[pre][post])) < 25:    # single fit
                    if 'exp' in connData['best_fit'][pre][post]:  # exponential            
                        prob = '%s*exp(-dist_2D/%s)*(dist_2D<%s)' % (a0mat[pre][post],lmat[pre][post],dfinal[pre][post])      
                    else: # gaussian
                        prob = '%s*exp(-(dist_2D-%s)**2/(2*%s**2))*(dist_2D<%s)' % (a0g[pre][post],x0g[pre][post],l0g[pre][post],dfinal[pre][post])   
                        
                else:
                    if 'expl' in connData['best_fit'][pre][post]:  # exponential + linear interpolation [25:d0]
                        if int(float(d0[pre][post])) == 25:    #d0==25 -> exponential fit when dist_2D>25, else prob[0um:25um] = pmat[12.5]
                            prob = '%s*exp(-dist_2D/%s)*(dist_2D<%s) if dist_2D > %s else %f' % (a0mat[pre][post],lmat[pre][post],dfinal[pre][post],d0[pre][post],float(pmat[12.5][pre][post]))
                        else:    #d0>25 -> exponential fit when dist_2D>d0, else prob[0um:d0] = linear interpolation [25:d0]
                            d01 = int(float(d0[pre][post]))
                            y1 = float(pmat[25][pre][post])
                            y2 = float(pmat[d01][pre][post])
                            x1 = 25
                            x2 = d01                   
                            angular = (y2 - y1)/(x2 - x1)
                            linear = y2 - x2*angular
                            prob = '%s*exp(-dist_2D/%s)*(dist_2D<%s) if dist_2D > %s else %f * dist_2D + %f' % (a0mat[pre][post],lmat[pre][post],dfinal[pre][post],d0[pre][post],angular,linear)
                   
                    elif 'exp' in connData['best_fit'][pre][post]:  # exponential     
                        if float(pmat[12.5][pre][post]) > float(pmat[25][pre][post]):
                            prob = '%s*exp(-dist_2D/%s)*(dist_2D<%s)' % (a0e[pre][post],l0e[pre][post],dfinal[pre][post])
                        else:  
                            prob = '%s*exp(-dist_2D/%s)*(dist_2D<%s) if dist_2D > %s else %f' % (a0e[pre][post],l0e[pre][post],dfinal[pre][post],d0e[pre][post],float(pmat[12.5][pre][post]))      
                    
                    else: # gaussian
                        prob = '%s*exp(-(dist_2D-%s)**2/(2*%s**2))*(dist_2D<%s)' % (a0g[pre][post],x0g[pre][post],l0g[pre][post],dfinal[pre][post])             
                        
                # ------------------------------------------------------------------------------    
                # I -> I
                # ------------------------------------------------------------------------------
                if pre in IpopsS1:
                    if post in IpopsS1:                             
                        connID = ConnTypes[pre][post][0]                        
                        synMechType = 'S1_II_STP_Det_' + str(connID)   
                        contA+= 1
                        netParams.connParams['II_' + pre + '_' + post] = { 
                                        'preConds': {'pop': cfg.popLabelElS1[pre]}, 
                                        'postConds': {'pop': cfg.popLabelElS1[post]},
                                        'synMech': synMechType,
                                        'probability': prob,
                                        'weight': parameters_syn['gsyn',connID],
                                        'synMechWeightFactor': synWeightFactor,
                                        'delay': 'defaultDelayS1+dist_3D/propVelocityS1',
                                        'synsPerConn': int(synperconnNumber[pre][post]+0.5),
                                        'sec': 'spiny'}        
                # ------------------------------------------------------------------------------
                #  I -> E  # with ME conn diversity
                # ------------------------------------------------------------------------------
                if pre in IpopsS1:
                    if post in EpopsS1:                                                       
                        cellpreList_A = []
                        cellpreList_B = []
                        cellpreList_C = []
                        connID_B = -1    
                        connID_C = -1                               
                        if 'SBC' in pre or 'LBC' in pre or 'NBC' in pre:                              
                            cellpost = cfg.popLabelElS1[post][0]   
                            for npre,cellpre in enumerate(cfg.popLabelElS1[pre]):   
                                premtype = pre[-3:]
                                preetype = cellpre[-3:]                                    
                                connID = connIEtype[premtype][preetype]                                     
                                if connID == ConnTypes[pre][post][0]:
                                    cellpreList_A.append(cellpre)    
                                elif connID == ConnTypes[pre][post][1]:
                                    cellpreList_B.append(cellpre)
                                    connID_B = ConnTypes[pre][post][1]
                                elif connID == ConnTypes[pre][post][2]:
                                    cellpreList_C.append(cellpre)
                                    connID_C = ConnTypes[pre][post][2]
                                else:
                                    print('ERROR')                                    
                        else:   
                            cellpreList_A = cfg.popLabelElS1[pre]                              
                            
                        connID = ConnTypes[pre][post][0]                            
                        synMechType = 'S1_IE_STP_Det_' + str(connID)
                        
                        contA+= 1                          
                        netParams.connParams['IE_'+pre+'_'+post] = { 
                                    'preConds': {'pop': cellpreList_A}, 
                                    'postConds': {'pop': cfg.popLabelElS1[post]},
                                    'synMech': synMechType,
                                    'probability': prob,
                                    'weight': parameters_syn['gsyn',connID],
                                    'synMechWeightFactor': synWeightFactor,
                                    'delay': 'defaultDelayS1+dist_3D/propVelocityS1',
                                    'synsPerConn': int(synperconnNumber[pre][post]+0.5),
                                    'sec': 'spiny'}  
                

                        if connID_B >= 0:          
                            connID = connID_B
                            synMechType = 'S1_IE_STP_Det_' + str(connID)         
                            netParams.connParams['IE_'+pre+'_'+post+'_B'] = { 
                                        'preConds': {'pop': cellpreList_B}, 
                                        'postConds': {'pop': cfg.popLabelElS1[post]},
                                        'synMech': synMechType,
                                        'probability': prob,
                                        'weight': parameters_syn['gsyn',connID],
                                        'synMechWeightFactor': synWeightFactor,
                                        'delay': 'defaultDelayS1+dist_3D/propVelocityS1',
                                        'synsPerConn': int(synperconnNumber[pre][post]+0.5),
                                        'sec': 'spiny'}                       
                
                                
                            if connID_C >= 0:          
                                connID = connID_C
                                synMechType = 'S1_IE_STP_Det_' + str(connID)         
                                netParams.connParams['IE_'+pre+'_'+post+'_C'] = { 
                                            'preConds': {'pop': cellpreList_C}, 
                                            'postConds': {'pop': cfg.popLabelElS1[post]},
                                            'synMech': synMechType,
                                            'probability': prob,
                                            'weight': parameters_syn['gsyn',connID],
                                            'synMechWeightFactor': synWeightFactor,
                                            'delay': 'defaultDelayS1+dist_3D/propVelocityS1',
                                            'synsPerConn': int(synperconnNumber[pre][post]+0.5),
                                            'sec': 'spiny'}                       
                                
                #------------------------------------------------------------------------------   
                # E -> E
                #------------------------------------------------------------------------------
                if pre in EpopsS1:
                    if post in EpopsS1:    
                        connID = ConnTypes[pre][post][0]                        
                        synMechType = 'S1_EE_STP_Det_' + str(connID)   
                        contA+= 1   
                        netParams.connParams['EE_'+pre+'_'+post] = { 
                            'preConds': {'pop': cfg.popLabelElS1[pre]}, 
                            'postConds': {'pop': cfg.popLabelElS1[post]},
                            'synMech': synMechType,
                            'probability': prob, 
                            'weight': parameters_syn['gsyn',connID],
                            'synMechWeightFactor': synWeightFactor,
                            'delay': 'defaultDelayS1+dist_3D/propVelocityS1',
                            'synsPerConn': int(synperconnNumber[pre][post]+0.5),
                            'sec': 'spinyEE'}    
    
                #------------------------------------------------------------------------------               
                #  E -> I  with ME conn diversity
                #------------------------------------------------------------------------------   
                if pre in EpopsS1:
                    if post in IpopsS1:                        
                        cellpostList_A = []
                        cellpostList_B = []
                        connID_B = -1                          
                        if ConnTypes[pre][post][0] == 131 or ConnTypes[pre][post][0] == 132: # EXCEPTIONS -> L6_IPC:L6_(DBC-LBC-NBC-SBC) and  L6_TPC_L:L6_(DBC-LBC-NBC-SBC)    
                            cellpostList_A = cfg.popLabelElS1[post]     
                        elif 'LBC' in post or 'NBC' in post or 'BP' in post or 'DBC' in post or 'BTC' in post:    
                            cellpre = cfg.popLabelElS1[pre][0]
                            for npost,cellpost in enumerate(cfg.popLabelElS1[post]):                                
                                postmtype = post[-3:]
                                postetype = cellpost[-3:]
                                if 'BP' in postmtype:
                                    postmtype = post[-2:]       
                                connID = connEItype[postmtype][postetype]                                
                                if connID == ConnTypes[pre][post][0]:
                                    cellpostList_A.append(cellpost)    
                                elif connID == ConnTypes[pre][post][1]:
                                    cellpostList_B.append(cellpost)
                                    connID_B = ConnTypes[pre][post][1]
                                else:
                                    print('ERROR')                                
                        else:                           
                            cellpostList_A = cfg.popLabelElS1[post]         
                             
                        connID = ConnTypes[pre][post][0]      
                        synMechType = 'S1_EI_STP_Det_' + str(connID)  
                        contA+= 1                                                              
                        netParams.connParams['EI_'+pre+'_'+post] = { 
                                        'preConds': {'pop': cfg.popLabelElS1[pre]}, 
                                        'postConds': {'pop': cellpostList_A},
                                        'synMech': synMechType,
                                        'probability': prob, 
                                        'weight': parameters_syn['gsyn',connID],
                                        'synMechWeightFactor': synWeightFactor,
                                        'delay': 'defaultDelayS1+dist_3D/propVelocityS1',
                                        'synsPerConn': int(synperconnNumber[pre][post]+0.5),
                                        'sec': 'spiny'}   

                        if connID_B >= 0:          
                            connID = connID_B
                            synMechType = 'S1_EI_STP_Det_' + str(connID)        
                            netParams.connParams['EI_'+pre+'_'+post+'_B'] = { 
                                            'preConds': {'pop': cfg.popLabelElS1[pre]}, 
                                            'postConds': {'pop': cellpostList_B},
                                            'synMech': synMechType,
                                            'probability': prob, 
                                            'weight': parameters_syn['gsyn',connID],
                                            'synMechWeightFactor': synWeightFactor,
                                            'delay': 'defaultDelayS1+dist_3D/propVelocityS1',
                                            'synsPerConn': int(synperconnNumber[pre][post]+0.5),
                                            'sec': 'spiny'}   

#------------------------------------------------------------------------------
# NetStim inputs to simulate Spontaneous synapses + background in S1 neurons - data from Rat
#------------------------------------------------------------------------------
SourcesNumber = 5 # for each post Mtype - sec distribution
synperNeuronStimI = connData['synperNeuronStimI']
synperNeuronStimE = connData['synperNeuronStimE']
GsynStimI = connData['GsynStimI']
GsynStimE = connData['GsynStimE']
   
if cfg.addStimSynS1:      
    for post in IpopsS1 + EpopsS1:

        synperNeuron = synperNeuronStimI[post]
        ratespontaneous = cfg.rateStimIS1
        for qSnum in range(SourcesNumber):
            ratesdifferentiation = (0.8 + 0.4*qSnum/(SourcesNumber-1)) * (synperNeuron*ratespontaneous)/SourcesNumber
            netParams.stimSourceParams['StimSynS1_S_all_INH->' + post + '_' + str(qSnum)] = {'type': 'NetStim', 'rate': ratesdifferentiation, 'noise': 1.0}

        synperNeuron = synperNeuronStimE[post]
        ratespontaneous = cfg.rateStimES1
        for qSnum in range(SourcesNumber):
            ratesdifferentiation = (0.8 + 0.4*qSnum/(SourcesNumber-1)) * (synperNeuron*ratespontaneous)/SourcesNumber
            netParams.stimSourceParams['StimSynS1_S_all_EXC->' + post + '_' + str(qSnum)] = {'type': 'NetStim', 'rate': ratesdifferentiation, 'noise': 1.0}
            
    #------------------------------------------------------------------------------
    for post in EpopsS1:
        for qSnum in range(SourcesNumber):
            netParams.stimTargetParams['StimSynS1_T_all_EXC->' + post + '_' + str(qSnum)] = {
                'source': 'StimSynS1_S_all_EXC->' + post + '_' + str(qSnum), 
                'conds': {'cellType': cfg.popLabelElS1[post]}, 
                'synMech': 'AMPA_S1', 
                'sec': 'spinyEE', 
                'weight': 0.001*GsynStimE[post],
                'delay': 0.1}

    for post in IpopsS1:
        for qSnum in range(SourcesNumber):
            netParams.stimTargetParams['StimSynS1_T_all_EXC->' + post + '_' + str(qSnum)] = {
                'source': 'StimSynS1_S_all_EXC->' + post + '_' + str(qSnum), 
                'synMech': 'AMPA_S1', 
                'conds': {'cellType': cfg.popLabelElS1[post]}, 
                'sec': 'spiny', 
                'weight': 0.001*GsynStimE[post],
                'delay': 0.1}

    for post in EpopsS1+IpopsS1:
        for qSnum in range(SourcesNumber):
            netParams.stimTargetParams['StimSynS1_T_all_INH->' + post + '_' + str(qSnum)] = {
                'source': 'StimSynS1_S_all_INH->' + post + '_' + str(qSnum), 
                'conds': {'cellType': cfg.popLabelElS1[post]}, 
                'synMech': 'GABAA_S1', 
                'sec': 'spiny', 
                'weight': 0.001*GsynStimI[post],
                'delay': 0.1}

#------------------------------------------------------------------------------
# Th->S1 connectivity parameters
#------------------------------------------------------------------------------
if cfg.connect_Th_S1:

    # mtype VPM_sTC POm_sTC_s1 nameref
    with open('../info/anatomy/convergence_Th_S1.txt') as mtype_file:
        mtype_content = mtype_file.read()       

    convergence_Th_S1 = {}
    convergence_Th_S1['VPM_sTC'] = {}
    convergence_Th_S1['VPL_sTC'] = {}
    convergence_Th_S1['POm_sTC_s1'] = {}

    for line in mtype_content.split('\n')[:-1]:
        mtype, preFO, preHO, nameref  = line.split()
        convergence_Th_S1['VPL_sTC'][mtype] = int(preFO) # First Order  
        convergence_Th_S1['VPM_sTC'][mtype] = int(preFO) # First Order
        convergence_Th_S1['POm_sTC_s1'][mtype] = int(preHO) # High Order 

    ## Connectivity rules
    radius_cilinder = netParams.sizeX/2.0
    synapsesperconnection_Th_S1 = 9.0
    radius2D_Th_S1 = 50.0

    for pre in ['VPL_sTC', 'VPM_sTC', 'POm_sTC_s1']:  #  
        if cfg.TC_S1[pre]:
            for post in EpopsS1+IpopsS1: 
                
                conn_convergence = np.ceil(convergence_Th_S1[pre][post]/synapsesperconnection_Th_S1)
                prob_conv = 1.0*(conn_convergence/cfg.popNumberS1[pre])*((radius_cilinder**2)/(radius2D_Th_S1**2)) # prob*(AreaS1/Area_Th_syn)  
                probability_rule = '%f if dist_2D < %f else 0.0' % (prob_conv,radius2D_Th_S1)

                netParams.connParams['thal_'+pre+'_'+post] = { 
                    'preConds': {'pop': pre}, 
                    'postConds': {'pop': cfg.popLabelElS1[post]},
                    'weight': 0.19,   # synaptic weight 
                    'sec': 'spinyEE', # target postsyn section
                    'delay': 'defaultDelay+dist_3D/propVelocity',
                    'synsPerConn': int(synapsesperconnection_Th_S1), 
                    'synMech': 'TC:S1'}  

                if pre=='POm_sTC_s1':
                    netParams.connParams['thal_'+pre+'_'+post]['convergence'] = conn_convergence # non-topographycal connectivity
                else:
                    netParams.connParams['thal_'+pre+'_'+post]['probability'] = probability_rule # FO (First Order)

#------------------------------------------------------------------------------
# S1-> Th
#------------------------------------------------------------------------------
if cfg.connect_S1_Th:

    ## load data from conn pre-processing file
    with open('connS1/conn_Th.pkl', 'rb') as fileObj: connData = pickle.load(fileObj)
    wmat = connData['wmat']
    cmat = connData['cmat']
    
    pops_TC     = ['VPL_sTC','VPM_sTC', 'POm_sTC_s1']
    pops_RTN    = ['ss_RTN_o', 'ss_RTN_m', 'ss_RTN_i']
    pops_FO     = ['VPL_sTC','VPM_sTC']
    pops_HO     = ['POm_sTC_s1'],

    pops_CT     = ['L5_TTPC2', 'L6_TPC_L4']

    radius2D_S1_RTN = 50.0
    radius2D_S1_TC = 50.0
    radius_cilinder = netParams.sizeX/2.0

    if cfg.connect_S1_RTN:
        for pre in pops_CT:
            for post in pops_RTN:

                syn = ['AMPA_Th'] # AMPA
                synWeightFactor = [1.0]

                conn_method = 'probability'

                conn_convergence = cfg.convergence_S1_RTN
                prob_conv = 1.0*(conn_convergence/cfg.popNumberS1[pre])*((radius_cilinder**2)/(radius2D_S1_RTN**2)) # prob*(AreaS1/Area_Th_syn)  
                prob_rule = '%f if dist_2D < %f else 0.0' % (prob_conv,radius2D_S1_RTN)

                netParams.connParams['thal_'+pre+'_'+post] = { 
                                'preConds': {'pop': cfg.popLabelElS1[pre]}, 
                                'postConds': {'pop': post},
                                'synMech': syn,
                                conn_method:  prob_rule,
                                'weight': cfg.connWeight_S1_RTN, 
                                'synMechWeightFactor': synWeightFactor,
                                'delay': 'defaultDelay+dist_3D/propVelocity',
                                'synsPerConn': 1,
                                'sec': 'soma'}

    if cfg.connect_S1_TC:
        for pre in pops_CT:
            for post in pops_TC:

                syn = ['AMPA_Th'] # AMPA
                synWeightFactor = [1.0]

                if post in pops_HO:
                    conn_method = 'divergence'
                    prob_rule = cfg.divergenceHO/2.0
                else: # topographycal connectivity
                    conn_method = 'probability'
                    conn_convergence = cfg.convergence_S1_TC
                    prob_conv = 1.0*(conn_convergence/cfg.popNumberS1[pre])*((radius_cilinder**2)/(radius2D_S1_TC**2)) # prob*(AreaS1/Area_Th_syn)  
                    prob_rule = '%f if dist_2D < %f else 0.0' % (prob_conv,radius2D_S1_TC)

                    netParams.connParams['thal_'+pre+'_'+post] = { 
                                'preConds': {'pop': cfg.popLabelElS1[pre]}, 
                                'postConds': {'pop': post},
                                'synMech': syn,
                                conn_method:  prob_rule,
                                'weight': cfg.connWeight_S1_TC, 
                                'synMechWeightFactor': synWeightFactor,
                                'delay': 'defaultDelay+dist_3D/propVelocity',
                                'synsPerConn': 1,
                                'sec': 'soma'}

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

# thalamictest = ['sTC_cell', 'sRE_cell']

# if cfg.singleCellPops:
#     for pop in netParams.popParams.values(): 
#         if pop['cellType'] not in thalamictest:
#             pop['numCells'] = 1
#             # print(pop['cellType'])

if cfg.singleCellPops:
    for pop in netParams.popParams.values(): 
        pop['numCells'] = 1

# Removes a specific model from the simulation (M1/S1/Thalamus)
if cfg.removeM1:
    for pop in netParams.popParams.keys(): 
        if pop in cfg.M1_pops:
            netParams.popParams[pop]['numCells'] = 0
            print("Removed ", pop)
if cfg.removeS1:
    for pop in netParams.popParams.keys(): 
        if pop in cfg.S1_pops:
            netParams.popParams[pop]['numCells'] = 0
            print("Removed ", pop)
if cfg.removeTh:
    for pop in netParams.popParams.keys(): 
        if pop in cfg.Th_pops:
            netParams.popParams[pop]['numCells'] = 0
            print("Removed ", pop)

#------------------------------------------------------------------------------
## Long-range input populations (VecStims)
if cfg.addLongConn:
    ## load experimentally based parameters for long range inputs
    with open('../conn/conn_long.pkl', 'rb') as fileObj: connLongData = pickle.load(fileObj)
    #ratesLong = connLongData['rates']

    numCells    = cfg.numCellsLong
    noise       = cfg.noiseLong
    start       = cfg.startLong

    longPops = ['S2', 'cM1', 'M2', 'OC']          # (2022_01_21 - removed S1 to add S1 model connections instead) (2021 - removed VL and POm because now they have separate connectivity
    # longPops = ['S1', 'S2', 'cM1', 'M2', 'OC']          # removed VL and POm because now they have separate connectivity
    # longPops = ['POm_sTC_m1', 'VL_sTC', 'S1', 'S2', 'cM1', 'M2', 'OC']
    ## create populations with fixed 
    for longPop in longPops:
        netParams.popParams[longPop] = {'cellModel': 'VecStim', 'numCells': numCells, 'rate': cfg.ratesLong[longPop], 
                                        'noise': noise, 'start': start, 'pulses': [], 'ynormRange': layer['long'+longPop]}
        if isinstance(cfg.ratesLong[longPop], str): # filename to load spikes from
            spikesFile = cfg.ratesLong[longPop]
            with open(spikesFile, 'r') as f: spks = json.load(f)
            netParams.popParams[longPop].pop('rate')
            netParams.popParams[longPop]['spkTimes'] = spks

#------------------------------------------------------------------------------
# Long range input pulses
#------------------------------------------------------------------------------
if cfg.addPulses:
    for key in [k for k in dir(cfg) if k.startswith('pulse')]:
        params = getattr(cfg, key, None)
        [pop, start, end, rate, noise] = [params[s] for s in ['pop', 'start', 'end', 'rate', 'noise']]
        if 'duration' in params and params['duration'] is not None and params['duration'] > 0:
            end = start + params['duration']

        if pop in netParams.popParams:
            if 'pulses' not in netParams.popParams[pop]: netParams.popParams[pop]['pulses'] = {}    
            netParams.popParams[pop]['pulses'].append({'start': start, 'end': end, 'rate': rate, 'noise': noise})

#------------------------------------------------------------------------------
# Current inputs (IClamp)
#------------------------------------------------------------------------------
if cfg.addIClamp:
    for key in [k for k in dir(cfg) if k.startswith('IClamp')]:
        params = getattr(cfg, key, None)
        [pop,sec,loc,start,dur,amp] = [params[s] for s in ['pop','sec','loc','start','dur','amp']]

        #cfg.analysis['plotTraces']['include'].append((pop,0))  # record that pop

        # add stim source
        netParams.stimSourceParams[key] = {'type': 'IClamp', 'delay': start, 'dur': dur, 'amp': amp}
        
        # connect stim source to target
        netParams.stimTargetParams[key+'_'+pop] =  {
            'source':   key, 
            'conds':    {'pop': pop},
            'sec':      sec, 
            'loc':      loc}

#------------------------------------------------------------------------------
# NetStim inputs
#------------------------------------------------------------------------------
if cfg.addNetStim:
    for key in [k for k in dir(cfg) if k.startswith('NetStim')]:
        params = getattr(cfg, key, None)
        [pop, ynorm, sec, loc, synMech, synMechWeightFactor, start, interval, noise, number, weight, delay] = \
        [params[s] for s in ['pop', 'ynorm', 'sec', 'loc', 'synMech', 'synMechWeightFactor', 'start', 'interval', 'noise', 'number', 'weight', 'delay']] 

        # cfg.analysis['plotTraces']['include'] = [(pop,0)]

        if synMech == ESynMech:
            wfrac = cfg.synWeightFractionEE
        elif synMech == SOMESynMech:
            wfrac = cfg.synWeightFractionSOME
        else:
            wfrac = [1.0]

        # add stim source
        netParams.stimSourceParams[key] = {'type': 'NetStim', 'start': start, 'interval': interval, 'noise': noise, 'number': number}

        # connect stim source to target
        # for i, syn in enumerate(synMech):
        netParams.stimTargetParams[key+'_'+pop] =  {
            'source':   key, 
            'conds':    {'pop': pop, 'ynorm': ynorm},
            'sec':      sec, 
            'loc':      loc,
            'synMech':  synMech,
            'weight':   weight,
            'synMechWeightFactor': synMechWeightFactor,
            'delay':    delay}

#------------------------------------------------------------------------------
# Local connectivity parameters
#------------------------------------------------------------------------------

if cfg.connectThalamusNetwork:

    ## load data from conn pre-processing file
    # with open('../conn/thal_conn_new.pkl', 'rb') as fileObj: connData = pickle.load(fileObj)
    with open('../conn/conn_thal_2021_12_14.pkl', 'rb') as fileObj: connData = pickle.load(fileObj)
    pmat = connData['pmat']
    lmat = connData['lmat']
    wmat = connData['wmat']
    cmat = connData['cmat']

    # bins = connData['bins']                       # removed in conn.py and conn.pkl
    # connDataSource = connData['connDataSource']   # removed in conn.py and conn.pkl

    #------------------------------------------------------------------------------
    # Thalamic connectivity parameters
    #------------------------------------------------------------------------------
    # Intrathalamic 
    Epops = [   'VL_sTC',       'VM_sTC_m1',
                'VPL_sTC',      'VPM_sTC',
                'POm_sTC_m1',   
                'POm_sTC_s1',
                ]
    Ipops = [   'mt_RTN', 'ss_RTN_o', 'ss_RTN_m', 'ss_RTN_i',
                ]    
    pops_TC     = [ 'VPL_sTC','VPM_sTC',
                    'POm_sTC_m1',
                    'POm_sTC_s1',
                    'VL_sTC','VM_sTC_m1',
                    ]
    pops_RTN    = [ 'mt_RTN',
                    'ss_RTN_o', 'ss_RTN_m', 'ss_RTN_i',
                    ]
    pops_FO     = [
                    'VPL_sTC','VPM_sTC',
                    'VL_sTC','VM_sTC_m1','VM_sTC_s1', # (Lam, 2015)
                    ]
    pops_HO     = ['POm_sTC_m1','POm_sTC_s1',]

    if cfg.addTopographicalConn:

        if cfg.connect_RTN_RTN:
            for pre in netParams.popParams.keys():
                if pre in pmat: # JV - added this fix so new pops don't have to be added in the conn file and can be declared in the netParams file
                    for post in netParams.popParams.keys():
                        if post in pmat[pre]:
                            if pre in pops_RTN and post in pops_RTN:

                                if cfg.connProb_RTN_RTN   is not None:  pmat[pre][post]=cfg.connProb_RTN_RTN
                                if cfg.connWeight_RTN_RTN is not None:  wmat[pre][post]=cfg.connWeight_RTN_RTN

                                l = cmat[pre][post]/4
                                syn = PVSynMech # only GABA A
                                synWeightFactor = [1.0]
                                netParams.connParams['thal_'+pre+'_'+post] = { 
                                    'preConds':     {'pop': pre}, 
                                    'postConds':    {'pop': post},
                                    'synMech':      syn,
                                    'probability':' %f * exp(-dist_3D/%f)\
                                                    *(dist_2D<%f)\
                                                    *(dist_y<(%f/%f))\
                                                    ' % (pmat[pre][post], l, 
                                                        cmat[pre][post], 
                                                        cmat[pre][post],cfg.yConnFactor),
                                    'weight':       wmat[pre][post] * cfg.intraThalamicGain, 
                                    'synMechWeightFactor': synWeightFactor,
                                    'delay':        'defaultDelay+dist_3D/propVelocity',
                                    'synsPerConn':  1,
                                    'sec':          'soma'}

        if cfg.connect_TC_RTN:
            for pre in netParams.popParams.keys():
                if pre in pops_TC:
                    for post in netParams.popParams.keys():
                        if post in pmat[pre]:

                            if cfg.connProb_TC_RTN   is not None:   pmat[pre][post]=cfg.connProb_TC_RTN
                            if cfg.connWeight_TC_RTN is not None:   wmat[pre][post]=cfg.connWeight_TC_RTN

                            l = cmat[pre][post]/4
                            y_thresh    = cmat[pre][post]/5

                            syn = ['AMPA'] # AMPA
                            synWeightFactor = [1.0]

                            if pre in pops_HO:
                                # non-topographycal connectivity - this rule works but still have to decide if need a more random and disperse connectivity
                                # prob_rule = '%f * exp(-dist_2D/500)\
                                #                 \
                                #                 ' % (pmat[pre][post]
                                #                     )
                                conn_method = 'divergence'
                                prob_rule = cfg.divergenceHO
                            else:
                                # topographycal connectivity
                                conn_method = 'probability'
                                prob_rule = '%f * exp(-dist_2D/%f)\
                                                *(dist_2D<%f)\
                                                *(abs(((((pre_y-%f)*(%f-%f))/(%f-%f))+%f)-post_y)<%f)\
                                                \
                                                ' % (pmat[pre][post], l, 
                                                    cmat[pre][post],
                                                    ymin[pre],ymax[post],ymin[post],ymax[pre],ymin[pre],ymin[post],y_thresh
                                                    )

                            netParams.connParams['thal_'+pre+'_'+post] = { 
                                'preConds':     {'pop': pre}, 
                                'postConds':    {'pop': post},
                                'synMech':      syn,
                                conn_method:    prob_rule,
                                'weight':       wmat[pre][post] * cfg.intraThalamicGain, 
                                'synMechWeightFactor': synWeightFactor,
                                'delay':        'defaultDelay+dist_3D/propVelocity',
                                'synsPerConn':  1,
                                'sec':          'soma'}

        if cfg.connect_RTN_TC:
            for pre in netParams.popParams.keys():
                if pre in pops_RTN:
                    for post in netParams.popParams.keys():
                        if (post in pmat[pre]) and (post in pops_TC):

                            if cfg.connProb_RTN_TC   is not None:   pmat[pre][post]=cfg.connProb_RTN_TC
                            if cfg.connWeight_RTN_TC is not None:   wmat[pre][post]=cfg.connWeight_RTN_TC

                            l = cmat[pre][post]/4
                            y_thresh    = cmat[pre][post]/5
                            if post in Epops:
                                if cfg.addBicuculline:  # Destexhe 1996 - Blocking GABAA (and also disconnecting RTN-RTN)
                                    syn = ['GABAB']     # GABA B
                                    synWeightFactor = [1.0]
                                    print('Bicuculline - Destexhe 1996')
                                else:                   # Regular RTN-TC mechanism (GABAa and GABAb)
                                    syn = NGFSynMech    # GABA A and GABA B
                                    synWeightFactor = [0.6,0.4]
                            elif post in Ipops:     # includes sTI and sRE
                                syn = PVSynMech     # only GABA A
                                synWeightFactor = [1.0]
                            
                            if post in pops_HO:
                                # non-topographycal connectivity - this rule works but still have to decide if need a more random and disperse connectivity
                                # prob_rule = '%f * exp(-dist_2D/500)\
                                #                 \
                                #                 ' % (pmat[pre][post]
                                #                     )
                                conn_method = 'divergence'
                                prob_rule = cfg.divergenceHO
                                # print(pre,post,'non-topographycal', prob_rule)
                            else:
                                # topographycal connectivity
                                conn_method = 'probability'
                                prob_rule = '%f * exp(-dist_2D/%f)\
                                                *(dist_2D<%f)\
                                                *(abs(((((pre_y-%f)*(%f-%f))/(%f-%f))+%f)-post_y)<%f)\
                                                \
                                                ' % (pmat[pre][post], l, 
                                                    cmat[pre][post],
                                                    ymin[pre],ymax[post],ymin[post],ymax[pre],ymin[pre],ymin[post],y_thresh
                                                    )
                                # print(pre,post,'topographycal', prob_rule)

                            netParams.connParams['thal_'+pre+'_'+post] = { 
                                'preConds':     {'pop': pre}, 
                                'postConds':    {'pop': post},
                                'synMech':      syn,
                                conn_method:    prob_rule,
                                'weight':       wmat[pre][post] * cfg.intraThalamicGain, 
                                'synMechWeightFactor': synWeightFactor,
                                'delay':        'defaultDelay+dist_3D/propVelocity',
                                'synsPerConn':  1,
                                'sec':          'soma'}

        # TODO ------------- next step: fix the output to the cortex

        if cfg.connect_TC_CTX:

            with open('../conn/conn_long.pkl', 'rb') as fileObj: connLongData = pickle.load(fileObj)

            # load load experimentally based parameters for long range inputs
            cmatLong    = connLongData['cmat']
            binsLong    = connLongData['bins']
            syns        = {'exc': ESynMech, 'inh': 'GABAA'}
            synFracs    = {'exc': cfg.synWeightFractionEE, 'inh': [1.0]}

                            #Pre pop        #post   #celltype   #proj       #Relative strength (conn_strength_shepherd)
            M1_thal_pops=[  # thalamus to M1
                            ('VL_sTC',      'PT',   'exc',      'core',     3),
                            ('VL_sTC',      'IT',   'exc',      'core',     3),
                            ('VL_sTC',      'CT',   'exc',      'core',     1),
                            
                            ('POm_sTC_m1',  'IT',   'exc',      'matrix',   4),]

            # bins=[[0.1, 0.195], [0.195, 0.29], [0.29, 0.37], [0.37, 0.47], [0.47, 0.58], [0.58, 0.6900000000000001], [0.6900000000000001, 0.8], [0.8, 0.9], [0.9, 1.0]]

            for (pre, post, synType, connType, conn_strength_shepherd) in M1_thal_pops:
                # print(pre, post, synType, connType, strength)
                for cellModel in cellModels:
                    for i, conv_key in enumerate(cmatLong.keys()):
                    # for i, binRange in enumerate(binsLong):
                        if conv_key == (pre, post, synType):
                            ruleLabel = pre+'_'+post+'_'+synType+'_'+cellModel+'_'+str(i)
                            # print(ruleLabel, conv_key, connLongData['cmat'][conv_key], connLongData['bins'][(pre,post)])
                            for j, conv in enumerate(connLongData['cmat'][conv_key]):
                                # print(conv_key, cellModel, connLongData['bins'][(pre,post)][j], conv)
                                netParams.connParams[ruleLabel+'_'+str(j)] = { 
                                        'preConds': {'pop': pre}, 
                                        'postConds': {
                                                        'cellModel': cellModel, 
                                                        'cellType': post, 
                                                        'ynorm': connLongData['bins'][(pre,post)][j],
                                                        },
                                        'synMech':      syns[synType],
                                        'convergence':  conv,
                                        'weight':       cfg.weightLong_thalM1 / cfg.synsperconn[cellModel]*conn_strength_shepherd, 
                                        # 'weight': cfg.weightLong_thalM1[pre] / cfg.synsperconn[cellModel]*conn_strength_shepherd, 
                                        'synMechWeightFactor': synFracs[synType],
                                        'delay':        'defaultDelay+dist_3D/propVelocity',
                                        'synsPerConn':  cfg.synsperconn[cellModel],
                                        'sec':          ['Adend1', 'Adend2', 'Adend3', 'Bdend', 'soma'],
                                        }

        if cfg.connect_CTX_TC:
            syns        = {'exc': ESynMech, 'inh': 'GABAA'}
            synFracs    = {'exc': cfg.synWeightFractionEE, 'inh': [1.0]}

                            #Pre pop #post          #celltype   #proj       #Relative strength
            thal_M1_pops=[  # M1 to thalamus
                            ('PT',  'VM_sTC_m1',    'exc',      'core',     4),
                            ('PT',  'POm_sTC_m1',   'exc',      'core',     4),
                            
                            ('CT',  'VM_sTC_m1',    'exc',      'core',     4),
                            ('CT',  'POm_sTC_m1',   'exc',      'core',     4),
                            
                            ('CT',  'VL_sTC',       'exc',      'core',     1),
                            
                            ('PT',  'mt_RTN',       'exc',      'core',     4),
                            ('CT',  'mt_RTN',       'exc',      'core',     4),]

            for (pre, post, synType, connType, conn_strength_shepherd) in thal_M1_pops:
                # print(pre, post, synType, connType, strength)
                for cellModel in cellModels:
                    for pop in netParams.popParams.keys():
                        if pop.startswith(pre):
                            # print(pop,post)
                            ruleLabel=pop+'_'+post+'_'+synType+'_'+cellModel
                            #print(ruleLabel,conn_strength_shepherd)
                            netParams.connParams[ruleLabel] = { 
                                                                'preConds':     {'pop':pop},
                                                                'postConds':    {'pop': post}, 
                                                                'synMech':      syns[synType],
                                                                'convergence':  5,                                                       #dummy value for now
                                                                'weight':       cfg.weightLong_M1thal / cfg.synsperconn[cellModel]*conn_strength_shepherd,      #dummy value for now
                                                                # 'weight': 0.5 / cfg.synsperconn[cellModel]*conn_strength_shepherd,      #dummy value for now
                                                                # 'weight': cfg.weightLong[pre] / cfg.synsperconn[cellModel], 
                                                                'synMechWeightFactor': synFracs[synType],
                                                                'delay':        'defaultDelay+dist_3D/propVelocity',
                                                                'synsPerConn':  cfg.synsperconn[cellModel],
                                                                'sec':          'soma',
                                                                }

#------------------------------------------------------------------------------
M1_allPops=['NGF1', 
            'IT2',  'SOM2',  'PV2',   'VIP2',  'NGF2', 
            'IT4',  'SOM4',  'PV4',   'VIP4',  'NGF4', 
            'IT5A', 'SOM5A', 'PV5A',  'VIP5A', 'NGF5A', 
            'IT5B', 'PT5B',  'SOM5B', 'PV5B',  'VIP5B', 'NGF5B', 
            'IT6',  'CT6',   'SOM6',  'PV6',   'VIP6',  'NGF6']

M1_pyrPops=[
            'IT2','IT4','IT5A','IT5B','PT5B','IT6','CT6'
            ]

M1_L1_pops  =['NGF1']
M1_L23_pops =['IT2',    'SOM2',     'PV2',      'VIP2',    'NGF2']
M1_L4_pops  =['IT4',    'SOM4',     'PV4',      'VIP4',    'NGF4']
M1_L5A_pops =['IT5A',   'SOM5A',    'PV5A',     'VIP5A',   'NGF5A']
M1_L5B_pops =['IT5B',   'PT5B',     'SOM5B',    'PV5B',    'VIP5B',     'NGF5B']
M1_L6_pops  =['IT6',    'CT6',      'SOM6',     'PV6',     'VIP6',      'NGF6']

Thal_pops   =['VPL_sTC','VPM_sTC',  'POm_sTC_s1', 'VL_sTC',  'VM_sTC_m1', 'POm_sTC_m1', 'mt_RTN', 'ss_RTN_o', 'ss_RTN_m', 'ss_RTN_i']
Long_pops   =['S1',     'S2',       'cM1',        'M2',      'OC']

#------------------------------------------------------------------------------
# M1->S1
'''
    References:
    Kinnischtzke, A. K., Simons, D. J., & Fanselow, E. E. (2014). Motor cortex broadly engages excitatory and inhibitory neurons in somatosensory barrel cortex. Cerebral Cortex , 24(8), 2237–2248.
    https://paperpile.com/app/p/ea0be1ed-6f8b-0046-9225-b447a8bdc8f3 

    Kinnischtzke, A. K., Fanselow, E. E., & Simons, D. J. (2016). Target-specific M1 inputs to infragranular S1 pyramidal neurons. Journal of Neurophysiology, 116(3), 1261–1274.
    https://paperpile.com/app/p/d623473e-41bd-07b3-93ce-44550d371fda 

    Yamawaki, N., Raineri Tapies, M. G., Stults, A., Smith, G. A., & Shepherd, G. M. (2021). Circuit organization of the excitatory sensorimotor loop through hand/forelimb S1 and M1. eLife, 10. https://doi.org/10.7554/eLife.66836
    https://paperpile.com/app/p/c91e1558-68fc-0bcc-a2f4-91c3f4cb0113 

    Mao, T., Kusefoglu, D., Hooks, B. M., Huber, D., Petreanu, L., & Svoboda, K. (2011). Long-range neuronal circuits underlying the interaction between sensory and motor cortex. Neuron, 72(1), 111–123.
    https://paperpile.com/app/p/9b23a6fc-70e7-08fc-af5a-29e1cccc19f1 
'''
if cfg.connect_M1_S1:
    
    synFracs=cfg.synWeightFractionEE
    synapsesperconnection_M1_S1 = 5

    # Extracted from (Kinnischtzke, 2014) (https://paperpile.com/app/p/ea0be1ed-6f8b-0046-9225-b447a8bdc8f3)
    M1_S1_conn=[  # M1->S1
                    # pre=IT???
                    #pre      post     type    syn         prob     mean    sd
                    ('IT',    'L23',   'RS',   ESynMech,   0.9  ,   6.731,  4.016),
                    ('IT',    'L4',    'RS',   ESynMech,   0.333,   3.727,  1.431),
                    ('IT',    'L5A',   'RS',   ESynMech,   0.703,   5.969,  3.038),
                    ('IT',    'L5B',   'RS',   ESynMech,   0.703,   5.969,  3.038),
                    ('IT',    'L6',    'RS',   ESynMech,   1.0  ,   11.756, 7.026),
                    
                    ('IT',    'L23',   'FS',   ESynMech,   1.0  ,   5.777,  1.614),
                    ('IT',    'L4',    'FS',   ESynMech,   0.143,   2.21,   None ),
                    ('IT',    'L5A',   'FS',   ESynMech,   0.75 ,   6.053,  0.974),
                    ('IT',    'L5B',   'FS',   ESynMech,   0.75 ,   6.053,  0.974),
                    ('IT',    'L6',    'FS',   ESynMech,   1.0  ,   8.177,  2.995),
                    
                    ('IT',    'L23',   'SOM',  ESynMech,   0.714,   6.147,  3.693),
                    ('IT',    'L4',    'SOM',  ESynMech,   0.2  ,   6.800,  None ),
                    ('IT',    'L5A',   'SOM',  ESynMech,   0.875,   4.974,  2.624),
                    ('IT',    'L5B',   'SOM',  ESynMech,   0.875,   4.974,  2.624),
                    # ('IT',    'L6',    'SOM',  ESynMech,   0.0  ,   None,   None),
                    ]
    
    for (prePop,postPop,cellSubtype,synMech,prob,mean,sd) in M1_S1_conn:
        for M1_pop in M1_allPops:
            if M1_pop.startswith(prePop):

                for post in EpopsS1: 
                    if post[0:2] == postPop[0:2] and cellSubtype == 'RS':
                        ruleLabel = 'M1_S1_'+M1_pop+'_'+post
                        netParams.connParams[ruleLabel] = {
                                                            'preConds':     {'pop': M1_pop},
                                                            'postConds':    {'pop': cfg.popLabelElS1[post]}, 
                                                            'probability':  0.1*prob,
                                                            'weight': 0.19,   # synaptic weight 
                                                            'sec': 'spinyEE', # target postsyn section
                                                            'delay': 'defaultDelay+dist_3D/propVelocity',
                                                            'synsPerConn': int(synapsesperconnection_M1_S1), 
                                                            'synMech': 'M1:S1'}                         
                        # print(ruleLabel,prob)

                for post in IpopsS1: 
                    if post[0:2] == postPop[0:2] and cellSubtype == 'FS' and '_MC' not in post:
                        ruleLabel = 'M1_S1_'+M1_pop+'_'+post
                        netParams.connParams[ruleLabel] = {
                                                            'preConds':     {'pop': M1_pop},
                                                            'postConds':    {'pop': cfg.popLabelElS1[post]}, 
                                                            'probability':  0.1*prob,
                                                            'weight': 0.19,   # synaptic weight 
                                                            'sec': 'spinyEE', # target postsyn section
                                                            'delay': 'defaultDelay+dist_3D/propVelocity',
                                                            'synsPerConn': int(synapsesperconnection_M1_S1), 
                                                            'synMech': 'M1:S1'}                         
                        # print(ruleLabel,prob)

                for post in IpopsS1: 
                    if post[0:2] == postPop[0:2] and cellSubtype == 'SOM' and '_MC' in post:
                        ruleLabel = 'M1_S1_'+M1_pop+'_'+post
                        netParams.connParams[ruleLabel] = {
                                                            'preConds':     {'pop': M1_pop},
                                                            'postConds':    {'pop': cfg.popLabelElS1[post]}, 
                                                            'probability':  0.1*prob,
                                                            'weight': 0.19,   # synaptic weight 
                                                            'sec': 'spinyEE', # target postsyn section
                                                            'delay': 'defaultDelay+dist_3D/propVelocity',
                                                            'synsPerConn': int(synapsesperconnection_M1_S1), 
                                                            'synMech': 'M1:S1'}                          
                        # print(ruleLabel,prob)

# S1->M1
if cfg.connect_S1_M1:

    synFracs=cfg.synWeightFractionEE

    # (Yamawaki, 2021 - Fig 6 C/G) (https://paperpile.com/app/p/c91e1558-68fc-0bcc-a2f4-91c3f4cb0113)
    ''' These results add key details about the excitatory connectivity in the last stage along the transcortical circuit leading to M1, showing that the main recipients of S1 corticocortical input are L2/3
    pyramidal neurons.'''
    # S1-L23 to M1                          S1-L5A to M1 
    # L23 = 3.183419112216958       |       2.7040816326530615
    # L4  = 0.4723746548962131      |       0.18707482993197289
    # L5  = 0.07361077973573238     |       0.06802721088435426
    # L6  = 0                       |       0
    S1_IT_23_5A = ['L23_PC_cAD', 'L5_UTPC_cAD']
    # Rescaling Projection weights to 0.5 mV
    S1_M1_conn  =[  #pre    post    syn         prob    input   rescaled input to 0.5
                    #S1-L2/3
                    ('L23', 'IT2',  ESynMech,   0.2,    3.183,  0.5),                    
                    ('L23', 'IT4',  ESynMech,   0.2,    0.472,  0.074),                    
                    ('L23', 'IT5A', ESynMech,   0.2,    0.074,  0.012),
                    ('L23', 'IT5B', ESynMech,   0.2,    0.074,  0.012),
                    ('L23', 'PT5B', ESynMech,   0.2,    0.074,  0.012),                    
                    # ('L23', 'IT6',  ESynMech,   0.2,    0,      0),
                    # ('L23', 'CT6',  ESynMech,   0.2,    0,      0),
                    #S1-L5A
                    ('L5A', 'IT2',  ESynMech,   0.2,    2.704,  0.425),                    
                    ('L5A', 'IT4',  ESynMech,   0.2,    0.187,  0.030),                    
                    ('L5A', 'IT5A', ESynMech,   0.2,    0.068,  0.011),
                    ('L5A', 'IT5B', ESynMech,   0.2,    0.068,  0.011),
                    ('L5A', 'PT5B', ESynMech,   0.2,    0.068,  0.011),                    
                    # ('L5A', 'IT6',  ESynMech,   0.2,    0,      0),
                    # ('L5A', 'CT6',  ESynMech,   0.2,    0,      0),
                    ]

    for (prePop,postPop,synMech,prob,input,rescaled_input) in S1_M1_conn:

        if prePop == 'L23':
            pre = 'L23_PC'
        if prePop == 'L5A':
            pre = 'L5_UTPC'
            
        ruleLabel = 'S1_M1_'+pre+'_'+postPop
        netParams.connParams[ruleLabel] = {
                                                'preConds':    {'pop': cfg.popLabelElS1[pre]}, 
                                                'postConds':    {'pop': postPop}, 
                                                'synMech':      synMech,
                                                'probability':  0.5*prob,
                                                'weight':       rescaled_input,
                                                'synMechWeightFactor': synFracs,
                                                'delay':        'defaultDelay+dist_3D/propVelocity',
                                                'sec':          'soma',   
                                                }
        # print(ruleLabel,prob,rescaled_input)

#------------------------------------------------------------------------------
# M1

with open('../conn/conn.pkl', 'rb') as fileObj: connData = pickle.load(fileObj)
pmat = connData['pmat']
wmat = connData['wmat']
bins = connData['bins']

#------------------------------------------------------------------------------
## E -> E
if cfg.addConn and cfg.EEGain > 0.0:
    labelsConns     = [ ('W+AS_norm', 'IT', 'L2/3,4'),  
                        ('W+AS_norm', 'IT', 'L5A,5B'), 
                        ('W+AS_norm', 'PT', 'L5B'),     
                        ('W+AS_norm', 'IT', 'L6'),      
                        ('W+AS_norm', 'CT', 'L6')]

    labelPostBins   = [ ('W+AS', 'IT', 'L2/3,4'),       
                        ('W+AS', 'IT', 'L5A,5B'),       
                        ('W+AS', 'PT', 'L5B'), 
                        ('W+AS', 'IT', 'L6'),           
                        ('W+AS', 'CT', 'L6')]

    labelPreBins    = [ 'W', 'AS', 'AS', 'W', 'W']
    preTypes        = [ ['IT'], ['IT'], ['IT', 'PT'], ['IT','CT'], ['IT','CT']] 
    postTypes       = [ 'IT', 'IT', 'PT', 'IT','CT']
    ESynMech        = [ 'AMPA','NMDA']

    for i,(label, preBinLabel, postBinLabel) in enumerate(zip(labelsConns,labelPreBins, labelPostBins)):
        for ipre, preBin in enumerate(bins[preBinLabel]):
            for ipost, postBin in enumerate(bins[postBinLabel]):
                for cellModel in cellModels:
                    ruleLabel = 'EE_'+cellModel+'_'+str(i)+'_'+str(ipre)+'_'+str(ipost)
                    netParams.connParams[ruleLabel] = { 
                        'preConds':     {'cellType': preTypes[i], 'ynorm': list(preBin)}, 
                        'postConds':    {'cellModel': cellModel, 'cellType': postTypes[i], 'ynorm': list(postBin)},
                        'synMech':      ESynMech,
                        'probability':  pmat[label][ipost,ipre],
                        'weight':       wmat[label][ipost,ipre] * cfg.EEGain / cfg.synsperconn[cellModel], 
                        'synMechWeightFactor': cfg.synWeightFractionEE,
                        'delay':        'defaultDelay+dist_3D/propVelocity',
                        'synsPerConn':  cfg.synsperconn[cellModel],
                        'sec':          'spiny'}

#------------------------------------------------------------------------------
## E -> I
if cfg.addConn and cfg.EIGain > 0.0:
    binsLabel   = 'inh'
    preTypes    = excTypes
    postTypes   = inhTypes
    ESynMech    = ['AMPA','NMDA']

    for i,postType in enumerate(postTypes):
        for ipre, preBin in enumerate(bins[binsLabel]):
            for ipost, postBin in enumerate(bins[binsLabel]):
                ruleLabel = 'EI_'+str(i)+'_'+str(ipre)+'_'+str(ipost)+'_'+str(postType)
                netParams.connParams[ruleLabel] = {
                    'preConds':     {'cellType': preTypes, 'ynorm': list(preBin)},
                    'postConds':    {'cellType': postType, 'ynorm': list(postBin)},
                    'synMech':      ESynMech,
                    'probability':  pmat[('E', postType)][ipost,ipre],
                    'weight':       wmat[('E', postType)][ipost,ipre] * cfg.EIGain * cfg.EICellTypeGain[postType],
                    'synMechWeightFactor': cfg.synWeightFractionEI,
                    'delay':        'defaultDelay+dist_3D/propVelocity',
                    'sec':          'soma'} # simple I cells used right now only have soma

#------------------------------------------------------------------------------
## I -> E
if cfg.addConn and cfg.IEGain > 0.0:

    binsLabel   = 'inh'
    preTypes    = inhTypes
    synMechs    = [PVSynMech, SOMESynMech, VIPSynMech, NGFSynMech] 
    weightFactors = [[1.0], cfg.synWeightFractionSOME, [1.0], cfg.synWeightFractionNGF] # Update VIP and NGF syns! 
    secs        = ['perisom', 'apicdend', 'apicdend', 'apicdend']
    postTypes   = excTypes

    for ipreType, (preType, synMech, weightFactor, sec) in enumerate(zip(preTypes, synMechs, weightFactors, secs)):
        for ipostType, postType in enumerate(postTypes):
            for ipreBin, preBin in enumerate(bins[binsLabel]):
                for ipostBin, postBin in enumerate(bins[binsLabel]):
                    for cellModel in ['HH_reduced', 'HH_full']:
                        ruleLabel = preType+'_'+postType+'_'+cellModel+'_'+str(ipreBin)+'_'+str(ipostBin)
                        netParams.connParams[ruleLabel] = {
                            'preConds':     {'cellType': preType, 'ynorm': list(preBin)},
                            'postConds':    {'cellModel': cellModel, 'cellType': postType, 'ynorm': list(postBin)},
                            'synMech':      synMech,
                            'probability':  '%f * exp(-dist_3D_border/probLambda)' % (pmat[(preType, 'E')][ipostBin,ipreBin]),
                            'weight':       cfg.IEweights[ipostBin] * cfg.IEGain/ cfg.synsperconn[cellModel],
                            'synMechWeightFactor': weightFactor,
                            'synsPerConn':  cfg.synsperconn[cellModel],
                            'delay':        'defaultDelay+dist_3D/propVelocity',
                            'sec':          sec} # simple I cells used right now only have soma

#------------------------------------------------------------------------------
## I -> I
if cfg.addConn and cfg.IIGain > 0.0:

    binsLabel   = 'inh'
    preTypes    = inhTypes
    synMechs    =  [PVSynMech, SOMESynMech, VIPSynMech, NGFSynMech]   
    sec         = 'perisom'
    postTypes   = inhTypes

    for ipre, (preType, synMech) in enumerate(zip(preTypes, synMechs)):
        for ipost, postType in enumerate(postTypes):
            for iBin, bin in enumerate(bins[binsLabel]):
                for cellModel in ['HH_reduced']:
                    ruleLabel = preType+'_'+postType+'_'+str(iBin)
                    netParams.connParams[ruleLabel] = {
                        'preConds':     {'cellType': preType, 'ynorm': bin},
                        'postConds':    {'cellModel': cellModel, 'cellType': postType, 'ynorm': bin},
                        'synMech':      synMech,
                        'probability':  '%f * exp(-dist_3D_border/probLambda)' % (pmat[(preType, postType)]),
                        'weight':       cfg.IIweights[iBin] * cfg.IIGain / cfg.synsperconn[cellModel],
                        'synsPerConn':  cfg.synsperconn[cellModel],
                        'delay':        'defaultDelay+dist_3D/propVelocity',
                        'sec':          sec} # simple I cells used right now only have soma

#------------------------------------------------------------------------------
# Long-range  connectivity parameters
#------------------------------------------------------------------------------
if cfg.addLongConn:

    # load load experimentally based parameters for long range inputs
    cmatLong = connLongData['cmat']
    binsLong = connLongData['bins']

    longPops    = ['S1', 'S2', 'cM1', 'M2', 'OC']          # removed VL and POm because now they have separate connectivity
    # longPops = ['POm_sTC_m1', 'VL_sTC', 'S1', 'S2', 'cM1', 'M2', 'OC']
    cellTypes   = ['IT', 'PT', 'CT', 'PV', 'SOM', 'VIP', 'NGF']
    EorI        = ['exc', 'inh']
    syns        = {'exc': ESynMech, 'inh': 'GABAA'}
    synFracs    = {'exc': cfg.synWeightFractionEE, 'inh': [1.0]}

    for longPop in longPops:
        for ct in cellTypes:
            for EorI in ['exc', 'inh']:
                for i, (binRange, convergence) in enumerate(zip(binsLong[(longPop, ct)], cmatLong[(longPop, ct, EorI)])):
                    for cellModel in cellModels:
                        ruleLabel = longPop+'_'+ct+'_'+EorI+'_'+cellModel+'_'+str(i)
                        netParams.connParams[ruleLabel] = { 
                            'preConds': {'pop': longPop}, 
                            'postConds': {'cellModel': cellModel, 'cellType': ct, 'ynorm': list(binRange)},
                            'synMech': syns[EorI],
                            'convergence': convergence,
                            'weight': cfg.weightLong[longPop] / cfg.synsperconn[cellModel], 
                            'synMechWeightFactor': synFracs[EorI],
                            #'synMechWeightFactor': synFracs,
                            'delay': 'defaultDelay+dist_3D/propVelocity',
                            'synsPerConn': cfg.synsperconn[cellModel],
                            'sec': 'spiny'}

#------------------------------------------------------------------------------
# Subcellular connectivity (synaptic distributions)
#------------------------------------------------------------------------------         
if cfg.addSubConn:
    with open('../conn/conn_dend_PT.json', 'r') as fileObj: connDendPTData = json.load(fileObj)
    with open('../conn/conn_dend_IT.json', 'r') as fileObj: connDendITData = json.load(fileObj)
    
    #------------------------------------------------------------------------------
    # L2/3,VM_sTC_m1,S2,cM1,M2 -> PT (Suter, 2015)
    lenY    = 30 
    spacing = 50
    gridY   = range(0, -spacing*lenY, -spacing)
    synDens, _, fixedSomaY = connDendPTData['synDens'], connDendPTData['gridY'], connDendPTData['fixedSomaY']

    for k in synDens.keys():
        prePop,postType = k.split('_')  # eg. split 'M2_PT'
        if prePop == 'L2': prePop = 'IT2'  # include conns from layer 2/3 and 4
        netParams.subConnParams[k] = {
        'preConds':     {'pop': prePop}, 
        'postConds':    {'cellType': postType},  
        'sec':          'spiny',
        'groupSynMechs': ESynMech, 
        'density':      {'type': '1Dmap', 'gridX': None, 'gridY': gridY, 'gridValues': synDens[k], 'fixedSomaY': fixedSomaY}} 

    #------------------------------------------------------------------------------
    # POm_sTC_m1, VM_sTC_m1, M2, OC  -> E (L2/3, L5A, L5B, L6) (Hooks 2013)
    lenY    = 26
    spacing = 50
    gridY   = range(0, -spacing*lenY, -spacing)
    synDens, _, fixedSomaY = connDendITData['synDens'], connDendITData['gridY'], connDendITData['fixedSomaY']

    for k in synDens.keys():
        prePop,post     = k.split('_')  # eg. split 'M2_L2'
        postCellTypes   = ['IT','PT','CT'] if prePop in ['OC','POm_sTC_m1'] else ['IT','CT']  # only OC,POm_sTC_m1 include PT cells
        postyRange      = list(layer[post.split('L')[1]]) # get layer yfrac range 
        if post == 'L2': postyRange[1] = layer['4'][1]  # apply L2 rule also to L4 
        netParams.subConnParams[k] = {
        'preConds':     {'pop': prePop}, 
        'postConds':    {'ynorm': postyRange , 'cellType': postCellTypes},  
        'sec':          'spiny',
        'groupSynMechs': ESynMech, 
        'density':      {'type': '1Dmap', 'gridX': None, 'gridY': gridY, 'gridValues': synDens[k], 'fixedSomaY': fixedSomaY}} 

    #------------------------------------------------------------------------------
    # S1, S2, cM1 -> E IT/CT; no data, assume uniform over spiny
    netParams.subConnParams['S1,S2,cM1->IT,CT'] = {
        'preConds':     {'pop': ['S1','S2','cM1']}, 
        'postConds':    {'cellType': ['IT','CT']},
        'sec':          'spiny',
        'groupSynMechs': ESynMech, 
        'density':      'uniform'} 

    #------------------------------------------------------------------------------
    # rest of local E->E (exclude IT2->PT); uniform distribution over spiny
    netParams.subConnParams['IT2->non-PT'] = {
        'preConds':     {'pop': ['IT2']}, 
        'postConds':    {'cellType': ['IT','CT']},
        'sec':          'spiny',
        'groupSynMechs': ESynMech, 
        'density':      'uniform'} 
        
    netParams.subConnParams['non-IT2->E'] = {
        'preConds':     {'pop': ['IT4','IT5A','IT5B','PT5B','IT6','CT6']}, 
        'postConds':    {'cellType': ['IT','PT','CT']},
        'sec':          'spiny',
        'groupSynMechs': ESynMech, 
        'density':      'uniform'} 

    #------------------------------------------------------------------------------
    # PV->E; perisomatic (no sCRACM)
    netParams.subConnParams['PV->E'] = {
        'preConds':     {'cellType': 'PV'}, 
        'postConds':    {'cellType': ['IT', 'CT', 'PT']},  
        'sec':          'perisom', 
        'density':      'uniform'} 

    #------------------------------------------------------------------------------
    # SOM->E; apical dendrites (no sCRACM)
    netParams.subConnParams['SOM->E'] = {
        'preConds':     {'cellType': 'SOM'}, 
        'postConds':    {'cellType': ['IT', 'CT', 'PT']},  
        'sec':          'apicdend',
        'groupSynMechs': SOMESynMech,
        'density':      'uniform'} 

    #------------------------------------------------------------------------------
    # VIP->E; apical dendrites (no sCRACM)
    netParams.subConnParams['VIP->E'] = {
        'preConds':     {'cellType': 'VIP'}, 
        'postConds':    {'cellType': ['IT', 'CT', 'PT']},  
        'sec':          'apicdend',
        'groupSynMechs': VIPSynMech,
        'density':      'uniform'} 

    #------------------------------------------------------------------------------
    # NGF->E; apical dendrites (no sCRACM)
    ## Add the following level of detail?
    # -- L1 NGF -> L2/3+L5 tuft
    # -- L2/3 NGF -> L2/3+L5 distal apical
    # -- L5 NGF -> L5 prox apical
    netParams.subConnParams['NGF->E'] = {
        'preConds':     {'cellType': 'NGF'}, 
        'postConds':    {'cellType': ['IT', 'CT', 'PT']},  
        'sec':          'apicdend',
        'groupSynMechs': NGFSynMech,
        'density':      'uniform'} 

    #------------------------------------------------------------------------------
    # All->I; apical dendrites (no sCRACM)
    netParams.subConnParams['All->I'] = {
        'preConds':     {'cellType': ['IT', 'CT', 'PT'] + inhTypes},# + longPops}, 
        'postConds':    {'cellType': inhTypes},  
        'sec':          'spiny',
        'groupSynMechs': ESynMech,
        'density':      'uniform'} 

#------------------------------------------------------------------------------
# Description
#------------------------------------------------------------------------------
netParams.description = """ 
- M1 net, 6 layers, 7 cell types 
- based in  M1 v103
- v1_batch0: S1+M1+Th
"""
