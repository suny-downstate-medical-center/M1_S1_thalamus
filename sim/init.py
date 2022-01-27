"""
init.py

Starting script to run NetPyNE-basedS1 model.

Usage:
    python init.py # Run simulation, optionally plot a raster

MPI usage:
    mpiexec -n 4 nrniv -python -mpi init.py

Contributors: salvadordura@gmail.com, fernandodasilvaborges@gmail.com, joaovvitor@gmail.com?
"""

import matplotlib; matplotlib.use('Agg')  # to avoid graphics error in servers
from netpyne import sim

# cfg, netParams = sim.readCmdLineArgs()

cfg, netParams = sim.readCmdLineArgs(simConfigDefault='cfg.py', netParamsDefault='netParams.py')

sim.initialize(
    simConfig = cfg, 	
    netParams = netParams)  				# create network object and set cfg and net params
sim.net.createPops()               			# instantiate network populations
sim.net.createCells()              			# instantiate network cells based on defined populations
sim.net.connectCells()            			# create connections between cells based on params
sim.net.addStims() 							# add network stimulation
sim.setupRecording()              			# setup variables to record for each cell (spikes, V traces, etc)
sim.runSim()                      			# run parallel Neuron simulation  
sim.gatherData()                  			# gather spiking data and cell info from each node
sim.saveData()                    			# save params, cell info and sim output to file (pickle,mat,txt,etc)#
sim.analysis.plotData()         			# plot spike raster etc

allpopsplotConn = [	'NGF1', 'IT2', 'SOM2', 'PV2', 'VIP2', 'NGF2', 
				'IT4', 'SOM4', 'PV4', 'VIP4', 'NGF4', 'IT5A', 
				'SOM5A', 'PV5A', 'VIP5A', 'NGF5A', 'IT5B', 
				'PT5B', 'SOM5B', 'PV5B', 'VIP5B', 'NGF5B', 
				'IT6', 'CT6', 'SOM6', 'PV6', 'VIP6', 'NGF6',	
				'VPL_sTC', 	'VPM_sTC', 		'POm_sTC_s1', 
				'VL_sTC', 	'VM_sTC_m1', 	'POm_sTC_m1',
				'mt_RTN', 	'ss_RTN_o', 	'ss_RTN_m', 	'ss_RTN_i',
				'L1_HAC_cNA', 'L23_LBC_dNA', 'L23_MC_cAC', 'L4_LBC_dNA', 'L4_MC_cAC', 'L5_LBC_dST', 'L5_MC_bAC', 'L6_LBC_bIR', 'L6_MC_bIR',
				'L23_PC_cAD', 'L4_PC_cAD', 'L4_SS_cAD', 'L4_SP_cAD', 
             	'L5_TTPC1_cAD', 'L5_TTPC2_cAD', 'L5_STPC_cAD', 'L5_UTPC_cAD',
             	'L6_TPC_L1_cAD', 'L6_TPC_L4_cAD', 'L6_BPC_cAD', 'L6_IPC_cAD', 'L6_UTPC_cAD']

features = ['numConns','convergence']
groups =['pop']
for feat in features:
   for group in groups:
       sim.analysis.plotConn(includePre=allpopsplotConn, includePost=allpopsplotConn, feature=feat, groupBy=group, figSize=(24,24), saveFig=True, orderBy='gid', graphType='matrix', fontSize=18, saveData='../data/'+cfg.simLabel[0:9]+'/'+cfg.simLabel + '_' + group + '_' + feat+ '_matrix.json')