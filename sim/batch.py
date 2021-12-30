"""
batch.py 

Batch simulation for M1+TH model using NetPyNE

Contributors: salvadordura@gmail.com, fernandodasilvaborges@gmail.com, joao
"""
from netpyne.batch import Batch
from netpyne import specs
import numpy as np

# ----------------------------------------------------------------------------------------------
# Custom
# ----------------------------------------------------------------------------------------------
def custom():
    params = specs.ODict()
    
    params['weightLong_thalM1']=[0.1,0.25,0.5,1.5]
    params['weightLong_M1thal']=[0.5]

    b = Batch(params=params, netParamsFile='netParams.py', cfgFile='cfg.py')

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
            'mpiCommand': 'mpirun -n 20 ./x86_64/special -mpi -python init.py', # --use-hwthread-cpus
            'skip': True}


    elif type=='hpc_slurm_gcp':
        b.runCfg = {'type': 'hpc_slurm', 
            'allocation': 'default',
            'walltime': '72:00:00', 
            'nodes': 1,
            'coresPerNode': 40,
            'email': 'fernandodasilvaborges@gmail.com',
            'folder': '/home/ext_fernandodasilvaborges_gmail_/S1_netpyne/sim/', 
            'script': 'init.py', 
            'mpiCommand': 'mpirun',
            'skipCustom': '_raster_gid.png'}

# ----------------------------------------------------------------------------------------------
# Main code
# ----------------------------------------------------------------------------------------------
if __name__ == '__main__': 
    b = custom() #

    simDate = '2021_12_30'
    simCode = simDate

    simType     = 'testingProjectionWeights_density100'

    b.date          = simDate
    b.dataFolder    = '../data/batch_sims_joao'
    b.batchLabel    = simDate+'fullModel_Th_M1_density'+'_'+simCode+'_'+simType
    b.saveFolder    = b.dataFolder + '/' +  b.batchLabel

    b.batchLabel = 'v100_batch0'  
    b.saveFolder = '../data/'+b.batchLabel
    b.method = 'grid'
    setRunCfg(b, 'mpi_direct')
    b.run() # run batch
