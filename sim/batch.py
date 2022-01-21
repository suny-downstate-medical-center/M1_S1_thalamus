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
    
    params['weightLong_thalM1']=[0.25,0.5,0.75,1.0,1.5]
    params['weightLong_M1thal']=[0.25,0.5,0.75,1.0,1.5]

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


    elif type=='mpi_direct2':
        b.runCfg = {'type': 'mpi_direct',
            'mpiCommand': 'mpirun -n 80 ./x86_64/special -mpi -python init.py', # --use-hwthread-cpus
            'skip': True}


    elif type=='mpi_direct':
        b.runCfg = {'type': 'mpi_direct',
            'cores': 40,
            'script': 'init.py',
            'mpiCommand': 'mpiexec',
            'skip': True}

    elif type=='hpc_slurm_gcp':
        b.runCfg = {'type': 'hpc_slurm', 
            'allocation': 'default',
            'walltime': '72:00:00', 
            'nodes': 1,
            'coresPerNode': 80,
            'email': 'fernandodasilvaborges@gmail.com',
            'folder': '/home/ext_fernandodasilvaborges_gmail_/ecas_model/sim/', 
            'script': 'init.py', 
            'mpiCommand': 'mpirun',
            'skipCustom': '_raster.png'}

# ----------------------------------------------------------------------------------------------
# Main code
# ----------------------------------------------------------------------------------------------
if __name__ == '__main__': 
    b = custom() #

    b.batchLabel = 'v0_batch2'  
    b.saveFolder = '../data/'+b.batchLabel
    b.method = 'grid'
    setRunCfg(b, 'hpc_slurm_gcp')
    b.run() # run batch
