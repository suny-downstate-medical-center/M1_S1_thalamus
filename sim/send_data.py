#!/usr/bin/python3

import sys
import os

    
#     '#--------------------------------------------------------------------------------------------------# \n')
# sys.stderr.write('#                                                                                                  # \n')
# sys.stderr.write('#                                      Invalid input argument                                      # \n')
# sys.stderr.write('#                                                                                                  # \n')
# sys.stderr.write('#
# input

sys.stderr.write('\n\n###   Insert path to file(s):   ###\n')
sys.stderr.write('      examples\n\n')
sys.stderr.write('M1/data/\n')
sys.stderr.write('M1/data/init_sims_joao\n')
sys.stderr.write('M1/data/batch_sims_joao\n')
# sys.stderr.write('insert path to file(s)\n')
sys.stderr.write('\nfilepath:')

projectPath_vm='/home/joao/Research/Models/NetPyNE/'
projectPath_local='/Users/joao/Research/Models/NetPyNE/M1/'
filePath=input()
filePath_vm = projectPath_vm+filePath
filePath_local = projectPath_local+filePath


sys.stderr.write('\n\n###   Insert path to bucket folder:   ###\n')
sys.stderr.write('      examples\n\n')
sys.stderr.write('gs://joao-bucket\n')
sys.stderr.write('gs://joao-bucket/thal_M1_joao_data\n')
sys.stderr.write('\nfilepath:')
bucketPath = input()
  
# # printing the sum in integer
# print(num1 + num2)



print('\n\ngsutil command: \n\n')
print('submit data:       gsutil -m cp -r '+filePath_vm+' '+bucketPath+'\n\n')
print('retrieve data:     gsutil -m cp -r '+bucketPath+' '+filePath_local+'\n\n')
