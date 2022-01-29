"""
script to load sim and plot
"""

from netpyne import sim
from matplotlib import pyplot as plt
import os
import IPython as ipy
import pickle as pkl


thalamicpops = ['ss_RTN_o', 'ss_RTN_m', 'ss_RTN_i', 'VPL_sTC', 'VPM_sTC', 'POm_sTC_s1', 'VL_sTC', 'VM_sTC_m1', 'POm_sTC_m1', 'mt_RTN' ]

if __name__ == '__main__':

    dataType = 'spont' #'speech' #'spont'

    if dataType == 'spont':
        timeRange = [0, 2500]
        filenames = ['../../VMdata/v1_batch6/v1_batch6_%d_%d_data.pkl' % (iseed, cseed) for iseed in [0] for cseed in [0]] 


    allData = []

    for filename in filenames:
        sim.load(filename, instantiate=False)

        # standardd plots
        sim.analysis.plotRaster(**{'include': thalamicpops, 'saveFig': filename[:-8]+'_raster_TH', 'showFig': False, 'popRates': 'minimal', 'orderInverse': True, 'timeRange': [0,2500], 'figSize': (24,12), 'lw': 0.3, 'markerSize': 3, 'marker': '.', 'dpi': 300})
 

