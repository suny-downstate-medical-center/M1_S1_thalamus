"""
script to load matrix
"""

from netpyne import sim
from matplotlib import pyplot as plt
import os
import IPython as ipy
import pickle

filename = '../data/v1_batch4/v1_batch4_0_0_data.pkl'
sim.load(filename, instantiate=False)

