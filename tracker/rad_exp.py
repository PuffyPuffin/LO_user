"""
Plot results of a particle tracking experiment.
"""
# Import packages - add mpl.set_backend() right after plt.close()?
import matplotlib
from matplotlib import pyplot as plt

import xarray as xr
import numpy as np
import csv

from lo_tools import Lfun
from lo_tools import plotting_functions as pfun
Ldir = Lfun.Lstart()
#

c = open('copper_RAD.csv','r')

data = csv.reader(c)
