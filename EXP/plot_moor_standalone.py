"""
Stand-alone code to plot a user-specified mooring extraction.
"""
from pathlib import Path
moor_fn = Path('/mnt/c/Users/kudou/Documents/PostBack_Work/LO'
    +'/WestPoint_2020.01.01_2020.12.31.nc')

import xarray as xr
import matplotlib
matplotlib.use('WebAgg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# load everything using xarray
xs = xr.load_dataset(moor_fn)
ot = xs.ocean_time.values
ot_dt = pd.to_datetime(ot)
t = (ot_dt - ot_dt[0]).total_seconds().to_numpy()
T = t/86400 # time in days from start
print('time step of mooring'.center(60,'-'))
print(t[1])
print('time limits'.center(60,'-'))
print('start ' + str(ot_dt[0]))
print('end   ' + str(ot_dt[-1]))
print('info'.center(60,'-'))
VN_list = []
for vn in xs.data_vars:
    print('%s %s' % (vn, xs[vn].shape))
    VN_list.append(vn)
    
# populate lists of variables to plot
vn2_list = ['zeta']
if 'Pair' in VN_list:
    vn2_list += ['shflux', 'swrad']
vn3_list = []
if 'salt' in VN_list:
    vn3_list += ['salt', 'temp']
if 'NO3' in VN_list:
    vn3_list += ['oxygen']

# plot time series using a pandas DataFrame
df = pd.DataFrame(index=ot)
for vn in vn2_list:
    df[vn] = xs[vn].values
for vn in vn3_list:
    # the -1 means surface values
    df[vn] = xs[vn][:, -1].values

plt.close('all')
dfx = df.plot(subplots=True, figsize=(16,10))
x = [3, 4]
for i in vn3_list:
    if i == 'salt':
        aa = np.mean(df[i].values)
        dfx[x[0]].axhline(y=aa, color = 'black', linestyle = '--')
    elif i == 'temp':
        aa = np.mean(df[i].values)
        dfx[x[1]].axhline(y=aa, color = 'black', linestyle = '--')        
    else:
        continue

plt.show()