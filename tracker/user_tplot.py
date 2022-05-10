"""
Plot results of a particle tracking experiment.
"""
# Import packages - add mpl.set_backend() right after plt.close()?
import matplotlib
from matplotlib import pyplot as plt

import xarray as xr
import numpy as np

from lo_tools import Lfun
from lo_tools import plotting_functions as pfun
Ldir = Lfun.Lstart()
#


# Choose an experiment and release to plot.
in_dir0 = Ldir['LOo'] / 'tracks'
exp_name = Lfun.choose_item(in_dir0, tag='', exclude_tag='.csv',
    itext='** Choose experiment from list **', last=False)
rel1 = Lfun.choose_item(in_dir0 / exp_name, tag='.nc', exclude_tag='grid',
    itext='** Choose item from list **', last=False)

rel2 = Lfun.choose_item(in_dir0 / exp_name, tag='.nc', exclude_tag='grid',
    itext='** Choose item from list **', last=False)

rel3 = Lfun.choose_item(in_dir0 / exp_name, tag='.nc', exclude_tag='grid',
    itext='** Choose item from list **', last=False)

in_dir02 = Ldir['LOo'] / 'tracks'
exp_name2 = Lfun.choose_item(in_dir02, tag='', exclude_tag='.csv',
    itext='** Choose experiment from list **', last=False)
rel18 = Lfun.choose_item(in_dir02 / exp_name2, tag='.nc', exclude_tag='grid',
    itext='** Choose item from list **', last=False)
rel28 = Lfun.choose_item(in_dir02 / exp_name2, tag='.nc', exclude_tag='grid',
    itext='** Choose item from list **', last=False)
rel38 = Lfun.choose_item(in_dir02 / exp_name2, tag='.nc', exclude_tag='grid',
    itext='** Choose item from list **', last=False)
#
dir = [in_dir0, in_dir02]
exp = [exp_name, exp_name]
rel0, rel8 = [rel1, rel2, rel3], [rel18, rel28, rel38]
rel = [rel0, rel8]
dsr, dsg = [[0]*2 for _ in range(3)], [[0]*2 for _ in range(3)]
n, m =2, 3
for i, j in n,m:
# get Datasets
    fn = dir[i] / exp[i] / rel[i][j]
    fng = dir[i] / exp[i] / 'grid.nc'
    dsr[i][j] = xr.open_dataset(fn, decode_times=False)
    dsg[i][j] = xr.open_dataset(fng)


for i, j in n, m:
    NT, NP = dsr.lon.shape
    NT2, NP2 = dsr2.lon.shape
#

# get a list of datetimes
ot_vec = dsr.ot.values
dt_list = [Lfun.modtime_to_datetime(ot) for ot in ot_vec]

ot_vec2 = dsr2.ot.values
dt_list2 = [Lfun.modtime_to_datetime(ot) for ot in ot_vec2]
#

# gather some fields, for convenience
lonp, latp = pfun.get_plon_plat(dsg.lon_rho.values, dsg.lat_rho.values)
hh = dsg.h.values
maskr = dsg.mask_rho.values

lonp2, latp2 = pfun.get_plon_plat(dsg2.lon_rho.values, dsg2.lat_rho.values)
hh2 = dsg2.h.values
maskr2 = dsg2.mask_rho.values
#

# subsample output for plotting
npmax = 1000 # max number of points to plot
step = max(1,int(np.floor(NP/npmax)))

step2 = max(1,int(np.floor(NP2/npmax)))

lon = dsr.lon.values[:,::step]
lat = dsr.lat.values[:,::step]

lon2 = dsr2.lon.values[:,::step2]
lat2 = dsr2.lat.values[:,::step2]
#

# make a mask that is False from the time a particle first leaves the domain
# and onwards
AA = [dsg.lon_rho.values[0,0], dsg.lon_rho.values[0,-1],
        dsg.lat_rho.values[0,0], dsg.lat_rho.values[-1,0]]
ib_mask = np.ones(lon.shape, dtype=bool)
ib_mask[lon < AA[0]] = False
ib_mask[lon > AA[1]] = False
ib_mask[lat < AA[2]] = False
ib_mask[lat > AA[3]] = False
NTS, NPS = lon.shape
for pp in range(NPS):
    tt = np.argwhere(ib_mask[:,pp]==False)
    if len(tt) > 0:
        ib_mask[tt[0][0]:, pp] = False

AA2 = [dsg2.lon_rho.values[0,0], dsg2.lon_rho.values[0,-1],
        dsg2.lat_rho.values[0,0], dsg2.lat_rho.values[-1,0]]
ib_mask2 = np.ones(lon2.shape, dtype=bool)
ib_mask2[lon2 < AA2[0]] = False
ib_mask2[lon2 > AA2[1]] = False
ib_mask2[lat2 < AA2[2]] = False
ib_mask2[lat2 > AA2[3]] = False
NTS2, NPS2 = lon2.shape
for pp in range(NPS2):
    tt2 = np.argwhere(ib_mask2[:,pp]==False)
    if len(tt2) > 0:
        ib_mask2[tt2[0][0]:, pp] = False
#

# and apply the mask to lon and lat
lon[~ib_mask] = np.nan
lat[~ib_mask] = np.nan

lon2[~ib_mask2] = np.nan
lat2[~ib_mask2] = np.nan
#

# PLOTTING
plt.close('all')
matplotlib.use('TkAgg')
pfun.start_plot(figsize=(10,5))
fig, axs = plt.subplots(1,2, dpi=100)
plt.tight_layout()
# MAP
# set domain limits
ok = [0,1]
grap = [121,122]
boun_lon = [lon, lon2]
boun_lat = [lat, lat2]
plon = [lonp, lonp2]
plat = [latp, latp2]
h = [hh, hh2]
title = ['Sinking Speed 8 m/day', 'No Sinking']
part = [lon[-1,:], lon2[-1,:]]
partlat = [lat[-1,:], lat2[-1,:]]

col = ['dodgerblue', 'forestgreen']
fig.suptitle('Fall Cruise 2-day Release: Final Particle Locations', fontsize=20)
for i in ok:
    pad = .02
    aa = [-123.3, -121.9, 47.1, 48.3]
    pfun.add_coast(axs[i])
    axs[i].axis(aa)
    pfun.dar(axs[i])
    axs[i].set_xlabel('Longitude')
    axs[i].set_ylabel('Latitude')
    axs[i].set_title(title[i])

    ate, = axs[i].plot(part[i], partlat[i], marker = 'o', lw =0, markerfacecolor = col[i], markeredgecolor = 'k', alpha=.25, markersize=4)


#fig, axs = plt.subplots(2, 2, layout='constrained')
#pc = axs[0, 0].pcolormesh(X, Y, Z, vmin=-1, vmax=1, cmap='RdBu_r')
#fig.colorbar(pc, ax=axs[0, 0])
#axs[0, 0].set_title('pcolormesh()')

    stat_num = [5, 6, 7, 9, 8, 1, 2, 3, 4]
    stat_lat = [47.6864,  47.6, 47.3937, 47.2828, 47.287, 47.8839, 47.925, 47.9835, 48.0024]
    stat_lon = [-122.423,  -122.365, -122.3601, -122.4437, -122.5387,  -122.3675, -122.4685, -122.6201, -122.6618]
    start, = axs[i].plot(lon[0,:], lat[0,:], 'r+', alpha=.05, markersize=5)
    stations, = axs[i].plot(stat_lon, stat_lat, marker = 'o', markerfacecolor = 'yellow', lw = 0, alpha=1, markersize=5, markeredgecolor = 'k')
    proxy1, = axs[i].plot([0], [0], marker = '+', markerfacecolor = 'red', lw = 0, alpha = 1, markersize = 5, label = 'Release location')
    proxy2, = axs[i].plot([0], [0], marker = 'o', markerfacecolor = 'yellow', lw = 0, alpha = 0.75, markersize = 5, markeredgecolor = 'k', label = 'Sampling stations')
    proxy4, = axs[i].plot([0], [0], marker = 'o', markerfacecolor = 'dodgerblue', lw = 0, alpha = 0.75, markersize = 4, label = '8 m/day sinking')
    proxy3, = axs[i].plot([0], [0], marker = 'o', markerfacecolor = 'forestgreen', lw = 0, alpha = 0.75, markersize = 4, label = 'No Sinking')
    prox = [proxy4, proxy3]
    labelz = ('Release location', 'Sampling stations', 'Final location (2d)')
    secnleg = axs[i].legend(handles = [proxy1, proxy2, prox[i]], labels = labelz, loc='upper right', fontsize = 'medium')

    for j, txt in enumerate(stat_num):
        axs[i].annotate(txt, (stat_lon[j], stat_lat[j]), textcoords="offset points", # how to position the text
        xytext=(0,10), size = 'small', # distance from text to points (x,y)
        ha='center', bbox=dict(fc='w',boxstyle='square',ec='k', color = 'k', alpha = 0.75))

plt.savefig('fall_twoday_no8.png')
plt.show()
pfun.end_plot()

# time series

plt.close('all')
matplotlib.use('TkAgg')
pfun.start_plot(figsize=(16,10))
fig = plt.figure()
fig.suptitle('Fall Cruise 2-day Release: Water Conditions over Time', fontsize=16)

td = (ot_vec - ot_vec[0])/86400
tv_list = ['z', 'salt', 'temp', 'oxygen']
tv_list2 = ['Depth (m)', 'Salt (PSU)','Temperature (°C)', 'Oxygen (mg/L)']
#tv_list = ['u', 'v', 'lon', 'lat']
ntv = len(tv_list)
for ii in range(ntv):
    tv = tv_list[ii]
    lbl = tv_list2[ii]
    NC = 1
    ax = fig.add_subplot(ntv,NC, (ii+1)*NC)
    v = dsr[tv].values[:,::step]
    v2 = dsr2[tv].values[:,::step2]
    v[~ib_mask] = np.nan
    v2[~ib_mask2] = np.nan
    ax.plot(td, v2, lw=.5, color ='forestgreen', alpha = .3, label = 'No Sinking')
    ax.plot(td, v, lw=.5, color = 'dodgerblue', linestyle ='dashed', alpha = .08, label = '8 m/day sinking')
    ax.text(.05, .05, lbl, fontweight='bold', transform=ax.transAxes)
    if ii == 0:
        handle1, = ax.plot([0], [0], lw=.5, color = 'forestgreen', linestyle = 'dashed', alpha = 0.9)
        handle2, = ax.plot([0], [0], lw=.5, color = 'dodgerblue', alpha = 1)
        ax.legend(handles = [handle2, handle1], labels = ('No Sinking', '8 m/day sinking'), bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left', ncol=2, mode="expand", borderaxespad=0.)
    if ii == ntv-1:
        ax.set_xlabel('Time (days)')

plt.savefig('fall_twoday_no8_var.png')
plt.show()
pfun.end_plot()


dsr.close()
dsg.close()
