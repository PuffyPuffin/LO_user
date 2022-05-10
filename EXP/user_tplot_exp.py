
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

exp_name2 = Lfun.choose_item(in_dir0, tag='', exclude_tag='.csv',
    itext='** Choose experiment from list **', last=False)
rel18 = Lfun.choose_item(in_dir0 / exp_name2, tag='.nc', exclude_tag='grid',
    itext='** Choose item from list **', last=False)
rel28 = Lfun.choose_item(in_dir0 / exp_name2, tag='.nc', exclude_tag='grid',
    itext='** Choose item from list **', last=False)

exp_name3 = Lfun.choose_item(in_dir0, tag='', exclude_tag='.csv',
    itext='** Choose experiment from list **', last=False)
rel140 = Lfun.choose_item(in_dir0 / exp_name2, tag='.nc', exclude_tag='grid',
    itext='** Choose item from list **', last=False)
rel240 = Lfun.choose_item(in_dir0 / exp_name2, tag='.nc', exclude_tag='grid',
    itext='** Choose item from list **', last=False)
#
exp = [exp_name, exp_name2]
rel0, rel8, rel40 = [rel1, rel2], [rel18, rel28], [rel140, rel240]
rel = [rel0, rel8, rel40]
dsr, dsg = [[0]*2 for _ in range(3)], [[0]*1 for _ in range(3)]
NT, NP, step = [[0]*2 for _ in range(3)], [[0]*2 for _ in range(3)], [[0]*2 for _ in range(3)]
ot_vec = [[0]*2 for _ in range(3)]
lon, lat =[[0]*2 for _ in range(3)], [[0]*2 for _ in range(3)]
AA, NTS, NPS = [[0]*1 for _ in range(3)], [[0]*2 for _ in range(3)], [[0]*2 for _ in range(3)]
ib_mask = [[0]*2 for _ in range(3)]
proxy = [[0]*2 for _ in range(3)]

n = list(range(3))

for i in n:
    for j in n:
# get Datasets
        fn = in_dir0 / exp[i] / rel[i][j]

        fng = in_dir0 / exp[i] / 'grid.nc'
        dsr[i][j] = xr.open_dataset(fn, decode_times=False)
        dsg[i] = xr.open_dataset(fng)
        ot_vec[i][j] = dsr[i][j].ot.values
        NT[i][j], NP[i][j] = dsr[i][j].lon.shape
        lonp, latp = pfun.get_plon_plat(dsg[i].lon_rho.values, dsg[i].lat_rho.values)
        npmax = 1500 # max number of points to plot
        step = max(1,int(np.floor(NP[i][j]/npmax)))

        lon[i][j] = dsr[i][j].lon.values[:,::step]
        lat[i][j] = dsr[i][j].lat.values[:,::step]
        AA[i] = [dsg[i].lon_rho.values[0,0], dsg[i].lon_rho.values[0,-1],
                dsg[i].lat_rho.values[0,0], dsg[i].lat_rho.values[-1,0]]
        ib_mask[i][j] = np.ones(lon[i][j].shape, dtype=bool)
        ib_mask[i][j][lon[i][j] < AA[i][0]] = False
        ib_mask[i][j][lon[i][j] > AA[i][1]] = False
        ib_mask[i][j][lat[i][j] < AA[i][2]] = False
        ib_mask[i][j][lat[i][j] > AA[i][3]] = False
        NTS[i][j], NPS[i][j] = lon[i][j].shape

        for pp in range(NPS[i][j]):
            tt = np.argwhere(ib_mask[i][j][:,pp]==False)
            if len(tt) > 0:
                ib_mask[i][j][tt[0][0]:, pp] = False
        lon[i][j][~ib_mask[i][j]] = np.nan
        lat[i][j][~ib_mask[i][j]] = np.nan
plt.close('all')
matplotlib.use('TkAgg')
pfun.start_plot(figsize=(13,7))
fig, axs = plt.subplots(1,2, sharey=True)

title = ['No Sinking', 'Sinking 8 m/day' 'Sinking 40 m/day']
col = [["yellow", "indigo"], ["orange", 'b']]
pad = .02
aa = [-123.5, -121.9, 47.1, 48.3]
x = [-123.3, -123.0, -122.7, -122.4, -122.1]

size = [3, 2]
alph = [0.75, 0.5]
mek = [["orangered", "mediumslateblue"], ["r", "dodgerblue"]]
mk = ['o', '*']
stat_num = [5, 6, 7, 9, 8, 1, 2, 3, 4]
stat_lat = [47.6864,  47.6, 47.3937, 47.2828, 47.287, 47.8839, 47.925, 47.9835, 48.0024]
stat_lon = [-122.423,  -122.365, -122.3601, -122.4437, -122.5387,  -122.3675, -122.4685, -122.6201, -122.6618]
prox = [[0]*3 for _ in range(2)]
labelz = ('Release location', 'Sampling stations', '2-Day release', '15-Day release')

tv_list = ['z', 'salt', 'temp', 'oxygen']
tv_list2 = ['Depth (m)', 'Salt (PSU)','Temperature (°C)', 'Oxygen (mg/L)']
#tv_list = ['u', 'v', 'lon', 'lat']
ntv = len(tv_list)
fig, axs = plt.subplots(ntv,len(n), sharex=True)
fig.suptitle('Summer Cruise Water Conditions', fontsize=16)
td = [[0]*2 for _ in range(1)]
v = [[0]*2 for _ in range(2)]
for i in m:
    for j in n:
        td[i][j] = (ot_vec[i][j] - ot_vec[j][j][0])/86400
        for ii in range(ntv):
            tv = tv_list[ii]
            lbl = tv_list2[ii]
            NC = 1
            v[i][j] = dsr[i][j][tv].values[:,::step]
            v[i][j][~ib_mask[i][j]] = np.nan
            axs[0][j].plot(td[0][ii], v[0][ii], lw=.5, color = col[0][j], linestyle ='dashed', alpha = .08)
            axs[i][j].text(.05, .05, lbl[i][j], fontweight='bold', transform=subfigure.transSubfigure)
            if ii == 0:
                handle3, = axs[i][j].plot([0], [0], color = 'peachpuff', lw = 0.5, linestyle ='dashed', alpha = .08)
                handle4, = axs[i][j].plot([0], [0], color = 'coral', lw =0.5, linestyle ='dashed', alpha = .08)
                handle5, = axs[i][j].plot([0], [0], color = 'darkred', lw = 0.5, linestyle ='dashed', alpha = .08)
                handle6, = axs[i][j].plot([0], [0], color = 'lavenderblush', lw = 0.5, linestyle ='dashed', alpha = .08)
                handle7, = axs[i][j].plot([0], [0], color = 'plum', lw = 0.5, linestyle ='dashed', alpha = .08)
                handle8, = axs[i][j].plot([0], [0], color = 'darkviolet', lw = 0.5, linestyle ='dashed', alpha = .08)
                handz = [[handle3, handle4, handle5], [handle6, handle7, handle8]]
#                axs[i][j].legend(handles = handz[i][j], labels = ('15 Day', '5 Day', '2 Day'), bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left', ncol=2, mode="expand", borderaxespad=0.)
            if ii == ntv-1:
                axs[i][j].set_xlabel('Time (days)')
