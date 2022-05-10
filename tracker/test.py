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

#rel2 = Lfun.choose_item(in_dir0 / exp_name, tag='.nc', exclude_tag='grid',
#    itext='** Choose item from list **', last=False)

rel3 = Lfun.choose_item(in_dir0 / exp_name, tag='.nc', exclude_tag='grid',
    itext='** Choose item from list **', last=False)

exp_name2 = Lfun.choose_item(in_dir0, tag='', exclude_tag='.csv',
    itext='** Choose experiment from list **', last=False)
rel18 = Lfun.choose_item(in_dir0 / exp_name2, tag='.nc', exclude_tag='grid',
    itext='** Choose item from list **', last=False)
#rel28 = Lfun.choose_item(in_dir0 / exp_name2, tag='.nc', exclude_tag='grid',
#    itext='** Choose item from list **', last=False)
rel38 = Lfun.choose_item(in_dir0 / exp_name2, tag='.nc', exclude_tag='grid',
    itext='** Choose item from list **', last=False)
#
exp = [exp_name, exp_name2]
rel0, rel8 = [rel1, rel3], [rel18, rel38]
rel = [rel0, rel8]
dsr, dsg = [[0]*2 for _ in range(2)], [[0]*1 for _ in range(2)]
NT, NP, step = [[0]*2 for _ in range(2)], [[0]*2 for _ in range(2)], [[0]*2 for _ in range(2)]
ot_vec = [[0]*2 for _ in range(2)]
lon, lat =[[0]*2 for _ in range(2)], [[0]*2 for _ in range(2)]
AA, NTS, NPS = [[0]*1 for _ in range(2)], [[0]*2 for _ in range(2)], [[0]*2 for _ in range(2)]
ib_mask = [[0]*2 for _ in range(2)]
proxy = [[0]*2 for _ in range(2)]

n = list(range(2))

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

title = ['No Sinking', 'Sinking 8 m/day']
col = [["yellow", "indigo"], ["orange", 'b']]
fig.suptitle('Fall Cruise Final Particle Locations', fontsize=20)
for i in n:
    for j in n:
        pad = .02
        aa = [-123.5, -121.9, 47.1, 48.3]
        x = [-123.3, -123.0, -122.7, -122.4, -122.1]
        axs[i].set_xticks(x)

        pfun.add_coast(axs[i])
        axs[i].axis(aa)
        pfun.dar(axs[i])
        axs[i].set_xlabel('Longitude')
        axs[0].set_ylabel('Latitude')

        axs[i].set_title(title[i])
        size = [3, 2]
        alph = [0.75, 0.5]
        mek = [["orangered", "mediumslateblue"], ["r", "dodgerblue"]]
        mk = ['o', '*']
        ate, = axs[i].plot(lon[i][j][-1,:], lat[i][j][-1,:], marker = 'o', lw =0, markerfacecolor = col[i][j], markeredgecolor = mek[i][j], alpha= 0.5, markersize= size[j])


#fig, axs = plt.subplots(2, 2, layout='constrained')
#pc = axs[0, 0].pcolormesh(X, Y, Z, vmin=-1, vmax=1, cmap='RdBu_r')
#fig.colorbar(pc, ax=axs[0, 0])
#axs[0, 0].set_title('pcolormesh()')

        stat_num = [5, 6, 7, 9, 8, 1, 2, 3, 4]
        stat_lat = [47.6864,  47.6, 47.3937, 47.2828, 47.287, 47.8839, 47.925, 47.9835, 48.0024]
        stat_lon = [-122.423,  -122.365, -122.3601, -122.4437, -122.5387,  -122.3675, -122.4685, -122.6201, -122.6618]
        prox = [[0]*3 for _ in range(2)]
        start, = axs[i].plot(lon[1][1][0,:], lat[1][1][0,:], marker= 'd', mfc= 'w', mec="limegreen", alpha=.1, markersize=3)
        stations, = axs[i].plot(stat_lon, stat_lat, marker = 'o', mfc= 'yellow', lw = 0, alpha=1, markersize=5, mec = 'k')
        proxy1, = axs[i].plot([0], [0], marker = 'd', mfc= 'w', mec="limegreen", lw = 0, alpha = 1, markersize = 3)
        proxy2, = axs[i].plot([0], [0], marker = 'o', mfc = 'yellow', lw = 0, alpha = 1, markersize = 5, mec = 'k')
        proxy5, = axs[i].plot([0], [0], marker = 'o', mfc = 'yellow', lw = 0, alpha = 1, markersize = 3, mec = "orangered")
#        proxy4, = axs[i].plot([0], [0], marker = 'o', markerfacecolor = 'coral', lw = 0, alpha = 0.75, markersize = 4, markeredgecolor = 'k')
        proxy3, = axs[i].plot([0], [0], marker = 'o', mfc = 'indigo', lw = 0, alpha = 1, markersize = 2, mec = "mediumslateblue")
        proxy8, = axs[i].plot([0], [0], marker = 'o', mfc = 'orange', lw = 0, alpha = 1, markersize = 3, mec = "r")
#        proxy7, = axs[i].plot([0], [0], marker = 'o', markerfacecolor = 'plum', lw = 0, alpha = 0.75, markersize = 4, markeredgecolor = 'k')
        proxy6, = axs[i].plot([0], [0], marker = 'o', mfc ='b' , lw = 0, alpha = 1, markersize = 2, mec = "dodgerblue")

        labelz = ('Release location', 'Sampling stations', '2-Day release', '15-Day release')
        secnleg = axs[0].legend(handles = [proxy1, proxy2, proxy3, proxy5], labels = labelz, loc='lower left', fontsize = 'medium')
        secnleg2 = axs[1].legend(handles = [proxy1, proxy2, proxy6, proxy8], labels = labelz, loc='lower left', fontsize = 'medium')

        for k, txt in enumerate(stat_num):
            axs[i].annotate(txt, (stat_lon[k], stat_lat[k]), textcoords="offset points", # how to position the text
            xytext=(0,10), size = 'x-small', # distance from text to points (x,y)
            ha='center', bbox=dict(fc='w',boxstyle='circle',ec='none', color = 'k', alpha = 0.5, pad = 0.2))
plt.savefig('fall_08_map.png')
plt.show()
pfun.end_plot()
