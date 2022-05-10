
import numpy as np
import sys, os
import xarray as xr

sys.path.append(os.path.abspath('../../LO/lo_tools/lo_tools'))

from lo_tools import Lfun, plotting_functions
# Import cas6 grid for finding coast
ds = xr.open_dataset('../LO_data/grids/cas6/grid.nc', decode_times=False)

bin = [0, 1]
# Set area to search for coastline (at least 1 grid cell larger than desired coast)
latbounds = [47.217061, 47.305420]
lonbounds = [-122.5, -122.3]
# lon_psi and lat_psi as the corners for pcolormesh
latp = ds.lat_psi.values
lonp = ds.lon_psi.values
# lon_rho and lat_rho (grid cell centers) for release points
lats = ds.lat_rho.values
lons = ds.lon_rho.values
maskr = ds.mask_rho.values
# h field for plotting
h = ds.h.values
h[maskr==0] = np.nan

# Boolean of lat points greater than lower bound

x = lats[:,0]>=latbounds[0]

y = lats[:,0]<=latbounds[1]
latin = np.argmin(~x)

latax = np.argmax(~y)

# Store index of first value greater than lower bound
x = lons[0,:]>=lonbounds[0]
y = lons[0,:]<=lonbounds[1]
print(lonbounds)
lonin = np.argmin(~x)
print(lonin)
# Boolean of lon points less than upper bound

# Store index of last value less than upper bound
lonax = np.argmax(~y)
print(lonax)
# Restrict mask_rho (binary water/coast) using indices from above
maskr = maskr[latin:latax, lonin:lonax]

print(maskr)
