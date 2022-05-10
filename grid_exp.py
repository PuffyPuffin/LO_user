
import numpy as np
import sys, os
import xarray as xr
import argparse
from pathlib import Path as pth

path0 = pth.cwd().parent.parent / 'LO_data' / 'grids' / 'cas6' / 'grid.nc'
path1 = pth.cwd().parent / 'LO_user'
paths = [path0, path1]
for i in paths:
    if str(paths) not in sys.path:
        sys.path.append(str(paths))

from lo_tools import Lfun, plotting_functions
# Import cas6 grid for finding coast
ds = xr.open_dataset(path0, decode_times=False)
parser = argparse.ArgumentParser()

parser.add_argument('latbounds', default=[47.217061, 47.305420], nargs ='+')
parser.add_argument('lonbounds', default=[-122.58, -122.331], nargs ='+')
args = parser.parse_args()
TR = args.__dict__
slat = TR['lonbounds']
slon = TR['latbounds']
def get_boundaries(slat, slon):
    # Set area to search for coastline (at least 1 grid cell larger than desired coast)
    slat = latbounds
    slon = lonbounds
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

    lonin = np.argmin(~x)

    # Boolean of lon points less than upper bound

    # Store index of last value less than upper bound
    lonax = np.argmax(~y)

    # Restrict mask_rho (binary water/coast) using indices from above
    maskr = maskr[latin:latax, lonin:lonax]



    #Print restricted maskr to comapre to later paper matrix
    # Restrict lat and lon using indices from above
    # plaid grid, take one row/column to prevent repeated values
    lat = lats[latin:latax, 0]
    lon = lons[0, lonin:lonax]
    # Create empty list to store coastal water lat indices
    il=[]
    # Create empty list to store coastal water lon indices
    jl=[]
    # Create empty matrix to show coastal values
    paper=np.zeros_like(maskr)
    # For loop finds all adjacent values to water values and stores a water value
    # if that value is next to any coast
    # 1=water, 0=land
    # Make sure for loop does not inclue first and last columns and rows, as they will
    # compare to non-adjacent values (index 0 -1 = index 14)
    for i in range(1,len(lat)-1):
        for j in range(1,len(lon)-1):
            if maskr[i, j]==1 and maskr[i-1, j]==0:
                paper[i,j]=2
                il.append(i)
                jl.append(j)
            elif maskr[i, j]==1 and maskr[i+1, j]==0:
                paper[i,j]=2
                il.append(i)
                jl.append(j)
            elif maskr[i, j]==1 and maskr[i, j+1]==0:
                paper[i,j]=2
                il.append(i)
                jl.append(j)
            elif maskr[i, j]==1 and maskr[i, j-1]==0:
                paper[i,j]=2
                il.append(i)
                jl.append(j)
            elif maskr[i, j]==1 and maskr[i-1, j-1]==0:
                paper[i,j]=2
                il.append(i)
                jl.append(j)
            elif maskr[i, j]==1 and maskr[i-1, j+1]==0:
                paper[i,j]=2
                il.append(i)
                jl.append(j)
            elif maskr[i, j]==1 and maskr[i+1, j-1]==0:
                paper[i,j]=2
                il.append(i)
                jl.append(j)
            elif maskr[i, j]==1 and maskr[i+1, j+1]==0:
                paper[i,j]=2
                il.append(i)
                jl.append(j)
            else:
                continue
            # Print matrix that shows ID'ed coast water
            # Create vectors of coastal waters within bounds
    np.savetxt('paper', paper)
    latf = lat[il]
    lonf = lon[jl]
    lonvec = lonf
    latvec = latf
    pcs_vec = np.zeros_like(lonf)
    ds.close()
    # plon00, plat00, pcs00 = ic_from_list(lonvec, latvec, pcs_vec)
    return latf, lonf
