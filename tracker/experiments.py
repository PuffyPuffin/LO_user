"""
This is where you set the run "gtagex" and the initial condition
based on an experiment name passed by the calling code.

"""

import numpy as np

def get_exp_info(exp_name):

    # Defaults
    if exp_name == 'ai0':
        gridname = 'ai0'; tag = 'v0'; ex_name = 'n0k'
    else:
        gridname = 'cas6'; tag = 'v3'; ex_name = 'lo8b'

    EI = {}
    EI['exp_name'] = exp_name # tracker experiment name
    # ROMS names
    EI['gridname'] = gridname
    EI['tag'] = tag
    EI['ex_name'] = ex_name
    EI['gtagex'] = gridname + '_' + tag + '_' + ex_name
    return EI

def get_ic(EI, fn00):
    # routines to set particle initial locations, all numpy arrays

    # NOTE: "pcs" refers to fractional depth, and goes linearly from -1 to 0
    # between the local bottom and free surface.  It is how we keep track of
    # vertical position, only converting to z-position when needed.

    exp_name = EI['exp_name']

    if exp_name == 'jdf0': # Mid-Juan de Fuca
        lonvec = np.linspace(-123.85, -123.6, 20)
        latvec = np.linspace(48.2, 48.4, 20)
        pcs_vec = np.array([0])
        plon00, plat00, pcs00 = ic_from_meshgrid(lonvec, latvec, pcs_vec)

    if exp_name == 'ai0': # Mid-Admiralty Inlet
        lonvec = np.array([-122.6])
        latvec = np.array([48])
        # These are: (Slope off JdF, Middle of JdF, Whidbey Basin)
        pcs_vec = np.linspace(-1,-0.9,num=1000)
        plon00, plat00, pcs00 = ic_from_list(lonvec, latvec, pcs_vec)

    elif exp_name == 'vmix': # three vertical profiles to test mixing
        # use with the new flag: -no_advection True, so a full command would be
        # python tracker.py -exp vmix -3d True -clb True -no_advection True
        lonvec = np.array([-125.35, -124.0, -122.581])
        latvec = np.array([47.847, 48.3, 48.244])
        # These are: (Slope off JdF, Middle of JdF, Whidbey Basin)
        pcs_vec = np.linspace(-1,0,num=4000)
        plon00, plat00, pcs00 = ic_from_list(lonvec, latvec, pcs_vec)
    elif exp_name == 'elb': # Elliot Bay
        # Import Parker's packages for mapping
        sys.path.append(os.path.abspath('../../LO/lo_tools/lo_tools'))
        import Lfun, plotting_functions

        # Import cas6 grid for finding coast
        ds = nc.Dataset('../../LiveOcean_data/grids/cas6/grid.nc')

        # Set area to search for coastline (at least 1 grid cell larger than desired coast)
        latbounds = [47.581, 47.65]
        lonbounds = [-122.441, -122.331]

        # lon_psi and lat_psi as the corners for pcolormesh
        latp = ds['lat_psi'][:].data
        lonp = ds['lon_psi'][:].data
        # lon_rho and lat_rho (grid cell centers) for release points
        lats = ds['lat_rho'][:].data
        lons = ds['lon_rho'][:].data
        maskr = ds['mask_rho'][:].data
        # h field for plotting
        h = ds['h'][:].data
        h[maskr==0] = np.nan

        # Boolean of lat points greater than lower bound
        rlat = lats[:,0]>=latbounds[0]
        # Store index of first value greater than lower bound
        latin = np.argmin(~rlat)

        # Boolean of lat points less than upper bound
        rlat = lats[:,0]<=latbounds[1]
        # Store index of last value less than upper bound
        latax = np.argmax(~rlat)

        # Boolean of lon points greater than lower bound
        rlon = lons[0,:]>=lonbounds[0]
        # Store index of first value greater than lower bound
        lonin = np.argmin(~rlon)

        # Boolean of lon points less than upper bound
        rlon = lons[0,:]<=lonbounds[1]
        # Store index of last value less than upper bound
        lonax = np.argmax(~rlon)

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
        for i in range(1,14):
            for j in range(1,15):
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
        latf = lat[il]
        lonf = lon[jl]
        lonvec = lonf
        latvec = latf
        pcs_vec = np.zeros_like(lonf)
        ds.close()
        plon00, plat00, pcs00 = ic_from_list(lonvec, latvec, pcs_vec)
    return plon00, plat00, pcs00

def ic_from_meshgrid(lonvec, latvec, pcs_vec):
    # First create three vectors of initial locations (as done in some cases above).
    # plat00 and plon00 should be the same length, and the length of pcs00 is
    # as many vertical positions you have at each lat, lon
    # (expressed as fraction of depth -1 < pcs < 0).
    # Then we create full output vectors (each has one value per point).
    # This code takes each lat, lon location and then assigns it to NSP points
    # corresponding to the vector of pcs values.
    lonmat, latmat = np.meshgrid(lonvec, latvec)
    plon_vec = lonmat.flatten()
    plat_vec = latmat.flatten()
    if len(plon_vec) != len(plat_vec):
        print('WARNING: Problem with length of initial lat, lon vectors')
    NSP = len(pcs_vec)
    NXYP = len(plon_vec)
    plon_arr = plon_vec.reshape(NXYP,1) * np.ones((NXYP,NSP))
    plat_arr = plat_vec.reshape(NXYP,1) * np.ones((NXYP,NSP))
    pcs_arr = np.ones((NXYP,NSP)) * pcs_vec.reshape(1,NSP)
    plon00 = plon_arr.flatten()
    plat00 = plat_arr.flatten()
    pcs00 = pcs_arr.flatten()
    return plon00, plat00, pcs00

def ic_from_list(lonvec, latvec, pcs_vec):
    # Like ic_from_meshgrid() but treats the lon, lat lists like lists of mooring locations.
    plon_vec = lonvec
    plat_vec = latvec
    if len(plon_vec) != len(plat_vec):
        print('WARNING: Problem with length of initial lat, lon lists')
    NSP = len(pcs_vec)
    NXYP = len(plon_vec)
    plon_arr = plon_vec.reshape(NXYP,1) * np.ones((NXYP,NSP))
    plat_arr = plat_vec.reshape(NXYP,1) * np.ones((NXYP,NSP))
    pcs_arr = np.ones((NXYP,NSP)) * pcs_vec.reshape(1,NSP)
    plon00 = plon_arr.flatten()
    plat00 = plat_arr.flatten()
    pcs00 = pcs_arr.flatten()
    return plon00, plat00, pcs00

def ic_from_TEFsegs(fn00, gridname, seg_list, DZ, NPmax=10000):
    import pickle
    import sys
    # select the indir
    from lo_tools import Lfun, zrfun
    Ldir = Lfun.Lstart()
    indir = Ldir['LOo'] / 'tef' / ('volumes_' + gridname)
    # load data
    j_dict = pickle.load(open(indir / 'j_dict.p', 'rb'))
    i_dict = pickle.load(open(indir / 'i_dict.p', 'rb'))
    G = zrfun.get_basic_info(fn00, only_G=True)
    h = G['h']
    xp = G['lon_rho']
    yp = G['lat_rho']
    plon_vec = np.array([])
    plat_vec = np.array([])
    hh_vec = np.array([])
    for seg_name in seg_list:
        jjj = j_dict[seg_name]
        iii = i_dict[seg_name]
        # untested 2021.10.05
        hh_vec = np.append(hh_vec, h[jjj,iii])
        plon_vec = np.append(plon_vec, xp[jjj,iii])
        plat_vec = np.append(plat_vec, yp[jjj,iii])
        # ji_seg = ji_dict[seg_name]
        # for ji in ji_seg:
        #     plon_vec = np.append(plon_vec, xp[ji])
        #     plat_vec = np.append(plat_vec, yp[ji])
        #     hh_vec = np.append(hh_vec, h[ji])
    plon00 = np.array([]); plat00 = np.array([]); pcs00 = np.array([])
    for ii in range(len(plon_vec)):
        x = plon_vec[ii]
        y = plat_vec[ii]
        hdz = DZ*np.floor(hh_vec[ii]/DZ) # depth to closest DZ m (above the bottom)
        if hdz >= DZ:
            zvec = np.arange(-hdz,DZ,DZ) # a vector that goes from -hdz to 0 in steps of DZ m
            svec = zvec/hh_vec[ii]
            ns = len(svec)
            if ns > 0:
                plon00 = np.append(plon00, x*np.ones(ns))
                plat00 = np.append(plat00, y*np.ones(ns))
                pcs00 = np.append(pcs00, svec)
    # subsample the I.C. vectors to around max length around NPmax
    NP = len(plon00)
    print(len(plon00))
    nstep = max(1,int(NP/NPmax))
    plon00 = plon00[::nstep]
    plat00 = plat00[::nstep]
    pcs00 = pcs00[::nstep]
    print(len(plon00))
    sys.stdout.flush()
    return plon00, plat00, pcs00
