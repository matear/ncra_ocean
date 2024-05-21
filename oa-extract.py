# # Notebook to extract the OA data for NCRA 

import xarray as xr
from dask_jobqueue import PBSCluster
from dask.distributed import Client
import numpy as np
#import xrft
import scipy
import matplotlib.pyplot as plt
import datetime
import pandas as pd
import matplotlib.dates as mdates
from matplotlib.dates import DateFormatter

client = Client()
client


# ## Read in data

def climatology(dsx,TIME1):
    clim = dsx.groupby(TIME1+'.month').mean(dim=TIME1)
    anom = dsx.groupby(TIME1+'.month') - clim
    season=dsx.groupby(TIME1+'.season').mean(dim=TIME1)
    return(clim,season,anom)


dir1='/scratch/xv83/rxm599/historical/'
# %time dbgc1 = xr.open_mfdataset(dir1+"files*.nc", parallel=True,chunks={"TIME41": 10})

dir2='/scratch/xv83/rxm599/future/'
# %time dbgc2 = xr.open_mfdataset(dir2+"files*.nc", parallel=True,chunks={"TIME41": 10})
    #.chunk(chunks={"TIME41":-1})
# #%time dbgc0 = xr.open_mfdataset(dir1+"ave*.nc" )

dir0='/scratch/xv83/rxm599/pi/'
# %time dbgc0 = xr.open_mfdataset(dir0+"files*.nc", parallel=True,chunks={"TIME41": 10})

dbgc0.PH

# ## Set GWL 

# +
# set GWL 
GW1p1_start='1995-01-01'
GW1p1_end='2014-12-31'
#GW1p1_end='2014-12-31'
GW1p5_start='2015-01-01'
GW1p5_end='2034-12-31'

GW2p0_start='2030-01-01'
GW2p0_end='2049-12-31'
GW3p0_start='2053-01-01'
GW3p0_end='2072-12-31'
GW4p0_start='2074-01-01'
GW4p0_end='2093-12-31'
# -

# ## Present-day

aragGW1p1=dbgc1.OAR.sel(TIME41=slice(GW1p1_start,GW1p1_end))
# %time arag=aragGW1p1.mean('TIME41').compute()
arag.to_netcdf('oa_current.nc')

phGW1p1=dbgc1.PH.sel(TIME41=slice(GW1p1_start,GW1p1_end))
# %time ph=phGW1p1.mean('TIME41').compute()
ph.to_netcdf('ph_current.nc')

# + jupyter={"source_hidden": true}
arag.plot(levels=30)
# -

# ## PI



aragGW1p1=dbgc0.OAR.sel(TIME41=slice(GW1p1_start,GW1p1_end))
# %time arag=aragGW1p1.mean('TIME41').compute()
arag.to_netcdf('oa_PI.nc')

phGW1p1=dbgc0.PH.sel(TIME41=slice(GW1p1_start,GW1p1_end))
# %time ph=phGW1p1.mean('TIME41').compute()
ph.to_netcdf('ph_PI.nc')

# ## Future

aragGW1p5=dbgc2.OAR.sel(TIME41=slice(GW1p5_start,GW1p5_end))
# %time arag=aragGW1p5.mean('TIME41').compute()

arag.to_netcdf('oa_GW1p5.nc')

aragGW2p0=dbgc2.OAR.sel(TIME41=slice(GW2p0_start,GW2p0_end))
# %time arag=aragGW2p0.mean('TIME41').compute()

arag.to_netcdf('oa_GW2p0.nc')

aragGW3p0=dbgc2.OAR.sel(TIME41=slice(GW3p0_start,GW3p0_end))
# %time arag=aragGW3p0.mean('TIME41').compute()

arag.to_netcdf('oa_GW3p0.nc')

aragGW4p0=dbgc2.OAR.sel(TIME41=slice(GW4p0_start,GW4p0_end))
# %time arag=aragGW4p0.mean('TIME41').compute()

arag.to_netcdf('oa_GW4p0.nc')

arag.plot(levels=[0,1,2,3,4,20])

# ## pH

phGW1p1=dbgc2.PH.sel(TIME41=slice(GW1p1_start,GW1p1_end))
# %time ph=phGW1p1.mean('TIME41').compute()
ph.to_netcdf('ph_GW1p1.nc')

phGW1p5=dbgc2.PH.sel(TIME41=slice(GW1p5_start,GW1p5_end))
# %time ph=phGW1p5.mean('TIME41').compute()
ph.to_netcdf('ph_GW1p5.nc')

phGW2p0=dbgc2.PH.sel(TIME41=slice(GW2p0_start,GW2p0_end))
# %time ph=phGW2p0.mean('TIME41').compute()
ph.to_netcdf('ph_GW2p0.nc')

phGW3p0=dbgc2.PH.sel(TIME41=slice(GW3p0_start,GW3p0_end))
# %time ph=phGW3p0.mean('TIME41').compute()
ph.to_netcdf('ph_GW3p0.nc')

phGW4p0=dbgc2.PH.sel(TIME41=slice(GW4p0_start,GW4p0_end))
# %time ph=phGW4p0.mean('TIME41').compute()
ph.to_netcdf('ph_GW4p0.nc')



