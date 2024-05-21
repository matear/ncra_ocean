# # Notebook to extract the Temp data for NCRA 

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


# +
dir1='/g/data/fp2/OFAM3/jra55_historical.1/surface/'

# %time dsst1 = xr.open_mfdataset(dir1+"ocean_temp_sfc_*", parallel=True,chunks={"Time": 10})
# -

dir2='/g/data/fp2/OFAM3/jra55_rcp8p5/surface/'
# %time dsst2 = xr.open_mfdataset(dir2+"ocean_temp_sfc*.nc", parallel=True,chunks={"Time": 10})
    #.chunk(chunks={"TIME41":-1})

dsst1.temp

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

# + jupyter={"outputs_hidden": true}
sstGW1p1=dsst1.temp.sel(Time=slice(GW1p1_start,GW1p1_end))
# %time sst=sstGW1p1.mean('Time').compute()
sst.to_netcdf('sst_mcurrent.nc')
# -

sstGW1p1

clim,seas,anom=climatology(sstGW1p1,'Time')

clim[0,:,:,:].plot()

# %time clim.to_netcdf('sst_current.nc')

# ## Future

# +

sstGW1p5=dsst2.temp.sel(Time=slice(GW1p5_start,GW1p5_end))
clim,seas,anom=climatology(sstGW1p5,'Time')
# %time clim.to_netcdf('sst_GW1p5.nc')
# -

sstGW2p0=dsst2.temp.sel(Time=slice(GW2p0_start,GW2p0_end))
clim,seas,anom=climatology(sstGW2p0,'Time')
# %time clim.to_netcdf('sst_GW2p0.nc')

sstGW3p0=dsst2.temp.sel(Time=slice(GW3p0_start,GW3p0_end))
clim,seas,anom=climatology(sstGW3p0,'Time')
# %time clim.to_netcdf('sst_GW3p0.nc')


sstGW4p0=dsst2.temp.sel(Time=slice(GW4p0_start,GW4p0_end))
clim,seas,anom=climatology(sstGW4p0,'Time')
# %time clim.to_netcdf('sst_GW4p0.nc')

# +

clim.sel(yt_ocean=-34,xt_ocean=165, method="nearest").plot()
# -

#dir2='/g/data/fp2/OFAM3/jra55_rcp8p5/surface/'
# %time dtmp = xr.open_mfdataset("sst_GW2p0.nc")
    #.chunk(chunks={"TIME41":-1})

dtmp.temp.sel(yt_ocean=-34,xt_ocean=165, method="nearest").plot()

dtmp.temp.sel(month=12).plot()

dtmp.temp


