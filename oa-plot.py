# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.1
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Notebook to plot OA data for NCRA 

# %%
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
# %config Completer.use_jedi = False

# %%
client = Client()
client

# %% [markdown]
# # Read in data

# %%
dir1='/home/599/rxm599/'
file2a='oa_GW2p0.nc'
file2b='ph_GW2p0.nc'
file0a='oa_PI.nc'
file0b='ph_PI.nc'

# %time dbgc2a = xr.open_dataset(dir1+file2a )
# %time dbgc2b = xr.open_dataset(dir1+file2b )
# %time dbgc0a = xr.open_dataset(dir1+file0a)
# %time dbgc0b = xr.open_dataset(dir1+file0b)
# #%time dbgc0 = xr.open_mfdataset(dir1+"ave*.nc", chunks={"TIME41": -1,"XT_OCEAN": 1,"YT_OCEAN": 1} )#.chunk(chunks={"TIME41":-1})

# #%time dbgc0 = xr.open_mfdataset(dir1+"files*.nc", chunks={"TIME41": -1,"XT_OCEAN": 1,"YT_OCEAN": 1} )#.chunk(chunks={"TIME41":-1})
#dbgc1 = xr.open_mfdataset(dir1+"files*.nc", chunks={"TIME41": -1,"XT_OCEAN":1,"YT_OCEAN":1} ).chunk(chunks={"TIME41":-1})

# %% [markdown]
# ## Plot

# %%
dbgc0b

# %%
ph_diff=dbgc2b.PH-dbgc0b.PH

fig, axs = plt.subplots(ncols=2,figsize=[10,5])
dbgc2b.PH.plot(levels=20,ax=axs[0])
plt.title('pH ')
plt.xlabel("Longitude")
plt.ylabel("Latitude")

ph_diff.plot(levels=20,ax=axs[1])
plt.xlabel("Longitude")
plt.ylabel("Latitude")
plt.tight_layout()
plt.show()

# %%

# %%
