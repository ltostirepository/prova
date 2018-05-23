import numpy as np
import matplotlib.pyplot as plt

from gammapy.data import DataStore
from gammapy.scripts import SpectrumAnalysisIACT

# Convenience classes to define analsys inputs
# At some point we'll add a convenience layer to run the analysis starting from a plain text config file.
from gammapy.utils.energy import EnergyBounds
from gammapy.image import SkyImage
from gammapy.spectrum import models
from regions import CircleSkyRegion
from astropy.coordinates import SkyCoord
import astropy.units as u


import logging
logging.basicConfig()
log = logging.getLogger('gammapy.spectrum')
log.setLevel(logging.ERROR)
print("\nlog created")

obs_cols = ['OBS_ID', 'GLON_PNT', 'GLAT_PNT']
#store_dir = '/Users/lucatosti/gammapy-extra/datasets/hess-crab4-hd-hap-prod2'
store_dir = '/Users/lucatosti/Desktop/GEMINGA/1dc/index/all'
data_store = DataStore.from_dir(store_dir)
data_store.info()


table = data_store.obs_table                                                                                      
pos_obs = SkyCoord(table['GLON_PNT'], table['GLAT_PNT'], frame='galactic', unit='deg')
crab_pos = SkyCoord.from_name( 'crab' , frame='galactic')                             
#crab_pos = SkyCoord( 0 , 0  ,unit='deg', frame='galactic')#184.5575 , -05.7844                                                             
offset = crab_pos.separation(pos_obs).deg
#print(offset)

#if(offset.any()<10):
 #   print("qualcuno")
mask = (offset<0.5)#(1 < offset) & (offset < 2)                                                                                

h=plt.hist(offset,100)
plt.show(h)

#mask =  (offset < 50)                                                                                
#print (mask)

table = table[mask]                                                                                               


#table.show_in_browser(jsviewer=True)  
#print("crab_prom_name")
print("posizione crab")
print(crab_pos)
#print("posizione osservazioni")
#print(pos_obs)
#if(mask.any()==True):
 #   print("almeno uno")
#else:
 #   print("nessuno")

#store_dir = '/Users/lucatosti/gammapy-extra/datasets/hess-crab4-hd-hap-prod2'
data_store = DataStore.from_dir(store_dir)
obs_id = data_store.obs_table['OBS_ID'].data
print(len(obs_id))
exit()
print("Use observations {}".format(obs_id))

obs_list = data_store.obs_list(obs_id)

#crab_pos = SkyCoord.from_name('crab')
on_region = CircleSkyRegion(crab_pos, 0.15 * u.deg)

model = models.LogParabola(
    alpha = 2.3,
    beta = 0,
    amplitude = 1e-11 * u.Unit('cm-2 s-1 TeV-1'),
    reference = 1 * u.TeV,
)

flux_point_binning = EnergyBounds.equal_log_spacing(0.7, 30, 5, u.TeV)# 0.7, 30, 5, u.TeV

exclusion_mask = SkyImage.read('/Users/lucatosti/gammapy-extra/datasets/exclusion_masks/tevcat_exclusion.fits')

config = dict(
    outdir = None,
    background = dict(
        on_region=on_region,
        exclusion_mask=exclusion_mask,
        min_distance = 0.1 * u.rad,
    ),
    extraction = dict(containment_correction=False),
    fit = dict(
        model=model,
        stat='wstat',
        forward_folded=True,
        fit_range = flux_point_binning[[0, -1]]
    ),
    fp_binning=flux_point_binning
)

######farà un sacco di warning innocui
ana = SpectrumAnalysisIACT(
    observations=obs_list,
    config=config,
)

ana.run() #run_extraction se non voglio fare il fit

print(ana.fit.result[0])

a = ana.spectrum_result.plot(
    energy_range=ana.fit.fit_range,
    energy_power=2,
    flux_unit='erg-1 cm-2 s-1',
    fig_kwargs=dict(figsize = (8,8)),
)


plt.show(a)


print("COMPLETE")


"""
fit = dict(
        model=model,
        stat='wstat',
        forward_folded=True,
        fit_range = flux_point_binning[[0, -1]]
    ),
"""
