import matplotlib.pyplot as plt

import gammapy
import numpy as np
import astropy
import regions
import sherpa

print('gammapy:', gammapy.__version__)
print('numpy:', np.__version__)
print('astropy', astropy.__version__)
print('regions', regions.__version__)
print('sherpa', sherpa.__version__)

import astropy.units as u
from astropy.coordinates import SkyCoord, Angle
from astropy.table import vstack as vstack_table
from regions import CircleSkyRegion
from gammapy.data import DataStore, ObservationList
from gammapy.data import ObservationStats, ObservationSummary
from gammapy.background.reflected import ReflectedRegionsBackgroundEstimator
from gammapy.utils.energy import EnergyBounds
from gammapy.spectrum import SpectrumExtraction, SpectrumObservation, SpectrumFit, SpectrumResult
from gammapy.spectrum.models import PowerLaw, ExponentialCutoffPowerLaw, LogParabola
from gammapy.spectrum import FluxPoints, SpectrumEnergyGroupMaker, FluxPointEstimator
from gammapy.image import SkyImage



import logging
logging.basicConfig()
log = logging.getLogger('gammapy.spectrum')
log.setLevel(logging.WARNING)


#DATA_DIR = '/Users/lucatosti/gammapy-extra/datasets/hess-crab4-hd-hap-prod2'#'$GAMMAPY_EXTRA/datasets/hess-crab4-hd-hap-prod2'
DATA_DIR = '/Users/lucatosti/Desktop/GEMINGA/1dc/index/all'
datastore = DataStore.from_dir(DATA_DIR)



table = datastore.obs_table
pos_obs = SkyCoord(table['GLON_PNT'], table['GLAT_PNT'], frame='galactic', unit='deg')
#crab_pos = SkyCoord.from_name( 'crab' , frame='galactic')
crab_pos = SkyCoord( 0 , 0  ,unit='deg', frame='galactic')#184.5575 , -05.7844                                                                
offset = crab_pos.separation(pos_obs).deg
mask = (0.5<offset) & (offset<2)#(1 < offset) & (offset < 2)      
h=plt.hist(offset,100)
plt.show(h)
table = table[mask]




datastore = DataStore.from_dir(DATA_DIR)
#obs_ids = [23523, 23526, 23559, 23592]
obs_ids = [510000, 510001, 510002, 510003]  
#obs_ids = datastore.obs_table['OBS_ID'].data
obs_list = datastore.obs_list(obs_ids)


target_position = SkyCoord.from_name('crab',  frame='galactic')
on_region_radius = Angle('0.11 deg')
on_region = CircleSkyRegion(center=target_position, radius=on_region_radius)


exclusion_region = CircleSkyRegion(
    center=SkyCoord.from_name('crab', frame='galactic'),
    radius=0.1 * u.deg,
)

image_center = target_position.galactic
exclusion_mask = SkyImage.empty(
    nxpix=150, nypix=150, binsz=0.05,
    xref=image_center.l.deg, yref=image_center.b.deg,
    proj='TAN', coordsys='GAL',
)

exclusion_mask = exclusion_mask.region_mask(exclusion_region)
exclusion_mask.data = 1. - exclusion_mask.data


background_estimator = ReflectedRegionsBackgroundEstimator(
    obs_list=obs_list,
    on_region=on_region,
    exclusion_mask = exclusion_mask)

background_estimator.run()

plt.figure(figsize=(8,8))
bkg_ex_plot=background_estimator.plot()
plt.show(bkg_ex_plot)

#######################
stats = []
for obs, bkg in zip(obs_list, background_estimator.result):
    stats.append(ObservationStats.from_obs(obs, bkg))

print(stats[1])

obs_summary = ObservationSummary(stats)
fig = plt.figure(figsize=(10,6))
ax1=fig.add_subplot(121)

aa = obs_summary.plot_excess_vs_livetime(ax=ax1)
ax2=fig.add_subplot(122)
bb = obs_summary.plot_significance_vs_livetime(ax=ax2)
plt.show(aa)
plt.show(bb)
########################

e_reco = EnergyBounds.equal_log_spacing(0.1, 40, 40, unit='TeV')
e_true = EnergyBounds.equal_log_spacing(0.05, 100., 200, unit='TeV')



ANALYSIS_DIR = 'crab_analysis'

extraction = SpectrumExtraction(
    obs_list=obs_list,
    bkg_estimate=background_estimator.result,
    containment_correction=False,
)
extraction.run()

# Add a condition on correct energy range in case it is not set by default
extraction.compute_energy_threshold(method_lo='area_max', area_percent_lo=10.0)

print(extraction.observations[0])
# Write output in the form of OGIP files: PHA, ARF, RMF, BKG
# extraction.run(obs_list=obs_list, bkg_estimate=background_estimator.result, outdir=ANALYSIS_DIR)



cc = extraction.observations[0].peek()
plt.show(cc)


#######################



model = PowerLaw(
    index=2 * u.Unit(''),
    amplitude=2e-11 * u.Unit('cm-2 s-1 TeV-1'),
    reference=1 * u.TeV,
)

joint_fit = SpectrumFit(obs_list=extraction.observations, model=model)

joint_fit.fit()
joint_fit.est_errors()
#fit.run(outdir = ANALYSIS_DIR)

joint_result = joint_fit.result


ax0, ax1 = joint_result[0].plot(figsize=(8,8))
ax0.set_ylim(0, 20)
dd = joint_result[0]
plt.show(dd)


#######################


ebounds = [0.3, 1.1, 3, 10.1, 30] * u.TeV

stacked_obs = extraction.observations.stack()

seg = SpectrumEnergyGroupMaker(obs=stacked_obs)
seg.compute_range_safe()
seg.compute_groups_fixed(ebounds=ebounds)

print(seg.groups)


fpe = FluxPointEstimator(
    obs=stacked_obs,
    groups=seg.groups,
    model=joint_result[0].model,
)
fpe.compute_points()


ee=fpe.flux_points.plot()
fpe.flux_points.table
plt.show(ee)

#####################

spectrum_result = SpectrumResult(
    points=fpe.flux_points,
    model=joint_result[0].model,
)
ax0, ax1 = spectrum_result.plot(
    energy_range=joint_fit.result[0].fit_range,
    energy_power=2, flux_unit='erg-1 cm-2 s-1',
    fig_kwargs=dict(figsize=(8,8)),
    point_kwargs=dict(color='navy')
)

ax0.set_xlim(0.4, 50)
plt.show(ax0)
plt.show(ax1)


stacked_obs = extraction.observations.stack()

stacked_fit = SpectrumFit(obs_list=stacked_obs, model=model)
stacked_fit.fit()
stacked_fit.est_errors()


stacked_result = stacked_fit.result
print(stacked_result[0])

stacked_table = stacked_result[0].to_table(format='.3g')
stacked_table['method'] = 'stacked'
joint_table = joint_result[0].to_table(format='.3g')
joint_table['method'] = 'joint'
total_table = vstack_table([stacked_table, joint_table])
print(total_table['method', 'index', 'index_err', 'amplitude', 'amplitude_err'])
