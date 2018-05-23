import matplotlib.pyplot as plt
import gammapy
import numpy as np
import astropy
import regions
import sherpa
import uncertainties
import photutils

import astropy.units as u
from astropy.coordinates import SkyCoord, Angle
from regions import CircleSkyRegion
from photutils.detection import find_peaks
from gammapy.data import DataStore
from gammapy.spectrum import (SpectrumExtraction,SpectrumFit,SpectrumResult,models,SpectrumEnergyGroupMaker,FluxPointEstimator)
from gammapy.image import SkyImage, IACTBasicImageEstimator
from gammapy.background import RingBackgroundEstimator, ReflectedRegionsBackgroundEstimator
from gammapy.utils.energy import EnergyBounds
from gammapy.detect import TSImageEstimator

from gammapy.catalog import SourceCatalog2FHL
from astropy.convolution import Ring2DKernel, Tophat2DKernel
from astropy.visualization import simple_norm

from gammapy.data import DataStore
from gammapy.image import SkyImage, SkyImageList
from gammapy.detect import KernelBackgroundEstimator as KBE

######## creo log per evitare il verbose   
import logging 
logging.basicConfig()
log = logging.getLogger('gammapy.spectrum')
log.setLevel(logging.ERROR)
print("\nlog created")

def show_image(image, radius=3, vmin=0, vmax=3,name='geminga'):
    """Little helper function to show the images for this application here."""
    image.smooth(radius=radius).show(vmin=vmin, vmax=vmax, add_cbar=True)
    image.cutout(position=SkyCoord.from_name(name,frame='galactic'),size=(2*u.deg, 3*u.deg)).smooth(radius=radius).show(vmin=vmin, vmax=vmax, add_cbar=True)



print("\n")
print("CHECK PACKAGE VERSION ########")
print('gammapy:', gammapy.__version__)
print('numpy:', np.__version__)
print('astropy:', astropy.__version__)
print('regions:', regions.__version__)
print('sherpa:', sherpa.__version__)
print('uncertainties:', uncertainties.__version__)
print('photutils:', photutils.__version__)
print("##############################\n")

###metti all al posto di gps o cambia all in extragalactic survey

data_store = DataStore.from_dir('/Users/lucatosti/Desktop/GEMINGA/1dc/index/gps') #### oggetto generale tipo database per coordinate
#data_store = DataStore.from_dir('$GAMMAPY_EXTRA/datasets/hess-crab4-hd-hap-prod2/')
data_store.info() 

show_it=True

# in generale selezione osservazioni---->>> funziona anche direttamente  SkyCoord.from_name('crab') --->>>metti geminga

# from astropy.coordinates import SkyCoord
# table = data_store.obs_table
# pos_obs = SkyCoord(table['GLON_PNT'], table['GLAT_PNT'], frame='galactic', unit='deg')
# pos_target = SkyCoord(0, 0, frame='galactic', unit='deg')
# offset = pos_target.separation(pos_obs).deg
# mask = (1 < offset) & (offset < 2)
# table = table[mask]
# table.show_in_browser(jsviewer=True)

"""
obs_id = [110380, 111140, 111159]
obs_list = data_store.obs_list(obs_id) ###modificare per creare lista completa
    ### applicare filtro -> solo osservazioni su geminga

#print('##########')
print(obs_list)
#print('##########')

obs_cols = ['OBS_ID', 'GLON_PNT', 'GLAT_PNT', 'LIVETIME']
#table_orig=data_store.obs_table.select_obs_id(obs_id)[obs_cols]

table = data_store.obs_table[obs_cols]

obs_list2 = data_store.obs_list(table['OBS_ID'])
#print(obs_list2)

pos_obs = SkyCoord(table['GLON_PNT'], table['GLAT_PNT'], frame='galactic', unit='deg')
#target_position = SkyCoord(0, 0, unit='deg', frame='galactic')
#pos_target = SkyCoord(0,0,unit='deg',frame='galactic')
pos_target = SkyCoord.from_name('geminga',frame='galactic')
print(pos_target)
offset = pos_target.separation(pos_obs).deg
bool1= (1 < offset) 
bool2= (offset < 1.001)
#bool3= (table['OBS_ID']==110380 | table['OBS_ID']== 111140 | table['OBS_ID']==111159)
mask= bool1 & bool2
table = table[mask]
table.show_in_browser(jsviewer=True)

#print(table_orig)
#print(table)

on_radius = 0.2 * u.deg
on_region = CircleSkyRegion(center=pos_target, radius=on_radius)


# Define reference image centered on the target
xref = pos_target.galactic.l.value
yref = pos_target.galactic.b.value
size = 10 * u.deg
binsz = 0.02 # degree per pixel
npix = int((size / binsz).value)

ref_image = SkyImage.empty(
    nxpix=800, nypix=600, binsz=0.02,
    xref=xref, yref=yref,
    proj='TAN', coordsys='GAL',
)

print(ref_image)


exclusion_mask = ref_image.region_mask(on_region) #### creazione maschera circolare per stime bkg
exclusion_mask.data = 1 - exclusion_mask.data
exclusion_mask.plot()####non vuole stamparlo

bkg_estimator = RingBackgroundEstimator(
    r_in=0.5 * u.deg,
    width=0.2 * u.deg,
)
image_estimator = IACTBasicImageEstimator(
    reference=ref_image,
    emin=100 * u.GeV,
    emax=100 * u.TeV,
    #offset_max=3 * u.deg,
    background_estimator=bkg_estimator,
    exclusion_mask=exclusion_mask,
)

"""
##c'è qualcosa che non va con le immagini -> scarta tutto il tutorial gammapy e prova : 

"""NON QUESTO PER CARITà
images = image_estimator.run(obs_list2)


images.names
show_image(images['counts'], radius=0, vmax=10)
show_image(images['counts'], vmax=5)
show_image(images['background'], vmax=4)
show_image(images['excess'], vmax=2)

"""
#110380, 111140, 111159

#source_pos = SkyCoord.from_name('geminga',frame='galactic')#(83.633083, 22.0145, unit='deg')
# If you have internet access, you could also use this to define the `source_pos`:
#source_pos = SkyCoord.from_name('crab')
source_pos = SkyCoord(0, 0, unit='deg',frame='galactic')
#source_pos = SkyCoord.from_name('crab')
print(source_pos)

ref_image = SkyImage.empty(
    nxpix=800, nypix=800, binsz=0.02,
    #xref=source_pos.ra.deg, yref=source_pos.dec.deg,
    xref=source_pos.l.deg , yref=source_pos.b.deg,
    coordsys='GAL', proj='TAN', #coordsys='CEL'
)
print("1")
#events = data_store.obs(obs_id=23523).events###110380
events = data_store.obs(obs_id=110380).events
counts_image = SkyImage.empty_like(ref_image)
counts_image.fill_events(events)


norm = simple_norm(counts_image.data, stretch='sqrt', min_cut=0, max_cut=0.3)#0.3
q1=counts_image.smooth(radius=0.1 * u.deg).plot(norm=norm, add_cbar=True)
q=counts_image.cutout(position=SkyCoord(-5, 0, unit='deg', frame='galactic'),size=(20*u.deg, 20*u.deg)).smooth(radius=0.1 * u.deg).plot(norm=norm, add_cbar=True)
#plt.show(counts_image)
if(show_it==False):
    plt.show(q)
    
#obs_ids = [23523, 23526]  ###[111140, 111159]
obs_ids = [110380,111140, 111159]
counts_image2 = SkyImage.empty_like(ref_image)
for obs_id in obs_ids:
    events = data_store.obs(obs_id=obs_id).events
    counts_image2.fill_events(events)

norm = simple_norm(counts_image2.data, stretch='sqrt', min_cut=0, max_cut=0.5)#0.5
counts_image2.smooth(radius=0.1 * u.deg).plot(norm=norm, add_cbar=True)


q2=counts_image2#.cutout(position=SkyCoord(0, 0, unit='deg', frame='galactic'),size=(10*u.deg, 10*u.deg)).plot(norm=norm, add_cbar=True)
if(show_it==False):
    plt.show(q2)


print("2")


################QUI C'è IL CONTO SUL BKG


source_kernel = Tophat2DKernel(radius=5)#5
source_kernel.normalize(mode='peak')
source_kernel = source_kernel.array

background_kernel = Ring2DKernel(radius_in=5, width=20)#20-10
background_kernel.normalize(mode='peak')
background_kernel = background_kernel.array


a = plt.imshow(source_kernel, interpolation='nearest', cmap='gray')
plt.colorbar()
plt.grid('off')
if(show_it==False):
    plt.show(a)

b = plt.imshow(background_kernel, interpolation='nearest', cmap='gray')
plt.colorbar()
plt.grid('off')
if(show_it==False):
    plt.show(b)
# To use the `KernelBackgroundEstimator` you first have to set
# up a source and background kernel and put the counts image input
# into a container `SkyImageList` class.
images = SkyImageList()
images['counts'] = counts_image2

print("3")

kbe = KBE(
    kernel_src=source_kernel,
    kernel_bkg=background_kernel,
    significance_threshold=5,
    mask_dilation_radius=0.06 * u.deg,
)
# This takes about 10 seconds on my machine
result2 = kbe.run(images)
print("results2")
print(result2.names)
####altro modo stima bkg

obs_list = data_store.obs_list(obs_ids)
target_position = SkyCoord(0, 0, unit='deg', frame='galactic')
on_radius = 0.3 * u.deg
on_region = CircleSkyRegion(center=target_position, radius=on_radius)

exclusion_mask = ref_image.region_mask(on_region)
exclusion_mask.data = 1 - exclusion_mask.data
bkg_estimator = RingBackgroundEstimator(
    r_in=0.5 * u.deg,
    width=0.2 * u.deg,
)
image_estimator = IACTBasicImageEstimator(
    reference=ref_image,
    emin=100 * u.GeV,
    emax=100 * u.TeV,
    #offset_max=3 * u.deg,
    background_estimator=bkg_estimator,
    exclusion_mask=exclusion_mask,
)
images_new = image_estimator.run(obs_list)
print("\n\n\n")
print(images_new.names)
print(result2.names)
print("\n\n\n")
# Let's have a look at the background image and the exclusion mask

# This doesn't work yet ... need to do SkyImage.plot fixes:
# fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(10, 3))
# background_image.plot(ax=axes[0])
# exclusion_image.plot(ax=axes[1])
# significance_image.plot(ax=axes[2])
print("bkg metodo 1")
background_image = result2['background']
background_image2 = images_new['background']
norm = simple_norm(background_image.data, stretch='sqrt', min_cut=0, max_cut=0.5)
background_image.plot(norm=norm, add_cbar=True)
if(show_it==False):
    plt.show(background_image.plot(norm=norm, add_cbar=True))
norm2 = simple_norm(background_image2.data, stretch='sqrt', min_cut=0, max_cut=0.5)
background_image2.plot(norm=norm, add_cbar=True)
print("bkg metodo 2")
if(show_it==False):
    plt.show(background_image2.plot(norm=norm, add_cbar=True))
    plt.show(background_image.cutout(position=SkyCoord(0, 0, unit='deg', frame='galactic'),size=(2*u.deg, 2*u.deg)))
result2['exclusion'].plot()
if(show_it==True):
    plt.show(result2['exclusion'].plot())
    plt.show(result2['background'].plot())
    plt.show(images_new['background'].plot())
    plt.show(result2['background'].cutout(position=SkyCoord(0, 0, unit='deg', frame='galactic'),size=(2*u.deg,2*u.deg)).plot())
    plt.show(images_new['background'].cutout(position=SkyCoord(0, 0, unit='deg', frame='galactic'),size=(2*u.deg, 2*u.deg)).plot())

significance_image = result2['significance']
c=significance_image.cutout(position=SkyCoord(0, 0, unit='deg', frame='galactic'),size=(2*u.deg, 5*u.deg)).plot(add_cbar=True, vmin=0, vmax=20)
if(show_it==False):
    plt.show(c)
print("4")
if(show_it==False):
    plt.show(images_new['psf'].plot())
    plt.show(images_new['excess'].cutout(position=SkyCoord(0, 0, unit='deg', frame='galactic'),size=(2*u.deg, 2*u.deg)).plot())
