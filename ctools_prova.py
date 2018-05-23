import gammalib
import ctools
import cscripts

env CTADATA=/Users/lucatosti/Desktop/GEMINGA/1dc
env CALDB=/Users/lucatosti/Desktop/GEMINGA/1dc/caldb 

obsselect = cscripts.csobsselect()
obsselect['inobs']     = '$CTADATA/obs/obs_gps_baseline.xml'
obsselect['pntselect'] = 'CIRCLE'
obsselect['coordsys']  = 'GAL'
obsselect['glon']      = 0.0
obsselect['glat']      = 0.0
obsselect['rad']       = 3.0
obsselect['tmin']      = 'NONE'
obsselect['tmax']      = 'NONE'
obsselect['outobs']    = 'obs.xml'
obsselect.execute()


select = ctools.ctselect()
select['inobs']   = 'obs.xml'
select['ra']      = 'NONE'
select['dec']     = 'NONE'
select['rad']     = 'NONE'
select['tmin']    = 'NONE'
select['tmax']    = 'NONE'
select['emin']    = 1.0
select['emax']    = 100.0
select['outobs']  = 'obs_selected.xml'
select.execute()

binning = ctools.ctbin()
binning['inobs']    = 'obs_selected.xml'
binning['xref']     = 0.0
binning['yref']     = 0.0
binning['coordsys'] = 'GAL'
binning['proj']     = 'CAR'
binning['binsz']    = 0.02
binning['nxpix']    = 300
binning['nypix']    = 300
binning['ebinalg']  = 'LOG'
binning['emin']     = 1.0
binning['emax']     = 100.0
binning['enumbins'] = 20
binning['outcube']  = 'cntcube.fits'
binning.execute()


expcube = ctools.ctexpcube()
expcube['inobs']   = 'obs_selected.xml'
expcube['incube']  = 'cntcube.fits'
expcube['outcube'] = 'expcube.fits'
expcube.execute()


psfcube = ctools.ctpsfcube()
psfcube['inobs']    = 'obs_selected.xml'
psfcube['incube']   = 'NONE'
psfcube['ebinalg']  = 'LOG'
psfcube['emin']     = 1.0
psfcube['emax']     = 100.0
psfcube['enumbins'] = 20
psfcube['nxpix']    = 10
psfcube['nypix']    = 10
psfcube['binsz']    = 1.0
psfcube['coordsys'] = 'GAL'
psfcube['proj']     = 'CAR'
psfcube['xref']     = 0.0
psfcube['yref']     = 0.0
psfcube['outcube']  = 'psfcube.fits'
psfcube.execute()


bkgcube = ctools.ctbkgcube()
bkgcube['inobs']    = 'obs_selected.xml'
bkgcube['inmodel']  = '$CTOOLS/share/models/bkg_irf.xml'
bkgcube['incube']   = 'cntcube.fits'
bkgcube['outcube']  = 'bkgcube.fits'
bkgcube['outmodel'] = 'bkgcube.xml'
bkgcube.execute()

