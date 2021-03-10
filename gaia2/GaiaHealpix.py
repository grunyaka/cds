from astropy.coordinates import SkyCoord, EarthLocation
from astropy import coordinates as coord
from astropy.coordinates.tests.utils import randomly_sample_sphere
from astropy.time import Time
from astropy import units as u
import matplotlib.pyplot as plt
import numpy as np
from astroquery.vizier import Vizier
from astropy.coordinates import Angle
from astroquery.gaia import Gaia
import os
import math
import sys
import healpy as hp
import pandas as pd


filename = os.path.join('./','GaiaSource-my_close_stars.csv')
df = pd.read_csv(filename)[['ra', 'ra_error', 'dec', 'dec_error','pmra', 'pmra_error', 'pmdec', 'pmdec_error','g_rp','ref_epoch']]
NSIDE = 1
NPIX = hp.nside2npix(NSIDE)

hp_numbers = hp.ang2pix(NSIDE, df['ra'].to_numpy(), df['dec'].to_numpy(), lonlat=True, nest=True)
m = np.zeros(NPIX)
atall = 0
for number in hp_numbers:
    m[number] += 1
    atall += 1
print(m)
print(atall)
hp.mollview(m, title="Mollview image", nest = True)
hp.graticule()

