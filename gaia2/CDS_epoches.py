from astropy.coordinates import SkyCoord, EarthLocation
from astropy import coordinates as coord
from astropy.coordinates.tests.utils import randomly_sample_sphere
from astropy.time import Time
from astropy import units as u
import numpy as np
from astroquery.vizier import Vizier
from astropy.coordinates import Angle
from astroquery.gaia import Gaia
import os
import math
import sys



names = ['CMC15', 'Pan-STARRS DR1', '2MASS All-Sky', 'SDSS DR12', 'URAT1', 'UCAC4']
Nnames = ['0.CMC__15', '2.PanSta1', '3.2massAS', '4.SDSS_12', '5.URAT__1', '6.UCAC__4']#это только чтобы логи на экран выводить, в теории можно убрать, если убрать print'ы
links = ['I/327/cmc15', 'II/349/ps1', 'II/246/out', 'V/147/sdss12', 'I/329/urat1', 'I/322A/out']
signs = ['RA_ICRS', 'DE_ICRS', 'e_RA_ICRS', 'e_DE_ICRS', 'MJD-51263', 'Jmag',
        'RAJ2000', 'DEJ2000', 'e_RAJ2000', 'e_DEJ2000', '_tab1_10', 'gmag',
        'RAJ2000', 'DEJ2000', 'errMaj', 'errMin', '_tab1_36', 'Jmag',  # _tab1_36 - JD
        'RA_ICRS', 'DE_ICRS', 'e_RA_ICRS', 'e_RA_ICRS', 'ObsDate', 'gmag',
        'RAJ2000', 'DEJ2000', 'sigs', 'sigm', 'Epoch', 'gmag',
        'RAJ2000', 'DEJ2000', 'e_RAJ2000', 'e_DEJ2000', 'EpRA', 'gmag']

filename = os.path.join('./','list.txt')
outfolder = os.path.join('./','CDS_epoches')

