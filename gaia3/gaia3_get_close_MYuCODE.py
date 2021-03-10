import ssl
from astroquery.gaia import Gaia

import math
import numpy as np

plx_lim = 50
mag_lim = 12.5
pmra_lim = 100
pmdec_lim = 100


ssl._create_default_https_context = ssl._create_unverified_context
job = Gaia.launch_job_async(
    'SELECT*FROM gaiaedr3.gaia_source WHERE parallax>{} AND phot_g_mean_mag>={} AND abs(pmra)>{} '
    'AND abs(pmdec)>{}'.format(plx_lim, mag_lim, pmra_lim, pmdec_lim))

gaiat = job.get_results()

print(gaiat)



