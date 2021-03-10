import numpy as np
from astroquery.vizier import Vizier
from astroquery.gaia import Gaia



ra0 = 65.40019016727416
dec0 = -27.046093160540174
limmag = 18
fov = 1

job = Gaia.launch_job_async("SELECT * \
        FROM gaiadr2.gaia_source \
        WHERE CONTAINS(POINT('ICRS',gaiadr2.gaia_source.ra,gaiadr2.gaia_source.dec),CIRCLE('ICRS',%f,%f,%f))=1\
                                   AND pmra IS NOT NULL AND abs(pmra)>0 \
        AND pmdec IS NOT NULL AND abs(pmdec)>0\
                                   AND  phot_g_mean_mag<%f;" % (ra0, dec0, fov, limmag) \
                                    , dump_to_file=False)
table = job.get_results()
print(table.columns)


ra_g = np.array(table['ra'])
dec_g = np.array(table['dec'])
era_g = np.array(table['ra_error'])
edec_g = np.array(table['dec_error'])
t_g = np.array(table['ref_epoch'])
pmra = np.array(table['pmra'])
pmdec = np.array(table['pmdec'])
gmag_g = np.array(table['phot_g_mean_mag'])
