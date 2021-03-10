from astroquery.gaia import Gaia
import numpy as np
import math

import ssl


def tangFromRADE(ra, de, RA, DE):
    ksi = math.cos(de)*math.sin(ra-RA)/(math.sin(de)*math.sin(DE)+math.cos(de)*math.cos(DE)*math.cos(ra-RA))
    eta = (math.sin(de)*math.cos(DE)-math.cos(de)*math.sin(DE)*math.cos(ra-RA))/(math.sin(de)*math.sin(DE)+math.cos(de)*math.cos(DE)*math.cos(ra-RA));
    return ksi,eta

def RADecFromTang(ksi, eta, RA, Dec):
    secRho = math.sqrt(1+ksi*ksi+eta*eta)
    t11 = -math.sin(RA)
    t12 = math.cos(RA)
    t13 = 0
    t21 = -math.cos(RA)*math.sin(Dec)
    t22 = -math.sin(RA)*math.sin(Dec)
    t23 = math.cos(Dec)
    t31 = math.cos(RA)*math.cos(Dec)
    t32 = math.sin(RA)*math.cos(Dec)
    t33 = math.sin(Dec)
    x = (ksi*t11+eta*t21+t31)/secRho
    y = (ksi*t12+eta*t22+t32)/secRho
    z = (ksi*t13+eta*t23+t33)/secRho
    ra = math.atan2(y,x)
    dec = math.atan2(z,math.sqrt(x*x+y*y))
    if(ra<0):
        ra+=2*math.pi
    return ra,dec


def getGaiaDR2(ra0,dec0,limgmag,fov,yr):
    ssl._create_default_https_context = ssl._create_unverified_context
    job = Gaia.launch_job_async("SELECT * \
    FROM gaiadr2.gaia_source \
    WHERE CONTAINS(POINT('ICRS',gaiadr2.gaia_source.ra,gaiadr2.gaia_source.dec),CIRCLE('ICRS',%f,%f,%f))=1\
                               AND pmra IS NOT NULL AND abs(pmra)>0 \
    AND pmdec IS NOT NULL AND abs(pmdec)>0\
                               AND  phot_g_mean_mag<%f;"%(ra0,dec0,fov,limgmag) \
                                , dump_to_file=False)

    gaiat = job.get_results()
    ra = np.array([gaiat['ra']]).transpose()
    dec = np.array([gaiat['dec']]).transpose()
    pmra = np.array([gaiat['pmra']]).transpose()
    pmdec = np.array([gaiat['pmdec']]).transpose()

    gmag = np.array([gaiat['phot_g_mean_mag']]).transpose()
    bpmag = np.array([gaiat['phot_bp_mean_mag']]).transpose()
    rpmag = np.array([gaiat['phot_rp_mean_mag']]).transpose()

    eqpost = np.empty((0, 2), np.float_)
    gmagt = np.empty(0, np.float_)
    bpmagt = np.empty(0, np.float_)
    rpmagt = np.empty(0, np.float_)

    for k in range(np.size(gmag)):
        gmagt = np.append(gmagt, [gmag[k, 0]], axis=0)
        bpmagt = np.append(bpmagt, [bpmag[k, 0]], axis=0)
        rpmagt = np.append(rpmagt, [rpmag[k, 0]], axis=0)
        alpha, delta = RADecFromTang((yr - 2015.5) * math.pi * pmra[k, 0] / 648000000.0, \
                                          (yr - 2015.5) * math.pi * pmdec[k, 0] / 648000000.0, \
                                          math.radians(ra[k, 0]), math.radians(dec[k, 0]))
        eqpost = np.append(eqpost, [[math.degrees(alpha), math.degrees(delta)]], axis=0)

    return eqpost, gmagt, bpmagt, rpmagt

def jobGaiaDR2(ra0,dec0,limgmag,fov,yr):
    ssl._create_default_https_context = ssl._create_unverified_context
    job = Gaia.launch_job_async("SELECT * \
    FROM gaiadr2.gaia_source \
    WHERE CONTAINS(POINT('ICRS',gaiadr2.gaia_source.ra,gaiadr2.gaia_source.dec),CIRCLE('ICRS',%f,%f,%f))=1\
                               AND pmra IS NOT NULL AND abs(pmra)>0 \
    AND pmdec IS NOT NULL AND abs(pmdec)>0\
                               AND  phot_g_mean_mag<%f;"%(ra0,dec0,fov,limgmag) \
                                , dump_to_file=False)

    gaiat = job.get_results()
    # print(gaiat)
    ra = np.array(gaiat['ra'])
    dec = np.array(gaiat['dec'])
    pmra = np.array(gaiat['pmra'])
    pmdec = np.array(gaiat['pmdec'])

    gmag = np.array(gaiat['phot_g_mean_mag'])
    bpmag = np.array(gaiat['phot_bp_mean_mag'])
    rpmag = np.array(gaiat['phot_rp_mean_mag'])

    des = np.array(gaiat['designation'])

    N = np.size(gmag)

    for k in range(N):
        alpha, delta = RADecFromTang((yr - 2015.5) * math.pi * pmra[k] / 648000000.0, \
                                          (yr - 2015.5) * math.pi * pmdec[k] / 648000000.0, \
                                          math.radians(ra[k]), math.radians(dec[k]))
        print(ra[k],dec[k],math.degrees(alpha),math.degrees(delta))
        ra[k],dec[k] = math.degrees(alpha),math.degrees(delta)

    return ra,dec, gmag, bpmag, rpmag, des

# eqpost, gmagt, bpmagt, rpmagt = getGaiaDR2(100,23,12.5,0.2,2020.3)


ra,dec, gmag, bpmag, rpmag, des = jobGaiaDR2(100,23,14.5,0.2,1960.3)


print(gmag)
