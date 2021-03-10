import pandas as pd
import numpy as np
import os
import healpy as hp
from astroquery.vizier import Vizier
from astroquery.gaia import Gaia
#import astropy.coordinates as coord
from astropy.coordinates import SkyCoord
from astropy.io import ascii
import astropy.units as u
from astropy.table import Table
from astropy.time import Time
import math


def GaiaTable(gaia_ra, gaia_dec, starname):
    global radsec
    global info
    global limmag

    fov = radsec/3600

    try:
        job = Gaia.launch_job_async("SELECT * FROM gaiadr2.gaia_source \
            WHERE CONTAINS(POINT('ICRS',gaiadr2.gaia_source.ra,gaiadr2.gaia_source.dec), \
            CIRCLE('ICRS',%f,%f,%f))=1 AND  phot_g_mean_mag<%f \
            AND pmra IS NOT NULL AND abs(pmra)>0 AND pmdec IS NOT NULL AND abs(pmdec)>0;" \
                                    % (gaia_ra, gaia_dec, fov, limmag), dump_to_file=False)

        gaiat = job.get_results()[
            ['source_id'] + list(info.columns)[:-2] + ['pmra', 'pmdec', 'pmra_error', 'pmdec_error']]
        gaiat['phot_g_mean_mag'].name = 'mag'

        df = gaiat.to_pandas()
        df = df.set_index('source_id')
        try:
            df = df.drop(index=starname)
            return df, True
        except:
            return df, True
    except Exception as e:
        print('Gaia \t' + ' ERROR:', e)
        return ['NaN', False]


def RenameColumns(table, catalog):
    global info

    old_columns = info.drop(['link', 'pecularity'], axis=1)
    new_columns = info.drop(['link', 'pecularity'], axis=1)
    new_columns = new_columns.rename(columns={'phot_g_mean_mag': 'mag'})

    old_columns = old_columns.loc[catalog].tolist()
    new_columns = new_columns.columns.tolist()

    for i in range(len(new_columns)):
        table = table.rename(columns={old_columns[i]: new_columns[i]})

    return table


def errors_2MASS(row):
    a = row['errMaj']
    b = row['errMin']
    angle = row['errPA']

    row['errMaj'] = abs(a * b) / (math.sqrt(
        (a ** 2) * (math.cos(math.radians(angle)) ** 2) +
        (b ** 2) * (math.sin(math.radians(angle)) ** 2)))
    row['errMin'] = abs(a * b) / (math.sqrt(
        (a ** 2) * (math.sin(math.radians(angle)) ** 2) +
        (b ** 2) * (math.cos(math.radians(angle)) ** 2)))

    return row


def epoch_UCAC4(row):
    a = row['EpRA']
    b = row['EpDE']

    row['EpRA'] = (a + b) / 2

    return row


def CatTable(cat_ra, cat_dec, catalog):
    global info
    global radsec

    try:
        link = info.loc[catalog]['link']
        viz = Vizier(columns=['all'])
        table = viz.query_region(SkyCoord(
            ra=cat_ra, dec=cat_dec, unit=(u.deg, u.deg), frame='icrs'),
            radius=radsec * u.arcsec, catalog=link)[0]

        df = pd.DataFrame.from_records(table.as_array())

        pecularities = info.pecularity[catalog]
        pecs = pecularities.split('|')

        if (not pecs[0] == 'mas'):

            column1 = info['ra_error'][catalog]
            column2 = info['dec_error'][catalog]

            df[column1] = df[column1].apply(lambda x: x * 1000)
            df[column2] = df[column2].apply(lambda x: x * 1000)

            if len(pecs[0].split(':')) == 2:
                df[['errMaj', 'errMin', 'errPA']] = \
                    df[['errMaj', 'errMin', 'errPA']].apply(errors_2MASS, axis=1)

        if (not pecs[1] == 'jy'):
            ps = pecs[1].split(':')

            if len(ps) == 2:
                column = info['ref_epoch'][catalog]

                df[column] = \
                    df[column].apply(lambda x: Time((x + int(ps[0])), format=ps[1]).jyear)

            else:
                df[['EpRA', 'EpDE']] = \
                    df[['EpRA', 'EpDE']].apply(epoch_UCAC4, axis=1)

        cols = info.loc[catalog].tolist()[:-2]
        df = df[['_r'] + cols]
        df = RenameColumns(df, catalog)
        df._r = round(df._r, 4)
        df.ref_epoch = round(df.ref_epoch, 1)
        df = df.sort_values(by=['_r'])
        #print(catalog)  # + ': success')

        return [df, True]

    except Exception as e:
        print(catalog + '\t' + ' ERROR:', e)
        return ['NaN', False]


def tangFromRADE(ra, de, RA, DE):
    ksi = math.cos(de)*math.sin(ra-RA)/(math.sin(de)*math.sin(DE)+math.cos(de)*math.cos(DE)*math.cos(ra-RA))
    eta = (math.sin(de)*math.cos(DE)-math.cos(de)*math.sin(DE)*math.cos(ra-RA))/(math.sin(de)*math.sin(DE)+math.cos(de)*math.cos(DE)*math.cos(ra-RA))
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


def calibration(x,y,ksi,eta,n,Nmin,rmax):
    try:
        Q = int((n + 1) * (n + 2) / 2)
        C = np.ones((np.size(x), Q))
        q = 0
        for i in range(n + 1):
            for j in range(n + 1):
                if (i + j <= n):
                    C[:, q] = (x ** i) * (y ** j)
                    q += 1
        We = np.diag(np.ones(np.size(x)))
        flag = 0
        it = 0
        while (flag == 0):
            Zx = np.dot(np.linalg.inv(np.dot(np.dot(np.transpose(C), We), C)), \
                        np.dot(np.dot(np.transpose(C), We), ksi))
            Zy = np.dot(np.linalg.inv(np.dot(np.dot(np.transpose(C), We), C)), \
                        np.dot(np.dot(np.transpose(C), We), eta))
            rx = np.dot(We, ksi - np.dot(C, Zx))
            ry = np.dot(We, eta - np.dot(C, Zy))
            r = np.sqrt(rx ** 2 + ry ** 2)
            kmax = np.argmax(r)
            flag = 1
            if (np.size(ksi) - it <= Nmin):
                break
            if (r[kmax] > rmax):
                We[kmax, kmax] = 0
                flag = 0
                it += 1
        uwex = np.dot(np.dot(np.transpose(rx), We), rx) / (np.size(x) - it - Q)
        uwey = np.dot(np.dot(np.transpose(ry), We), ry) / (np.size(y) - it - Q)
        return Zx, Zy, math.sqrt(uwex), math.sqrt(uwey), np.size(x) - it, True
    except Exception as e:
        print('Calibration \t' + ' ERROR:', e)
        return 0, 0, 0, 0, 0, False


def transform(x,y,Z):
    Q = np.size(Z)
    n = int(np.max(np.roots(np.array([1,3,-2*(Q-1)]))))
    C = np.ones((np.size(x), Q))
    q = 0
    for i in range(n+1):
        for j in range(n+1):
            if (i + j <= n):
                C[:, q] = (x ** i) * (y ** j)
                q += 1
    return np.dot(C,Z)


def sqrerror(e0, e, tg, tc):
    e2 = e / (tg - tc)
    e2 = e2**2 + e0**2
    #e - error from calibration, e0 - pm error from gaia, tg - gaia epoch, tc - cat epoch
    return e2


def GudeTableGaiaCat(star, gaiadf, catdf):
    gudedf = pd.DataFrame(columns=['kg', 'eg', 'kc', 'ec'])
    # tg -- t_new
    # t -- t_
    tg = 2015.5
    for gaia_id in list(gaiadf.index):
        gstar = gaiadf.loc[gaia_id]
        pmr = gstar.pmra
        pmd = gstar.pmdec
        bias = 10000
        found = False
        for cat_id in list(catdf.index):
            cstar = catdf.loc[cat_id]
            t = cstar.ref_epoch
            k_c, e_c = tangFromRADE(cstar.ra * math.pi / 180., cstar.dec * math.pi / 180.,
                                    star.gaia_ra * math.pi / 180., star.gaia_dec * math.pi / 180.)  # [rad]
            # to take into account pm
            alpha, delta = RADecFromTang((t - tg) * math.pi * pmr / 648000000.0,
                                         (t - tg) * math.pi * pmd / 648000000.0,
                                         math.radians(gstar.ra), math.radians(gstar.dec))
            # rr, dd = math.degrees(alpha), math.degrees(delta)
            k_g, e_g = tangFromRADE(alpha, delta, star.gaia_ra * math.pi / 180.0,
                                    star.gaia_dec * math.pi / 180.0)  # [radians]
            biasK = abs(k_c - k_g) * 648000000.0 / math.pi  # [mas]
            biasE = abs(e_c - e_g) * 648000000.0 / math.pi  # [mas]
            current_bias = math.sqrt(biasE ** 2 + biasK ** 2)
            if current_bias < bias:
                bias = current_bias
                found = True
                drop_id = cat_id
                gude_star = [k_g * 648000000.0 / math.pi, e_g * 648000000.0 / math.pi,
                             k_c * 648000000.0 / math.pi, e_c * 648000000.0 / math.pi]
        if found:
            gudedf.loc[gaia_id] = gude_star
            catdf = catdf.drop(index=drop_id)
    return gudedf


def MU(ra_new, ra_old, dec_new, dec_old, t_new, t_old, Zk, Ze):
    #ra_new, ra_old, dec_new, dec_old, t_new, t_old
    '''New proper motion calculation.'''
    Kc, Ec = tangFromRADE(ra_old * math.pi / 180., dec_old * math.pi / 180.,
                          ra_new * math.pi / 180., dec_new * math.pi / 180.)  # [rad] Catalog tang coords.
    Kg = transform(0, 0, Zk)[0]  # Gaia tang coords
    Eg = transform(0, 0, Ze)[0]

    return (Kg - Kc) / (t_new - t_old), (Eg - Ec) / (t_new - t_old)


def F (star, muk, mue, ek, ee):
    '''Calcilation of statistical parameter F'''
    dmuk = (muk - star.pmra)
    dmue = (mue - star.pmdec)
    ek2 = sqrerror(star.gaia_ra_error, ek, 2015.5, star.cat_epoch)
    ee2 = sqrerror(star.gaia_dec_error, ee, 2015.5, star.cat_epoch)
    return math.sqrt(dmuk ** 2 / ek2 + dmue ** 2 / ee2)



'''Global variables'''
radsec = 600
min_delta_mag = 10
limmag = 18
n = 1  # order of reduction he
Nmin = 10  # minimal number of reference stars
rmax = 100  # maximal residual in [mas]
info = pd.read_csv('cat_info.csv').fillna('NaN')
outdir = 'gaia3F_results'
indir = 'catalog_positions_mag'

files = os.listdir(indir)
print(files)
#['Pan-STARRS DR1.csv', 'URAT1.csv', '2MASS All-Sky.csv', 'SDSS DR12.csv', 'CMC15.csv', 'UCAC4.csv']
#source_id,gaia_ra,gaia_dec,gaia_ra_error,gaia_dec_error,pmra,pmdec,pmra_error,pmdec_error,gaia_epoch,cat_ra,cat_dec,cat_epoch,delta_mag,dist_sec

'''Temporary variables for debbuging'''

file = files[5]
catalog = file.split('.')[0]
inpath = os.path.join(indir, file)
df = pd.read_csv(inpath, index_col = 0)
outpath = os.path.join(outdir,catalog+'.csv')

check = os.path.exists(outpath)

if check:
    result = pd.read_csv(outpath, index_col=0)
    done_indexes = list(result.index)
    del result
else:
    result = pd.DataFrame(columns=['gaia_ra', 'gaia_dec', 'cat_ra', 'cat_dec', 'dt', 'dm', 'dist_sec', 'gaia_pmra', 'gaia_pmdec', 'new_pmksi', 'new_pmeta', 'F'])
    result.to_csv(outpath)
    done_indexes = []
    del result


i = 0
for starid in df.index:
    i+=1

    if starid in done_indexes:
        continue

    print(i, 'of', len(df.index))
    print(starid, 'in prosses')

    star = df.loc[starid]  # the object information with spherical coords

    '''Getting a DataFrame with Gaia spherical coords of gude stars.'''
    gaiadf, check = GaiaTable(star.gaia_ra, star.gaia_dec, star.name)
    if not check:
        continue

    '''Getting a DataFrame with the Catalog spherical coords of gude stars.'''
    catdf, check = CatTable(star.cat_ra, star.cat_dec, catalog)
    if not check:
        continue

    '''Compilation of Catalog gude stars DataFrame and Gaia gude stars Data frame 
    to Calculate a common DataFrame with gude stars tangentional coords from both cats.'''
    gudedf = GudeTableGaiaCat(star, gaiadf, catdf)
    kg = gudedf.kg.to_numpy()# ksi_gaia - first tangentional coords of gude stars from Gaia
    eg = gudedf.eg.to_numpy()# eta_gaia - second tangentional coords of gude stars from Gaia
    kc = gudedf.kc.to_numpy()# ksi_catalog - first tangentional coords of gude stars from Catalog
    ec = gudedf.ec.to_numpy()# eta_catalog - second tangentional coords of gude stars from Catalog

    '''Calibration inputs: Gaia and Catalog tang coors of gude stars, see 288-290 lines.
    Calibration outputs: 
    Zk, Ze - calibration vectors to calculate calibrated gaia tang coors
    ek, ee - residuals, Ns - number of gude stars after calibration
    check - to check if calibration was successful.'''
    Zk, Ze, ek, ee, Ns, check = calibration(kg, eg, kc, ec, n, Nmin, rmax)  # [mas]
    if not check:
        continue

    print('success!')
    '''Calculating of new proper motions by the object information (star) and by the calibration vectors (Zk, Ze).'''
    muk, mue = MU(star.gaia_ra, star.cat_ra,
                  star.gaia_dec, star.cat_dec,
                  2015.5, star.cat_epoch, Zk, Ze)

    f = F(star, muk, mue, ek, ee)

    star_line = [star.gaia_ra, star.gaia_dec, star.cat_ra, star.cat_dec, 2015.5 - star.cat_epoch,
                 star.delta_mag, star.dist_sec, star.pmra, star.pmdec, muk, mue, f]

    result = pd.read_csv(outpath, index_col=0)
    result.loc[star.name] = star_line
    result.to_csv(outpath)
    del result
    #break# for debugging
print('all prosseses are finally finished')




