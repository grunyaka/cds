# To calculate coordinates of stars from different catalogues in Gaia DR2 system of reference stars.
# Used cats: ucac4, urat1, cmc15, pan-starrs dr1, 2mass all-sky, sdss dr12.
# In case of Gaia ADQL is used. For other cats CDS interface (Vizier) is used. All of this realised with astroquerty.
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

def createFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print ('Error: Creating directory. ' +  directory)

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
    Q = int((n+1)*(n+2)/2)
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
        Zx = np.dot(np.linalg.inv(np.dot(np.dot(np.transpose(C), We), C)),\
                    np.dot(np.dot(np.transpose(C), We), ksi))
        Zy = np.dot(np.linalg.inv(np.dot(np.dot(np.transpose(C), We), C)),\
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
    return Zx,Zy,math.sqrt(uwex),math.sqrt(uwey),np.size(x)-it


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

#####################################################################################
foldout = os.path.join('./','out')
#filename = os.path.join('/home','agrippina','coding','current','source','list.txt')
filename = os.path.join('./','list.txt')
createFolder(foldout)
logfile = 'cat-to-gaia.txt'
new = True
new = False
lastRA = 254.39300256751923

#field size in arcdeg
fovlarge = 10./60.#10 arcmin
fovsmall = 10./3600.#10 arcsec
fov = fovlarge

limgmag = 18.0

with open(filename, 'rt') as f:
    file_content = f.read()
lines = file_content.split('\n')
print('stars in list.txt:', len(lines))

names = ['CMC15', 'Pan-STARRS DR1', '2MASS All-Sky', 'SDSS DR12', 'URAT1', 'UCAC4']
Nnames = ['1.CMC__15', '2.PanSta1', '3.2massAS', '4.SDSS_12', '5.URAT__1', '6.UCAC__4']
links = ['I/327/cmc15', 'II/349/ps1', 'II/246/out', 'V/147/sdss12', 'I/329/urat1', 'I/322A/out']
signs = ['RA_ICRS', 'DE_ICRS', 'e_RA_ICRS', 'e_DE_ICRS', 'MJD-51263', 'Jmag',
        'RAJ2000', 'DEJ2000', 'e_RAJ2000', 'e_DEJ2000', '_tab1_10', 'gmag',
        'RAJ2000', 'DEJ2000', 'errMaj', 'errMin', '_tab1_36', 'Jmag',  # ]
        'RA_ICRS', 'DE_ICRS', 'e_RA_ICRS', 'e_RA_ICRS', 'ObsDate', 'gmag',
        'RAJ2000', 'DEJ2000', 'sigs', 'sigm', 'Epoch', 'gmag',
        'RAJ2000', 'DEJ2000', 'e_RAJ2000', 'e_DEJ2000', 'EpRA', 'gmag']
logstring = 'RA [deg]'+'\t'+'Decl [deg]'+'\t'+'gmag [mag]'+'\t'+'PM_RA [mas/yr]'+'\t'+'PM_Decl [mas/yr]'

if new:
    with open(logfile, 'w') as f:  # 'a'
        f.write('# Gaia DR2' + '\n')
        for i in range(len(names)):
            f.write('# ' + str(names[i]) + '\n')
            logstring = logstring + '\t' + '|' + signs[i * 6] + ' [deg]' + '\t' + signs[i * 6 + 1] + ' [deg]' + '\t' + \
                        signs[i * 6 + 5] + ' [mag]' + '\t' + 'PM Ksi [mas]' + '\t' + 'PM Eta [mas]' + '\t' + 'F'
        f.write(logstring + '\n')

skip = True
for line in lines:
    idata = line.split(',')
    if (len(idata) > 1):
        bulstar = False
        ra0 = float(idata[5])
        er0 = float(idata[6])
        dec0 = float(idata[7])
        ed0 = float(idata[8])
        pmra0 = float(idata[12])
        epr0 = float(idata[13])
        pmdec0 = float(idata[14])
        epd0 = float(idata[15])
        gmag0 = float(idata[50])
        t0 = float(idata[4])
        logstring = str(ra0)+'\t'+str(dec0)+'\t'+str(gmag0)+'\t'+str(pmra0)+'\t'+str(pmdec0)
        outf = str(round(ra0,4))+str(round(dec0,4))
        #print(outf)
        outfile = os.path.join(foldout, outf+'.txt')
        with open(outfile,'w') as f:#'a'
            f.write( str(ra0)+' '+str(dec0)+' '+str(gmag0)+ '\n')
        #break
        #print('begin')
        if not new:
            if (ra0 == lastRA):
                skip = False
                continue
            if skip:
                continue
        print(ra0, dec0, gmag0)
        job = Gaia.launch_job_async("SELECT * \
        FROM gaiadr2.gaia_source \
        WHERE CONTAINS(POINT('ICRS',gaiadr2.gaia_source.ra,gaiadr2.gaia_source.dec),CIRCLE('ICRS',%f,%f,%f))=1\
                                   AND pmra IS NOT NULL AND abs(pmra)>0 \
        AND pmdec IS NOT NULL AND abs(pmdec)>0\
                                   AND  phot_g_mean_mag<%f;" % (ra0, dec0, fov, limgmag) \
                                    , dump_to_file=False)
        gaiat = job.get_results()
        ra_g = np.array(gaiat['ra'])
        dec_g = np.array(gaiat['dec'])
        era_g = np.array(gaiat['ra_error'])
        edec_g = np.array(gaiat['dec_error'])
        t_g = np.array(gaiat['ref_epoch'])
        pmra = np.array(gaiat['pmra'])
        pmdec = np.array(gaiat['pmdec'])
        gmag_g = np.array(gaiat['phot_g_mean_mag'])

        #print('finding the object in Gaia2 data array begin')
        Ng = 0
        for i in range(len(ra_g)):  # to find the object in Gaia2 data array
            r = ra_g[i]
            d = dec_g[i]
            g = gmag_g[i]
            if (ra0 == r) & (dec0 == d):
                #print(i, r, d, ra0, dec0)
                Ng = i
        #print('finding the object in Gaia2 data array finished')

       # Put the object to the begining of the list
        ra_g = np.insert(ra_g, 0, ra_g[Ng])
        ra_g = np.delete(ra_g, Ng+1)

        dec_g = np.insert(dec_g, 0, dec_g[Ng])
        dec_g = np.delete(dec_g, Ng + 1)

        era_g = np.insert(era_g, 0, era_g[Ng])
        era_g = np.delete(era_g, Ng + 1)

        edec_g = np.insert(edec_g, 0, edec_g[Ng])
        edec_g = np.delete(edec_g, Ng + 1)

        t_g = np.insert(t_g, 0, t_g[Ng])
        t_g = np.delete(t_g, Ng + 1)

        pmra = np.insert(pmra, 0, pmra[Ng])
        pmra = np.delete(pmra, Ng + 1)

        pmdec = np.insert(pmdec, 0, pmdec[Ng])
        pmdec = np.delete(pmdec, Ng + 1)

        gmag_g = np.insert(gmag_g, 0, gmag_g[Ng])
        gmag_g = np.delete(gmag_g, Ng + 1)

        Ng = 0 #we put the object to the beggining, so its index = 0 now
        with open(outfile, 'a') as f:
            f.write('Gaia index of the object: ' + str(Ng) + '\n')
            f.write('-----------------------------------------------------------------------------------------------' + '\n')
        #all the cats without ucac4 (it will be after all ##############################################################
        # --------------------------------------------------------------------------------------------------------------
        for icat in range(len(names)):
            object = True
            Ncatname = Nnames[icat]
            catname = names[icat]
            link = links[icat]
            #print(catname)
            test = True
            try:
                table = Vizier(columns=['all']).query_region(
                    coord.SkyCoord(ra=str(ra0), dec=str(dec0), unit=(u.deg, u.deg), frame='icrs'), Angle(fov, "deg"),
                    catalog=[link])[0]
            except IndexError:
                test = False
                print(Ncatname + ': connection attempt FAILED')
                logstring = logstring + '\t' + '|' + 'nan' + '\t' + 'nan' + '\t' + 'nan' + '\t' + 'nan' + '\t' + 'nan' + '\t' + 'nan'
                with open(outfile, 'a') as f:
                    f.write(catname + ': connection attempt FAILED' + '\n')
            #break
            if test:
                print(Ncatname + ': SUCCESS!  PREVIOUS NUMBER OF STARS IN FOV = ' + str(fov*3600) + ' arcsec: ', len(table))
                # break
                #if (catname == '2MASS All-Sky'):
                #    print(table.columns)
                # print(table.columns)
                # for column in table.columns:
                #    print(column)
                with open(outfile, 'a') as f:
                    f.write(catname + '\n')
                    f.write(
                        'Ncat' + '\t' + 'KSIcat [mas]   ' + '\t' + 'ETAcat [mas]   ' + '\t' + 'Ngaia' + '\t' + 'KSIgaia [mas]   ' + '\t' + 'ETAgaia [mas]   ' + '\t' + 'bias [mas]          ' + '\t' + 'mgaia' + '\t' + 'mcat' + '\t' + 'tgaia' + '\t' + 'tcat' + '\n')

                ra_c = np.array(table[signs[6 * icat]])  # dec
                dec_c = np.array(table[signs[6 * icat + 1]])  # dec
                era_c = np.array(table[signs[6 * icat + 2]], dtype=np.float64)  # mas
                edec_c = np.array(table[signs[6 * icat + 3]], dtype=np.float64)  # mas
                t_c = np.array(table[signs[6 * icat + 4]])  # mjd-51263 or mjd or jd or jy
                mag_c = np.array(table[signs[6 * icat + 5]])  # mag in branch g or j
                #era_c = np.array(era_c, dtype=np.float64)
                #edec_c = np.array(edec_c, dtype=np.float64)

                todelete = []
                if (catname == 'CMC15'):
                    #t_c = t_c + 51263  # mjd
                    for index in range(len(era_c)):
                        era_c[index] = era_c[index]*1000 #from arcsec to mas
                        edec_c[index] = edec_c[index] * 1000  # from arcsec to mas
                        try:
                            t_c[index] = Time((t_c[index]+ 51263), format='mjd').jyear
                        except:
                            todelete.append(index)

                if (catname == 'Pan-STARRS DR1'):
                    for index in range(len(era_c)):
                        era_c[index] = era_c[index] * 1000  # from arcsec to mas
                        edec_c[index] = edec_c[index] * 1000  # from arcsec to mas
                        try:
                            t_c[index] = Time(t_c[index], format='mjd').jyear
                        except:
                            todelete.append(index)

                if (catname == '2MASS All-Sky'):
                    # print(Time(t_c[0], format='jd').iso)
                    # print(np.array(table['Date'])[0])
                    errPA = np.array(table['errPA'])
                    for index in range(len(era_c)):
                        era_c[index] = era_c[index] * 1000  # from arcsec to mas
                        edec_c[index] = edec_c[index] * 1000  # from arcsec to mas
                        a = era_c[index]
                        b = edec_c[index]
                        era_c[index] = abs(a * b) / (math.sqrt(
                            (a ** 2) * (math.cos(math.radians(errPA[index])) ** 2) +
                            (b ** 2) * (math.sin(math.radians(errPA[index])) ** 2)))
                        edec_c[index] = abs(a * b) / (math.sqrt(
                            (a ** 2) * (math.sin(math.radians(errPA[index])) ** 2) +
                            (b ** 2) * (math.cos(math.radians(errPA[index])) ** 2)))
                        try:
                            t_c[index] = Time(t_c[index], format='jd').jyear
                        except:
                            todelete.append(index)

                if (catname == 'SDSS DR12'):
                    for index in range(len(era_c)):
                        era_c[index] = era_c[index] * 1000  # from arcsec to mas
                        edec_c[index] = edec_c[index] * 1000  # from arcsec to mas

                if (catname == 'URAT1'):
                    #era_c = np.array(era_c, dtype=np.float64)
                    #edec_c = np.array(edec_c, dtype=np.float64)
                    #print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@", type(era_c))
                    for index in range(len(era_c)):
                        # print('URAT1 pos errors:', era_c[index], edec_c[index])
                        try:
                            era_c[index] = math.sqrt(era_c[index] * era_c[index] + edec_c[index] * edec_c[index])
                            edec_c[index] = era_c[index]
                            # print('calc  pos errors:', era_c[index], edec_c[index])
                        except:
                            print('calc urat errors FAILED', era_c[index], edec_c[index])

                if (catname == 'UCAC4'):
                    #print('tra_c', 'tdec_c', 't_c')
                    tra_c = t_c
                    tdec_c = np.array(table['EpDE'])  # yr
                    t_c = (tra_c + tdec_c) / 2.
                    #for index in range(len(t_c)):
                    #    print(tra_c[index], tdec_c[index], t_c[index])

                np.delete(ra_c, todelete)
                np.delete(dec_c, todelete)
                np.delete(era_c, todelete)
                np.delete(edec_c, todelete)
                np.delete(t_c, todelete)
                np.delete(mag_c, todelete)

                odata = []
                odatanum = 0  # number of ref.stars
                indexdone = []  # to remember indexes of cat stars which have benn already found,
                            # to remove them from next serachings
                #print(era_c)
                #print(np.square(era_c))
                R3sigma = 3 * math.sqrt((np.sum(np.square(era_c)))/len(era_c))
                #print(edec_c)
                #print(np.square(edec_c))
                D3sigma = 3 * math.sqrt((np.sum(np.square(edec_c))) / len(edec_c))
                print('3sigmas:', R3sigma, D3sigma)

                #sigma3 = 3 * math.sqrt((er ** 2 + ed ** 2 + erg ** 2 + edg ** 2) / 4.)

                for j in range(len(ra_g)):
                    rg = ra_g[j] # [arcdeg]
                    dg = dec_g[j] # [arcdeg]
                    erg = era_g[j] # [mas]
                    edg = edec_g[j] # [mas]
                    tg = t_g[j] # [yr]
                    mg = gmag_g[j] # [mag]
                    pmr = pmra[j]  # [mas/yr]
                    pmd = pmdec[j] # [mas/yr]
                    catind = [i for i in list(range(len(ra_c))) if not i in indexdone]
                    star = 'nan'
                    dmmin = 100
                    indexj = 0
                    for i in catind:
                        r = ra_c[i]  # [arcdeg]
                        d = dec_c[i]  # [arcdeg]
                        er = era_c[i] # [mas]
                        ed = edec_c[i] # [mas]
                        t = t_c[i]  # [yr]
                        m = mag_c[i]  # [mag]
                        #sigma3 = 3 * math.sqrt((er ** 2 + ed ** 2 + erg ** 2 + edg ** 2) / 4.)

                        try:
                            dm = abs(m - mg)
                        except:
                            dm = 99

                        k_c, e_c = tangFromRADE(r * math.pi / 180., d * math.pi / 180.,
                                                ra0 * math.pi / 180., dec0 * math.pi / 180.)  # [rad]
                        # to take into account pm
                        alpha, delta = RADecFromTang((t - tg) * math.pi * pmr / 648000000.0, \
                                                     (t - tg) * math.pi * pmd / 648000000.0, \
                                                     math.radians(rg), math.radians(dg))
                        #rr, dd = math.degrees(alpha), math.degrees(delta)

                        k_g, e_g = tangFromRADE(alpha, delta, ra0 * math.pi / 180.0,
                                                dec0 * math.pi / 180.0)  # [radians]
                        biasK = abs(k_c - k_g)* 648000000.0 / math.pi  # [mas]
                        biasE = abs(e_c - e_g) * 648000000.0 / math.pi  # [mas]
                        #if ((biasK < R3sigma)&(biasE < D3sigma)):
                        if ((biasK < 10000) & (biasE < 10000)):
                            if star == 'nan':
                                star = [j, k_g * 648000000.0 / math.pi, e_g * 648000000.0 / math.pi,
                                        i, k_c * 648000000.0 / math.pi, e_c * 648000000.0 / math.pi,
                                        mg, m, tg, t, r, d]
                                indexi = i
                            elif (dm < dmmin):
                                dmmin = dm
                                star = [j, k_g * 648000000.0 / math.pi, e_g * 648000000.0 / math.pi,
                                        i, k_c * 648000000.0 / math.pi, e_c * 648000000.0 / math.pi,
                                        mg, m, tg, t, r, d]
                                indexi = i
                    if star == 'nan':
                        if (j == 0):
                            object = False
                            #print(Ncatname + ': ERROR. The object is not found.')
                            logstring = logstring + '\t' + '|' + 'nan' + '\t' + 'nan' + '\t' + 'nan' + '\t' + 'nan' + '\t' + 'nan' +'\t' + 'nan'
                            with open(outfile, 'a') as f:
                                f.write('ERROR. The object is not found.')
                            break
                    else:
                        indexdone.append(indexi)
                        odata.append(star)
                        odatanum += 1
                        string = ''
                        for e in star:
                            string = string + str(e) + '\t'
                        with open(outfile, 'a') as f:
                            f.write(string + '\n')
                if object:
                    #print('number of stars before calibration: ', odatanum)
                    with open(outfile, 'a') as f:
                        f.write('Number of stars before calibration: ' + str(odatanum) + '\n')
                    #Kg = odata[0][1]
                    #Eg = odata[0][2]
                    #Kg, Eg = tangFromRADE(ra0, ra0, ra0 * math.pi / 180.0,
                    #                        dec0 * math.pi / 180.0)
                    Kc = odata[0][4]
                    Ec = odata[0][5]
                    t = odata[0][8]
                    tg = odata[0][9]
                    #KG = 0
                    #EG = 0

                    n = 1  # order of reduction model
                    Nmin = 10  # minimal number of reference stars
                    rmax = 100  # maximal residual in [mas]
                    with open(outfile, 'a') as f:
                        f.write('Settings of the model:' + '\n' + 'order of the model (n) = ' + str(
                            n) + '\n' + 'min num of ref stars (Nmin) = ' + str(
                            Nmin) + '\n' + 'max residual (rmax) = ' + str(
                            rmax) + ' [mas]' + '\n')

                    kc = np.empty(odatanum - 1)  # minus 1 because we should exclude the object
                    ec = np.empty(odatanum - 1)
                    kg = np.empty(odatanum - 1)
                    eg = np.empty(odatanum - 1)
                    i = 0
                    for star in odata:
                        if not (i == 0):
                            kc[i-1] = star[1]
                            ec[i-1] = star[2]
                            kg[i-1] = star[4]
                            eg[i-1] = star[5]
                        i += 1

                    try:
                        Zk, Ze, ek, ee, Ns = calibration(kg, eg, kc, ec, n, Nmin, rmax)  # [mas]
                    except:
                        with open(outfile, 'a') as f:
                            f.write('CALIBRATION FAILED' + '\n')
                        logstring = logstring + '\t' + '|' + 'nan' + '\t' + 'nan' + '\t' + 'nan' + '\t' + 'nan' + '\t' + 'nan' + '\t' + 'nan'
                        continue
                    bulstar = True
                    #print('number of stars after calibration:', Ns)
                    with open(outfile, 'a') as f:
                        f.write('Number of stars after calibration: ' + str(Ns) + '\n')
                    print(Zk, Ze, ek, ee, Ns)
                    Zkstring = ''
                    for e in Zk:
                        Zkstring = Zkstring + str(e) + '\t'
                    Zestring = ''
                    for e in Ze:
                        Zestring = Zestring + str(e) + '\t'
                    with open(outfile, 'a') as f:
                        f.write(
                            '\n' + 'Parameters of the model:' + '\n' + Zkstring + '\n' + Zestring + '\n' + 'Errors:' + str(
                                ek) + ' [mas], ' + str(ee) + ' [mas]' + '\n')
                    KG = transform(0, 0, Zk)
                    EG = transform(0, 0, Ze)
                    #print(KG, EG)
                    muk = (KG - Kc) / (tg - t)
                    mue = (EG - Ec) / (tg - t)
                    print('PM:')
                    print(muk, pmra0)
                    print(mue, pmdec0)
                    dmuk = (muk - pmra0)
                    dmue = (mue - pmdec0)

                    ek2 = sqrerror(er0, ek, tg, t)
                    ee2 = sqrerror(ed0, ee, tg, t)
                    print(math.sqrt(ek2), math.sqrt(ee2))

                    F = math.sqrt(dmuk ** 2 / ek2 + dmue ** 2 / ee2)

                    logstring = logstring + '\t' + '|' + str(odata[0][10]) + '\t' + str(odata[0][11]) + '\t' + str(
                        odata[0][7]) + '\t' + str(muk) + '\t' + str(mue) + '\t' + str(F)

                    #print(ra0, dec0, 'UCAC4')
                    #print('Gaia PM errors: ', er0, ed0)
                    #print('Calc residuals: ', ek, ee)
                    #print('Epoch diff: ', tg - t)
                    #print('PM erroros: ', math.sqrt(ek2), math.sqrt(ee2))
                    #print('Gaia2 PM [mas]: ', pmra0, pmdec0)
                    #print('Calc  PM [mas]: ', muk, mue)
                    #print('F = ', F)
                    with open(outfile, 'a') as f:
                        f.write('Gaia PM errors: ' + str(er0) + ' ' + str(ed0) + '\n')
                        f.write('Calc residuals: ' + str(ek) + ' ' + str(ee) + '\n')
                        f.write('Epoch diff: ' + str(tg - t) + '\n')
                        f.write('PM erroros: ' + str(math.sqrt(ek2)) + ' ' + str(math.sqrt(ee2)) + '\n')
                        f.write('Gaia2 PM [mas]: ' + str(pmra0) + ' ' + str(pmdec0) + '\n')
                        f.write('Calc PM [mas]: ' + str(muk) + ' ' + str(mue) + '\n')
                        f.write('F: ' + str(F) + '\n')
                with open(outfile, 'a') as f:
                    f.write(
                        '-----------------------------------------------------------------------------------------------' + '\n')
        if bulstar:
            with open(logfile, 'a') as f:
                f.write(logstring + '\n')
            #break
        #break
        ################################################################################################################
    #break
else:
    print('The Gaia cycle is finally completed.')