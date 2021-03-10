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
import sys

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


#####################################################################################
filename = os.path.join('./','list.txt')
logfile = 'findObject.txt'
log = 'findObject_log.txt'

new = True
#new = False # Это на случай, если программа почему-то прервется, тогда тут нужно будет раскомментить,
lastRA = 292.81970856588407#а тут написать Gaia RA, на котором прервалось

#field size
fov = 1./60.#15 arcmin в arcdeg
epsilon = 10000 #mas, это кружок, чтобы понять, ту ли мы нашли звезду или нет,
                # я сперва пробовала 3сигму считать, но она какая-то сильно маленькая получается, неск сот mas
limgmag = 21.0



names = ['CMC15', 'Pan-STARRS DR1', '2MASS All-Sky', 'SDSS DR12', 'URAT1', 'UCAC4']
Nnames = ['0.CMC__15', '2.PanSta1', '3.2massAS', '4.SDSS_12', '5.URAT__1', '6.UCAC__4']#это только чтобы логи на экран выводить, в теории можно убрать, если убрать print'ы
links = ['I/327/cmc15', 'II/349/ps1', 'II/246/out', 'V/147/sdss12', 'I/329/urat1', 'I/322A/out']
signs = ['RA_ICRS', 'DE_ICRS', 'e_RA_ICRS', 'e_DE_ICRS', 'MJD-51263', 'Jmag',
        'RAJ2000', 'DEJ2000', 'e_RAJ2000', 'e_DEJ2000', '_tab1_10', 'gmag',
        'RAJ2000', 'DEJ2000', 'errMaj', 'errMin', '_tab1_36', 'Jmag',  # _tab1_36 - JD
        'RA_ICRS', 'DE_ICRS', 'e_RA_ICRS', 'e_RA_ICRS', 'ObsDate', 'gmag',
        'RAJ2000', 'DEJ2000', 'sigs', 'sigm', 'Epoch', 'gmag',
        'RAJ2000', 'DEJ2000', 'e_RAJ2000', 'e_DEJ2000', 'EpRA', 'gmag']

if new:#делаю шапки для выходных файлов
    with open (log, 'w') as f:
        f.write('RA[deg]          \tDecl [deg]          \t|CMC15  |PanSt1 |2MASSAS|SDSS12 |URAT1  |UCAC4  |total\n')
    with open(logfile, 'w') as f:  # 'a'
        #f.write('Gaia' + 14 * '\t' + '|CMC15' + 8 * '\t' + '|Pan-Starrs DR1' + 7 * '\t' + '|2MASS All-Sky' + 7 * '\t' + '|SDSS DR12' + 7 * '\t' + '|URAT1' + 8 * '\t' + '|UCAC4' + 8 * '\t' + '\n')
        f.write('Gaia\t\t\t\t\t\t\t\t\t\t\t\t\t\t|CMC15\t\t\t\t\t\t\t\t|Pan-Starrs DR1\t\t\t\t\t\t\t|2MASS All-Sky\t\t\t\t\t\t\t|SDSS DR12\t\t\t\t\t\t\t|URAT1\t\t\t\t\t\t\t\t|UCAC4\n')
        logstring = 'RA[deg]          \tDecl [deg]          \tmag [mag]\tPM_RA [mas/yr]     \tPM_Decl [mas/yr]\t'
        for i in range(len(names)):
            logstring = logstring+'|RA[deg]          \tDecl [deg]          \tmag [mag]\t'
        f.write(logstring + '\n')
skip = True#это чтобы пропускать уже пройденные зезды, если мы запускаем прерванную программу, пока не дойдём до той, на которой остановились

with open(filename, 'rt') as f:#это читаем входной файл по каждой звезде
    file_content = f.read()
lines = file_content.split('\n')
#на период отладки сделаем список поменьше, и из тех, что точно есть в ucac4

#print('lines in list.txt:', len(lines))
for line in lines:#берем по одной звезде из списка Gaia
    idata = line.split(',')
    if (len(idata) > 1):
        object = False#есть ли звезда хотя бы в одном каталоге (кроме Gaia)
        cf = 0#счетчик для подсчета, в скольких каталогах найден объект
        ra0 = float(idata[5])
        er0 = float(idata[6])
        dec0 = float(idata[7])
        ed0 = float(idata[8])
        pmra0 = float(idata[12])
        pmdec0 = float(idata[14])
        gmag0 = float(idata[50])
        t0 = float(idata[4])
        logstring = '{:19.15f}\t{:19.15f}\t{:9.6f}\t{:21.15f}\t{:21.15f}\t'.format(ra0, dec0, gmag0, pmra0, pmdec0)#начало строки для Файла findOblect.txt
        logs = '{:19.15f}\t{:19.15f}\t{:9.6f}\t'.format(ra0, dec0, gmag0)#начало строки для файла findObject_log.txt
        if not new:#если запускаем прерванную программу, нужно пропустить уже посчитанные звезды
            if (ra0 == lastRA):
                skip = False
                #continue
            if skip:
                continue
        print(ra0, dec0, gmag0)
        for icat, catname in enumerate(names):
            #print(icat, catname)
            object_cat = False#есть ли звезда в этом каталоге
            Ncatname = Nnames[icat]
            link = links[icat]
            test = False#это если в каталоге не будет нужной области неба
            disconnect = True#это если сервер страсбурга будет недоступен
            while disconnect:
                disconnect = False
                fail = 0
                try:
                    test = True
                    table = Vizier(columns=['all']).query_region(
                        coord.SkyCoord(ra=str(ra0), dec=str(dec0), unit=(u.deg, u.deg), frame='icrs'), Angle(fov, "deg"),
                        catalog=[link])[0]
                except IndexError:
                    #print(IndexError)
                    test = False
                    fail = 1
                except:# not IndexError:
                    test = False
                    print("Unexpected error:", sys.exc_info()[0])
                    #print(IndexError)
                    #print ('.', end='')
                    disconnect = True
            if not test:
                print(Ncatname + ': connection attempt FAILED')
                logstring = logstring + '|' + 'nan' + 14 * ' ' + '\t' + 'nan' + 14 * ' ' + '\t' + 'nan' + 5 * ' ' + '\t'
                logs = logs + '| fail{} '.format(fail)
            else:
                print(Ncatname + ': CONNECTED!')

                if (catname == '2MASS All-Sky'):# это нужно, так как в 2mass иногда меняется название столбца с датой
                    for index, column in enumerate(table.columns):
                        if index == 36:
                            signs[16] = column

                ra_c = np.array(table[signs[6 * icat]])  # dec
                dec_c = np.array(table[signs[6 * icat + 1]])  # dec
                t_c = np.array(table[signs[6 * icat + 4]])  # mjd-51263 or mjd or jd or jy
                mag_c = np.array(table[signs[6 * icat + 5]])  # mag in branch g or j

                todelete = []#это для отбраковки, если с датой наблюдения что-то не то
                if (catname == 'CMC15'):
                    for index in range(len(t_c)):
                        try:
                            t_c[index] = Time((t_c[index]+ 51263), format='mjd').jyear
                        except:
                            todelete.append(index)

                if (catname == 'Pan-STARRS DR1'):
                    for index in range(len(t_c)):
                        try:
                            t_c[index] = Time(t_c[index], format='mjd').jyear
                        except:
                            todelete.append(index)

                if (catname == '2MASS All-Sky'):
                    for index in range(len(t_c)):
                        try:
                            t_c[index] = Time(t_c[index], format='jd').jyear
                        except:
                            todelete.append(index)

                if (catname == 'UCAC4'):#здесь моменты наблюдения для ra и dec даны отдельно, так что я просто беру среднее
                    tra_c = t_c
                    tdec_c = np.array(table['EpDE'])  # yr
                    t_c = (tra_c + tdec_c) / 2.

                np.delete(ra_c, todelete)
                np.delete(dec_c, todelete)
                np.delete(t_c, todelete)
                np.delete(mag_c, todelete)

                rg = ra0  # [arcdeg]
                dg = dec0  # [arcdeg]
                erg = er0  # [mas]
                edg = ed0  # [mas]
                tg = t0  # [yr]
                mg = gmag0  # [mag]
                pmr = pmra0  # [mas/yr]
                pmd = pmdec0  # [mas/yr]
                star = 'nan'
                dmmin = 100
                for i, r in enumerate(ra_c):
                    d = dec_c[i]  # [arcdeg]
                    t = t_c[i]  # [yr]
                    m = mag_c[i]  # [mag]
                    #считаем разности зв величин, если в малом кружке epsilon будет более 1 звезды, иногда mag = 'nan'
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
                    k_g, e_g = tangFromRADE(alpha, delta, ra0 * math.pi / 180.0,
                                            dec0 * math.pi / 180.0)  # [radians]
                    # считаем разность положений, чтобы понять, может ли звезда из каталога быть Объектом
                    biasK = abs(k_c - k_g) * 648000000.0 / math.pi  # [mas]
                    biasE = abs(e_c - e_g) * 648000000.0 / math.pi  # [mas]
                    if ((biasK < epsilon) & (biasE < epsilon)):#epsilon = 10000 mas = 10 arcsec
                        object = True#звезда есть хотя бы в одном каталоге
                        object_cat = True# звезда есть в этом каталоге
                        if star == 'nan':
                            star = [r, d, m]
                            #indexi = i
                        elif (dm < dmmin):
                            dmmin = dm
                            star = [r, d, m]
                            #indexi = i
                if object_cat:#дополняем строчки для выходных файлов
                    logstring = logstring + '|{}\t{:19.15f}\t{:19.15f}\t'.format(star[0], star[1], star[2])# для findfObject.txt
                    logs = logs+'|success'# для findfObject_log.txt
                    cf+=1
                else:# это в случае если ни одной звезды из каталога не было в маленьком epsilon-кружке
                    logstring = logstring + '|' + 'nan' + 14 * ' ' + '\t' + 'nan' + 14 * ' ' + '\t' + 'nan' + 5 * ' ' + '\t'
                    logs = logs + '|  NF   '
        if object:#если звезда нашлась хотя бы в одном ещё каталоге, то записываем её в файл findfObject.txt
            with open(logfile, 'a') as f:
                f.write(logstring + '\n')
        with open(log, 'a') as f:# в файл findObject_log.txt информацию о звезде записываем в любом случае
            f.write(logs + '|{}\n'.format(cf))
print('The Gaia cycle is finally completed.')