from astropy import coordinates as coord
from astropy import units as u
from astroquery.vizier import Vizier
from astropy.coordinates import Angle
import os


filename = os.path.join('./','list.txt')
outfile = 'find_gaiastar_incats.txt'
fovmin = 10#field size in arcmin
fovsec = 120#field size in arcsec
fov = fovmin/60.#field size in arcdeg
fov = fovsec/3600.
fovlarge = 10./60.
fovsmall = 10./3600.
fov = fovsmall
limgmag = 18.0
with open(filename, 'rt') as f:
    file_content = f.read()
lines = file_content.split('\n')
skip = True
for line in lines:
    idata = line.split(',')
    if (len(idata) > 1):
        r0 = float(idata[5])
        d0 = float(idata[7])
        gmag0 = float(idata[50])
        t0 = float(idata[4])
        skip = False
        if (r0 == 290.80842276227105):
            skip = False
        if skip:
            continue
        print( r0, d0, gmag0)
        with open(outfile, 'a') as f:  # 'a'
            f.write(str(r0) + ' ' + str(d0) + ' ' + str(gmag0) + '\n')
        names = ['CMC15', 'Pan-STARRS DR1', '2MASS All-Sky', 'SDSS DR12', 'URAT1', 'UCAC4']
        Nnames = ['1.CMC__15', '2.PanSta1', '3.2massAS', '4.SDSS_12', '5.URAT__1', '6.UCAC__4']
        links = ['I/327/cmc15', 'II/349/ps1', 'II/246/out', 'V/147/sdss12', 'I/329/urat1', 'I/322A/out']
        signs = ['RA_ICRS', 'DE_ICRS', 'e_RA_ICRS', 'e_DE_ICRS', 'MJD-51263', 'Jmag',
                 'RAJ2000', 'DEJ2000', 'e_RAJ2000', 'e_DEJ2000', '_tab1_10', 'gmag',
                 'RAJ2000', 'DEJ2000', 'errMaj', 'errMin', '_tab1_36', 'Jmag',  # ]
                 'RA_ICRS', 'DE_ICRS', 'e_RA_ICRS', 'e_RA_ICRS', 'ObsDate', 'gmag',
                 'RAJ2000', 'DEJ2000', 'sigs', 'sigm', 'Epoch', 'gmag',
                 'RAJ2000', 'DEJ2000', 'e_RAJ2000', 'e_DEJ2000', 'EpRA', 'gmag']
        for icat in range(len(names)):
            Ncatname = Nnames[icat]
            catname = names[icat]
            link = links[icat]
            test = True
            try:
                table = Vizier(columns=['all']).query_region(
                    coord.SkyCoord(ra=str(r0), dec=str(d0), unit=(u.deg, u.deg), frame='icrs'), Angle(fov, "deg"),
                    catalog=[link])[0]
                print(Ncatname + ': SUCCESS!  NUMBER OF STARS IN FOV=' + str(fov*3600) + ' arcsec: ', str(len(table)))

                with open(outfile, 'a') as f:  # 'a'
                    f.write(Ncatname + ': SUCCESS!  NUMBER OF STARS IN FOV=' + str(fov*3600) + ' arcsec: ' + str(len(table)) + '\n')
                    for star in table:
                        ra = star[signs[6 * icat]]
                        dec = star[signs[6 * icat + 1]]
                        mag = star[signs[6 * icat + 5]]
                        f.write(str(ra) + ' ' + str(dec) + ' ' + str(mag) + '\n')

            except IndexError:
                test = False
                print(Ncatname + ': connection attempt FAILED')
                with open(outfile, 'a') as f:  # 'a'
                    f.write(Ncatname + ': connection attempt FAILED' + '\n')