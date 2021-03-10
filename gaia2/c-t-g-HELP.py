# i=0
# for g in gaiat.columns:
#    print (i, g)
#    i+=1
# print(gaiat.columns)

#to transform "Gaia2 equatorial coords -pm"  to "Gaia2 tangential coords", the optical center is coords of the Gaia2 object
print('Gaia RA:', ra_g)
print('PM RA:',pmra)
dt = tra_c - t_g
print('dt:', dt)
ksi_g = (ra_g*3600000 + pmra * (tra_c - t_g))/3600000
print('Gaia RA-PM', ksi_g)
exit()
eta_g = dec_g + pmdec * (tdec_c - t_g)/3600000
i=0
for r in ra_g:
    d = dec_g[i]
    ksi_g[i], eta_g[i]=tangFromRADE(r,d,ra0,dec0)
    #print(i,') ',ksi_g[i],',',eta_g[i],';',r,',',d)
    i+=1
exit()
i = 0
for r in ra_g:  # to find the object in Gaia2 data array
    d = dec_g[i]
    g = gmag_g[i]
    if (ra0 == r) & (dec0 == d):
print(i, r, d, ra0, dec0)
Ng = i
break
    i += 1

# to transform Gaia2 equatorial coords to Gaia2 tangential coords, the optical center is coords of the Gaia2 object
ksi_g = ra_g
eta_g = dec_g
i = 0
for r in ra_g:
    d = dec_g[i]
    ksi_g[i], eta_g[i] = tangFromRADE(r, d, ra0, dec0)
    # print(i,') ',ksi_g[i],',',eta_g[i],';',r,',',d)
    i += 1
exit()

Ns = []
catbiases = []
biaslim = 5#arcsec
print('proverka N: ', len(Ns))
i=0
cat = False
print('GAIA: ', ra0, dec0, gmag0)
for k in ksi_c:
    e = eta_c[i]
    tr = tra_c[i]
    td = tdec_c[i]
    gmag = gmag_c[i]
    biasksi = abs(pmra0 * (tr - t0) - 3600000 * k)
    biaseta = abs(pmdec0 * (td - t0) - 3600000 * e)
    bias = math.sqrt(biasksi*biasksi+biaseta*biaseta)/1000.
    if (bias < biaslim):
        Ns.append(str(i))
        catbiases.append(str(bias))
        catbias = bias
        Nc[N] = i
print(catname,': ', r, d, gmag)
print('bias [ arcsec]: ', catbias)
    i+=1
print('proverka N: ', len(Ns))
if(len(Ns)>1):
    mindmag = 10
    j = 0
    for i in int(Ns):
dmag = abs(gmag0 - gmag_c[i])
if (mindmag > dmag):
    nimgmag = dmag
    Nc[N] = i
    catbias = float(catbiases[j])
j+=1

if (catbias > 0.6):
    print('Error of current object is to large in this cat')
    print('GAIA: ', ra0, dec0, gmag0)
    print(catname, ': ', r, d, gmag)
    print('bias [ arcsec]: ', catbias)
    print('PM [mas]:', pmra0, pmdec0)
    print('epoch differences:', t0-tr, t0-td)

else:
    print()
exit()
g_u4_k = []
g_u4_e = []
u4_g_k = []
u4_g_e = []
u4_tk = []
u4_te = []

i = 0
for r in ra_g:
    gmag = gmag_g[i]
    d = dec_g[i]
    pmr = pmra[i]
    pmd = pmdec[i]
    t = t_g[i]
    print(r, d, gmag)
    if not (r == ra0) & (d == dec0) & (gmag0 == gmag):
        j = 0
        for rr in ra_u4:
            dd = dec_u4[j]
            ttrr = tra_u4[j]
            ttdd = tdec_u4[j]
            err = era_u4[j]
            edd = edec_u4[j]
    print('ucac4 data:')
    print(rr, dd)
    print(err, edd)
    if (abs((r - rr) * 3600000 + pmr * (ttrr - t)) < 3 * err):
        if (abs((d - dd) * 3600000 + pmd * (ttdd - t)) < 3 * edd):
            if not (r == ra0) & (d == dec0):
                g_u4_ra.append(str(r*3600000+pmr(ttrr-t)))
                g_u4_dec.append(str(d*3600000+pmd(ttdd-t)))
                u4_g_ra.append(str(rr*3600000))
                u4_g_dec.append(str(dd*3600000))
                u4_tr.append(str(ttrr))
                u4_td.append(str(ttdd))
                print('urah')
                print(r, d)
                print(rr, dd)
    else:
print('object is found')
g_ra0 = r * 3600000
g_dec0 = d * 3600000
u4_ra0 = rr * 3600000
u4_dec0 = dd * 3600000
    j += 1
    else:
print('OBJECT FOUND')
    i+=1

# coord = SkyCoord(ra=ra0, dec=dec0, unit=(u.degree, u.degree), frame='icrs')
# print(coord)
# width = u.Quantity(0.01, u.deg)
# height = u.Quantity(0.01, u.deg)
# object = Gaia.query_object_async(coordinate=coord, width=width, height=height)
# object.pprint()
# print(object.columns)
# print(object['phot_g_mean_mag'])
# exit()
