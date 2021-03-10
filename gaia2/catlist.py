from astroquery.vizier import Vizier
from astropy import coordinates as coord
from astropy.coordinates import Angle
from astropy import units as u


catalog_list = Vizier.find_catalogs('Kang W51')
print({k:v.description for k,v in catalog_list.items()})
print(catalog_list)


catalogs = Vizier.get_catalogs(catalog_list.keys())
print(catalogs)


fovsec = 10# field size in arcsec
fov = fovsec/3600.# field size in deg

table = Vizier(columns=['all']).query_region(coord.SkyCoord(ra=str(r), dec=str(d), unit=(u.deg, u.deg), frame='icrs'), Angle(fov, "deg"), catalog=[link])[0]