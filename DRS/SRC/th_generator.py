from datetime import datetime, timedelta
import pyfits

temps = datetime.utcnow()
name = "../Brut/simu" + "00"+ temps.strftime("_%Y%m%d_%H%M%S_") + "th0.fts"

hdulist = pyfits.open("../DATA/230164c.fits")
header = hdulist[0].header
header["DATE"] = temps.strftime("%Y-%m-%dT%H:%M:%S.%f")
header["TIMEEND"] = (temps + timedelta(minutes=1)).strftime("%H:%M:%S.%f")

hdulist.writeto(name)
hdulist.close()