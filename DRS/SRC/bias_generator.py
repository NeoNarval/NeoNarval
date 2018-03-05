from datetime import datetime, timedelta
import pyfits
import sys

try:
    date_image = sys.argv[1]
    temps = datetime.strptime(date_image, "%Y-%m-%dT%H:%M:%S.%f")
except:
    temps = datetime.utcnow()


name = "../Brut/simu" + "00"+ temps.strftime("_%Y%m%d_%H%M%S_") + "bia.fts"

hdulist = pyfits.open("../DATA/223110b.fits")
header = hdulist[0].header
header["DATE"] = temps.strftime("%Y-%m-%dT%H:%M:%S.%f")
header["TIMEEND"] = (temps + timedelta(minutes=1)).strftime("%H:%M:%S.%f")

hdulist.writeto(name)
hdulist.close()
