import numpy as np
import pickle
import pyfits

h=pyfits.open("../DATA/luna_calibered.fits")
wave=h[1].data.field(0)
h.close()

h=pyfits.open("../TEMP/th_spectrum.fits")
inti=h[1].data.field(1)
h.close()

# We create the new fits file
col1 = pyfits.Column(name='wavelength_lane1', format='E', array=wave)
col2 = pyfits.Column(name='intensity_lane1', format='E', array=inti)
# col3 = pyfits.Column(name='wavelength_lane2', format='E', array=tableau[2])
# col4 = pyfits.Column(name='intensity_lane2', format='E', array=tableau[3])

cols = pyfits.ColDefs([col1, col2]) # , col3, col4])
tbhdu2 = pyfits.BinTableHDU.from_columns(cols)

tbhdu1 = pyfits.PrimaryHDU(header=pyfits.open("../DATA/luna_calibered.fits")[0].header)
#tbhdu1.header["FP_FILE"] = cfg.get_param("Fabry Perot fts file").split("/").pop()
#tbhdu1.header["STAR_FILE"] = cfg.get_param("Test fts file").split("/").pop()

final_HDU = pyfits.HDUList([tbhdu1, tbhdu2])

final_HDU.writeto('../DATA/th_calibered.fits', clobber=True) # changer le nom du fichier pour mettre le nouveau

print "Reussi !"
