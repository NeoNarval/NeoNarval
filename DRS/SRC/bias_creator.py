import numpy as np
import methods as meth

num_ligne = 4096
num_colonnes = 4096

image = 300 * np.ones((num_ligne, num_colonnes))
meth.ecriture_FITS(image, "00", "bia")

