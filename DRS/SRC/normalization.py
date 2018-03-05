import cPickle
import numpy as np
import sys
import os.path
# ARGUMENT 1 = path du fichier a normaliser
# ARRGUMENT 2 = path du fichier flat

path_fts_file = os.path.abspath(sys.argv[1])
path_flat_file = os.path.abspath(sys.argv[2])

fts_name = os.path.basename(path_fts_file)
flat_name = os.path.basename(path_flat_file)
NBR_ORDRE = 40
NBR_VOIES = 2

for ordre in xrange(NBR_ORDRE):
    for voie in xrange(NBR_VOIES):
        spectre_fts = "../TEMP/SP_" + fts_name[:-4] + "_OR{}_LA{}.p".format(ordre, voie + 1)
        spectre_flat = "../TEMP/SP_" + flat_name[:-4] + "_OR{}_LA{}.p".format(ordre, voie + 1)
        
        fichier_fts = cPickle.load(open(spectre_fts, 'rb'))
        fichier_flat = cPickle.load(open(spectre_flat, 'rb'))
        
        fichier_flat = fichier_flat / abs(np.mean(fichier_flat))
        # fichier_flat = fichier_flat / np.max(fichier_flat)
        
        star_lane = np.zeros(len(fichier_flat))
        for i, value in enumerate(fichier_flat):
            if value > 0:
                star_lane[i] = fichier_fts[i] / value
            else:
                star_lane[i] = -1
                
        new_path_fts = (os.path.dirname(spectre_fts)
                        + "/SN"
                        + os.path.basename(spectre_fts)[2:])

        cPickle.dump(star_lane, open(new_path_fts, 'wb'))
