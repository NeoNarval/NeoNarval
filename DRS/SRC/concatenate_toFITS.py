import cPickle
import pyfits
import sys
import utils
from datetime import datetime
from scipy import signal


def filter_spectrum(spectrum):

    normal_cutoff = 0.53
    b, a = signal.butter(3, normal_cutoff, btype='low', analog=False)
    y = signal.lfilter(b, a, spectrum)

    return y


# ARGUMENT 1 : path du fichier traite comme il est dans FILES
# ARGUMENT 2 : path du flat comme il est dans FILES
# ARGUMENT 3 : normalise

path = sys.argv[1]
# path_flat_comp = sys.argv[2]
tableau = []
cfg = utils.CfgFile("../DATA/Amatrix_DATA.txt")

for lane in [1, 2]:
    list_lane = []
    for ordre in xrange(40):

        try:
            normalized = int(sys.argv[3])
        except:
            normalized = 1
        if normalized:
            path_comp = "../TEMP/SN_" + path.split("/")[-1][:-4] \
                        + "_OR" + str(ordre) + "_LA" + str(lane) + ".p"  # chemin image
        else:
            path_comp = "../TEMP/SP_" + path.split("/")[-1][:-4] \
                        + "_OR" + str(ordre) + "_LA" + str(lane) + ".p"  # chemin image
        
        try:
            star_lane = cPickle.load(open(path_comp, 'rb'))
            star_lane_new = [0] * len(star_lane)
            star_lane_new = star_lane
            star_lane_new_filtered = star_lane_new  # filter_spectrum(star_lane_new) A DECOMMENTER POUR FILTRER
            list_lane += list(star_lane_new_filtered)
        except IOError:
            print(path_comp)
            list_lane += [-2] * 4612  # Quelle taille prendre ?..
            print "concatenate: probleme avec l'ordre %d et la voie %d" % (ordre, lane)
    long_lane = range(len(list_lane))
    try:
        
        true_tabl = pyfits.open('../DATA/th_calibered.fits')[1].data.field(0)
        #long_lane = true_tabl  #Spectre en arturo si commentee
    except IOError:
        print "concatenate: spectre du thorium manquant !, spectre en Arturo"
    tableau.append(long_lane)
    tableau.append(list_lane)

col1 = pyfits.Column(name='wavelength_lane1', format='E', array=tableau[0])
col2 = pyfits.Column(name='intensity_lane1', format='E', array=tableau[1])
col3 = pyfits.Column(name='wavelength_lane2', format='E', array=tableau[2])
col4 = pyfits.Column(name='intensity_lane2', format='E', array=tableau[3])

cols = pyfits.ColDefs([col1, col2, col3, col4])
tbhdu2 = pyfits.BinTableHDU.from_columns(cols)

tbhdu1 = pyfits.PrimaryHDU(header=pyfits.open(path)[0].header)
tbhdu1.header["FP_FILE"] = cfg.get_param("Fabry Perot fts file").split("/").pop()
star_file = cfg.get_param("Test fts file").split("/").pop()
tbhdu1.header["STARFILE"] = star_file
final_HDU = pyfits.HDUList([tbhdu1, tbhdu2])

temps = datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%S.%f")

path = '../TEMP/spectrum_' + star_file[:-4] + '.fits'
final_HDU.writeto(path, clobber=True)
path = '/home/main/Documents/DRS/REDUCED/' + star_file[:-4] + '.fits'
final_HDU.writeto(path, clobber=True)
print "concatenate: SPECTRE ECRIT A " + path
