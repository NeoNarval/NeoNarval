from datetime import datetime
import pyfits
import sys
import os.path as path
# Script Python pour renommer les flats en fonction
# de leur date reelle dans le header
nom_fichier = path.abspath(sys.argv[1])


hdulist = pyfits.open(nom_fichier)
header = hdulist[0].header

date = datetime.strptime(header["DATE"], "%Y-%m-%dT%H:%M:%S.%f")

# chemin_split = chemin.split("/")
# old_name = chemin_split.pop().split("_")

type_laurent = path.basename(nom_fichier).split(".")[0][-1]
true_type = "null"

if type_laurent == "b":
    true_type = "bia"
elif type_laurent == "f":
    true_type = "fla"
elif type_laurent == "c":
    true_type = "th0"
elif type_laurent == "o":
    true_type = "st0"
elif type_laurent == "a":
    true_type = "fp0"

new_name = "/Narval" + date.strftime("_%Y%m%d_%H%M%S_") + true_type + ".fts"
new_path = path.dirname(nom_fichier) + new_name

hdulist.writeto(new_path)
hdulist.close()

print new_path
