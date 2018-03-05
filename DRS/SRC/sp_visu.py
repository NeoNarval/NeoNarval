import cPickle
import matplotlib.pyplot as plt
import subprocess
import pyfits

cmd = 'zenity --file-selection --multiple --filename=/home/main/Documents/DRS/TEMP/ --file-filter="S*.p *.fits" 2>/dev/null'

file = subprocess.check_output(cmd, shell=True)[:-1]
files = file.split("|")
for file in files:
    if file[-1] == "p":
        spectre = cPickle.load(open(file, 'rb'))
        print file.split("/").pop()
        plt.figure()
        plt.plot(spectre)
    elif file[-1] == "s":
        fichier = pyfits.open(file)[1].data
        X = fichier.field(0)
        Y = fichier.field(1)
        plt.figure()
        plt.plot(X, Y)

plt.show()
