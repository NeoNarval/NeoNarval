#!/bin/bash

################################################################################
######### Script servant à effectuer le prétraitement d'une nuit ###############
# Paramètres :                                                                 #
# - v : nombre de voies                                                        #
# - f : nombre de flat nécessaires                                             #
# - b : nombre de bias nécessaires                                             #
# - F : nombre de fabry-pérot nécessaires                                      #
################################################################################

. functions.sh 
NB_FLAT=10
NB_BIAS=1
NB_FP=1

FLAT_DONE=0
BIAS_DONE=0
FP_DONE=0

files="../DATA/test.txt"

while getopts v:f:b:F: OPTIONS; do
    case $OPTIONS in
        v)
            VOIES=$OPTARG;;
        f)
            NB_FLAT=$OPTARG;;
        b)
            NB_BIAS=$OPTARG;;
        F)
            NB_FP=$OPTARG;;
    esac
done

rm -f .reussi

if [ -e ".lock" ]
then
    echo "Fichier lock présent, le script s'arrête."
    exit 0
fi

echo "lock file" >> .lock

trap "rm .lock" EXIT		# Au cas où le script se termine, on supprime .lock

echo "\
################################################################################
################# ATTENTE DES FICHIERS DE PRETRAITEMENT ########################
################################################################################"


while [[ $FLAT_DONE != $NB_FLAT || $BIAS_DONE != $NB_BIAS || $FP_DONE != $NB_FP ]]
do
    FICHIER=$(inotifywait -q -e create -e moved_to --format %f ../Brut/.)

    if [[ $FICHIER =~ fla\.fts$ ]]; then
	# add_to_list FLAT "$FICHIER" $files (add_to_list pas necessaires)
	(( FLAT_DONE++ ))
	echo "Flats available !" $FLAT_DONE
    fi
    if [[ $FICHIER =~ bia\.fts$ ]]; then
	echo "Bias available !"
       	# add_to_list BIAS "$FICHIER" $files
        (( BIAS_DONE++ ))
    fi
    if [[ $FICHIER =~ fp0\.fts$ ]]; then
	echo "Fabry-Perot available"
	# add_to_list FABRY "$FICHIER" $files 
        (( FP_DONE++ ))
    fi
done
rm .lock



echo "\
########################################
########## PRETRAITEMENT ###############
########################################"

sleep 5s
python Search_bias.py

echo "Bias handled by arthur2"


sleep 5s
python Search_flat.py
echo "Flat handled by arthur2"


sleep 5s
python Search_fp.py
echo "Fabry-Perot handled by arthur2"


sleep 5s
python envelope.py $VOIES



set_param "Lanes per order" $VOIES ../DATA/Amatrix_DATA.txt
set_param "Lanes per order" $VOIES ../DATA/find_fp_slits_DATA.txt 

sleep 5s
python find_fp_slits.py 2>/dev/null

echo "Boundaries and FP lines found by arthur2"

sleep 5s


echo "\
########################################
####### CREATION DES MATRICES A1 #######
########################################"

set_param "Evolution level" 1 ../DATA/Amatrix_DATA.txt
fp_file=$(get_param "Fabry Perot fts file" ../DATA/Amatrix_DATA.txt)

create_A_matrix $fp_file
echo "Fichiers Amatrix_v1_*_*.npz generees dans TEMP"
echo "\
########################################
####### CREATION DES MATRICES A2 #######
########################################"

echo "~ Création des spectres du Fabry-Pérot avec les matrices A1 ~"
get_spectrum $fp_file
echo "Fichiers SP_ONE_*_*.p generees dans TEMP"




echo "~ Création des spectres de l'image one ~"
one_file="../DATA/one.fts"
flat_file=$(get_param "Flat fts file" ../DATA/Amatrix_DATA.txt)
get_spectrum $one_file

echo "~~~ Normalisation des spectre du FP ~~~"

python normalization.py $fp_file $one_file
if [[ $? == 0 ]]; then
    echo "Normalisation finie"
fi

echo "~~~~~~~~~ Création du FP2 ~~~~~~~~~~~~~"
fp_file=$(get_param "Fabry Perot fts file" ../DATA/Amatrix_DATA.txt)
set_param "Test fts file" $fp_file ../DATA/Amatrix_DATA.txt
fp2_file=$(echo $fp_file | sed "s/fp1/fp2/")
python cosmic_auto.py $fp2_file 1

echo "~~~~~ Création de la matrice A2 ~~~~~~~"
set_param "Evolution level" 2 ../DATA/Amatrix_DATA.txt
create_A_matrix $fp2_file 


echo "\
########################################
###### CALCUL DU SPECTRE DU FLAT #######
########################################"

./stars3.sh fla	


# if [[ $MATRIX_DONE == $(($NBR_ORDRE_TOTAL * $VOIES)) ]]; then
echo reussi >> .reussi
# fi
