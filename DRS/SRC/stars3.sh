#!/bin/bash

################################################################################
#### Ce script s'occupe d'une étoile ou d'un thorium-argon en particulier ######
################################################################################

. functions.sh
DATA_FILE="../DATA/Amatrix_DATA.txt"
################################################################################
###################### GESTION DES FICHIERS LOCK ###############################
################################################################################

# Afin de gérer le lancement de plusieurs étoiles en même temps, mais de tout de même
# les traiter les unes après les autres, il y a un système de plusieurs fichiers lock :
# chaque étoile a son fichier lock, qu'elle supprime une fois finie

last="./.lock0"			   # On initialise last, qui représente le fichier lock de l'étoile d'avant
last=$(ls -Art .star* 2>/dev/null | tail -n 1) # last = dernier fichier de type .star* du dossier
newfile=./.star$(printf "%02d" $(( 10#${last:5:2} + 1 )))  # creation du nouveau fichier incrementé (max = 99)
echo lock >> $newfile

trap "rm $newfile" EXIT # pour supprimer le lock au cas ou il est killed

if [ -e $last ]; then
    echo "\
~~~~~~~~~~~~~~~~~~~~~~~~ ETOILE $1 EN ATTENTE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    inotifywait -qq -e delete $last  # quand l'étoile précédente est terminée, on part
fi
echo "\
################################################################################
################ PRETRAITEMENT ETOILE OU THORIUM ARGON #########################
################################################################################"

flat_file=$(get_param "Flat fts file" $DATA_FILE)

if [[ $1 == "st" ]]; then
    type="st"
    shift
    stars=${@}
    python Search_star.py $stars
    
elif [[ $1 == "th" ]]; then
    type="th"
    shift
    th=${@}
    echo "Readying $th"
    python Search_thar.py $th
    
elif [[ $1 == "fla" ]]; then
    type="fla"
    set_param "Test fts file" $flat_file $DATA_FILE
fi

echo "\
################################################################################
########### CALCUL DU SPECTRE 1, ORDRE PAR ORDRE, VOIE PAR VOIE ################
################################################################################"

DOSSIER="..\/TEMP"
fts_file=$(get_param "Test fts file" $DATA_FILE)

if [[ $type == "st" ]]; then
    sp2_file=$(echo $fts_file | sed "s/st1/st2/")
elif [[ $type == "th" ]]; then
    sp2_file=$(echo $fts_file | sed "s/th1/th2/")
elif [[ $type == "fla" ]]; then
    sp2_file=$(echo $fts_file | sed "s/10f/11f/")
fi
echo "Extracting $fts_file"
get_spectrum $fts_file

norm=0

if [[ $type != "fla" ]]; then
    echo "~~~ Normalisation du spectre  avec le flat $flat_file ~~~"
    python normalization.py $fts_file $flat_file
    if [[ $? == 0 ]]; then
	echo "Normalisation finie"
    echo "Files in TEMP have changed name:  SP -> SN"
	norm=1
    fi
fi


echo "~~~~~~~~~ Création de SP2 ~~~~~~~~~~~~~"

echo norm = $norm
python cosmic_auto.py $sp2_file $norm

get_spectrum $sp2_file

# if [[ $SPECTRE_DONE == $(($NBR_ORDRE_TOTAL * $VOIES)) ]]; then
#     echo "Toutes les matrices sont réussies !"
# fi

if [[ $type != "fla" ]]; then
    echo "~~~ Normalisation des spectre du SP ~~~"
    python normalization.py $sp2_file $flat_file
    if [[ $? == 0 ]]; then
	echo "Normalisation finie"
	norm=1
    fi
fi

# CONCATENATION DES DIFFERENTS SPECTRES OBTENUS

if [[ $1 != "fla" ]]; then
   #ALApython concatenate_toFITS.py $sp2_file $flat_file
python concatenate_toFITS.py $fts_file $flat_file 1
#ALA At this point $flat_file is of no use. The 3rd argument is 0 if un_normalized. This is only used to identify the rifht file name
fi


echo "ETOILE TERMINEE"

rm $newfile

